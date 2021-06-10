/*
OFDM modem encoder

Copyright 2021 Ahmet Inan <inan@aicodix.de>
*/

#include <iostream>
#include <cassert>
#include <cmath>
#include "complex.hh"
#include "utils.hh"
#include "decibel.hh"
#include "fft.hh"
#include "wav.hh"
#include "pcm.hh"
#include "mls.hh"
#include "crc.hh"
#include "psk.hh"
#include "galois_field.hh"
#include "bose_chaudhuri_hocquenghem_encoder.hh"

template <typename value, typename cmplx, int rate>
struct Encoder
{
	typedef PhaseShiftKeying<8, cmplx, value> Mod;
	static const int symbol_len = (1280 * rate) / 8000;
	static const int guard_len = symbol_len / 8;
	static const int data_bits = 64800;
	static const int data_cols = 432;
	static const int data_rows = data_bits / data_cols / Mod::BITS;
	static const int mls0_len = 127;
	static const int mls0_poly = 0b10001001;
	static const int mls1_len = 255;
	static const int mls1_poly = 0b100101011;
	static const int mls2_poly = 0b100101010001;
	DSP::WritePCM<value> *pcm;
	DSP::FastFourierTransform<symbol_len, cmplx, 1> bwd;
	CODE::CRC<uint16_t> crc;
	CODE::BoseChaudhuriHocquenghemEncoder<255, 71> bchenc;
	cmplx fdom[symbol_len];
	cmplx tdom[symbol_len];
	cmplx guard[guard_len];
	cmplx papr_min, papr_max;
	int data_off;
	int mls0_off;
	int mls1_off;

	static int bin(int carrier)
	{
		return (carrier + symbol_len) % symbol_len;
	}
	void symbol()
	{
		bwd(tdom, fdom);
		for (int i = 0; i < symbol_len; ++i)
			tdom[i] /= sqrt(value(8 * symbol_len));
		for (int i = 0; i < guard_len; ++i) {
			value x = value(i) / value(guard_len - 1);
			x = value(0.5) * (value(1) - std::cos(DSP::Const<value>::Pi() * x));
			guard[i] = DSP::lerp(guard[i], tdom[i+symbol_len-guard_len], x);
		}
		cmplx peak, mean;
		for (int i = 0; i < symbol_len; ++i) {
			cmplx power(tdom[i].real() * tdom[i].real(), tdom[i].imag() * tdom[i].imag());
			peak = cmplx(std::max(peak.real(), power.real()), std::max(peak.imag(), power.imag()));
			mean += power;
		}
		if (mean.real() > 0 && mean.imag() > 0) {
			cmplx papr(peak.real() / mean.real(), peak.imag() / mean.imag());
			papr *= value(symbol_len);
			papr_min = cmplx(std::min(papr_min.real(), papr.real()), std::min(papr_min.imag(), papr.imag()));
			papr_max = cmplx(std::max(papr_max.real(), papr.real()), std::max(papr_max.imag(), papr.imag()));
		}
		pcm->write(reinterpret_cast<value *>(guard), guard_len, 2);
		pcm->write(reinterpret_cast<value *>(tdom), symbol_len, 2);
		for (int i = 0; i < guard_len; ++i)
			guard[i] = tdom[i];
	}
	void pilot_block()
	{
		CODE::MLS seq2(mls2_poly);
		value data_fac = sqrt(value(symbol_len) / value(data_cols));
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] = 0;
		for (int i = data_off; i < data_off + data_cols; ++i) {
			value tmp[Mod::BITS];
			for (int k = 0; k < Mod::BITS; ++k)
				tmp[k] = 1 - 2 * seq2();
			fdom[bin(i)] = data_fac * Mod::map(tmp);
		}
		symbol();
	}
	void schmidl_cox()
	{
		CODE::MLS seq0(mls0_poly);
		value mls0_fac = sqrt(value(symbol_len) / value(mls0_len));
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] = 0;
		fdom[bin(mls0_off-2)] = mls0_fac;
		for (int i = 0; i < mls0_len; ++i)
			fdom[bin(2*i+mls0_off)] = (1 - 2 * seq0());
		for (int i = 0; i < mls0_len; ++i)
			fdom[bin(2*i+mls0_off)] *= fdom[bin(2*(i-1)+mls0_off)];
		symbol();
	}
	void meta_data(uint64_t md)
	{
		uint8_t data[9] = { 0 }, parity[23] = { 0 };
		for (int i = 0; i < 55; ++i)
			CODE::set_be_bit(data, i, (md>>i)&1);
		crc.reset();
		uint16_t cs = crc(md << 9);
		for (int i = 0; i < 16; ++i)
			CODE::set_be_bit(data, i+55, (cs>>i)&1);
		bchenc(data, parity);
		CODE::MLS seq4(mls1_poly);
		value mls1_fac = sqrt(value(symbol_len) / value(mls1_len));
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] = 0;
		fdom[bin(mls1_off-1)] = mls1_fac;
		for (int i = 0; i < 71; ++i)
			fdom[bin(i+mls1_off)] = (1 - 2 * CODE::get_be_bit(data, i));
		for (int i = 71; i < mls1_len; ++i)
			fdom[bin(i+mls1_off)] = (1 - 2 * CODE::get_be_bit(parity, i-71));
		for (int i = 0; i < mls1_len; ++i)
			fdom[bin(i+mls1_off)] *= fdom[bin(i-1+mls1_off)];
		for (int i = 0; i < mls1_len; ++i)
			fdom[bin(i+mls1_off)] *= (1 - 2 * seq4());
		symbol();
	}
	Encoder(DSP::WritePCM<value> *pcm, uint8_t *inp, int freq_off, uint64_t call_sign) :
		pcm(pcm), crc(0xA8F4), bchenc({
			0b100011101, 0b101110111, 0b111110011, 0b101101001,
			0b110111101, 0b111100111, 0b100101011, 0b111010111,
			0b000010011, 0b101100101, 0b110001011, 0b101100011,
			0b100011011, 0b100111111, 0b110001101, 0b100101101,
			0b101011111, 0b111111001, 0b111000011, 0b100111001,
			0b110101001, 0b000011111, 0b110000111, 0b110110001})
	{
		data_off = (freq_off * symbol_len) / rate - data_cols / 2;
		mls0_off = data_off + 90;
		mls1_off = data_off + 89;
		papr_min = cmplx(1000, 1000), papr_max = cmplx(-1000, -1000);
		pilot_block();
		schmidl_cox();
		meta_data((call_sign << 8) | 2);
		pilot_block();
		for (int j = 0; j < data_rows; ++j) {
			for (int i = 0; i < data_cols; ++i) {
				value tmp[Mod::BITS];
				for (int k = 0; k < Mod::BITS; ++k) {
					int l = Mod::BITS * (data_cols * j + i) + k;
					tmp[k] = 1 - 2 * ((inp[l/8] >> (l % 8)) & 1);
				}
				fdom[bin(i+data_off)] *= Mod::map(tmp);
			}
			symbol();
		}
		schmidl_cox();
		meta_data((call_sign << 8) | 2);
		pilot_block();
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] = 0;
		symbol();
		std::cerr << "real PAPR: " << DSP::decibel(papr_min.real()) << " .. " << DSP::decibel(papr_max.real()) << " dB" << std::endl;
		if (pcm->channels() == 2)
			std::cerr << "imag PAPR: " << DSP::decibel(papr_min.imag()) << " .. " << DSP::decibel(papr_max.imag()) << " dB" << std::endl;
	}
};

long long int base37_encoder(const char *str)
{
	long long int acc = 0;
	for (char c = *str++; c; c = *str++) {
		acc *= 37;
		if (c >= '0' && c <= '9')
			acc += c - '0' + 1;
		else if (c >= 'a' && c <= 'z')
			acc += c - 'a' + 11;
		else if (c >= 'A' && c <= 'Z')
			acc += c - 'A' + 11;
		else if (c != ' ')
			return -1;
	}
	return acc;
}

int main(int argc, char **argv)
{
	if (argc < 6 || argc > 8) {
		std::cerr << "usage: " << argv[0] << " OUTPUT RATE BITS CHANNELS INPUT [OFFSET] [CALLSIGN]" << std::endl;
		return 1;
	}

	const char *output_name = argv[1];
	int output_rate = std::atoi(argv[2]);
	int output_bits = std::atoi(argv[3]);
	int output_chan = std::atoi(argv[4]);
	const char *input_name = argv[5];

	int freq_off = output_chan == 1 ? 2000 : 0;
	if (argc >= 7)
		freq_off = std::atoi(argv[6]);

	if ((output_chan == 1 && freq_off < 1350) || freq_off < 1350 - output_rate / 2 || freq_off > output_rate / 2 - 1350) {
		std::cerr << "Unsupported frequency offset." << std::endl;
		return 1;
	}

	long long int call_sign = base37_encoder("ANONYMOUS");
	if (argc >= 8)
		call_sign = base37_encoder(argv[7]);

	if (call_sign <= 0 || call_sign >= 129961739795077L) {
		std::cerr << "Unsupported call sign." << std::endl;
		return 1;
	}

	typedef float value;
	typedef DSP::Complex<value> cmplx;

	std::ifstream input_file(input_name, std::ios::binary);
	if (input_file.bad()) {
		std::cerr << "Couldn't open file \"" << input_name << "\" for reading." << std::endl;
		return 1;
	}
	const int length = 64800 / 8;
	uint8_t *input_data = new uint8_t[length];
	for (int i = 0; i < length; ++i)
		input_data[i] = input_file.get();

	DSP::WriteWAV<value> output_file(output_name, output_rate, output_bits, output_chan);
	output_file.silence(output_rate);
	switch (output_rate) {
	case 8000:
		delete new Encoder<value, cmplx, 8000>(&output_file, input_data, freq_off, call_sign);
		break;
	case 16000:
		delete new Encoder<value, cmplx, 16000>(&output_file, input_data, freq_off, call_sign);
		break;
	case 44100:
		delete new Encoder<value, cmplx, 44100>(&output_file, input_data, freq_off, call_sign);
		break;
	case 48000:
		delete new Encoder<value, cmplx, 48000>(&output_file, input_data, freq_off, call_sign);
		break;
	default:
		std::cerr << "Unsupported sample rate." << std::endl;
		return 1;
	}
	output_file.silence(output_rate);
	delete []input_data;

	return 0;
}


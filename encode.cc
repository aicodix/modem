/*
OFDM modem encoder

Copyright 2021 Ahmet Inan <inan@aicodix.de>
*/

#include <iostream>
#include <cassert>
#include <cmath>
#include "xorshift.hh"
#include "complex.hh"
#include "utils.hh"
#include "bitman.hh"
#include "decibel.hh"
#include "fft.hh"
#include "wav.hh"
#include "pcm.hh"
#include "mls.hh"
#include "crc.hh"
#include "psk.hh"
#include "qam.hh"
#include "polar_tables.hh"
#include "polar_encoder.hh"
#include "bose_chaudhuri_hocquenghem_encoder.hh"

template <typename value, typename cmplx, int rate>
struct Encoder
{
	typedef int8_t code_type;
	static const int symbol_len = (1280 * rate) / 8000;
	static const int guard_len = symbol_len / 8;
	static const int bits_max = 65536;
	static const int data_max = 1024;
	static const int cols_max = 273 + 32;
	static const int mls0_len = 320;
	static const int mls0_poly = 0b1100110001;
	static const int mls0_seed = 214;
	static const int mls1_len = 255;
	static const int mls1_poly = 0b100101011;
	static const int mls2_poly = 0b100101010001;
	DSP::WritePCM<value> *pcm;
	DSP::FastFourierTransform<symbol_len, cmplx, -1> fwd;
	DSP::FastFourierTransform<symbol_len, cmplx, 1> bwd;
	CODE::CRC<uint16_t> crc0;
	CODE::CRC<uint32_t> crc1;
	CODE::BoseChaudhuriHocquenghemEncoder<255, 71> bchenc;
	CODE::PolarEncoder<code_type> polarenc;
	uint8_t input_data[data_max];
	code_type code[bits_max], perm[bits_max], mesg[bits_max];
	cmplx fdom[symbol_len];
	cmplx tdom[symbol_len];
	cmplx temp[symbol_len];
	cmplx kern[symbol_len];
	cmplx guard[guard_len];
	cmplx prev[cols_max];
	value weight[guard_len];
	value papr_min, papr_max;
	int mod_bits;
	int oper_mode;
	int code_order;
	int code_off;
	int cons_cols;
	int cons_rows;
	int mls0_off;
	int mls1_off;

	static int bin(int carrier)
	{
		return (carrier + symbol_len) % symbol_len;
	}
	static int nrz(bool bit)
	{
		return 1 - 2 * bit;
	}
	void clipping_and_filtering(value scale, bool limit)
	{
		for (int i = 0; i < symbol_len; ++i) {
			value pwr = norm(tdom[i]);
			if (pwr > value(1))
				tdom[i] /= sqrt(pwr);
		}
		fwd(temp, tdom);
		for (int i = 0; i < symbol_len; ++i) {
			if (norm(fdom[i])) {
				temp[i] *= scale / std::sqrt(value(symbol_len));
				cmplx err = temp[i] - fdom[i];
				value mag = abs(err);
				value lim = 0.1 * mod_distance();
				if (limit && mag > lim)
					temp[i] -= ((mag - lim) / mag) * err;
			} else {
				temp[i] = 0;
			}
		}
		bwd(tdom, temp);
		for (int i = 0; i < symbol_len; ++i)
			tdom[i] /= scale * std::sqrt(value(symbol_len));
	}
	void tone_reservation()
	{
		for (int n = 0; n < 100; ++n) {
			int peak = 0;
			for (int i = 1; i < symbol_len; ++i)
				if (norm(tdom[peak]) < norm(tdom[i]))
					peak = i;
			cmplx orig = tdom[peak];
			if (norm(orig) <= value(1))
				break;
			for (int i = 0; i < symbol_len; ++i)
				tdom[i] -= orig * kern[(symbol_len-peak+i)%symbol_len];
		}
	}
	void symbol(bool papr_reduction = true, bool guard_interval = true)
	{
		bwd(tdom, fdom);
		value scale = 2;
		for (int i = 0; i < symbol_len; ++i)
			tdom[i] /= scale * std::sqrt(value(symbol_len));
		clipping_and_filtering(scale, oper_mode > 0 && papr_reduction);
		if (oper_mode > 0 && papr_reduction)
			tone_reservation();
		auto clamp = [](value v){ return v < value(-1) ? value(-1) : v > value(1) ? value(1) : v; };
		for (int i = 0; i < symbol_len; ++i)
			tdom[i] = cmplx(clamp(tdom[i].real()), clamp(tdom[i].imag()));
		for (int i = 0; i < guard_len; ++i)
			guard[i] = DSP::lerp(guard[i], tdom[i+symbol_len-guard_len], weight[i]);
		value peak = 0, mean = 0;
		for (int i = 0; i < symbol_len; ++i) {
			value power(norm(tdom[i]));
			peak = std::max(peak, power);
			mean += power;
		}
		mean /= symbol_len;
		if (mean > 0) {
			value papr(peak / mean);
			papr_min = std::min(papr_min, papr);
			papr_max = std::max(papr_max, papr);
		}
		if (guard_interval)
			pcm->write(reinterpret_cast<value *>(guard), guard_len, 2);
		pcm->write(reinterpret_cast<value *>(tdom), symbol_len, 2);
		for (int i = 0; i < guard_len; ++i)
			guard[i] = tdom[i];
	}
	void pilot_block()
	{
		CODE::MLS seq2(mls2_poly);
		value code_fac = std::sqrt(value(symbol_len) / value(cons_cols));
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] = 0;
		for (int i = code_off; i < code_off + cons_cols; ++i)
			fdom[bin(i)] = code_fac * nrz(seq2());
		symbol();
	}
	void schmidl_cox()
	{
		CODE::MLS seq0(mls0_poly, mls0_seed);
		value mls0_fac = std::sqrt(value(symbol_len) / value(mls0_len));
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] = 0;
		for (int i = 0; i < mls0_len; ++i)
			fdom[bin(i+mls0_off)] = mls0_fac * nrz(seq0());
		symbol(false);
		symbol(false, false);
	}
	void meta_data(uint64_t md)
	{
		uint8_t data[9] = { 0 }, parity[23] = { 0 };
		for (int i = 0; i < 55; ++i)
			CODE::set_be_bit(data, i, (md>>i)&1);
		crc0.reset();
		uint16_t cs = crc0(md << 9);
		for (int i = 0; i < 16; ++i)
			CODE::set_be_bit(data, i+55, (cs>>i)&1);
		bchenc(data, parity);
		CODE::MLS seq1(mls1_poly);
		value cons_fac = std::sqrt(value(symbol_len) / value(cons_cols));
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] = 0;
		fdom[bin(mls1_off-1)] = cons_fac;
		for (int i = 0; i < 71; ++i)
			fdom[bin(i+mls1_off)] = nrz(CODE::get_be_bit(data, i));
		for (int i = 71; i < mls1_len; ++i)
			fdom[bin(i+mls1_off)] = nrz(CODE::get_be_bit(parity, i-71));
		for (int i = 0; i < mls1_len; ++i)
			fdom[bin(i+mls1_off)] *= fdom[bin(i-1+mls1_off)];
		for (int i = 0; i < mls1_len; ++i)
			fdom[bin(i+mls1_off)] *= nrz(seq1());
		if (oper_mode > 0) {
			for (int i = code_off; i < code_off + cons_cols; ++i) {
				if (i == mls1_off-1)
					i += mls1_len + 1;
				fdom[bin(i)] = cons_fac * nrz(seq1());
			}
		}
		symbol();
	}
	cmplx mod_map(code_type *b)
	{
		switch (mod_bits) {
		case 2:
			return PhaseShiftKeying<4, cmplx, code_type>::map(b);
		case 4:
			return QuadratureAmplitudeModulation<16, cmplx, code_type>::map(b);
		case 6:
			return QuadratureAmplitudeModulation<64, cmplx, code_type>::map(b);
		}
		return 0;
	}
	value mod_distance()
	{
		switch (mod_bits) {
		case 2:
			return PhaseShiftKeying<4, cmplx, code_type>::DIST;
		case 4:
			return QuadratureAmplitudeModulation<16, cmplx, code_type>::DIST;
		case 6:
			return QuadratureAmplitudeModulation<64, cmplx, code_type>::DIST;
		}
		return 2;
	}
	void shuffle(code_type *dest, const code_type *src)
	{
		if (code_order == 11) {
			CODE::XorShiftMask<int, 11, 1, 3, 4, 1> seq;
			dest[0] = src[0];
			for (int i = 1; i < 2048; ++i)
				dest[i] = src[seq()];
		} else if (code_order == 12) {
			CODE::XorShiftMask<int, 12, 1, 1, 4, 1> seq;
			dest[0] = src[0];
			for (int i = 1; i < 4096; ++i)
				dest[i] = src[seq()];
		} else if (code_order == 13) {
			CODE::XorShiftMask<int, 13, 1, 1, 9, 1> seq;
			dest[0] = src[0];
			for (int i = 1; i < 8192; ++i)
				dest[i] = src[seq()];
		} else if (code_order == 14) {
			CODE::XorShiftMask<int, 14, 1, 5, 10, 1> seq;
			dest[0] = src[0];
			for (int i = 1; i < 16384; ++i)
				dest[i] = src[seq()];
		} else if (code_order == 15) {
			CODE::XorShiftMask<int, 15, 1, 1, 3, 1> seq;
			dest[0] = src[0];
			for (int i = 1; i < 32768; ++i)
				dest[i] = src[seq()];
		} else if (code_order == 16) {
			CODE::XorShiftMask<int, 16, 1, 1, 14, 1> seq;
			dest[0] = src[0];
			for (int i = 1; i < 65536; ++i)
				dest[i] = src[seq()];
		}
	}
	Encoder(DSP::WritePCM<value> *pcm, const char *const *input_names, int input_count, int freq_off, uint64_t call_sign, int oper_mode) :
		pcm(pcm), crc0(0xA8F4), crc1(0x8F6E37A0), bchenc({
			0b100011101, 0b101110111, 0b111110011, 0b101101001,
			0b110111101, 0b111100111, 0b100101011, 0b111010111,
			0b000010011, 0b101100101, 0b110001011, 0b101100011,
			0b100011011, 0b100111111, 0b110001101, 0b100101101,
			0b101011111, 0b111111001, 0b111000011, 0b100111001,
			0b110101001, 0b000011111, 0b110000111, 0b110110001}),
			oper_mode(oper_mode)
	{
		const uint32_t *frozen_bits = nullptr;
		int code_cols = 0;
		int comb_cols = 32;
		int comb_dist = 1;
		int comb_off = 1;
		int data_bits = 0;
		int reserved_tones = 0;
		switch (oper_mode) {
		case 0:
			code_cols = 256;
			comb_cols = 0;
			break;
		case 23:
			mod_bits = 2;
			cons_rows = 8;
			code_order = 12;
			code_cols = 256;
			data_bits = 2048;
			reserved_tones = 32;
			frozen_bits = frozen_4096_2080;
			break;
		case 24:
			mod_bits = 2;
			cons_rows = 16;
			code_order = 13;
			code_cols = 256;
			data_bits = 4096;
			reserved_tones = 32;
			frozen_bits = frozen_8192_4128;
			break;
		case 25:
			mod_bits = 2;
			cons_rows = 32;
			code_order = 14;
			code_cols = 256;
			data_bits = 8192;
			reserved_tones = 32;
			frozen_bits = frozen_16384_8224;
			break;
		case 26:
			mod_bits = 4;
			cons_rows = 4;
			code_order = 12;
			code_cols = 256;
			data_bits = 2048;
			reserved_tones = 32;
			frozen_bits = frozen_4096_2080;
			break;
		case 27:
			mod_bits = 4;
			cons_rows = 8;
			code_order = 13;
			code_cols = 256;
			data_bits = 4096;
			reserved_tones = 32;
			frozen_bits = frozen_8192_4128;
			break;
		case 28:
			mod_bits = 4;
			cons_rows = 16;
			code_order = 14;
			code_cols = 256;
			data_bits = 8192;
			reserved_tones = 32;
			frozen_bits = frozen_16384_8224;
			break;
		case 29:
			mod_bits = 6;
			cons_rows = 5;
			code_order = 13;
			code_cols = 273;
			data_bits = 4096;
			reserved_tones = 15;
			frozen_bits = frozen_8192_4128;
			break;
		case 30:
			mod_bits = 6;
			cons_rows = 10;
			code_order = 14;
			code_cols = 273;
			data_bits = 8192;
			reserved_tones = 15;
			frozen_bits = frozen_16384_8224;
			break;
		default:
			return;
		}
		int data_bytes = data_bits / 8;
		int offset = (freq_off * symbol_len) / rate;
		mls0_off = offset - mls0_len / 2;
		mls1_off = offset - mls1_len / 2;
		cons_cols = code_cols + comb_cols;
		code_off = offset - cons_cols / 2;
		if (oper_mode > 0) {
			comb_dist = comb_cols ? cons_cols / comb_cols : 1;
			comb_off = comb_cols ? comb_dist / 2 : 1;
			if (reserved_tones) {
				value kern_fac = 1 / value(10 * reserved_tones);
				for (int i = 0, j = code_off - reserved_tones / 2; i < reserved_tones; ++i, ++j) {
					if (j == code_off)
						j += cons_cols;
					fdom[bin(j)] = kern_fac;
				}
				bwd(kern, fdom);
			}
		}
		for (int i = 0; i < guard_len / 4; ++i)
			weight[i] = 0;
		for (int i = guard_len / 4; i < guard_len / 4 + guard_len / 2; ++i) {
			value x = value(i - guard_len / 4) / value(guard_len / 2 - 1);
			weight[i] = value(0.5) * (value(1) - std::cos(DSP::Const<value>::Pi() * x));
		}
		for (int i = guard_len / 4 + guard_len / 2; i < guard_len; ++i)
			weight[i] = 1;
		papr_min = 1000, papr_max = -1000;
		pilot_block();
		if (!oper_mode) {
			schmidl_cox();
			meta_data(call_sign << 8);
		}
		for (int input_index = 0; input_index < input_count; ++input_index) {
			const char *input_name = input_names[input_index];
			if (input_count == 1 && input_name[0] == '-' && input_name[1] == 0)
				input_name = "/dev/stdin";
			std::ifstream input_file(input_name, std::ios::binary);
			if (input_file.bad()) {
				std::cerr << "Couldn't open file \"" << input_name << "\" for reading." << std::endl;
				continue;
			}
			for (int i = 0; i < data_bytes; ++i)
				input_data[i] = std::max(input_file.get(), 0);
			CODE::Xorshift32 scrambler;
			for (int i = 0; i < data_bytes; ++i)
				input_data[i] ^= scrambler();
			schmidl_cox();
			meta_data((call_sign << 8) | oper_mode);
			for (int i = 0; i < data_bits; ++i)
				mesg[i] = nrz(CODE::get_le_bit(input_data, i));
			crc1.reset();
			for (int i = 0; i < data_bytes; ++i)
				crc1(input_data[i]);
			for (int i = 0; i < 32; ++i)
				mesg[i+data_bits] = nrz((crc1()>>i)&1);
			polarenc(code, mesg, frozen_bits, code_order);
			shuffle(perm, code);
			for (int i = 0; i < cons_cols; ++i)
				prev[i] = fdom[bin(i+code_off)];
			CODE::MLS seq0(mls0_poly);
			for (int j = 0, k = 0; j < cons_rows; ++j) {
				for (int i = 0; i < cons_cols; ++i) {
					if (i % comb_dist == comb_off) {
						prev[i] *= nrz(seq0());
						fdom[bin(i+code_off)] = prev[i];
					} else {
						fdom[bin(i+code_off)] = prev[i] * mod_map(perm+k);
						k += mod_bits;
					}
				}
				symbol();
			}
		}
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] = 0;
		symbol();
		std::cerr << "PAPR: " << DSP::decibel(papr_min) << " .. " << DSP::decibel(papr_max) << " dB" << std::endl;
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
	if (argc < 8) {
		std::cerr << "usage: " << argv[0] << " OUTPUT RATE BITS CHANNELS OFFSET MODE CALLSIGN INPUT.." << std::endl;
		return 1;
	}

	const char *output_name = argv[1];
	if (output_name[0] == '-' && output_name[1] == 0)
		output_name = "/dev/stdout";
	int output_rate = std::atoi(argv[2]);
	int output_bits = std::atoi(argv[3]);
	int output_chan = std::atoi(argv[4]);

	int freq_off = std::atoi(argv[5]);
	if (freq_off % 50) {
		std::cerr << "Frequency offset must be divisible by 50." << std::endl;
		return 1;
	}
	int input_count = argc - 8;
	int oper_mode = std::atoi(argv[6]);
	if (!oper_mode != !input_count) {
		std::cerr << "Using operation mode " << oper_mode << " but " << input_count << " input file" << (input_count == 1 ? "" : "s") << " provided." << std::endl;
		return 1;
	}
	long long int call_sign = base37_encoder(argv[7]);
	if (call_sign <= 0 || call_sign >= 129961739795077L) {
		std::cerr << "Unsupported call sign." << std::endl;
		return 1;
	}

	if (oper_mode && (oper_mode < 23 || oper_mode > 30)) {
		std::cerr << "Unsupported operation mode." << std::endl;
		return 1;
	}
	int band_width = 2000;
	if ((output_chan == 1 && freq_off < band_width / 2) || freq_off < band_width / 2 - output_rate / 2 || freq_off > output_rate / 2 - band_width / 2) {
		std::cerr << "Unsupported frequency offset." << std::endl;
		return 1;
	}
	typedef float value;
	typedef DSP::Complex<value> cmplx;
	DSP::WriteWAV<value> output_file(output_name, output_rate, output_bits, output_chan);
	output_file.silence(output_rate);
	switch (output_rate) {
	case 8000:
		delete new Encoder<value, cmplx, 8000>(&output_file, argv+8, input_count, freq_off, call_sign, oper_mode);
		break;
	case 16000:
		delete new Encoder<value, cmplx, 16000>(&output_file, argv+8, input_count, freq_off, call_sign, oper_mode);
		break;
	case 44100:
		delete new Encoder<value, cmplx, 44100>(&output_file, argv+8, input_count, freq_off, call_sign, oper_mode);
		break;
	case 48000:
		delete new Encoder<value, cmplx, 48000>(&output_file, argv+8, input_count, freq_off, call_sign, oper_mode);
		break;
	default:
		std::cerr << "Unsupported sample rate." << std::endl;
		return 1;
	}
	output_file.silence(output_rate);

	return 0;
}


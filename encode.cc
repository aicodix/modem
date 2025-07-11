/*
OFDM modem encoder

Copyright 2021 Ahmet Inan <inan@aicodix.de>
*/

#include <iostream>
#include <cassert>
#include <cmath>
#include "common.hh"
#include "xorshift.hh"
#include "complex.hh"
#include "utils.hh"
#include "bitman.hh"
#include "decibel.hh"
#include "fft.hh"
#include "wav.hh"
#include "pcm.hh"
#include "mls.hh"
#include "psk.hh"
#include "qam.hh"
#include "polar_encoder.hh"

template <typename value, typename cmplx, int rate>
struct Encoder : public Common
{
	typedef int8_t code_type;
	static const int guard_len = rate / 100;
	static const int symbol_len = guard_len * 16;
	DSP::WritePCM<value> *pcm;
	DSP::FastFourierTransform<symbol_len, cmplx, -1> fwd;
	DSP::FastFourierTransform<symbol_len, cmplx, 1> bwd;
	CODE::PolarEncoder<code_type> polar_encoder;
	code_type code[bits_max], perm[bits_max], mesg[bits_max];
	cmplx fdom[symbol_len];
	cmplx tdom[symbol_len];
	cmplx best[symbol_len];
	cmplx kern[symbol_len];
	cmplx guard[guard_len];
	cmplx tone[tone_count];
	cmplx prev[tone_count];
	cmplx temp[tone_count];
	value weight[guard_len];

	static int bin(int carrier)
	{
		return (carrier + symbol_len) % symbol_len;
	}
	static int nrz(bool bit)
	{
		return 1 - 2 * bit;
	}
	void clipping_and_filtering(value scale)
	{
		for (int i = 0; i < symbol_len; ++i) {
			value pwr = norm(tdom[i]);
			if (pwr > value(1))
				tdom[i] /= sqrt(pwr);
		}
		fwd(fdom, tdom);
		for (int i = 0; i < symbol_len; ++i) {
			int j = bin(i + tone_off);
			if (i >= tone_count)
				fdom[j] = 0;
			else
				fdom[j] *= 1 / (scale * symbol_len);
		}
		bwd(tdom, fdom);
		for (int i = 0; i < symbol_len; ++i)
			tdom[i] *= scale;
		auto clamp = [](value v){ return v < value(-1) ? value(-1) : v > value(1) ? value(1) : v; };
		for (int i = 0; i < symbol_len; ++i)
			tdom[i] = cmplx(clamp(tdom[i].real()), clamp(tdom[i].imag()));
	}
	void symbol(int symbol_number)
	{
		value scale = value(0.5) / std::sqrt(value(tone_count));
		for (int i = 0; differential && symbol_number > 0 && i < tone_count; ++i)
			tone[i] *= prev[i];
		int trials = symbol_number ? 128 : 4;
		value best_papr = 1000;
		for (int trial = trials - 1; trial >= 0; --trial) {
			for (int i = 0; i < tone_count; ++i)
				temp[i] = tone[i];
			if (symbol_number >= 0) {
				int meta_data = symbol_number ? trial : (oper_mode | (trial << 5));
				hadamard_encoder(meta, meta_data);
				CODE::MLS seq(0x163, meta_data);
				for (int i = 0, m = 0; i < tone_count; ++i)
					if (i % block_length == pilot_off)
						temp[i] *= meta[m++];
					else
						temp[i] *= nrz(seq());
			}
			for (int i = 0; i < symbol_len; ++i)
				fdom[i] = 0;
			for (int i = 0; i < tone_count; ++i)
				fdom[bin(i+tone_off)] = temp[i];
			bwd(tdom, fdom);
			for (int i = 0; i < symbol_len; ++i)
				tdom[i] *= scale;
			if (symbol_number < 0)
				break;
			value peak = 0, mean = 0;
			for (int i = 0; i < symbol_len; ++i) {
				value power(norm(tdom[i]));
				peak = std::max(peak, power);
				mean += power;
			}
			mean /= symbol_len;
			value papr(peak / mean);
			if (papr < best_papr) {
				best_papr = papr;
				for (int i = 0; differential && symbol_number >= 0 && i < tone_count; ++i)
					prev[i] = temp[i];
				for (int i = 0; i < symbol_len; ++i)
					best[i] = tdom[i];
			}
		}
		if (symbol_number >= 0) {
			for (int i = 0; i < symbol_len; ++i)
				tdom[i] = best[i];
			std::cerr << " " << DSP::decibel(best_papr);
		}
		clipping_and_filtering(scale);
		if (symbol_number != -1) {
			for (int i = 0; i < guard_len; ++i)
				guard[i] = DSP::lerp(guard[i], tdom[i+symbol_len-guard_len], weight[i]);
			pcm->write(reinterpret_cast<value *>(guard), guard_len, 2);
		}
		for (int i = 0; i < guard_len; ++i)
			guard[i] = tdom[i];
		pcm->write(reinterpret_cast<value *>(tdom), symbol_len, 2);
	}
	void finish()
	{
		for (int i = 0; i < guard_len; ++i)
			guard[i] *= 1 - weight[i];
		pcm->write(reinterpret_cast<value *>(guard), guard_len, 2);
		for (int i = 0; i < guard_len; ++i)
			guard[i] = 0;
	}
	void leading_noise(int num = 1)
	{
		CODE::MLS noise(0x163);
		for (int j = 0; j < num; ++j) {
			for (int i = 0; i < tone_count; ++i)
				tone[i] = nrz(noise());
			symbol(-3);
		}
	}
	void schmidl_cox()
	{
		CODE::MLS seq0(mls0_poly, mls0_seed);
		for (int i = 0; i < tone_count; ++i)
			tone[i] = nrz(seq0());
		symbol(-2);
		symbol(-1);
	}
	cmplx map_bits(code_type *b, int bits)
	{
		switch (bits) {
		case 1:
			return PhaseShiftKeying<2, cmplx, code_type>::map(b);
		case 2:
			return PhaseShiftKeying<4, cmplx, code_type>::map(b);
		case 3:
			return PhaseShiftKeying<8, cmplx, code_type>::map(b);
		case 4:
			return QuadratureAmplitudeModulation<16, cmplx, code_type>::map(b);
		case 6:
			return QuadratureAmplitudeModulation<64, cmplx, code_type>::map(b);
		case 8:
			return QuadratureAmplitudeModulation<256, cmplx, code_type>::map(b);
		}
		return 0;
	}
	value mod_distance()
	{
		switch (mod_bits) {
		case 1:
			return PhaseShiftKeying<2, cmplx, code_type>::DIST;
		case 2:
			return PhaseShiftKeying<4, cmplx, code_type>::DIST;
		case 3:
			return PhaseShiftKeying<8, cmplx, code_type>::DIST;
		case 4:
			return QuadratureAmplitudeModulation<16, cmplx, code_type>::DIST;
		case 6:
			return QuadratureAmplitudeModulation<64, cmplx, code_type>::DIST;
		case 8:
			return QuadratureAmplitudeModulation<256, cmplx, code_type>::DIST;
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
	void guard_interval_weights()
	{
		for (int i = 0; i < guard_len / 4; ++i)
			weight[i] = 0;
		for (int i = guard_len / 4; i < guard_len / 4 + guard_len / 2; ++i) {
			value x = value(i - guard_len / 4) / value(guard_len / 2 - 1);
			weight[i] = value(0.5) * (value(1) - std::cos(DSP::Const<value>::Pi() * x));
		}
		for (int i = guard_len / 4 + guard_len / 2; i < guard_len; ++i)
			weight[i] = 1;
	}
	Encoder(DSP::WritePCM<value> *pcm, const char *const *input_names, int input_count, int freq_off, int oper_mode) : pcm(pcm)
	{
		setup(oper_mode);
		int offset = (freq_off * symbol_len) / rate;
		tone_off = offset - tone_count / 2;
		guard_interval_weights();
		leading_noise();
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
				data[i] = std::max(input_file.get(), 0);
			CODE::Xorshift32 scrambler;
			for (int i = 0; i < data_bytes; ++i)
				data[i] ^= scrambler();
			schmidl_cox();
			for (int i = 0; i < data_bits; ++i)
				mesg[i] = nrz(CODE::get_le_bit(data, i));
			crc0.reset();
			for (int i = 0; i < data_bytes; ++i)
				crc0(data[i]);
			for (int i = 0; i < 32; ++i)
				mesg[i+data_bits] = nrz((crc0()>>i)&1);
			polar_encoder(code, mesg, frozen_bits, code_order);
			shuffle(perm, code);
			CODE::MLS seq1(mls1_poly);
			std::cerr << "PAPR (dB):";
			for (int j = 0, k = 0; j < symbol_count; ++j) {
				pilot_off = (block_skew * j + first_pilot) % block_length;
				for (int i = 0; i < tone_count; ++i) {
					if (i % block_length == pilot_off) {
						tone[i] = nrz(seq1());
					} else {
						int bits = mod_bits;
						if (oper_mode >= 7 && oper_mode <= 9 && k % 32 == 30)
							bits = 2;
						if (oper_mode >= 21 && oper_mode <= 23 && k % 64 == 60)
							bits = 4;
						tone[i] = map_bits(perm+k, bits);
						k += bits;
					}
				}
				symbol(j);
			}
			std::cerr << std::endl;
		}
		finish();
	}
};

int main(int argc, char **argv)
{
	if (argc < 8) {
		std::cerr << "usage: " << argv[0] << " OUTPUT RATE BITS CHANNELS OFFSET MODE INPUT.." << std::endl;
		return 1;
	}

	const char *output_name = argv[1];
	if (output_name[0] == '-' && output_name[1] == 0)
		output_name = "/dev/stdout";
	int output_rate = std::atoi(argv[2]);
	int output_bits = std::atoi(argv[3]);
	int output_chan = std::atoi(argv[4]);

	int freq_off = std::atoi(argv[5]);
	if (freq_off % 100) {
		std::cerr << "Frequency offset must be divisible by 100." << std::endl;
		return 1;
	}
	int input_count = argc - 7;
	int oper_mode = std::atoi(argv[6]);
	if (oper_mode < 0 || oper_mode > 27) {
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
		delete new Encoder<value, cmplx, 8000>(&output_file, argv+7, input_count, freq_off, oper_mode);
		break;
	case 16000:
		delete new Encoder<value, cmplx, 16000>(&output_file, argv+7, input_count, freq_off, oper_mode);
		break;
	case 44100:
		delete new Encoder<value, cmplx, 44100>(&output_file, argv+7, input_count, freq_off, oper_mode);
		break;
	case 48000:
		delete new Encoder<value, cmplx, 48000>(&output_file, argv+7, input_count, freq_off, oper_mode);
		break;
	default:
		std::cerr << "Unsupported sample rate." << std::endl;
		return 1;
	}
	output_file.silence(output_rate);

	return 0;
}


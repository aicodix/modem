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
#include "hadamard_encoder.hh"

template <typename value, typename cmplx, int rate>
struct Encoder
{
	typedef int8_t code_type;
	static const int guard_len = rate / 100;
	static const int symbol_len = guard_len * 16;
	static const int bits_max = 65536;
	static const int data_max = 4096;
	static const int mls0_poly = 0x331;
	static const int mls0_seed = 214;
	static const int mls1_poly = 0x25;
	static const int data_tones = 256;
	static const int pilot_tones = 32;
	static const int reserved_tones = 32;
	static const int tone_count = data_tones + pilot_tones + reserved_tones;
	static const int block_length = 10;
	static const int block_skew = 3;
	static const int pilot_offset = 4;
	static const int reserved_offset = 9;
	DSP::WritePCM<value> *pcm;
	DSP::FastFourierTransform<symbol_len, cmplx, -1> fwd;
	DSP::FastFourierTransform<symbol_len, cmplx, 1> bwd;
	CODE::CRC<uint32_t> crc0;
	CODE::HadamardEncoder<6> hadamardenc;
	CODE::PolarEncoder<code_type> polarenc;
	uint8_t input_data[data_max];
	code_type code[bits_max], perm[bits_max], mesg[bits_max];
	int8_t mode[32];
	cmplx fdom[symbol_len];
	cmplx tdom[symbol_len];
	cmplx temp[symbol_len];
	cmplx kern[symbol_len];
	cmplx guard[guard_len];
	cmplx tone[tone_count];
	value weight[guard_len];
	value papr_min, papr_max;
	const uint32_t *frozen_bits;
	int mod_bits;
	int data_bits;
	int data_bytes;
	int code_order;
	int tone_off;
	int symbol_count;

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
				temp[i] *= 1 / (scale * symbol_len);
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
			tdom[i] *= scale;
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
				tdom[i] -= orig * kern[bin(i-peak)];
		}
	}
	void symbol(bool papr_reduction = true, bool guard_interval = true)
	{
		for (int i = 0; i < tone_count; ++i)
			fdom[bin(i+tone_off)] = tone[i];
		bwd(tdom, fdom);
		value scale = value(0.5) / std::sqrt(value(tone_count));
		for (int i = 0; i < symbol_len; ++i)
			tdom[i] *= scale;
		clipping_and_filtering(scale, papr_reduction);
		if (papr_reduction)
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
	void leading_noise(int num = 1)
	{
		CODE::MLS noise(0x163);
		for (int j = 0; j < num; ++j) {
			for (int i = 0; i < tone_count; ++i)
				tone[i] = nrz(noise());
			symbol(false);
		}
	}
	void schmidl_cox()
	{
		CODE::MLS seq0(mls0_poly, mls0_seed);
		for (int i = 0; i < tone_count; ++i)
			tone[i] = nrz(seq0());
		symbol(false);
		symbol(false, false);
	}
	cmplx map_bits(code_type *b, int bits)
	{
		switch (bits) {
		case 2:
			return PhaseShiftKeying<4, cmplx, code_type>::map(b);
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
		case 2:
			return PhaseShiftKeying<4, cmplx, code_type>::DIST;
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
	void setup(int oper_mode, int freq_off)
	{
		switch (oper_mode) {
		case 0:
			symbol_count = 1;
			break;
		case 1:
			mod_bits = 2;
			symbol_count = 4;
			code_order = 11;
			data_bits = 1024;
			frozen_bits = frozen_2048_1056;
			break;
		case 2:
			mod_bits = 2;
			symbol_count = 8;
			code_order = 12;
			data_bits = 2048;
			frozen_bits = frozen_4096_2080;
			break;
		case 3:
			mod_bits = 2;
			symbol_count = 16;
			code_order = 13;
			data_bits = 4096;
			frozen_bits = frozen_8192_4128;
			break;
		case 4:
			mod_bits = 2;
			symbol_count = 32;
			code_order = 14;
			data_bits = 8192;
			frozen_bits = frozen_16384_8224;
			break;
		case 5:
			mod_bits = 4;
			symbol_count = 4;
			code_order = 12;
			data_bits = 2048;
			frozen_bits = frozen_4096_2080;
			break;
		case 6:
			mod_bits = 4;
			symbol_count = 8;
			code_order = 13;
			data_bits = 4096;
			frozen_bits = frozen_8192_4128;
			break;
		case 7:
			mod_bits = 4;
			symbol_count = 16;
			code_order = 14;
			data_bits = 8192;
			frozen_bits = frozen_16384_8224;
			break;
		case 8:
			mod_bits = 4;
			symbol_count = 32;
			code_order = 15;
			data_bits = 16384;
			frozen_bits = frozen_32768_16416;
			break;
		case 9:
			mod_bits = 6;
			symbol_count = 11;
			code_order = 14;
			data_bits = 8192;
			frozen_bits = frozen_16384_8224;
			break;
		case 10:
			mod_bits = 6;
			symbol_count = 22;
			code_order = 15;
			data_bits = 16384;
			frozen_bits = frozen_32768_16416;
			break;
		case 11:
			mod_bits = 6;
			symbol_count = 44;
			code_order = 16;
			data_bits = 32768;
			frozen_bits = frozen_65536_32800;
			break;
		case 12:
			mod_bits = 8;
			symbol_count = 4;
			code_order = 13;
			data_bits = 4096;
			frozen_bits = frozen_8192_4128;
			break;
		case 13:
			mod_bits = 8;
			symbol_count = 8;
			code_order = 14;
			data_bits = 8192;
			frozen_bits = frozen_16384_8224;
			break;
		case 14:
			mod_bits = 8;
			symbol_count = 16;
			code_order = 15;
			data_bits = 16384;
			frozen_bits = frozen_32768_16416;
			break;
		case 15:
			mod_bits = 8;
			symbol_count = 32;
			code_order = 16;
			data_bits = 32768;
			frozen_bits = frozen_65536_32800;
			break;
		default:
			return;
		}
		data_bytes = data_bits / 8;
		int offset = (freq_off * symbol_len) / rate;
		tone_off = offset - tone_count / 2;
	}
	void tone_reservation_kernel(int roff)
	{
		value mag = 1 / value(10 * reserved_tones);
		for (int i = 0; i < symbol_len; ++i)
			temp[i] = 0;
		for (int i = 0; i < reserved_tones; ++i)
			temp[bin(i*block_length+tone_off+roff)] = mag;
		bwd(kern, temp);
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
	Encoder(DSP::WritePCM<value> *pcm, const char *const *input_names, int input_count, int freq_off, int oper_mode) :
		pcm(pcm), crc0(0x8F6E37A0)
	{
		setup(oper_mode, freq_off);
		guard_interval_weights();
		papr_min = 1000, papr_max = -1000;
		leading_noise();
		hadamardenc(mode, oper_mode);
		if (!oper_mode) {
			schmidl_cox();
			CODE::MLS seq1(mls1_poly);
			for (int i = 0, m = 0; i < tone_count; ++i) {
				if (i % block_length == pilot_offset) {
					tone[i] = nrz(seq1()) * mode[m++];
				} else {
					tone[i] = 0;
				}
			}
			symbol(false);
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
			for (int i = 0; i < data_bits; ++i)
				mesg[i] = nrz(CODE::get_le_bit(input_data, i));
			crc0.reset();
			for (int i = 0; i < data_bytes; ++i)
				crc0(input_data[i]);
			for (int i = 0; i < 32; ++i)
				mesg[i+data_bits] = nrz((crc0()>>i)&1);
			polarenc(code, mesg, frozen_bits, code_order);
			shuffle(perm, code);
			CODE::MLS seq1(mls1_poly);
			for (int j = 0, k = 0; j < symbol_count; ++j) {
				int poff = (block_skew * j + pilot_offset) % block_length;
				int roff = (block_skew * j + reserved_offset) % block_length;
				for (int i = 0, m = 0; i < tone_count; ++i) {
					if (i % block_length == poff) {
						tone[i] = nrz(seq1()) * mode[m++];
					} else if (i % block_length == roff) {
						tone[i] = 0;
					} else {
						int bits = mod_bits;
						if (oper_mode >= 9 && oper_mode <= 11 && k % 64 == 60)
							bits = 4;
						tone[i] = map_bits(perm+k, bits);
						k += bits;
					}
				}
				tone_reservation_kernel(roff);
				symbol();
			}
		}
		for (int i = 0; i < tone_count; ++i)
			tone[i] = 0;
		symbol(false);
		std::cerr << "PAPR: " << DSP::decibel(papr_min) << " .. " << DSP::decibel(papr_max) << " dB" << std::endl;
	}
};

int main(int argc, char **argv)
{
	if (argc < 7) {
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
	if (!oper_mode != !input_count) {
		std::cerr << "Using operation mode " << oper_mode << " but " << input_count << " input file" << (input_count == 1 ? "" : "s") << " provided." << std::endl;
		return 1;
	}
	if (oper_mode < 0 || oper_mode > 15) {
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


/*
OFDM modem common bits

Copyright 2025 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

#include <cstdint>
#include <cstring>
#include <cassert>
#include <cmath>
#include <algorithm>
namespace DSP { using std::abs; using std::min; using std::nearbyint; using std::copysign; }
#include "sincos.hh"
#include "atan2.hh"
#include "xorshift.hh"
#include "complex.hh"
#include "decibel.hh"
#include "bitman.hh"
#include "quick.hh"
#include "fft.hh"
#include "mls.hh"
#include "psk.hh"
#include "qam.hh"
#include "crc.hh"
#include "polar_tables.hh"
#include "hadamard_encoder.hh"

struct Common
{
	static const int mod_max = 8;
	static const int code_max = 16;
	static const int bits_max = 1 << code_max;
	static const int data_max = 8192;
	// the preamble symbol and the payload symbols of the longest frames, 8PSK and QAM64 normal
	static const int symbols_max = 1 + 22;
	static const int mls0_poly = 0x2d1;
	static const int mls0_seed = 214;
	static const int mls1_poly = 0x43;
	static const int mls2_poly = 0x163;
	static const int data_tones = 256;
	static const int seed_tones = 64;
	static const int tone_count = data_tones + seed_tones;
	static const int block_length = 5;
	static const int block_skew = 3;
	static const int first_seed = 4;
	CODE::CRC<uint16_t> crc0;
	CODE::CRC<uint32_t> crc1;
	CODE::HadamardEncoder<7> hadamard_encoder;
	int8_t seed[seed_tones];
	uint8_t data[data_max];
	const uint32_t *frozen_bits;
	int mod_bits;
	int data_bits;
	int data_bytes;
	int code_order;
	int oper_mode;
	int seed_off;
	int symbol_count;

	Common() : crc0(0xA8F4), crc1(0x8F6E37A0) {}

	template <typename cmplx, typename code_type>
	static typename cmplx::value_type mod_distance(int bits)
	{
		switch (bits) {
		case 1:
			return PhaseShiftKeying<1, cmplx, code_type>::DIST;
		case 2:
			return PhaseShiftKeying<2, cmplx, code_type>::DIST;
		case 3:
			return PhaseShiftKeying<3, cmplx, code_type>::DIST;
		case 4:
			return QuadratureAmplitudeModulation<4, cmplx, code_type>::DIST;
		case 6:
			return QuadratureAmplitudeModulation<6, cmplx, code_type>::DIST;
		case 8:
			return QuadratureAmplitudeModulation<8, cmplx, code_type>::DIST;
		}
		return 2;
	}

	// some tones carry fewer bits, so the code length fits
	int bits_at(int k)
	{
		if (mod_bits == 3 && k % 32 == 30)
			return 2;
		if (mod_bits == 6 && k % 64 == 60)
			return 4;
		return mod_bits;
	}

	template <typename SEQ, typename TYPE>
	static void permute(TYPE *dest, const TYPE *src, int size, bool forward)
	{
		SEQ seq;
		dest[0] = src[0];
		for (int i = 1; i < size; ++i) {
			int j = seq();
			if (forward)
				dest[i] = src[j];
			else
				dest[j] = src[i];
		}
	}

	// forward interleaves the code bits like the encoder, backward undoes it like the decoder
	template <typename TYPE>
	static void shuffle(TYPE *dest, const TYPE *src, int order, bool forward)
	{
		if (order == 8)
			permute<CODE::XorShiftMask<int, 8, 1, 1, 2, 1>>(dest, src, 256, forward);
		else if (order == 11)
			permute<CODE::XorShiftMask<int, 11, 1, 3, 4, 1>>(dest, src, 2048, forward);
		else if (order == 12)
			permute<CODE::XorShiftMask<int, 12, 1, 1, 4, 1>>(dest, src, 4096, forward);
		else if (order == 13)
			permute<CODE::XorShiftMask<int, 13, 1, 1, 9, 1>>(dest, src, 8192, forward);
		else if (order == 14)
			permute<CODE::XorShiftMask<int, 14, 1, 5, 10, 1>>(dest, src, 16384, forward);
		else if (order == 15)
			permute<CODE::XorShiftMask<int, 15, 1, 1, 3, 1>>(dest, src, 32768, forward);
		else if (order == 16)
			permute<CODE::XorShiftMask<int, 16, 1, 1, 14, 1>>(dest, src, 65536, forward);
	}

	bool setup(int mode)
	{
		if (mode & 128)
			return false;
		int modulation = (mode >> 4) & 7;
		switch (modulation) {
		case 0:
			mod_bits = 1;
			symbol_count = 8;
			code_order = 11;
			break;
		case 1:
			mod_bits = 2;
			symbol_count = 4;
			code_order = 11;
			break;
		case 2:
			mod_bits = 3;
			symbol_count = 11;
			code_order = 13;
			break;
		case 3:
			mod_bits = 4;
			symbol_count = 4;
			code_order = 12;
			break;
		case 4:
			mod_bits = 6;
			symbol_count = 11;
			code_order = 14;
			break;
		case 5:
			mod_bits = 8;
			symbol_count = 8;
			code_order = 14;
			break;
		default:
			return false;
		}
		bool frame_size = mode & 1;
		if (frame_size) {
			if (symbol_count == 4) {
				symbol_count *= 4;
				code_order += 2;
			} else {
				symbol_count *= 2;
				++code_order;
			}
		}
		int code_rate = (mode >> 1) & 7;
		if (code_rate == 0) {
			switch (code_order) {
			case 11:
				data_bits = 1024;
				frozen_bits = frozen_2048_1056;
				break;
			case 12:
				data_bits = 2048;
				frozen_bits = frozen_4096_2080;
				break;
			case 13:
				data_bits = 4096;
				frozen_bits = frozen_8192_4128;
				break;
			case 14:
				data_bits = 8192;
				frozen_bits = frozen_16384_8224;
				break;
			case 15:
				data_bits = 16384;
				frozen_bits = frozen_32768_16416;
				break;
			case 16:
				data_bits = 32768;
				frozen_bits = frozen_65536_32800;
				break;
			default:
				return false;
			}
		} else if (code_rate == 1) {
			switch (code_order) {
			case 11:
				data_bits = 1368;
				frozen_bits = frozen_2048_1400;
				break;
			case 12:
				data_bits = 2736;
				frozen_bits = frozen_4096_2768;
				break;
			case 13:
				data_bits = 5472;
				frozen_bits = frozen_8192_5504;
				break;
			case 14:
				data_bits = 10944;
				frozen_bits = frozen_16384_10976;
				break;
			case 15:
				data_bits = 21888;
				frozen_bits = frozen_32768_21920;
				break;
			case 16:
				data_bits = 43776;
				frozen_bits = frozen_65536_43808;
				break;
			default:
				return false;
			}
		} else if (code_rate == 2) {
			switch (code_order) {
			case 11:
				data_bits = 1536;
				frozen_bits = frozen_2048_1568;
				break;
			case 12:
				data_bits = 3072;
				frozen_bits = frozen_4096_3104;
				break;
			case 13:
				data_bits = 6144;
				frozen_bits = frozen_8192_6176;
				break;
			case 14:
				data_bits = 12288;
				frozen_bits = frozen_16384_12320;
				break;
			case 15:
				data_bits = 24576;
				frozen_bits = frozen_32768_24608;
				break;
			case 16:
				data_bits = 49152;
				frozen_bits = frozen_65536_49184;
				break;
			default:
				return false;
			}
		} else if (code_rate == 3) {
			switch (code_order) {
			case 11:
				data_bits = 1704;
				frozen_bits = frozen_2048_1736;
				break;
			case 12:
				data_bits = 3408;
				frozen_bits = frozen_4096_3440;
				break;
			case 13:
				data_bits = 6816;
				frozen_bits = frozen_8192_6848;
				break;
			case 14:
				data_bits = 13632;
				frozen_bits = frozen_16384_13664;
				break;
			case 15:
				data_bits = 27264;
				frozen_bits = frozen_32768_27296;
				break;
			case 16:
				data_bits = 54528;
				frozen_bits = frozen_65536_54560;
				break;
			default:
				return false;
			}
		} else {
			return false;
		}
		oper_mode = mode;
		data_bytes = data_bits / 8;
		return true;
	}
};


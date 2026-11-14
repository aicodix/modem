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
	// the code bits of the longest frame, in codewords of up to bits_max bits
	static const int frame_max = 256 * 256 * mod_max;
	static const int data_max = frame_max / 8;
	// the preamble symbol and the payload symbols of the longest frame
	static const int symbols_max = 1 + 256;
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
	uint32_t frozen_bits[bits_max / 32];
	int mod_bits;
	int data_bits;
	int data_bytes;
	int code_order;
	int oper_mode;
	int seed_off;
	int symbol_count;
	// the code bits of a frame, in codewords of 1 << code_order bits each
	int frame_bits;
	int codewords;

	Common() : crc0(0xA8F4), crc1(0x8F6E37A0) {}

	template <typename cmplx, typename code_type>
	static typename cmplx::value_type mod_distance(int bits)
	{
		switch (bits) {
		case 1:
			return PhaseShiftKeying<1, cmplx, code_type>::DIST;
		case 2:
			return PhaseShiftKeying<2, cmplx, code_type>::DIST;
		case 4:
			return QuadratureAmplitudeModulation<4, cmplx, code_type>::DIST;
		case 6:
			return QuadratureAmplitudeModulation<6, cmplx, code_type>::DIST;
		case 8:
			return QuadratureAmplitudeModulation<8, cmplx, code_type>::DIST;
		}
		return 2;
	}


	template <typename SEQ, typename TYPE>
	static void permute(TYPE *dest, const TYPE *src, int size, bool forward)
	{
		SEQ seq;
		dest[0] = src[0];
		for (int i = 1; i < size; ++i) {
			// below a power of two, skip the indices beyond the size
			int j;
			do
				j = seq();
			while (j >= size);
			if (forward)
				dest[i] = src[j];
			else
				dest[j] = src[i];
		}
	}

	// forward interleaves the given number of code bits like the encoder, backward undoes it like the decoder
	// the xorshift triples have the fewest linear dependencies between the low bits, and between the top bits,
	// of an index and of the indices up to eight steps later, as ranked by tools/rank.cc of code
	template <typename TYPE>
	static void shuffle(TYPE *dest, const TYPE *src, int bits, bool forward)
	{
		int order = 0;
		while (1 << order < bits)
			++order;
		if (order == 8)
			permute<CODE::XorShiftMask<int, 8, 3, 5, 5, 1>>(dest, src, bits, forward);
		else if (order == 10)
			permute<CODE::XorShiftMask<int, 10, 3, 7, 8, 1>>(dest, src, bits, forward);
		else if (order == 11)
			permute<CODE::XorShiftMask<int, 11, 7, 5, 3, 1>>(dest, src, bits, forward);
		else if (order == 12)
			permute<CODE::XorShiftMask<int, 12, 5, 7, 11, 1>>(dest, src, bits, forward);
		else if (order == 13)
			permute<CODE::XorShiftMask<int, 13, 4, 7, 5, 1>>(dest, src, bits, forward);
		else if (order == 14)
			permute<CODE::XorShiftMask<int, 14, 3, 7, 9, 1>>(dest, src, bits, forward);
		else if (order == 15)
			permute<CODE::XorShiftMask<int, 15, 13, 8, 5, 1>>(dest, src, bits, forward);
		else if (order == 16)
			permute<CODE::XorShiftMask<int, 16, 9, 7, 13, 1>>(dest, src, bits, forward);
		else if (order == 17)
			permute<CODE::XorShiftMask<int, 17, 8, 9, 14, 1>>(dest, src, bits, forward);
		else if (order == 18)
			permute<CODE::XorShiftMask<int, 18, 3, 7, 4, 1>>(dest, src, bits, forward);
		else if (order == 19)
			permute<CODE::XorShiftMask<int, 19, 7, 8, 1, 1>>(dest, src, bits, forward);
	}

	// the bits per tone and k of the code rate k / 16 of each efficiency index. Each modulation starts at
	// rate 1 / 2 and ends below where the next one starts, rates below 1 / 2 only with BPSK, spaced from
	// 1 / 16 to 7.5 bits per tone as evenly as the rates allow
	static constexpr int efficiency[32][2] = {
		{ 1, 1 }, { 1, 2 }, { 1, 3 }, { 1, 4 }, { 1, 5 }, { 1, 6 }, { 1, 7 }, { 1, 8 },
		{ 1, 9 }, { 1, 10 }, { 1, 11 }, { 1, 12 }, { 1, 14 }, { 2, 8 }, { 2, 9 }, { 2, 10 },
		{ 2, 11 }, { 2, 12 }, { 2, 14 }, { 4, 8 }, { 4, 9 }, { 4, 10 }, { 4, 11 }, { 6, 8 },
		{ 6, 9 }, { 6, 10 }, { 8, 8 }, { 8, 9 }, { 8, 10 }, { 8, 11 }, { 8, 13 }, { 8, 15 },
	};
	// the mode carries 2 << (mode >> 5) payload symbols, from 4 to 256, and the efficiency index mode & 31,
	// or no payload symbols at all for mode >> 5 == 0, a ping with the value mode & 31.
	// The code bits of a frame are a power of two, or three times one, and longer frames have several
	// codewords of the longest code, interleaved across the whole frame
	bool setup(int mode)
	{
		if (mode < 0 || mode > 255)
			return false;
		if (!(mode >> 5)) {
			mod_bits = 1;
			symbol_count = 0;
			frame_bits = 0;
			codewords = 0;
			code_order = 0;
			data_bits = 0;
			oper_mode = mode;
			data_bytes = 0;
			return true;
		}
		mod_bits = efficiency[mode & 31][0];
		symbol_count = 2 << (mode >> 5);
		frame_bits = symbol_count * data_tones * mod_bits;
		code_order = 0;
		while (code_order < code_max && !((frame_bits >> code_order) & 1))
			++code_order;
		codewords = frame_bits >> code_order;
		const uint8_t *nested;
		switch (code_order) {
		case 10:
			nested = nested_1024_32;
			break;
		case 11:
			nested = nested_2048_32;
			break;
		case 12:
			nested = nested_4096_32;
			break;
		case 13:
			nested = nested_8192_32;
			break;
		case 14:
			nested = nested_16384_32;
			break;
		case 15:
			nested = nested_32768_32;
			break;
		case 16:
			nested = nested_65536_32;
			break;
		default:
			return false;
		}
		// the payload is k / 16 of each codeword, with its 32 bit CRC on top
		int k = efficiency[mode & 31][1];
		data_bits = codewords * (k << (code_order - 4));
		// bit i is frozen, if its entry is at least k
		for (int i = 0; i < 1 << (code_order - 5); ++i)
			frozen_bits[i] = 0;
		for (int i = 0; i < 1 << code_order; ++i)
			if (((nested[i/2] >> (i%2*4)) & 15) >= k)
				frozen_bits[i/32] |= uint32_t(1) << (i%32);
		oper_mode = mode;
		data_bytes = data_bits / 8;
		return true;
	}
	// the duration of a frame in seconds is 41 / 300 times the payload symbols plus three,
	// pick the payload symbols, as 2 << code, whose frame lasts closest to the given seconds
	static int symbols_code(int seconds)
	{
		auto off = [seconds](int code){ int d = 41 * (3 + (2 << code)) - 300 * seconds; return d < 0 ? -d : d; };
		int best = 1;
		for (int code = 2; code < 8; ++code)
			if (off(code) < off(best))
				best = code;
		return best;
	}
	// the most robust mode, the lowest efficiency index, that carries length bytes in about the given seconds,
	// set up for it, or -1 if none does
	int pick_mode(int seconds, int length)
	{
		int code = symbols_code(seconds);
		for (int index = 0; index < 32; ++index) {
			int mode = code << 5 | index;
			if (setup(mode) && data_bytes >= length)
				return mode;
		}
		return -1;
	}
};


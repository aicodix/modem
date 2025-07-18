/*
OFDM modem common bits

Copyright 2025 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

#include "crc.hh"
#include "polar_tables.hh"
#include "hadamard_encoder.hh"

struct Common
{
	static const int mod_max = 8;
	static const int code_max = 16;
	static const int bits_max = 1 << code_max;
	static const int data_max = 4096;
	static const int symbols_max = 26 + 1;
	static const int mls0_poly = 0x331;
	static const int mls0_seed = 214;
	static const int mls1_poly = 0x43;
	static const int mls2_poly = 0x163;
	static const int data_tones = 256;
	static const int head_tones = 32;
	static const int tail_tones = 32;
	static const int tone_count = data_tones + head_tones + tail_tones;
	static const int block_length = 10;
	static const int block_skew = 3;
	static const int first_head = 4;
	static const int first_tail = 9;
	static constexpr int slm_poly[16] = {
			0x11d, 0x12b, 0x12d, 0x14d, 0x15f, 0x163, 0x165, 0x169,
			0x171, 0x187, 0x18d, 0x1a9, 0x1c3, 0x1cf, 0x1e7, 0x1f5 };
	CODE::CRC<uint16_t> crc0;
	CODE::CRC<uint32_t> crc1;
	CODE::HadamardEncoder<6> hadamard_encoder;
	int8_t head[32];
	int8_t tail[32];
	uint8_t data[data_max];
	const uint32_t *frozen_bits;
	int mod_bits;
	int data_bits;
	int data_bytes;
	int code_order;
	int oper_mode;
	int tone_off;
	int head_off;
	int tail_off;
	int symbol_count;
	bool differential;

	Common() : crc0(0xA8F4), crc1(0x8F6E37A0) {}

	bool setup(int mode)
	{
		bool analog_mode = mode & 128;
		if (analog_mode) {
			std::cerr << "analog mode not supported yet" << std::endl;
			return false;
		}
		int modulation = (mode >> 4) & 7;
		switch (modulation) {
		case 0:
			mod_bits = 1;
			symbol_count = 8;
			differential = true;
			code_order = 11;
			break;
		case 1:
			mod_bits = 2;
			symbol_count = 4;
			differential = true;
			code_order = 11;
			break;
		case 2:
			mod_bits = 3;
			symbol_count = 11;
			differential = true;
			code_order = 13;
			break;
		case 3:
			mod_bits = 4;
			symbol_count = 4;
			differential = false;
			code_order = 12;
			break;
		case 4:
			mod_bits = 6;
			symbol_count = 11;
			differential = false;
			code_order = 14;
			break;
		case 5:
			mod_bits = 8;
			symbol_count = 8;
			differential = false;
			code_order = 14;
			break;
		case 6:
			mod_bits = 10;
			symbol_count = 13;
			differential = false;
			code_order = 15;
			break;
		case 7:
			mod_bits = 12;
			symbol_count = 11;
			differential = false;
			code_order = 15;
			break;
		default:
			return false;
		}
		bool frame_length = mode & 1;
		if (frame_length) {
			symbol_count *= 2;
			++code_order;
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
		} else {
			std::cerr << "code rate " << code_rate << " not supported yet" << std::endl;
			return false;
		}
		oper_mode = mode;
		data_bytes = data_bits / 8;
		return true;
	}
};


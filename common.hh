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
	static const int mls0_poly = 0x331;
	static const int mls0_seed = 214;
	static const int mls1_poly = 0x43;
	static const int data_tones = 256;
	static const int pilot_tones = 64;
	static const int tone_count = data_tones + pilot_tones;
	static const int block_length = 5;
	static const int block_skew = 3;
	static const int first_pilot = 4;
	CODE::CRC<uint32_t> crc0;
	CODE::HadamardEncoder<7> hadamard_encoder;
	int8_t meta[64];
	uint8_t data[data_max];
	const uint32_t *frozen_bits;
	int mod_bits;
	int data_bits;
	int data_bytes;
	int code_order;
	int oper_mode;
	int tone_off;
	int pilot_off;
	int reserved_off;
	int symbol_count;
	bool differential;

	Common() : crc0(0x8F6E37A0) {}

	void setup(int mode)
	{
		switch (mode) {
		case 0:
			mod_bits = 1;
			symbol_count = 8;
			differential = true;
			code_order = 11;
			data_bits = 1024;
			frozen_bits = frozen_2048_1056;
			break;
		case 1:
			mod_bits = 1;
			symbol_count = 16;
			differential = true;
			code_order = 12;
			data_bits = 2048;
			frozen_bits = frozen_4096_2080;
			break;
		case 2:
			mod_bits = 1;
			symbol_count = 32;
			differential = true;
			code_order = 13;
			data_bits = 4096;
			frozen_bits = frozen_8192_4128;
			break;
		case 3:
			mod_bits = 2;
			symbol_count = 4;
			differential = true;
			code_order = 11;
			data_bits = 1024;
			frozen_bits = frozen_2048_1056;
			break;
		case 4:
			mod_bits = 2;
			symbol_count = 8;
			differential = true;
			code_order = 12;
			data_bits = 2048;
			frozen_bits = frozen_4096_2080;
			break;
		case 5:
			mod_bits = 2;
			symbol_count = 16;
			differential = true;
			code_order = 13;
			data_bits = 4096;
			frozen_bits = frozen_8192_4128;
			break;
		case 6:
			mod_bits = 2;
			symbol_count = 32;
			differential = true;
			code_order = 14;
			data_bits = 8192;
			frozen_bits = frozen_16384_8224;
			break;
		case 7:
			mod_bits = 3;
			symbol_count = 11;
			differential = true;
			code_order = 13;
			data_bits = 4096;
			frozen_bits = frozen_8192_4128;
			break;
		case 8:
			mod_bits = 3;
			symbol_count = 22;
			differential = true;
			code_order = 14;
			data_bits = 8192;
			frozen_bits = frozen_16384_8224;
			break;
		case 9:
			mod_bits = 3;
			symbol_count = 44;
			differential = true;
			code_order = 15;
			data_bits = 16384;
			frozen_bits = frozen_32768_16416;
			break;
		case 10:
			mod_bits = 1;
			symbol_count = 8;
			differential = false;
			code_order = 11;
			data_bits = 1024;
			frozen_bits = frozen_2048_1056;
			break;
		case 11:
			mod_bits = 1;
			symbol_count = 16;
			differential = false;
			code_order = 12;
			data_bits = 2048;
			frozen_bits = frozen_4096_2080;
			break;
		case 12:
			mod_bits = 1;
			symbol_count = 32;
			differential = false;
			code_order = 13;
			data_bits = 4096;
			frozen_bits = frozen_8192_4128;
			break;
		case 13:
			mod_bits = 2;
			symbol_count = 4;
			differential = false;
			code_order = 11;
			data_bits = 1024;
			frozen_bits = frozen_2048_1056;
			break;
		case 14:
			mod_bits = 2;
			symbol_count = 8;
			differential = false;
			code_order = 12;
			data_bits = 2048;
			frozen_bits = frozen_4096_2080;
			break;
		case 15:
			mod_bits = 2;
			symbol_count = 16;
			differential = false;
			code_order = 13;
			data_bits = 4096;
			frozen_bits = frozen_8192_4128;
			break;
		case 16:
			mod_bits = 2;
			symbol_count = 32;
			differential = false;
			code_order = 14;
			data_bits = 8192;
			frozen_bits = frozen_16384_8224;
			break;
		case 17:
			mod_bits = 4;
			symbol_count = 4;
			differential = false;
			code_order = 12;
			data_bits = 2048;
			frozen_bits = frozen_4096_2080;
			break;
		case 18:
			mod_bits = 4;
			symbol_count = 8;
			differential = false;
			code_order = 13;
			data_bits = 4096;
			frozen_bits = frozen_8192_4128;
			break;
		case 19:
			mod_bits = 4;
			symbol_count = 16;
			differential = false;
			code_order = 14;
			data_bits = 8192;
			frozen_bits = frozen_16384_8224;
			break;
		case 20:
			mod_bits = 4;
			symbol_count = 32;
			differential = false;
			code_order = 15;
			data_bits = 16384;
			frozen_bits = frozen_32768_16416;
			break;
		case 21:
			mod_bits = 6;
			symbol_count = 11;
			differential = false;
			code_order = 14;
			data_bits = 8192;
			frozen_bits = frozen_16384_8224;
			break;
		case 22:
			mod_bits = 6;
			symbol_count = 22;
			differential = false;
			code_order = 15;
			data_bits = 16384;
			frozen_bits = frozen_32768_16416;
			break;
		case 23:
			mod_bits = 6;
			symbol_count = 44;
			differential = false;
			code_order = 16;
			data_bits = 32768;
			frozen_bits = frozen_65536_32800;
			break;
		case 24:
			mod_bits = 8;
			symbol_count = 4;
			differential = false;
			code_order = 13;
			data_bits = 4096;
			frozen_bits = frozen_8192_4128;
			break;
		case 25:
			mod_bits = 8;
			symbol_count = 8;
			differential = false;
			code_order = 14;
			data_bits = 8192;
			frozen_bits = frozen_16384_8224;
			break;
		case 26:
			mod_bits = 8;
			symbol_count = 16;
			differential = false;
			code_order = 15;
			data_bits = 16384;
			frozen_bits = frozen_32768_16416;
			break;
		case 27:
			mod_bits = 8;
			symbol_count = 32;
			differential = false;
			code_order = 16;
			data_bits = 32768;
			frozen_bits = frozen_65536_32800;
			break;
		default:
			return;
		}
		oper_mode = mode;
		data_bytes = data_bits / 8;
	}
};


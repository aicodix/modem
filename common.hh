/*
OFDM modem common bits

Copyright 2025 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

#include <iomanip>
#include <iostream>
#include <cstdint>
#include <cstring>
#include <cassert>
#include <cmath>
namespace DSP { using std::abs; using std::min; using std::cos; using std::sin; }
#include "xorshift.hh"
#include "complex.hh"
#include "decibel.hh"
#include "bitman.hh"
#include "quick.hh"
#include "wav.hh"
#include "pcm.hh"
#include "fft.hh"
#include "mls.hh"
#include "psk.hh"
#include "qam.hh"
#include "crc.hh"
#include "polar_tables.hh"
#include "hadamard_encoder.hh"

struct Common
{
	static const int mod_max = 12;
	static const int code_max = 16;
	static const int bits_max = 1 << code_max;
	static const int data_max = 8192;
	static const int symbols_max = 26 + 1;
	static const int mls0_poly = 0x331;
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
	int tone_off;
	int seed_off;
	int symbol_count;

	Common() : crc0(0xA8F4), crc1(0x8F6E37A0) {}

	bool setup(int mode)
	{
		bool analog_mode = mode & 128;
		if (analog_mode) {
			std::cerr << "analog mode not supported yet" << std::endl;
			return false;
		}
		std::cerr << "modulation: ";
		int modulation = (mode >> 4) & 7;
		switch (modulation) {
		case 0:
			mod_bits = 1;
			symbol_count = 8;
			code_order = 11;
			std::cerr << "BPSK";
			break;
		case 1:
			mod_bits = 2;
			symbol_count = 4;
			code_order = 11;
			std::cerr << "QPSK";
			break;
		case 2:
			mod_bits = 3;
			symbol_count = 11;
			code_order = 13;
			std::cerr << "8PSK";
			break;
		case 3:
			mod_bits = 4;
			symbol_count = 4;
			code_order = 12;
			std::cerr << "QAM16";
			break;
		case 4:
			mod_bits = 6;
			symbol_count = 11;
			code_order = 14;
			std::cerr << "QAM64";
			break;
		case 5:
			mod_bits = 8;
			symbol_count = 8;
			code_order = 14;
			std::cerr << "QAM256";
			break;
		case 6:
			mod_bits = 10;
			symbol_count = 13;
			code_order = 15;
			std::cerr << "QAM1024";
			break;
		case 7:
			mod_bits = 12;
			symbol_count = 11;
			code_order = 15;
			std::cerr << "QAM4096";
			break;
		default:
			return false;
		}
		std::cerr << std::endl;
		bool frame_size = mode & 1;
		std::cerr << "frame size: " << (frame_size ? "normal" : "short") << std::endl;
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
		std::cerr << "code rate: ";
		if (code_rate == 0) {
			std::cerr << "1/2";
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
			std::cerr << "2/3";
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
			std::cerr << "3/4";
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
			std::cerr << "5/6";
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
			std::cerr << "unsupported" << std::endl;
			return false;
		}
		std::cerr << std::endl;
		oper_mode = mode;
		data_bytes = data_bits / 8;
		float duration = 41. / 300. * (3 + symbol_count);
		std::cerr << "duration: " << duration << "s" << std::endl;
		std::cerr << "payload: " << data_bytes << "B" << std::endl;
		std::cerr << "bitrate: " << data_bits / duration / 1000. << "kb/s" << std::endl;
		return true;
	}
};


/*
OFDM modem helpers for the native programs

Copyright 2026 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

#include <iomanip>
#include <iostream>
#include "wav.hh"
#include "pcm.hh"
#include "common.hh"

static void print_setup(const Common &common, int mode)
{
	// a mode that failed to set up leaves oper_mode alone
	if (common.oper_mode != mode) {
		std::cerr << "oper mode: unsupported" << std::endl;
		return;
	}
	if (!common.symbol_count) {
		std::cerr << "ping: " << (mode & 31) << std::endl;
		std::cerr << "duration: " << 41. / 300. * 3 << "s" << std::endl;
		return;
	}
	const char *modulation = common.mod_bits == 1 ? "BPSK" : common.mod_bits == 2 ? "QPSK" : common.mod_bits == 4 ? "QAM16" : common.mod_bits == 6 ? "QAM64" : "QAM256";
	std::cerr << "symbols: " << common.symbol_count << std::endl;
	std::cerr << "efficiency index: " << (mode & 31) << std::endl;
	std::cerr << "modulation: " << modulation << std::endl;
	std::cerr << "code rate: " << (common.data_bits / common.codewords >> (common.code_order - 4)) << "/16" << std::endl;
	std::cerr << "codewords: " << common.codewords << " of " << (1 << common.code_order) << " bits" << std::endl;
	float duration = 41. / 300. * (3 + common.symbol_count);
	std::cerr << "duration: " << duration << "s" << std::endl;
	std::cerr << "payload: " << common.data_bytes << "B" << std::endl;
	std::cerr << "bitrate: " << common.data_bits / duration / 1000. << "kb/s" << std::endl;
}

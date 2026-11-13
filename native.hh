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
	static const char *const modulations[6] = { "BPSK", "QPSK", nullptr, "QAM16", "QAM64", "QAM256" };
	if (mode & 128 || ((mode >> 4) & 7) > 5 || !modulations[(mode >> 4) & 7]) {
		std::cerr << "oper mode: unsupported" << std::endl;
		return;
	}
	static const char *const code_rates[4] = { "1/2", "11/16", "3/4", "13/16" };
	std::cerr << "modulation: " << modulations[(mode >> 4) & 7] << std::endl;
	std::cerr << "frame size: " << (mode & 1 ? "normal" : "short") << std::endl;
	int code_rate = (mode >> 1) & 7;
	if (code_rate > 3) {
		std::cerr << "code rate: unsupported" << std::endl;
		return;
	}
	std::cerr << "code rate: " << code_rates[code_rate] << std::endl;
	float duration = 41. / 300. * (3 + common.symbol_count);
	std::cerr << "duration: " << duration << "s" << std::endl;
	std::cerr << "payload: " << common.data_bytes << "B" << std::endl;
	std::cerr << "bitrate: " << common.data_bits / duration / 1000. << "kb/s" << std::endl;
}

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

// a guard interval, the copies of the sync symbol, the metadata and payload symbols with their guard intervals and the last guard interval
static float frame_duration(const Common &common)
{
	return (2 + 40 * Common::sync_copies + 41 * (common.symbol_count + 1)) / 300.f;
}

static void print_setup(const Common &common, int mode)
{
	if (mode & 128) {
		float duration = frame_duration(common);
		std::cerr << "analog mode" << std::endl;
		std::cerr << "symbols: " << common.symbol_count << std::endl;
		std::cerr << "duration: " << duration << "s" << std::endl;
		std::cerr << "payload: " << common.symbol_count * common.data_tones << " complex values" << std::endl;
		return;
	}
	static const char *const modulations[8] = { "BPSK", "QPSK", "8PSK", "QAM16", "QAM64", "QAM256", "QAM1024", "QAM4096" };
	static const char *const code_rates[4] = { "1/2", "2/3", "3/4", "5/6" };
	std::cerr << "modulation: " << modulations[(mode >> 4) & 7] << std::endl;
	std::cerr << "frame size: " << (mode & 1 ? "normal" : "short") << std::endl;
	int code_rate = (mode >> 1) & 7;
	if (code_rate > 3) {
		std::cerr << "code rate: unsupported" << std::endl;
		return;
	}
	std::cerr << "code rate: " << code_rates[code_rate] << std::endl;
	float duration = frame_duration(common);
	std::cerr << "duration: " << duration << "s" << std::endl;
	std::cerr << "payload: " << common.data_bytes << "B" << std::endl;
	std::cerr << "bitrate: " << common.data_bits / duration / 1000. << "kb/s" << std::endl;
}

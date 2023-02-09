/*
OFDM modem encoder

Copyright 2023 Ahmet Inan <inan@aicodix.de>
*/

#include <iostream>
#include <cassert>
#include <cmath>
#include "xorshift.hh"
#include "complex.hh"
#include "phasor.hh"
#include "bitman.hh"
#include "utils.hh"
#include "fft.hh"
#include "wav.hh"
#include "pcm.hh"
#include "mls.hh"
#include "crc.hh"
#include "psk.hh"
#include "polar_tables.hh"
#include "polar_helper.hh"
#include "polar_encoder.hh"

struct Encoder
{
	typedef float value;
	typedef DSP::Complex<value> cmplx;
	typedef int8_t code_type;
	static const int mod_bits = 2;
	static const int code_order = 12;
	static const int code_len = 1 << code_order;
	static const int symbol_len = 256;
	static const int guard_len = symbol_len / 8;
	static const int data_bits = 2048;
	static const int mesg_bits = data_bits + 32;
	static const int subcarrier_count = 64;
	static const int payload_symbols = 32;
	static const int first_subcarrier = 16;
	DSP::WritePCM<value> *pcm;
	DSP::FastFourierTransform<symbol_len, cmplx, 1> bwd;
	CODE::CRC<uint32_t> crc;
	CODE::PolarSysEnc<code_type> polarenc;
	code_type code[code_len], mesg[mesg_bits];
	cmplx fdom[symbol_len], tdom[symbol_len], guard[guard_len];

	static int nrz(bool bit)
	{
		return 1 - 2 * bit;
	}
	void symbol()
	{
		bwd(tdom, fdom);
		for (int i = 0; i < symbol_len; ++i)
			tdom[i] /= std::sqrt(value(8*symbol_len));
		for (int i = 0; i < guard_len; ++i) {
			value x = value(i) / value(guard_len - 1);
			x = value(0.5) * (value(1) - std::cos(DSP::Const<value>::Pi() * x));
			guard[i] = DSP::lerp(guard[i], tdom[i+symbol_len-guard_len], x);
		}
		pcm->write(reinterpret_cast<value *>(guard), guard_len, 2);
		pcm->write(reinterpret_cast<value *>(tdom), symbol_len, 2);
		for (int i = 0; i < guard_len; ++i)
			guard[i] = tdom[i];
	}
	void schmidl_cox()
	{
		CODE::MLS seq(0b1100111);
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] = 0;
#if 0
		fdom[first_subcarrier+1] = std::sqrt(value(2 * symbol_len) / value(subcarrier_count));
		for (int i = first_subcarrier + 2; i < first_subcarrier + subcarrier_count + 1; ++i)
			fdom[i] = fdom[i-1] * cmplx(nrz(seq()));
		symbol();
		for (int i = first_subcarrier; i < first_subcarrier + subcarrier_count; ++i)
			fdom[i] = fdom[i+1];
		fdom[first_subcarrier+subcarrier_count] = 0;
		symbol();
#else
		fdom[first_subcarrier] = std::sqrt(value(2 * symbol_len) / value(subcarrier_count));
		for (int i = first_subcarrier + 1; i < first_subcarrier + subcarrier_count; ++i)
			fdom[i] = fdom[i-1] * cmplx(nrz(seq()));
		symbol();
		symbol();
#endif
	}
	cmplx mod_map(code_type *b)
	{
		return PhaseShiftKeying<4, cmplx, code_type>::map(b);
	}
	Encoder(DSP::WritePCM<value> *pcm, const uint8_t *inp) : pcm(pcm), crc(0x8F6E37A0)
	{
		schmidl_cox();
		for (int i = 0; i < data_bits; ++i)
			mesg[i] = nrz(CODE::get_le_bit(inp, i));
		crc.reset();
		for (int i = 0; i < data_bits / 8; ++i)
			crc(inp[i]);
		for (int i = 0; i < 32; ++i)
			mesg[i+data_bits] = nrz((crc()>>i)&1);
		polarenc(code, mesg, frozen_4096_2080, code_order);
		for (int j = 0; j < payload_symbols; ++j) {
			for (int i = 0; i < subcarrier_count; ++i)
				fdom[first_subcarrier+i] *=
					mod_map(code+mod_bits*(subcarrier_count*j+i));
			symbol();
		}
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] = 0;
		symbol();
	}
};

int main(int argc, char **argv)
{
	if (argc != 3) {
		std::cerr << "usage: " << argv[0] << " OUTPUT INPUT" << std::endl;
		return 1;
	}

	const char *output_name = argv[1];
	if (output_name[0] == '-' && output_name[1] == 0)
		output_name = "/dev/stdout";
	const char *input_name = argv[2];
	if (input_name[0] == '-' && input_name[1] == 0)
		input_name = "/dev/stdin";

	std::ifstream input_file(input_name, std::ios::binary);
	if (input_file.bad()) {
		std::cerr << "Couldn't open file \"" << input_name << "\" for reading." << std::endl;
		return 1;
	}
	const int data_len = 2048 / 8;
	uint8_t *input_data = new uint8_t[data_len];
	for (int i = 0; i < data_len; ++i)
		input_data[i] = std::max(input_file.get(), 0);

	CODE::Xorshift32 scrambler;
	for (int i = 0; i < data_len; ++i)
		input_data[i] ^= scrambler();

	const int output_rate = 8000;
	const int output_bits = 16;
	const int output_chan = 1;
	DSP::WriteWAV<float> output_file(output_name, output_rate, output_bits, output_chan);
	output_file.silence(output_rate);
	delete new Encoder(&output_file, input_data);
	output_file.silence(output_rate);
	delete []input_data;

	return 0;
}


/*
OFDM modem decoder

Copyright 2021 Ahmet Inan <inan@aicodix.de>
*/

#include <iostream>
#include <cassert>
#include <cmath>
namespace DSP { using std::abs; using std::min; using std::cos; using std::sin; }
#include "schmidl_cox.hh"
#include "bip_buffer.hh"
#include "theil_sen.hh"
#include "xorshift.hh"
#include "complex.hh"
#include "decibel.hh"
#include "blockdc.hh"
#include "hilbert.hh"
#include "phasor.hh"
#include "bitman.hh"
#include "delay.hh"
#include "wav.hh"
#include "pcm.hh"
#include "fft.hh"
#include "mls.hh"
#include "crc.hh"
#include "psk.hh"
#include "qam.hh"
#include "polar_tables.hh"
#include "polar_list_decoder.hh"
#include "hadamard_encoder.hh"
#include "hadamard_decoder.hh"

template <typename value, typename cmplx, int rate>
struct Decoder
{
	typedef int8_t code_type;
#ifdef __AVX2__
	typedef SIMD<code_type, 32 / sizeof(code_type)> mesg_type;
#else
	typedef SIMD<code_type, 16 / sizeof(code_type)> mesg_type;
#endif
	typedef DSP::Const<value> Const;
	static const int symbol_len = (1280 * rate) / 8000;
	static const int filter_len = (((21 * rate) / 8000) & ~3) | 1;
	static const int guard_len = symbol_len / 8;
	static const int extended_len = symbol_len + guard_len;
	static const int mod_max = 6;
	static const int code_max = 16;
	static const int bits_max = 1 << code_max;
	static const int data_max = 1024;
	static const int cols_max = 256 + 32;
	static const int rows_max = 32;
	static const int cons_max = cols_max * rows_max;
	static const int mls0_len = 320;
	static const int mls0_poly = 0b1100110001;
	static const int mls0_seed = 214;
	static const int mls0_off = - mls0_len / 2;
	static const int mls1_poly = 0b100101011;
	static const int buffer_len = 5 * extended_len;
	static const int search_pos = extended_len;
	DSP::ReadPCM<value> *pcm;
	DSP::FastFourierTransform<symbol_len, cmplx, -1> fwd;
	DSP::BlockDC<value, value> blockdc;
	DSP::Hilbert<cmplx, filter_len> hilbert;
	DSP::BipBuffer<cmplx, buffer_len> input_hist;
	DSP::TheilSenEstimator<value, cols_max> tse;
	SchmidlCox<value, cmplx, search_pos, symbol_len, guard_len> correlator;
	CODE::CRC<uint32_t> crc1;
	CODE::HadamardEncoder<6> hadamardenc;
	CODE::HadamardDecoder<6> hadamarddec;
	CODE::PolarListDecoder<mesg_type, code_max> polardec;
	uint8_t output_data[data_max];
	mesg_type mesg[bits_max];
	code_type code[bits_max], perm[bits_max];
	int8_t mode[32];
	cmplx cons[cons_max], chan[cols_max];
	cmplx fdom[symbol_len], tdom[symbol_len];
	value index[cols_max], phase[cols_max];
	value cfo_rad, sfo_rad;
	const uint32_t *frozen_bits;
	int mod_bits;
	int code_order;
	int symbol_pos;
	int oper_mode;
	int crc_bits;

	static int bin(int carrier)
	{
		return (carrier + symbol_len) % symbol_len;
	}
	static value nrz(bool bit)
	{
		return 1 - 2 * bit;
	}
	static cmplx demod_or_erase(cmplx curr, cmplx prev)
	{
		if (!(norm(prev) > 0))
			return 0;
		cmplx cons = curr / prev;
		if (!(norm(cons) <= 4))
			return 0;
		return cons;
	}
	const cmplx *mls0_seq()
	{
		CODE::MLS seq0(mls0_poly, mls0_seed);
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] = 0;
		value cur = 0, prv = 0;
		for (int i = 0; i < mls0_len; ++i, prv = cur)
			fdom[bin(i+mls0_off)] = prv * (cur = nrz(seq0()));
		return fdom;
	}
	cmplx mod_map(code_type *b)
	{
		switch (mod_bits) {
		case 2:
			return PhaseShiftKeying<4, cmplx, code_type>::map(b);
		case 4:
			return QuadratureAmplitudeModulation<16, cmplx, code_type>::map(b);
		case 6:
			return QuadratureAmplitudeModulation<64, cmplx, code_type>::map(b);
		}
		return 0;
	}
	void mod_hard(code_type *b, cmplx c)
	{
		switch (mod_bits) {
		case 2:
			return PhaseShiftKeying<4, cmplx, code_type>::hard(b, c);
		case 4:
			return QuadratureAmplitudeModulation<16, cmplx, code_type>::hard(b, c);
		case 6:
			return QuadratureAmplitudeModulation<64, cmplx, code_type>::hard(b, c);
		}
	}
	void mod_soft(code_type *b, cmplx c, value precision)
	{
		switch (mod_bits) {
		case 2:
			return PhaseShiftKeying<4, cmplx, code_type>::soft(b, c, precision);
		case 4:
			return QuadratureAmplitudeModulation<16, cmplx, code_type>::soft(b, c, precision);
		case 6:
			return QuadratureAmplitudeModulation<64, cmplx, code_type>::soft(b, c, precision);
		}
	}
	void shuffle(code_type *dest, const code_type *src)
	{
		if (code_order == 11) {
			CODE::XorShiftMask<int, 11, 1, 3, 4, 1> seq;
			dest[0] = src[0];
			for (int i = 1; i < 2048; ++i)
				dest[seq()] = src[i];
		} else if (code_order == 12) {
			CODE::XorShiftMask<int, 12, 1, 1, 4, 1> seq;
			dest[0] = src[0];
			for (int i = 1; i < 4096; ++i)
				dest[seq()] = src[i];
		} else if (code_order == 13) {
			CODE::XorShiftMask<int, 13, 1, 1, 9, 1> seq;
			dest[0] = src[0];
			for (int i = 1; i < 8192; ++i)
				dest[seq()] = src[i];
		} else if (code_order == 14) {
			CODE::XorShiftMask<int, 14, 1, 5, 10, 1> seq;
			dest[0] = src[0];
			for (int i = 1; i < 16384; ++i)
				dest[seq()] = src[i];
		} else if (code_order == 15) {
			CODE::XorShiftMask<int, 15, 1, 1, 3, 1> seq;
			dest[0] = src[0];
			for (int i = 1; i < 32768; ++i)
				dest[seq()] = src[i];
		} else if (code_order == 16) {
			CODE::XorShiftMask<int, 16, 1, 1, 14, 1> seq;
			dest[0] = src[0];
			for (int i = 1; i < 65536; ++i)
				dest[seq()] = src[i];
		}
	}
	const cmplx *next_sample()
	{
		cmplx tmp;
		pcm->read(reinterpret_cast<value *>(&tmp), 1);
		if (pcm->channels() == 1)
			tmp = hilbert(blockdc(tmp.real()));
		return input_hist(tmp);
	}
	Decoder(DSP::ReadPCM<value> *pcm, const char *const *output_names, int output_count) :
		pcm(pcm), correlator(mls0_seq()), crc1(0x8F6E37A0)
	{
		blockdc.samples(filter_len);
		DSP::Phasor<cmplx> osc;
		const cmplx *buf;
		int output_index = 0;
		int comb_cols = 32;
		int code_cols = 256;
		int comb_dist = 9;
		int comb_off = 4;
		int cons_cols = code_cols + comb_cols;
		int code_off = - cons_cols / 2;
		while (output_index < output_count) {
			do {
				if (!pcm->good())
					return;
				buf = next_sample();
			} while (!correlator(buf));

			symbol_pos = correlator.symbol_pos;
			cfo_rad = correlator.cfo_rad;
			std::cerr << "symbol pos: " << symbol_pos << std::endl;
			std::cerr << "coarse cfo: " << cfo_rad * (rate / Const::TwoPi()) << " Hz " << std::endl;

			osc.omega(-cfo_rad);
			for (int i = 0; i < symbol_len; ++i)
				tdom[i] = buf[i+symbol_pos+symbol_len] * osc();
			for (int i = 0; i < guard_len; ++i)
				osc();
			fwd(fdom, tdom);
			CODE::MLS seq0(mls0_poly, mls0_seed);
			for (int i = 0; i < mls0_len; ++i)
				fdom[bin(i+mls0_off)] *= nrz(seq0());
			for (int i = 0; i < cons_cols; ++i)
				chan[i] = fdom[bin(i+code_off)];
			for (int i = 0; i < symbol_len; ++i)
				tdom[i] = buf[i+symbol_pos+symbol_len+extended_len] * osc();
			for (int i = 0; i < guard_len; ++i)
				osc();
			fwd(fdom, tdom);
			CODE::MLS seq1(mls1_poly);
			auto clamp = [](int v){ return v < -127 ? -127 : v > 127 ? 127 : v; };
			for (int i = 0; i < comb_cols; ++i)
				mode[i] = clamp(std::nearbyint(127 * demod_or_erase(fdom[bin(i*comb_dist+comb_off+code_off)], chan[i]).real() * nrz(seq1())));
			oper_mode = hadamarddec(mode);
			if (oper_mode < 0 || oper_mode > 8) {
				std::cerr << "operation mode " << oper_mode << " unsupported." << std::endl;
				continue;
			}
			std::cerr << "oper mode: " << oper_mode << std::endl;
			if (!oper_mode)
				continue;

			int data_bits = 0;
			int cons_rows = 0;
			switch (oper_mode) {
			case 1:
				mod_bits = 2;
				cons_rows = 8;
				code_order = 12;
				data_bits = 2048;
				frozen_bits = frozen_4096_2080;
				break;
			case 2:
				mod_bits = 2;
				cons_rows = 16;
				code_order = 13;
				data_bits = 4096;
				frozen_bits = frozen_8192_4128;
				break;
			case 3:
				mod_bits = 2;
				cons_rows = 32;
				code_order = 14;
				data_bits = 8192;
				frozen_bits = frozen_16384_8224;
				break;
			case 4:
				mod_bits = 4;
				cons_rows = 4;
				code_order = 12;
				data_bits = 2048;
				frozen_bits = frozen_4096_2080;
				break;
			case 5:
				mod_bits = 4;
				cons_rows = 8;
				code_order = 13;
				data_bits = 4096;
				frozen_bits = frozen_8192_4128;
				break;
			case 6:
				mod_bits = 4;
				cons_rows = 16;
				code_order = 14;
				data_bits = 8192;
				frozen_bits = frozen_16384_8224;
				break;
			case 7:
				mod_bits = 6;
				cons_rows = 6;
				code_order = 13;
				data_bits = 4096;
				frozen_bits = frozen_8192_4128;
				break;
			case 8:
				mod_bits = 6;
				cons_rows = 11;
				code_order = 14;
				data_bits = 8192;
				frozen_bits = frozen_16384_8224;
				break;
			default:
				return;
			}
			int data_bytes = data_bits / 8;

			std::cerr << "demod ";
			for (int j = 0; j < cons_rows; ++j) {
				if (j) {
					for (int i = 0; i < extended_len; ++i)
						correlator(buf = next_sample());
					for (int i = 0; i < symbol_len; ++i)
						tdom[i] = buf[i] * osc();
					for (int i = 0; i < guard_len; ++i)
						osc();
					fwd(fdom, tdom);
				} else {
					for (int i = 0; i < symbol_pos+symbol_len+extended_len; ++i)
						correlator(buf = next_sample());
					seq1.reset();
					hadamardenc(mode, oper_mode);
				}
				for (int i = 0; i < comb_cols; ++i)
					fdom[bin(comb_dist*i+comb_off+code_off)] *= nrz(seq1()) * mode[i];
				for (int i = 0; i < cons_cols; ++i)
					cons[cons_cols*j+i] = demod_or_erase(fdom[bin(i+code_off)], chan[i]);
				for (int i = 0; i < comb_cols; ++i) {
					index[i] = code_off + comb_dist * i + comb_off;
					phase[i] = arg(cons[cons_cols*j+comb_dist*i+comb_off]);
				}
				tse.compute(index, phase, comb_cols);
				//std::cerr << "Theil-Sen slope = " << tse.slope() << std::endl;
				//std::cerr << "Theil-Sen yint = " << tse.yint() << std::endl;
				for (int i = 0; i < cons_cols; ++i)
					cons[cons_cols*j+i] *= DSP::polar<value>(1, -tse(i+code_off));
				for (int i = 0; i < cons_cols; ++i)
					if (i % comb_dist == comb_off)
						chan[i] = fdom[bin(i+code_off)];
					else
						chan[i] *= DSP::polar<value>(1, tse(i+code_off));
				std::cerr << ".";
			}
			std::cerr << " done" << std::endl;
			std::cerr << "Es/N0 (dB):";
			value sp = 0, np = 0;
			for (int j = 0, k = 0; j < cons_rows; ++j) {
				for (int i = 0; i < comb_cols; ++i) {
					cmplx hard(1, 0);
					cmplx error = cons[cons_cols*j+comb_dist*i+comb_off] - hard;
					sp += norm(hard);
					np += norm(error);
				}
				value precision = sp / np;
				// precision = 8;
				value snr = DSP::decibel(precision);
				std::cerr << " " << snr;
				if (std::is_same<code_type, int8_t>::value && precision > 32)
					precision = 32;
				for (int i = 0; i < cons_cols; ++i) {
					if (i % comb_dist == comb_off)
						continue;
					mod_soft(perm+k, cons[cons_cols*j+i], precision);
					k += mod_bits;
				}
			}
			std::cerr << std::endl;
			crc_bits = data_bits + 32;
			for (int i = code_cols * cons_rows * mod_bits; i < bits_max; ++i)
				perm[i] = 0;
			shuffle(code, perm);
			polardec(nullptr, mesg, code, frozen_bits, code_order);
			int best = -1;
			for (int k = 0; k < mesg_type::SIZE; ++k) {
				crc1.reset();
				for (int i = 0; i < crc_bits; ++i)
					crc1(mesg[i].v[k] < 0);
				if (crc1() == 0) {
					best = k;
					break;
				}
			}
			if (best < 0) {
				std::cerr << "payload decoding error." << std::endl;
				continue;
			}
			for (int i = 0; i < data_bits; ++i)
				CODE::set_le_bit(output_data, i, mesg[i].v[best] < 0);

			const char *output_name = output_names[output_index++];
			if (output_count == 1 && output_name[0] == '-' && output_name[1] == 0)
				output_name = "/dev/stdout";
			std::ofstream output_file(output_name, std::ios::binary | std::ios::trunc);
			if (output_file.bad()) {
				std::cerr << "Couldn't open file \"" << output_name << "\" for writing." << std::endl;
				continue;
			}
			CODE::Xorshift32 scrambler;
			for (int i = 0; i < data_bytes; ++i)
				output_data[i] ^= scrambler();
			for (int i = 0; i < data_bytes; ++i)
				output_file.put(output_data[i]);
		}
	}
};

int main(int argc, char **argv)
{
	if (argc < 3) {
		std::cerr << "usage: " << argv[0] << " INPUT OUTPUT.." << std::endl;
		return 1;
	}

	typedef float value;
	typedef DSP::Complex<value> cmplx;

	const char *input_name = argv[1];
	if (input_name[0] == '-' && input_name[1] == 0)
		input_name = "/dev/stdin";

	DSP::ReadWAV<value> input_file(input_name);

	if (input_file.channels() < 1 || input_file.channels() > 2) {
		std::cerr << "Only real or analytic signal (one or two channels) supported." << std::endl;
		return 1;
	}

	int output_count = argc - 2;
	switch (input_file.rate()) {
	case 8000:
		delete new Decoder<value, cmplx, 8000>(&input_file, argv+2, output_count);
		break;
	case 16000:
		delete new Decoder<value, cmplx, 16000>(&input_file, argv+2, output_count);
		break;
	case 44100:
		delete new Decoder<value, cmplx, 44100>(&input_file, argv+2, output_count);
		break;
	case 48000:
		delete new Decoder<value, cmplx, 48000>(&input_file, argv+2, output_count);
		break;
	default:
		std::cerr << "Unsupported sample rate." << std::endl;
		return 1;
	}
	return 0;
}


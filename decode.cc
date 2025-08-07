/*
OFDM modem decoder

Copyright 2021 Ahmet Inan <inan@aicodix.de>
*/

#include <iomanip>
#include <iostream>
#include <cstdint>
#include <cassert>
#include <cmath>
namespace DSP { using std::abs; using std::min; using std::cos; using std::sin; }
#include "common.hh"
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
#include "psk.hh"
#include "qam.hh"
#include "polar_list_decoder.hh"

template <typename value, typename cmplx, int rate>
struct Decoder : Common
{
	typedef int16_t code_type;
	typedef SIMD<code_type, 32> mesg_type;
	typedef DSP::Const<value> Const;
	static const int guard_len = rate / 300;
	static const int symbol_len = guard_len * 40;
	static const int filter_len = 129;
	static const int extended_len = symbol_len + guard_len;
	static const int buffer_len = 5 * extended_len;
	static const int search_pos = extended_len;
	static const int tone_off = - tone_count / 2;
	DSP::ReadPCM<value> *pcm;
	DSP::FastFourierTransform<symbol_len, cmplx, -1> fwd;
	DSP::BlockDC<value, value> blockdc;
	DSP::Hilbert<cmplx, filter_len> hilbert;
	DSP::BipBuffer<cmplx, buffer_len> input_hist;
	DSP::TheilSenEstimator<value, tone_count> tse;
	SchmidlCox<value, cmplx, search_pos, symbol_len, guard_len> correlator;
	CODE::PolarListDecoder<mesg_type, code_max> polar_decoder;
	mesg_type mesg[bits_max];
	code_type code[bits_max], perm[bits_max];
	cmplx demod[tone_count], chan[tone_count], tone[tone_count];
	cmplx fdom[symbol_len], tdom[symbol_len];
	value index[tone_count], phase[tone_count];
	value snr[symbols_max];
	value cfo_rad, sfo_rad;
	int symbol_pos;
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
		if (norm(prev) > 0) {
			cmplx demod = curr / prev;
			if (norm(demod) < 4)
				return demod;
		}
		return 0;
	}
	static void base37_decoder(char *str, int64_t val, int len)
	{
		for (int i = len-1; i >= 0; --i, val /= 37)
			str[i] = " 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"[val%37];
	}
	const cmplx *mls0_seq()
	{
		CODE::MLS seq0(mls0_poly, mls0_seed);
		value cur = 0, prv = 0;
		for (int i = 0; i < tone_count; ++i, prv = cur)
			fdom[bin(i+tone_off)] = prv * (cur = nrz(seq0()));
		return fdom;
	}
	cmplx map_bits(code_type *b, int bits)
	{
		switch (bits) {
		case 1:
			return PhaseShiftKeying<2, cmplx, code_type>::map(b);
		case 2:
			return PhaseShiftKeying<4, cmplx, code_type>::map(b);
		case 3:
			return PhaseShiftKeying<8, cmplx, code_type>::map(b);
		case 4:
			return QuadratureAmplitudeModulation<16, cmplx, code_type>::map(b);
		case 6:
			return QuadratureAmplitudeModulation<64, cmplx, code_type>::map(b);
		case 8:
			return QuadratureAmplitudeModulation<256, cmplx, code_type>::map(b);
		case 10:
			return QuadratureAmplitudeModulation<1024, cmplx, code_type>::map(b);
		case 12:
			return QuadratureAmplitudeModulation<4096, cmplx, code_type>::map(b);
		}
		return 0;
	}
	void demap_soft(code_type *b, cmplx c, value precision, int bits)
	{
		switch (bits) {
		case 1:
			return PhaseShiftKeying<2, cmplx, code_type>::soft(b, c, precision);
		case 2:
			return PhaseShiftKeying<4, cmplx, code_type>::soft(b, c, precision);
		case 3:
			return PhaseShiftKeying<8, cmplx, code_type>::soft(b, c, precision);
		case 4:
			return QuadratureAmplitudeModulation<16, cmplx, code_type>::soft(b, c, precision);
		case 6:
			return QuadratureAmplitudeModulation<64, cmplx, code_type>::soft(b, c, precision);
		case 8:
			return QuadratureAmplitudeModulation<256, cmplx, code_type>::soft(b, c, precision);
		case 10:
			return QuadratureAmplitudeModulation<1024, cmplx, code_type>::soft(b, c, precision);
		case 12:
			return QuadratureAmplitudeModulation<4096, cmplx, code_type>::soft(b, c, precision);
		}
	}
	void demap_hard(code_type *b, cmplx c, int bits)
	{
		switch (bits) {
		case 1:
			return PhaseShiftKeying<2, cmplx, code_type>::hard(b, c);
		case 2:
			return PhaseShiftKeying<4, cmplx, code_type>::hard(b, c);
		case 3:
			return PhaseShiftKeying<8, cmplx, code_type>::hard(b, c);
		case 4:
			return QuadratureAmplitudeModulation<16, cmplx, code_type>::hard(b, c);
		case 6:
			return QuadratureAmplitudeModulation<64, cmplx, code_type>::hard(b, c);
		case 8:
			return QuadratureAmplitudeModulation<256, cmplx, code_type>::hard(b, c);
		case 10:
			return QuadratureAmplitudeModulation<1024, cmplx, code_type>::hard(b, c);
		case 12:
			return QuadratureAmplitudeModulation<4096, cmplx, code_type>::hard(b, c);
		}
	}
	void shuffle(code_type *dest, const code_type *src, int order)
	{
		if (order == 8) {
			CODE::XorShiftMask<int, 8, 1, 1, 2, 1> seq;
			dest[0] = src[0];
			for (int i = 1; i < 256; ++i)
				dest[seq()] = src[i];
		} else if (order == 11) {
			CODE::XorShiftMask<int, 11, 1, 3, 4, 1> seq;
			dest[0] = src[0];
			for (int i = 1; i < 2048; ++i)
				dest[seq()] = src[i];
		} else if (order == 12) {
			CODE::XorShiftMask<int, 12, 1, 1, 4, 1> seq;
			dest[0] = src[0];
			for (int i = 1; i < 4096; ++i)
				dest[seq()] = src[i];
		} else if (order == 13) {
			CODE::XorShiftMask<int, 13, 1, 1, 9, 1> seq;
			dest[0] = src[0];
			for (int i = 1; i < 8192; ++i)
				dest[seq()] = src[i];
		} else if (order == 14) {
			CODE::XorShiftMask<int, 14, 1, 5, 10, 1> seq;
			dest[0] = src[0];
			for (int i = 1; i < 16384; ++i)
				dest[seq()] = src[i];
		} else if (order == 15) {
			CODE::XorShiftMask<int, 15, 1, 1, 3, 1> seq;
			dest[0] = src[0];
			for (int i = 1; i < 32768; ++i)
				dest[seq()] = src[i];
		} else if (order == 16) {
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
	int64_t meta_data()
	{
		shuffle(code, perm, 8);
		polar_decoder(nullptr, mesg, code, frozen_256_71, 8);
		int best = -1;
		for (int k = 0; k < mesg_type::SIZE; ++k) {
			crc0.reset();
			for (int i = 0; i < 71; ++i)
				crc0(mesg[i].v[k] < 0);
			if (crc0() == 0) {
				best = k;
				break;
			}
		}
		if (best < 0)
			return -1;
		uint64_t md = 0;
		for (int i = 0; i < 55; ++i)
			md |= uint64_t(mesg[i].v[best] < 0) << i;
		return md;
	}
	Decoder(DSP::ReadPCM<value> *pcm, const char *const *output_names, int output_count) : pcm(pcm), correlator(mls0_seq())
	{
		blockdc.samples(filter_len);
		DSP::Phasor<cmplx> osc;
		const cmplx *buf;
		int output_index = 0;
		int sample_count = 0;
		while (output_index < output_count) {
			do {
				if (!pcm->good())
					return;
				buf = next_sample();
				++sample_count;
			} while (!correlator(buf));

			symbol_pos = correlator.symbol_pos;
			cfo_rad = correlator.cfo_rad;
			std::cerr << "symbol pos: " << sample_count - buffer_len + symbol_pos << std::endl;
			std::cerr << "coarse cfo: " << cfo_rad * (rate / Const::TwoPi()) << " Hz " << std::endl;

			osc.omega(-cfo_rad);
			for (int i = 0; i < symbol_len; ++i)
				tdom[i] = buf[i+symbol_pos] * osc();
			fwd(fdom, tdom);
			for (int i = 0; i < tone_count; ++i)
				tone[i] = fdom[bin(i+tone_off)];
			for (int i = 0; i < symbol_len; ++i)
				tdom[i] = buf[i+symbol_pos+symbol_len] * osc();
			for (int i = 0; i < guard_len; ++i)
				osc();
			fwd(fdom, tdom);
			for (int i = 0; i < tone_count; ++i)
				chan[i] = fdom[bin(i+tone_off)];
			for (int i = 0; i < tone_count; ++i) {
				index[i] = tone_off + i;
				phase[i] = arg(demod_or_erase(chan[i], tone[i]));
			}
			tse.compute(index, phase, tone_count);
			std::cerr << "coarse sfo: " << -1000000 * tse.slope() / Const::TwoPi() << " ppm" << std::endl;
			std::cerr << "residual cfo: " << 1000 * tse.yint() * rate / (Const::TwoPi() * symbol_len) << " mHz" << std::endl;
			for (int i = 0; i < tone_count; ++i)
				tone[i] *= DSP::polar<value>(1, tse(i+tone_off));
			for (int i = 0; i < tone_count; ++i)
				chan[i] = DSP::lerp(chan[i], tone[i], value(0.5));
			CODE::MLS seq0(mls0_poly, mls0_seed);
			for (int i = 0; i < tone_count; ++i)
				chan[i] *= nrz(seq0());
			for (int i = 0; i < symbol_len; ++i)
				tdom[i] = buf[i+symbol_pos+symbol_len+extended_len] * osc();
			for (int i = 0; i < guard_len; ++i)
				osc();
			fwd(fdom, tdom);
			CODE::MLS seq1(mls1_poly);
			mod_bits = 1;
			oper_mode = -1;
			symbol_count = 0;
			for (int j = 0, k = 0; j < symbol_count + 1; ++j) {
				pilot_off = (block_skew * j + first_pilot) % block_length;
				if (j) {
					for (int i = 0; i < extended_len; ++i)
						correlator(buf = next_sample());
					for (int i = 0; i < symbol_len; ++i)
						tdom[i] = buf[i] * osc();
					for (int i = 0; i < guard_len; ++i)
						osc();
					fwd(fdom, tdom);
				}
				for (int i = 0; i < tone_count; ++i)
					tone[i] = fdom[bin(i+tone_off)];
				for (int i = pilot_off; i < tone_count; i += block_length)
					tone[i] *= nrz(seq1());
				for (int i = 0; i < tone_count; ++i)
					demod[i] = demod_or_erase(tone[i], chan[i]);
				auto rot = [](cmplx p){ return std::abs(p.real()) < std::abs(p.imag()) ? cmplx(0, p.imag() < 0 ? 1 : -1) : cmplx(p.real() < 0 ? -1 : 1); };
				cmplx ptsa_sum;
				for (int i = 0; i < tone_count; i += 8)
					if (i % block_length == pilot_off)
						ptsa_sum += demod[i];
				cmplx ptsa_phase = rot(ptsa_sum);
				if (abs(ptsa_sum) < pilot_tones / (8 * 2)) {
					std::cerr << "pilot phase damaged" << std::endl;
					oper_mode = -1;
					break;
				}
				for (int i = 0; i < tone_count; i += 8) {
					tone[i] *= ptsa_phase;
					demod[i] *= ptsa_phase;
				}
				cmplx ptsb_sum;
				for (int i = 1; i < tone_count; i += 8)
					if (i % block_length == pilot_off)
						ptsb_sum += demod[i];
				cmplx ptsb_phase = rot(ptsb_sum);
				if (abs(ptsb_sum) < pilot_tones / (8 * 2)) {
					std::cerr << "pilot phase damaged" << std::endl;
					oper_mode = -1;
					break;
				}
				for (int i = 1; i < tone_count; i += 8) {
					tone[i] *= ptsb_phase;
					demod[i] *= ptsb_phase;
				}
				cmplx ptsc_sum;
				for (int i = 2; i < tone_count; i += 8)
					if (i % block_length == pilot_off)
						ptsc_sum += demod[i];
				cmplx ptsc_phase = rot(ptsc_sum);
				if (abs(ptsc_sum) < pilot_tones / (8 * 2)) {
					std::cerr << "pilot phase damaged" << std::endl;
					oper_mode = -1;
					break;
				}
				for (int i = 2; i < tone_count; i += 8) {
					tone[i] *= ptsc_phase;
					demod[i] *= ptsc_phase;
				}
				cmplx ptsd_sum;
				for (int i = 3; i < tone_count; i += 8)
					if (i % block_length == pilot_off)
						ptsd_sum += demod[i];
				cmplx ptsd_phase = rot(ptsd_sum);
				if (abs(ptsd_sum) < pilot_tones / (8 * 2)) {
					std::cerr << "pilot phase damaged" << std::endl;
					oper_mode = -1;
					break;
				}
				for (int i = 3; i < tone_count; i += 8) {
					tone[i] *= ptsd_phase;
					demod[i] *= ptsd_phase;
				}
				cmplx ptse_sum;
				for (int i = 4; i < tone_count; i += 8)
					if (i % block_length == pilot_off)
						ptse_sum += demod[i];
				cmplx ptse_phase = rot(ptse_sum);
				if (abs(ptse_sum) < pilot_tones / (8 * 2)) {
					std::cerr << "pilot phase damaged" << std::endl;
					oper_mode = -1;
					break;
				}
				for (int i = 4; i < tone_count; i += 8) {
					tone[i] *= ptse_phase;
					demod[i] *= ptse_phase;
				}
				cmplx ptsf_sum;
				for (int i = 5; i < tone_count; i += 8)
					if (i % block_length == pilot_off)
						ptsf_sum += demod[i];
				cmplx ptsf_phase = rot(ptsf_sum);
				if (abs(ptsf_sum) < pilot_tones / (8 * 2)) {
					std::cerr << "pilot phase damaged" << std::endl;
					oper_mode = -1;
					break;
				}
				for (int i = 5; i < tone_count; i += 8) {
					tone[i] *= ptsf_phase;
					demod[i] *= ptsf_phase;
				}
				cmplx ptsg_sum;
				for (int i = 6; i < tone_count; i += 8)
					if (i % block_length == pilot_off)
						ptsg_sum += demod[i];
				cmplx ptsg_phase = rot(ptsg_sum);
				if (abs(ptsg_sum) < pilot_tones / (8 * 2)) {
					std::cerr << "pilot phase damaged" << std::endl;
					oper_mode = -1;
					break;
				}
				for (int i = 6; i < tone_count; i += 8) {
					tone[i] *= ptsg_phase;
					demod[i] *= ptsg_phase;
				}
				cmplx ptsh_sum;
				for (int i = 7; i < tone_count; i += 8)
					if (i % block_length == pilot_off)
						ptsh_sum += demod[i];
				cmplx ptsh_phase = rot(ptsh_sum);
				if (abs(ptsh_sum) < pilot_tones / (8 * 2)) {
					std::cerr << "pilot phase damaged" << std::endl;
					oper_mode = -1;
					break;
				}
				for (int i = 7; i < tone_count; i += 8) {
					tone[i] *= ptsh_phase;
					demod[i] *= ptsh_phase;
				}
				for (int i = 0; i < pilot_tones; ++i) {
					index[i] = tone_off + block_length * i + pilot_off;
					phase[i] = arg(demod[block_length*i+pilot_off]);
				}
				tse.compute(index, phase, pilot_tones);
				//std::cerr << "Theil-Sen slope = " << tse.slope() << std::endl;
				//std::cerr << "Theil-Sen yint = " << tse.yint() << std::endl;
				for (int i = 0; i < tone_count; ++i)
					demod[i] *= DSP::polar<value>(1, -tse(i+tone_off));
				for (int i = 0; i < tone_count; ++i)
					chan[i] *= DSP::polar<value>(1, tse(i+tone_off));
				value sp = 0, np = 0;
				for (int i = 0, l = k; i < tone_count; ++i) {
					cmplx hard(1, 0);
					if (i % block_length != pilot_off) {
						int bits = mod_bits;
						if (mod_bits == 3 && l % 32 == 30)
							bits = 2;
						if (mod_bits == 6 && l % 64 == 60)
							bits = 4;
						if (mod_bits == 10 && l % 128 == 120)
							bits = 8;
						if (mod_bits == 12 && l % 128 == 120)
							bits = 8;
						demap_hard(perm+l, demod[i], bits);
						hard = map_bits(perm+l, bits);
						l += bits;
					}
					cmplx error = demod[i] - hard;
					sp += norm(hard);
					np += norm(error);
				}
				value precision = sp / np;
				snr[j] = precision;
				precision = std::min(precision, value(1023));
				for (int i = 0; i < tone_count; ++i) {
					if (i % block_length != pilot_off) {
						int bits = mod_bits;
						if (mod_bits == 3 && k % 32 == 30)
							bits = 2;
						if (mod_bits == 6 && k % 64 == 60)
							bits = 4;
						if (mod_bits == 10 && k % 128 == 120)
							bits = 8;
						if (mod_bits == 12 && k % 128 == 120)
							bits = 8;
						demap_soft(perm+k, demod[i], precision, bits);
						k += bits;
					}
				}
				if (!j) {
					int64_t meta_info = meta_data();
					if (meta_info < 0) {
						std::cerr << "preamble decoding error." << std::endl;
						break;
					}
					int64_t call = meta_info >> 8;
					if (call == 0 || call >= 129961739795077L) {
						std::cerr << "call sign unsupported." << std::endl;
						break;
					}
					char call_sign[10];
					base37_decoder(call_sign, call, 9);
					call_sign[9] = 0;
					std::cerr << "call sign: " << call_sign << std::endl;
					int mode = meta_info & 255;
					if (!setup(mode))
						break;
					k = 0;
					for (int i = 0; i < symbol_pos+symbol_len+extended_len; ++i)
						correlator(buf = next_sample());
					std::cerr << "oper mode: " << oper_mode << std::endl;
				}
				for (int i = pilot_off; i < tone_count; i += block_length)
					chan[i] = DSP::lerp(chan[i], tone[i], value(0.5));
			}
			if (oper_mode < 0)
				continue;
			DSP::quick_sort(snr, symbol_count + 1);
			std::cerr << "Es/N0 (dB): " << DSP::decibel(snr[0]) << " .. " << DSP::decibel(snr[symbol_count/2]) << " .. " << DSP::decibel(snr[symbol_count]) << std::endl;
			crc_bits = data_bits + 32;
			shuffle(code, perm, code_order);
			polar_decoder(nullptr, mesg, code, frozen_bits, code_order);
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
				CODE::set_le_bit(data, i, mesg[i].v[best] < 0);

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
				data[i] ^= scrambler();
			for (int i = 0; i < data_bytes; ++i)
				output_file.put(data[i]);
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

	std::cerr << std::fixed << std::setprecision(1);
	int output_count = argc - 2;
	switch (input_file.rate()) {
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


/*
OFDM modem decoder

Copyright 2021 Ahmet Inan <inan@aicodix.de>
*/

#include <algorithm>
#include <iostream>
#include <cassert>
#include <cmath>
namespace DSP { using std::abs; using std::min; using std::cos; using std::sin; }
#include "bip_buffer.hh"
#include "theil_sen.hh"
#include "xorshift.hh"
#include "trigger.hh"
#include "complex.hh"
#include "permute.hh"
#include "decibel.hh"
#include "blockdc.hh"
#include "hilbert.hh"
#include "phasor.hh"
#include "bitman.hh"
#include "delay.hh"
#include "sma.hh"
#include "wav.hh"
#include "pcm.hh"
#include "fft.hh"
#include "mls.hh"
#include "crc.hh"
#include "osd.hh"
#include "psk.hh"
#include "qam.hh"
#include "polar_tables.hh"
#include "polar_helper.hh"
#include "polar_encoder.hh"
#include "polar_list_decoder.hh"

template <typename value, typename cmplx, int search_pos, int symbol_len, int guard_len>
struct SchmidlCox
{
	typedef DSP::Const<value> Const;
	static const int match_len = guard_len | 1;
	static const int match_del = (match_len - 1) / 2;
	DSP::FastFourierTransform<symbol_len, cmplx, -1> fwd;
	DSP::FastFourierTransform<symbol_len, cmplx, 1> bwd;
	DSP::SMA4<cmplx, value, symbol_len, false> cor;
	DSP::SMA4<value, value, 2*symbol_len, false> pwr;
	DSP::SMA4<value, value, match_len, false> match;
	DSP::Delay<value, match_del> align;
	DSP::SchmittTrigger<value> threshold;
	DSP::FallingEdgeTrigger falling;
	cmplx tmp0[symbol_len], tmp1[symbol_len];
	cmplx kern[symbol_len];
	cmplx cmplx_shift = 0;
	value timing_max = 0;
	value phase_max = 0;
	int index_max = 0;

	static int bin(int carrier)
	{
		return (carrier + symbol_len) % symbol_len;
	}
	static cmplx demod_or_erase(cmplx curr, cmplx prev, value pwr)
	{
		if (!(norm(curr) > pwr))
			return 0;
		if (!(norm(prev) > pwr))
			return 0;
		cmplx cons = curr / prev;
		if (!(norm(cons) <= 4))
			return 0;
		return cons;
	}
public:
	int symbol_pos = 0;
	value cfo_rad = 0;
	value frac_cfo = 0;

	SchmidlCox(const cmplx *sequence) : threshold(value(0.17*match_len), value(0.19*match_len))
	{
		fwd(kern, sequence);
		for (int i = 0; i < symbol_len; ++i)
			kern[i] = conj(kern[i]) / value(symbol_len);
	}
	bool operator()(const cmplx *samples)
	{
		cmplx P = cor(samples[search_pos+symbol_len] * conj(samples[search_pos+2*symbol_len]));
		value R = value(0.5) * pwr(norm(samples[search_pos+2*symbol_len]));
		value min_R = 0.00001 * symbol_len;
		R = std::max(R, min_R);
		value timing = match(norm(P) / (R * R));
		value phase = align(arg(P));

		bool collect = threshold(timing);
		bool process = falling(collect);

		if (!collect && !process)
			return false;

		if (timing_max < timing) {
			timing_max = timing;
			phase_max = phase;
			index_max = match_del;
		} else if (index_max < symbol_len + guard_len + match_del) {
			++index_max;
		}

		if (!process)
			return false;

		frac_cfo = phase_max / value(symbol_len);

		DSP::Phasor<cmplx> osc;
		osc.omega(frac_cfo);
		symbol_pos = search_pos - index_max;
		index_max = 0;
		timing_max = 0;
		for (int i = 0; i < symbol_len; ++i)
			tmp1[i] = samples[i+symbol_pos+symbol_len] * osc();
		fwd(tmp0, tmp1);
		value min_pwr = 0;
		for (int i = 0; i < symbol_len; ++i)
			min_pwr += norm(tmp0[i]);
		min_pwr /= symbol_len;
		for (int i = 0; i < symbol_len; ++i)
			tmp1[i] = demod_or_erase(tmp0[i], tmp0[bin(i-1)], min_pwr);
		fwd(tmp0, tmp1);
		for (int i = 0; i < symbol_len; ++i)
			tmp0[i] *= kern[i];
		bwd(tmp1, tmp0);

		int shift = 0;
		value peak = 0;
		value next = 0;
		for (int i = 0; i < symbol_len; ++i) {
			value power = norm(tmp1[i]);
			if (power > peak) {
				next = peak;
				peak = power;
				shift = i;
			} else if (power > next) {
				next = power;
			}
		}
		if (peak <= next * 4)
			return false;

		int pos_err = std::nearbyint(arg(tmp1[shift]) * symbol_len / Const::TwoPi());
		if (abs(pos_err) > guard_len / 2)
			return false;
		symbol_pos -= pos_err;

		cfo_rad = shift * (Const::TwoPi() / symbol_len) - frac_cfo;
		if (cfo_rad >= Const::Pi())
			cfo_rad -= Const::TwoPi();
		return true;
	}
};

void base37_decoder(char *str, long long int val, int len)
{
	for (int i = len-1; i >= 0; --i, val /= 37)
		str[i] = " 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"[val%37];
}

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
	static const int code_max = 14;
	static const int bits_max = 1 << code_max;
	static const int cols_max = 273 + 16;
	static const int rows_max = 32;
	static const int cons_max = cols_max * rows_max;
	static const int mls0_len = 127;
	static const int mls0_off = - mls0_len + 1;
	static const int mls0_poly = 0b10001001;
	static const int mls1_len = 255;
	static const int mls1_off = - mls1_len / 2;
	static const int mls1_poly = 0b100101011;
	static const int buffer_len = 4 * extended_len;
	static const int search_pos = extended_len;
	DSP::ReadPCM<value> *pcm;
	DSP::FastFourierTransform<symbol_len, cmplx, -1> fwd;
	DSP::BlockDC<value, value> blockdc;
	DSP::Hilbert<cmplx, filter_len> hilbert;
	DSP::BipBuffer<cmplx, buffer_len> input_hist;
	DSP::TheilSenEstimator<value, cols_max> tse;
	SchmidlCox<value, cmplx, search_pos, symbol_len/2, guard_len> correlator;
	CODE::CRC<uint16_t> crc0;
	CODE::CRC<uint32_t> crc1;
	CODE::OrderedStatisticsDecoder<255, 71, 4> osddec;
	CODE::PolarEncoder<mesg_type> polarenc;
	CODE::PolarListDecoder<mesg_type, code_max> polardec;
	CODE::ReverseFisherYatesShuffle<4096> shuffle_4096;
	CODE::ReverseFisherYatesShuffle<8192> shuffle_8192;
	CODE::ReverseFisherYatesShuffle<16384> shuffle_16384;
	int8_t genmat[255*71];
	mesg_type mesg[bits_max], mess[bits_max];
	code_type code[bits_max];
	cmplx cons[cons_max], prev[cols_max];
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
		CODE::MLS seq0(mls0_poly);
		for (int i = 0; i < symbol_len/2; ++i)
			fdom[i] = 0;
		for (int i = 0; i < mls0_len; ++i)
			fdom[(i+mls0_off/2+symbol_len/2)%(symbol_len/2)] = nrz(seq0());
		return fdom;
	}
	void systematic()
	{
		polarenc(mess, mesg, frozen_bits, code_order);
		int code_bits = 1 << code_order;
		for (int i = 0, j = 0; i < code_bits && j < crc_bits; ++i)
			if (!((frozen_bits[i/32] >> (i%32)) & 1))
				mesg[j++] = mess[i];
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
	void shuffle(code_type *c)
	{
		switch (code_order) {
		case 12:
			shuffle_4096(c);
			break;
		case 13:
			shuffle_8192(c);
			break;
		case 14:
			shuffle_16384(c);
			break;
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
	Decoder(uint8_t *out, int *len, DSP::ReadPCM<value> *pcm, int skip_count) :
		pcm(pcm), correlator(mls0_seq()), crc0(0xA8F4), crc1(0x8F6E37A0)
	{
		CODE::BoseChaudhuriHocquenghemGenerator<255, 71>::matrix(genmat, true, {
			0b100011101, 0b101110111, 0b111110011, 0b101101001,
			0b110111101, 0b111100111, 0b100101011, 0b111010111,
			0b000010011, 0b101100101, 0b110001011, 0b101100011,
			0b100011011, 0b100111111, 0b110001101, 0b100101101,
			0b101011111, 0b111111001, 0b111000011, 0b100111001,
			0b110101001, 0b000011111, 0b110000111, 0b110110001});

		blockdc.samples(filter_len);
		DSP::Phasor<cmplx> osc;
		const cmplx *buf;
		bool okay;
		do {
			okay = false;
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
				tdom[i] = buf[i+symbol_pos+extended_len] * osc();
			fwd(fdom, tdom);
			CODE::MLS seq1(mls1_poly);
			for (int i = 0; i < mls1_len; ++i)
				fdom[bin(i+mls1_off)] *= nrz(seq1());
			int8_t soft[mls1_len];
			uint8_t data[(mls1_len+7)/8];
			for (int i = 0; i < mls1_len; ++i)
				soft[i] = std::min<value>(std::max<value>(
					std::nearbyint(127 * demod_or_erase(
					fdom[bin(i+mls1_off)], fdom[bin(i-1+mls1_off)]).real()),
					-127), 127);
			bool unique = osddec(data, soft, genmat);
			if (!unique) {
				std::cerr << "OSD error." << std::endl;
				continue;
			}
			uint64_t md = 0;
			for (int i = 0; i < 55; ++i)
				md |= (uint64_t)CODE::get_be_bit(data, i) << i;
			uint16_t cs = 0;
			for (int i = 0; i < 16; ++i)
				cs |= (uint16_t)CODE::get_be_bit(data, i+55) << i;
			crc0.reset();
			if (crc0(md<<9) != cs) {
				std::cerr << "header CRC error." << std::endl;
				continue;
			}
			oper_mode = md & 255;
			if (oper_mode && (oper_mode < 23 || oper_mode > 30)) {
				std::cerr << "operation mode " << oper_mode << " unsupported." << std::endl;
				continue;
			}
			std::cerr << "oper mode: " << oper_mode << std::endl;
			if ((md>>8) == 0 || (md>>8) >= 129961739795077L) {
				std::cerr << "call sign unsupported." << std::endl;
				continue;
			}
			char call_sign[10];
			base37_decoder(call_sign, md>>8, 9);
			call_sign[9] = 0;
			std::cerr << "call sign: " << call_sign << std::endl;
			okay = true;
		} while (skip_count--);

		*len = 0;
		if (!okay || !oper_mode)
			return;

		int data_bits = 0;
		int cons_rows = 0;
		int comb_cols = 0;
		int code_cols = 0;
		switch (oper_mode) {
		case 23:
			mod_bits = 2;
			cons_rows = 8;
			comb_cols = 0;
			code_order = 12;
			code_cols = 256;
			data_bits = 2048;
			frozen_bits = frozen_4096_2080;
			break;
		case 24:
			mod_bits = 2;
			cons_rows = 16;
			comb_cols = 0;
			code_order = 13;
			code_cols = 256;
			data_bits = 4096;
			frozen_bits = frozen_8192_4128;
			break;
		case 25:
			mod_bits = 2;
			cons_rows = 32;
			comb_cols = 0;
			code_order = 14;
			code_cols = 256;
			data_bits = 8192;
			frozen_bits = frozen_16384_8224;
			break;
		case 26:
			mod_bits = 4;
			cons_rows = 4;
			comb_cols = 8;
			code_order = 12;
			code_cols = 256;
			data_bits = 2048;
			frozen_bits = frozen_4096_2080;
			break;
		case 27:
			mod_bits = 4;
			cons_rows = 8;
			comb_cols = 8;
			code_order = 13;
			code_cols = 256;
			data_bits = 4096;
			frozen_bits = frozen_8192_4128;
			break;
		case 28:
			mod_bits = 4;
			cons_rows = 16;
			comb_cols = 8;
			code_order = 14;
			code_cols = 256;
			data_bits = 8192;
			frozen_bits = frozen_16384_8224;
			break;
		case 29:
			mod_bits = 6;
			cons_rows = 5;
			comb_cols = 16;
			code_order = 13;
			code_cols = 273;
			data_bits = 4096;
			frozen_bits = frozen_8192_4128;
			break;
		case 30:
			mod_bits = 6;
			cons_rows = 10;
			comb_cols = 16;
			code_order = 14;
			code_cols = 273;
			data_bits = 8192;
			frozen_bits = frozen_16384_8224;
			break;
		default:
			return;
		}
		int cons_cols = code_cols + comb_cols;
		int comb_dist = comb_cols ? cons_cols / comb_cols : 1;
		int comb_off = comb_cols ? comb_dist / 2 : 1;
		int code_off = - cons_cols / 2;

		for (int i = 0; i < symbol_pos+extended_len; ++i)
			buf = next_sample();
		for (int i = 0; i < symbol_len; ++i)
			tdom[i] = buf[i] * osc();
		for (int i = 0; i < guard_len; ++i)
			osc();
		fwd(fdom, tdom);
		for (int i = 0; i < cons_cols; ++i)
			prev[i] = fdom[bin(i+code_off)];
		std::cerr << "demod ";
		CODE::MLS seq0(mls0_poly);
		for (int j = 0; j < cons_rows; ++j) {
			for (int i = 0; i < extended_len; ++i)
				buf = next_sample();
			for (int i = 0; i < symbol_len; ++i)
				tdom[i] = buf[i] * osc();
			for (int i = 0; i < guard_len; ++i)
				osc();
			fwd(fdom, tdom);
			for (int i = 0; i < cons_cols; ++i)
				cons[cons_cols*j+i] = demod_or_erase(fdom[bin(i+code_off)], prev[i]);
			if (oper_mode > 25) {
				for (int i = 0; i < comb_cols; ++i)
					cons[cons_cols*j+comb_dist*i+comb_off] *= nrz(seq0());
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
						prev[i] = fdom[bin(i+code_off)];
					else
						prev[i] *= DSP::polar<value>(1, tse(i+code_off));
				for (int i = 0; i < cons_cols; ++i) {
					index[i] = code_off + i;
					if (i % comb_dist == comb_off) {
						phase[i] = arg(cons[cons_cols*j+i]);
					} else {
						code_type tmp[mod_bits];
						mod_hard(tmp, cons[cons_cols*j+i]);
						phase[i] = arg(cons[cons_cols*j+i] * conj(mod_map(tmp)));
					}
				}
			}
			tse.compute(index, phase, cons_cols);
			//std::cerr << "Theil-Sen slope = " << tse.slope() << std::endl;
			//std::cerr << "Theil-Sen yint = " << tse.yint() << std::endl;
			for (int i = 0; i < cons_cols; ++i)
				cons[cons_cols*j+i] *= DSP::polar<value>(1, -tse(i+code_off));
			if (oper_mode > 25) {
				for (int i = 0; i < cons_cols; ++i)
					if (i % comb_dist != comb_off)
						prev[i] *= DSP::polar<value>(1, tse(i+code_off));
			} else {
				for (int i = 0; i < cons_cols; ++i)
					prev[i] = fdom[bin(i+code_off)];
			}
			std::cerr << ".";
		}
		std::cerr << " done" << std::endl;
		std::cerr << "Es/N0 (dB):";
		value sp = 0, np = 0;
		for (int j = 0, k = 0; j < cons_rows; ++j) {
			if (oper_mode > 25) {
				for (int i = 0; i < comb_cols; ++i) {
					cmplx hard(1, 0);
					cmplx error = cons[cons_cols*j+comb_dist*i+comb_off] - hard;
					sp += norm(hard);
					np += norm(error);
				}
			} else {
				for (int i = 0; i < cons_cols; ++i) {
					code_type tmp[mod_bits];
					mod_hard(tmp, cons[cons_cols*j+i]);
					cmplx hard = mod_map(tmp);
					cmplx error = cons[cons_cols*j+i] - hard;
					sp += norm(hard);
					np += norm(error);
				}
			}
			value precision = sp / np;
			// precision = 8;
			value snr = DSP::decibel(precision);
			std::cerr << " " << snr;
			for (int i = 0; i < cons_cols; ++i) {
				if (oper_mode > 25 && i % comb_dist == comb_off)
					continue;
				mod_soft(code+k, cons[cons_cols*j+i], precision);
				k += mod_bits;
			}
		}
		std::cerr << std::endl;
		*len = data_bits / 8;
		crc_bits = data_bits + 32;
		CODE::PolarHelper<mesg_type>::PATH metric[mesg_type::SIZE];
		for (int i = code_cols * cons_rows * mod_bits; i < bits_max; ++i)
			code[i] = 0;
		shuffle(code);
		polardec(metric, mesg, code, frozen_bits, code_order);
		systematic();
		int order[mesg_type::SIZE];
		for (int k = 0; k < mesg_type::SIZE; ++k)
			order[k] = k;
		std::sort(order, order+mesg_type::SIZE, [metric](int a, int b){ return metric[a] < metric[b]; });
		int best = -1;
		for (int k = 0; k < mesg_type::SIZE; ++k) {
			crc1.reset();
			for (int i = 0; i < crc_bits; ++i)
				crc1(mesg[i].v[order[k]] < 0);
			if (crc1() == 0) {
				best = order[k];
				break;
			}
		}
		if (best < 0) {
			std::cerr << "payload decoding error." << std::endl;
			*len = 0;
			return;
		}
		int flips = 0;
		for (int i = 0, j = 0; i < data_bits; ++i, ++j) {
			while ((frozen_bits[j / 32] >> (j % 32)) & 1)
				++j;
			bool received = code[j] < 0;
			bool decoded = mesg[i].v[best] < 0;
			flips += received != decoded;
			CODE::set_le_bit(out, i, decoded);
		}
		std::cerr << "bit flips: " << flips << std::endl;
	}
};

int main(int argc, char **argv)
{
	if (argc < 3 || argc > 4) {
		std::cerr << "usage: " << argv[0] << " OUTPUT INPUT [SKIP]" << std::endl;
		return 1;
	}

	typedef float value;
	typedef DSP::Complex<value> cmplx;

	const char *output_name = argv[1];
	if (output_name[0] == '-' && output_name[1] == 0)
		output_name = "/dev/stdout";
	const char *input_name = argv[2];
	if (input_name[0] == '-' && input_name[1] == 0)
		input_name = "/dev/stdin";

	DSP::ReadWAV<value> input_file(input_name);

	if (input_file.channels() < 1 || input_file.channels() > 2) {
		std::cerr << "Only real or analytic signal (one or two channels) supported." << std::endl;
		return 1;
	}

	int skip_count = 0;
	if (argc > 3)
		skip_count = std::atoi(argv[3]);

	const int data_max = 1024;
	uint8_t *output_data = new uint8_t[data_max];
	int data_len = 0;

	switch (input_file.rate()) {
	case 8000:
		delete new Decoder<value, cmplx, 8000>(output_data, &data_len, &input_file, skip_count);
		break;
	case 16000:
		delete new Decoder<value, cmplx, 16000>(output_data, &data_len, &input_file, skip_count);
		break;
	case 44100:
		delete new Decoder<value, cmplx, 44100>(output_data, &data_len, &input_file, skip_count);
		break;
	case 48000:
		delete new Decoder<value, cmplx, 48000>(output_data, &data_len, &input_file, skip_count);
		break;
	default:
		std::cerr << "Unsupported sample rate." << std::endl;
		return 1;
	}

	std::ofstream output_file(output_name, std::ios::binary | std::ios::trunc);
	if (output_file.bad()) {
		std::cerr << "Couldn't open file \"" << output_name << "\" for writing." << std::endl;
		return 1;
	}
	CODE::Xorshift32 scrambler;
	for (int i = 0; i < data_len; ++i)
		output_data[i] ^= scrambler();
	for (int i = 0; i < data_len; ++i)
		output_file.put(output_data[i]);
	delete []output_data;
	return 0;
}


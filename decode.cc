/*
OFDM modem decoder

Copyright 2023 Ahmet Inan <inan@aicodix.de>
*/

#include <algorithm>
#include <iostream>
#include <cassert>
#include <cmath>
namespace DSP { using std::abs; using std::min; using std::cos; using std::sin; }
#include "bip_buffer.hh"
#include "xorshift.hh"
#include "trigger.hh"
#include "complex.hh"
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
#include "psk.hh"
#include "hadamard_decoder.hh"
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
	DSP::SMA4<value, value, symbol_len, false> pwr;
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
	static cmplx demod_or_erase(cmplx curr, cmplx prev)
	{
		if (!(norm(prev) > 0))
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

	SchmidlCox(const cmplx *sequence) : threshold(value(0.2*match_len), value(0.3*match_len))
	{
		fwd(kern, sequence);
		for (int i = 0; i < symbol_len; ++i)
			kern[i] = conj(kern[i]) / value(symbol_len);
	}
	bool operator()(const cmplx *samples)
	{
		cmplx P = cor(samples[search_pos] * conj(samples[search_pos+symbol_len]));
		value R = value(0.5) * pwr(norm(samples[search_pos]) + norm(samples[search_pos+symbol_len]));
		value min_R = 0.00001 * symbol_len;
		R = std::max(R, min_R);
		value timing = match(norm(P) / (R * R));
		value phase = align(arg(P));

		bool collect = threshold(timing);
		bool process = falling(collect);

		// std::cout << timing << " " << process << " " << index_max << std::endl;
		// plot "< arecord -r 8000 -c 1 -f S16_LE | ./decode /dev/null - 1 | tee /tmp/data.txt" u ($1/33) w l, "/tmp/data.txt" u ($0-$3-1+16):2 w i, 0.2, 0.3

		if (!collect && !process)
			return false;

		if (timing_max < timing) {
			timing_max = timing;
			phase_max = phase;
			index_max = match_del;
		} else if (index_max < symbol_len + guard_len + match_del) {
			++index_max;
		} else if (process) {
			index_max = 0;
			timing_max = 0;
			return false;
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
			tmp1[i] = samples[i+symbol_pos] * osc();
		fwd(tmp0, tmp1);
		for (int i = 0; i < symbol_len; ++i)
			tmp1[i] = demod_or_erase(tmp0[i], tmp0[bin(i-1)]);
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

struct Decoder
{
	typedef float value;
	typedef DSP::Complex<value> cmplx;
	typedef int8_t code_type;
#ifdef __AVX2__
	typedef SIMD<code_type, 32 / sizeof(code_type)> mesg_type;
#else
	typedef SIMD<code_type, 16 / sizeof(code_type)> mesg_type;
#endif
	typedef DSP::Const<value> Const;
	static const int sample_rate = 8000;
	static const int code_order = 12;
	static const int mod_bits = 2;
	static const int code_len = 1 << code_order;
	static const int symbol_len = 256;
	static const int filter_len = 33;
	static const int guard_len = symbol_len / 8;
	static const int data_bits = 2048;
	static const int mesg_bits = data_bits + 32;
	static const int subcarrier_count = 64;
	static const int payload_symbols = 32;
	static const int first_subcarrier = 16;
	static const int cons_total = payload_symbols * subcarrier_count;
	static const int buffer_len = 2 * symbol_len + guard_len;
	static const int search_pos = symbol_len;
	DSP::ReadPCM<value> *pcm;
	DSP::FastFourierTransform<symbol_len, cmplx, -1> fwd;
	DSP::BlockDC<value, value> blockdc;
	DSP::Hilbert<cmplx, filter_len> hilbert;
	DSP::BipBuffer<cmplx, buffer_len> input_hist;
	SchmidlCox<value, cmplx, search_pos, symbol_len, guard_len> correlator;
	CODE::CRC<uint32_t> crc;
	CODE::HadamardDecoder<8> hadamard;
	CODE::PolarEncoder<mesg_type> polarenc;
	CODE::PolarListDecoder<mesg_type, code_order> polardec;
	mesg_type mesg[mesg_bits], mess[code_len];
	code_type code[code_len];
	cmplx cons[cons_total], prev[subcarrier_count];
	cmplx fdom[symbol_len], tdom[symbol_len];
	value cfo_rad, sfo_rad;
	const uint32_t *frozen_bits;
	int symbol_pos;

	static int nrz(bool bit)
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
	const cmplx *sync_seq()
	{
		CODE::MLS seq(0b1100111);
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] = 0;
		for (int i = first_subcarrier + 1; i < first_subcarrier + subcarrier_count; ++i)
			fdom[i] = nrz(seq());
		return fdom;
	}
	void systematic()
	{
		polarenc(mess, mesg, frozen_bits, code_order);
		int code_bits = 1 << code_order;
		for (int i = 0, j = 0; i < code_bits && j < mesg_bits; ++i)
			if (!((frozen_bits[i/32] >> (i%32)) & 1))
				mesg[j++] = mess[i];
	}
	cmplx mod_map(code_type *b)
	{
		return PhaseShiftKeying<4, cmplx, code_type>::map(b);
	}
	void mod_hard(code_type *b, cmplx c)
	{
		PhaseShiftKeying<4, cmplx, code_type>::hard(b, c);
	}
	void mod_soft(code_type *b, cmplx c, value precision)
	{
		PhaseShiftKeying<4, cmplx, code_type>::soft(b, c, precision);
	}
	const cmplx *next_sample()
	{
		value real;
		pcm->read(&real, 1);
		return input_hist(hilbert(blockdc(real)));
	}
	Decoder(uint8_t *out, int *len, DSP::ReadPCM<value> *pcm, int skip_count) :
		pcm(pcm), correlator(sync_seq()), crc(0x8F6E37A0)
	{
		frozen_bits = frozen_4096_2080;
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
			std::cerr << "coarse cfo: " << cfo_rad * (sample_rate / Const::TwoPi()) << " Hz " << std::endl;

			osc.omega(-cfo_rad);
			for (int i = 0; i < symbol_pos; ++i)
				buf = next_sample();
			for (int i = 0; i < symbol_len; ++i)
				tdom[i] = buf[i] * osc();
			for (int i = 0; i < guard_len; ++i)
				osc();
			fwd(fdom, tdom);
			for (int i = 0; i < subcarrier_count; ++i)
				prev[i] = fdom[first_subcarrier+i];
			for (int i = 0; i < symbol_len+guard_len; ++i)
				buf = next_sample();
			for (int i = 0; i < symbol_len; ++i)
				tdom[i] = buf[i] * osc();
			for (int i = 0; i < guard_len; ++i)
				osc();
			fwd(fdom, tdom);
			for (int i = 0; i < subcarrier_count; ++i)
				cons[i] = demod_or_erase(fdom[first_subcarrier+i], prev[i]);
			for (int i = 0; i < subcarrier_count; ++i)
				prev[i] = fdom[first_subcarrier+i];
			for (int i = 0; i < subcarrier_count; ++i)
				mod_soft(code+mod_bits*i, cons[i], 8);
			int oper_mode = hadamard(code);
			if (oper_mode != 1) {
				std::cerr << "operation mode " << oper_mode << " unsupported." << std::endl;
				continue;
			}
			std::cerr << "oper mode: " << oper_mode << std::endl;
			okay = true;
		} while (skip_count--);

		if (!okay)
			return;

		std::cerr << "demod ";
		for (int j = 0; j < payload_symbols; ++j) {
			for (int i = 0; i < symbol_len+guard_len; ++i)
				buf = next_sample();
			for (int i = 0; i < symbol_len; ++i)
				tdom[i] = buf[i] * osc();
			for (int i = 0; i < guard_len; ++i)
				osc();
			fwd(fdom, tdom);
			for (int i = 0; i < subcarrier_count; ++i)
				cons[subcarrier_count*j+i] = demod_or_erase(fdom[first_subcarrier+i], prev[i]);
			for (int i = 0; i < subcarrier_count; ++i)
				prev[i] = fdom[first_subcarrier+i];
			std::cerr << ".";
		}
		std::cerr << " done" << std::endl;
		if (1) {
			std::cerr << "Es/N0 (dB):";
			value sp = 0, np = 0;
			for (int j = 0; j < payload_symbols; ++j) {
				for (int i = 0; i < subcarrier_count; ++i) {
					code_type tmp[mod_bits];
					mod_hard(tmp, cons[subcarrier_count*j+i]);
					cmplx hard = mod_map(tmp);
					cmplx error = cons[subcarrier_count*j+i] - hard;
					sp += norm(hard);
					np += norm(error);
				}
				value precision = sp / np;
				value snr = DSP::decibel(precision);
				std::cerr << " " << snr;
				for (int i = 0; i < subcarrier_count; ++i)
					mod_soft(code+mod_bits*(subcarrier_count*j+i), cons[subcarrier_count*j+i], precision);
			}
			std::cerr << std::endl;
		} else {
			value precision = 8;
			for (int i = 0; i < cons_total; ++i)
				mod_soft(code+mod_bits*i, cons[i], precision);
		}
		CODE::PolarHelper<mesg_type>::PATH metric[mesg_type::SIZE];
		polardec(metric, mesg, code, frozen_bits, code_order);
		systematic();
		int order[mesg_type::SIZE];
		for (int k = 0; k < mesg_type::SIZE; ++k)
			order[k] = k;
		std::sort(order, order+mesg_type::SIZE, [metric](int a, int b){ return metric[a] < metric[b]; });
		int best = -1;
		for (int k = 0; k < mesg_type::SIZE; ++k) {
			crc.reset();
			for (int i = 0; i < mesg_bits; ++i)
				crc(mesg[i].v[order[k]] < 0);
			if (crc() == 0) {
				best = order[k];
				break;
			}
		}
		if (best < 0) {
			std::cerr << "payload decoding error." << std::endl;
			return;
		}
		*len = data_bits / 8;
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

	const char *output_name = argv[1];
	if (output_name[0] == '-' && output_name[1] == 0)
		output_name = "/dev/stdout";
	const char *input_name = argv[2];
	if (input_name[0] == '-' && input_name[1] == 0)
		input_name = "/dev/stdin";

	DSP::ReadWAV<float> input_file(input_name);

	if (input_file.channels() != 1) {
		std::cerr << "Only real signal (one channel) supported." << std::endl;
		return 1;
	}

	int skip_count = 0;
	if (argc > 3)
		skip_count = std::atoi(argv[3]);

	const int data_max = 2048 / 8;
	uint8_t *output_data = new uint8_t[data_max];
	int data_len = 0;

	if (input_file.rate() != 8000) {
		std::cerr << "Unsupported sample rate." << std::endl;
		return 1;
	}

	delete new Decoder(output_data, &data_len, &input_file, skip_count);

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


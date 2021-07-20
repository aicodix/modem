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
#include "resampler.hh"
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
#include "osd.hh"
#include "psk.hh"
#include "ldpc_tables.hh"
#include "ldpc_decoder.hh"
#include "galois_field.hh"
#include "bose_chaudhuri_hocquenghem_decoder.hh"

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
	DSP::Delay<value, match_del> delay;
	DSP::SchmittTrigger<value> threshold;
	DSP::FallingEdgeTrigger falling;
	cmplx tmp0[symbol_len], tmp1[symbol_len], tmp2[symbol_len];
	cmplx seq[symbol_len], kern[symbol_len];
	cmplx cmplx_shift = 0;
	value timing_max = 0;
	value phase_max = 0;
	int index_max = 0;

	static int bin(int carrier)
	{
		return (carrier + symbol_len) % symbol_len;
	}
public:
	int symbol_pos = 0;
	value cfo_rad = 0;
	value frac_cfo = 0;

	SchmidlCox(const cmplx *sequence) : threshold(value(0.17*match_len), value(0.19*match_len))
	{
		for (int i = 0; i < symbol_len; ++i)
			seq[i] = sequence[i];
		fwd(kern, sequence);
		for (int i = 0; i < symbol_len; ++i)
			kern[i] = conj(kern[i]) / value(symbol_len);
	}
	bool operator()(const cmplx *samples)
	{
		cmplx P = cor(samples[search_pos+symbol_len] * conj(samples[search_pos+2*symbol_len]));
		value R = value(0.5) * pwr(norm(samples[search_pos+2*symbol_len]));
		value min_R = 0.0001 * symbol_len;
		R = std::max(R, min_R);
		value timing = match(norm(P) / (R * R));
		value phase = delay(arg(P));

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
		for (int i = 0; i < symbol_len; ++i)
			tmp1[i] = 0;
		for (int i = 0; i < symbol_len; ++i)
			if (norm(tmp0[bin(i-1)]) > 0 &&
				std::min(norm(tmp0[i]), norm(tmp0[bin(i-1)])) * 2 >
				std::max(norm(tmp0[i]), norm(tmp0[bin(i-1)])))
					tmp1[i] = tmp0[i] / tmp0[bin(i-1)];
		fwd(tmp0, tmp1);
		for (int i = 0; i < symbol_len; ++i)
			tmp0[i] *= kern[i];
		bwd(tmp2, tmp0);

		int shift = 0;
		value peak = 0;
		value next = 0;
		for (int i = 0; i < symbol_len; ++i) {
			value power = norm(tmp2[i]);
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

		int pos_err = std::nearbyint(arg(tmp2[shift]) * symbol_len / Const::TwoPi());
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
	typedef DSP::Const<value> Const;
	static const int symbol_len = (1280 * rate) / 8000;
	static const int filter_len = (((21 * rate) / 8000) & ~3) | 1;
	static const int guard_len = symbol_len / 8;
	static const int ldpc_bits = 64800;
	static const int bch_bits = ldpc_bits - 21600;
	static const int data_bits = bch_bits - 12 * 16;
	static const int mod_min = 2;
	static const int mod_max = 3;
	static const int cons_max = ldpc_bits / mod_min;
	static const int cols_min = 360;
	static const int rows_max = cons_max / cols_min;
	static const int mls0_len = 127;
	static const int mls0_off = - mls0_len + 1;
	static const int mls0_poly = 0b10001001;
	static const int mls1_len = 255;
	static const int mls1_off = - mls1_len / 2;
	static const int mls1_poly = 0b100101011;
	static const int buffer_len = (rows_max + 8) * (symbol_len + guard_len);
	static const int search_pos = buffer_len - 4 * (symbol_len + guard_len);
	DSP::ReadPCM<value> *pcm;
	DSP::FastFourierTransform<symbol_len, cmplx, -1> fwd;
	DSP::FastFourierTransform<symbol_len, cmplx, 1> bwd;
	DSP::BlockDC<value, value> blockdc;
	DSP::Hilbert<cmplx, filter_len> hilbert;
	DSP::Resampler<value, filter_len, 3> resample;
	DSP::BipBuffer<cmplx, buffer_len> input_hist;
	SchmidlCox<value, cmplx, search_pos, symbol_len/2, guard_len> correlator;
	CODE::CRC<uint16_t> crc0;
	typedef CODE::GaloisField<16, 0b10000000000101101, uint16_t> GF;
	GF gf;
	CODE::BoseChaudhuriHocquenghemDecoder<24, 1, 65343, GF> bchdec1;
	CODE::OrderedStatisticsDecoder<255, 71, 4> osddec;
	CODE::LDPCDecoder<DVB_T2_TABLE_A3, 1> ldpcdec;
	int8_t genmat[255*71];
	int8_t code[ldpc_bits], bint[ldpc_bits];
	uint16_t erasures[24];
	cmplx head[symbol_len], tail[symbol_len], cons[cons_max];
	cmplx fdom[symbol_len], tdom[buffer_len], resam[buffer_len];
	value cfo_rad, sfo_rad;
	int symbol_pos;
	int oper_mode;
	int mod_bits;
	int cons_cnt;

	static int bin(int carrier)
	{
		return (carrier + symbol_len) % symbol_len;
	}
	static int nrz(bool bit)
	{
		return 1 - 2 * bit;
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
	int displacement(const cmplx *sym0, const cmplx *sym1)
	{
		fwd(head, sym0);
		fwd(tail, sym1);
		for (int i = 0; i < symbol_len; ++i)
			head[i] *= conj(tail[i]);
		bwd(tail, head);
		int idx = 0;
		for (int i = 0; i < symbol_len; ++i)
			if (norm(tail[i]) > norm(tail[idx]))
				idx = i;
		if (idx > symbol_len / 2)
			idx -= symbol_len;
		return -idx;
	}
	value frac_cfo(const cmplx *samples)
	{
		cmplx sum;
		for (int i = 0; i < symbol_len/2; ++i)
			sum += samples[i] * conj(samples[i+symbol_len/2]);
		return arg(sum) / (symbol_len/2);
	}
	void deinterleave()
	{
		for (int i = 0; i < cons_cnt; ++i)
			for (int k = 0; k < mod_bits; ++k)
				code[cons_cnt*k+i] = bint[mod_bits*i+k];
	}
	void interleave()
	{
		for (int i = 0; i < cons_cnt; ++i)
			for (int k = 0; k < mod_bits; ++k)
				bint[mod_bits*i+k] = code[cons_cnt*k+i];
	}
	cmplx mod_map(int8_t *b)
	{
		switch (oper_mode) {
		case 2:
		case 3:
			return PhaseShiftKeying<8, cmplx, int8_t>::map(b);
		case 4:
		case 5:
			return PhaseShiftKeying<4, cmplx, int8_t>::map(b);
		}
		return 0;
	}
	void mod_hard(int8_t *b, cmplx c)
	{
		switch (oper_mode) {
		case 2:
		case 3:
			PhaseShiftKeying<8, cmplx, int8_t>::hard(b, c);
			break;
		case 4:
		case 5:
			PhaseShiftKeying<4, cmplx, int8_t>::hard(b, c);
			break;
		}
	}
	void mod_soft(int8_t *b, cmplx c, value precision)
	{
		switch (oper_mode) {
		case 2:
		case 3:
			PhaseShiftKeying<8, cmplx, int8_t>::soft(b, c, precision);
			break;
		case 4:
		case 5:
			PhaseShiftKeying<4, cmplx, int8_t>::soft(b, c, precision);
			break;
		}
	}
	Decoder(uint8_t *out, DSP::ReadPCM<value> *pcm, int skip_count) :
		pcm(pcm), resample(rate, (rate * 19) / 40, 2), correlator(mls0_seq()), crc0(0xA8F4)
	{
		CODE::BoseChaudhuriHocquenghemGenerator<255, 71>::matrix(genmat, true, {
			0b100011101, 0b101110111, 0b111110011, 0b101101001,
			0b110111101, 0b111100111, 0b100101011, 0b111010111,
			0b000010011, 0b101100101, 0b110001011, 0b101100011,
			0b100011011, 0b100111111, 0b110001101, 0b100101101,
			0b101011111, 0b111111001, 0b111000011, 0b100111001,
			0b110101001, 0b000011111, 0b110000111, 0b110110001});

		bool real = pcm->channels() == 1;
		blockdc.samples(2*(symbol_len+guard_len));
		const cmplx *buf;
		do {
			do {
				if (!pcm->good())
					return;
				cmplx tmp;
				pcm->read(reinterpret_cast<value *>(&tmp), 1);
				if (real)
					tmp = hilbert(blockdc(tmp.real()));
				buf = input_hist(tmp);
			} while (!correlator(buf));
		} while (skip_count--);

		symbol_pos = correlator.symbol_pos;
		cfo_rad = correlator.cfo_rad;
		std::cerr << "symbol pos: " << symbol_pos << std::endl;
		std::cerr << "coarse cfo: " << cfo_rad * (rate / Const::TwoPi()) << " Hz " << std::endl;

		DSP::Phasor<cmplx> osc;
		osc.omega(-cfo_rad);
		for (int i = 0; i < symbol_len; ++i)
			tdom[i] = buf[i+symbol_pos+(symbol_len+guard_len)] * osc();
		fwd(fdom, tdom);
		CODE::MLS seq1(mls1_poly);
		for (int i = 0; i < mls1_len; ++i)
			fdom[bin(i+mls1_off)] *= nrz(seq1());
		int8_t soft[mls1_len];
		uint8_t data[(mls1_len+7)/8];
		for (int i = 0; i < mls1_len; ++i)
			soft[i] = std::min<value>(std::max<value>(
				std::nearbyint(127 * (fdom[bin(i+mls1_off)] /
				fdom[bin(i-1+mls1_off)]).real()), -128), 127);
		bool unique = osddec(data, soft, genmat);
		if (!unique) {
			std::cerr << "OSD error." << std::endl;
			return;
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
			return;
		}
		oper_mode = md & 255;
		int code_cols;
		switch (oper_mode) {
		case 2:
			code_cols = 432;
			mod_bits = 3;
			break;
		case 3:
			code_cols = 400;
			mod_bits = 3;
			break;
		case 4:
			code_cols = 400;
			mod_bits = 2;
			break;
		case 5:
			code_cols = 360;
			mod_bits = 2;
			break;
		default:
			std::cerr << "operation mode " << oper_mode << " unsupported." << std::endl;
			return;
		}
		cons_cnt = ldpc_bits / mod_bits;
		std::cerr << "oper mode: " << oper_mode << std::endl;
		if ((md>>8) == 0 || (md>>8) >= 129961739795077L) {
			std::cerr << "call sign unsupported." << std::endl;
			return;
		}
		char call_sign[10];
		base37_decoder(call_sign, md>>8, 9);
		call_sign[9] = 0;
		std::cerr << "call sign: " << call_sign << std::endl;

		int code_rows = cons_cnt / code_cols;
		int code_off = - code_cols / 2;

		int dis = displacement(buf+symbol_pos-(code_rows+1)*(symbol_len+guard_len), buf+symbol_pos+2*(symbol_len+guard_len));
		sfo_rad = (dis * Const::TwoPi()) / ((code_rows+3)*(symbol_len+guard_len));
		std::cerr << "coarse sfo: " << 1000000 * sfo_rad / Const::TwoPi() << " ppm" << std::endl;
		if (dis) {
			value diff = sfo_rad * (rate / Const::TwoPi());
			resample(resam, buf, -diff, buffer_len);
			symbol_pos = std::nearbyint(correlator.symbol_pos * (1 - sfo_rad / Const::TwoPi()));
			std::cerr << "resam pos: " << symbol_pos << std::endl;
		} else {
			for (int i = 0; i < buffer_len; ++i)
				resam[i] = buf[i];
		}
		cfo_rad = correlator.cfo_rad + correlator.frac_cfo - frac_cfo(resam+symbol_pos);
		std::cerr << "finer cfo: " << cfo_rad * (rate / Const::TwoPi()) << " Hz " << std::endl;

		osc.omega(-cfo_rad);
		for (int i = 0; i < buffer_len; ++i)
			tdom[i] = resam[i] * osc();

		cmplx *cur = tdom + symbol_pos - (code_rows + 1) * (symbol_len + guard_len);
		fwd(fdom, cur);
		for (int j = 0; j < code_rows; ++j) {
			for (int i = 0; i < code_cols; ++i)
				head[bin(i+code_off)] = fdom[bin(i+code_off)];
			fwd(fdom, cur += symbol_len+guard_len);
			for (int i = 0; i < code_cols; ++i)
				cons[code_cols*j+i] = fdom[bin(i+code_off)] / head[bin(i+code_off)];
		}
		if (1) {
			value sum = 0;
			for (int i = 0; i < cons_cnt; ++i) {
				int8_t tmp[mod_max];
				mod_hard(tmp, cons[i]);
				sum += arg(cons[i] * conj(mod_map(tmp)));
			}
			value avg = sum / cons_cnt;
			cfo_rad += avg / (symbol_len+guard_len);
			std::cerr << "finer cfo: " << cfo_rad * (rate / Const::TwoPi()) << " Hz " << std::endl;
			cmplx comp = DSP::polar<value>(1, -avg);
			for (int i = 0; i < cons_cnt; ++i)
				cons[i] *= comp;
		}
		value precision = 16;
		if (1) {
			value sp = 0, np = 0;
			for (int i = 0; i < cons_cnt; ++i) {
				int8_t tmp[mod_max];
				mod_hard(tmp, cons[i]);
				cmplx hard = mod_map(tmp);
				cmplx error = cons[i] - hard;
				sp += norm(hard);
				np += norm(error);
			}
			value snr = DSP::decibel(sp / np);
			std::cerr << "init Es/N0: " << snr << " dB" << std::endl;
			// $LLR=log(\frac{p(x=+1|y)}{p(x=-1|y)})$
			// $p(x|\mu,\sigma)=\frac{1}{\sqrt{2\pi}\sigma}}e^{-\frac{(x-\mu)^2}{2\sigma^2}}$
			value sigma = std::sqrt(np / (2 * sp));
			precision = 1 / (sigma * sigma);
		}
		for (int i = 0; i < cons_cnt; ++i)
			mod_soft(bint+mod_bits*i, cons[i], precision);
		deinterleave();
		int count = ldpcdec(code, code + bch_bits);
		if (count < 0)
			std::cerr << "payload LDPC decoding did not converge." << std::endl;
		if (1) {
			interleave();
			value sp = 0, np = 0;
			for (int i = 0; i < cons_cnt; ++i) {
				int8_t tmp[mod_max];
				for (int k = 0; k < mod_bits; ++k)
					tmp[k] = nrz(bint[mod_bits*i+k] < 0);
				cmplx hard = mod_map(tmp);
				cmplx error = cons[i] - hard;
				sp += norm(hard);
				np += norm(error);
			}
			value snr = DSP::decibel(sp / np);
			std::cerr << "corr Es/N0: " << snr << " dB" << std::endl;
			// $LLR=log(\frac{p(x=+1|y)}{p(x=-1|y)})$
			// $p(x|\mu,\sigma)=\frac{1}{\sqrt{2\pi}\sigma}}e^{-\frac{(x-\mu)^2}{2\sigma^2}}$
			value sigma = std::sqrt(np / (2 * sp));
			precision = 1 / (sigma * sigma);
		}
		for (int i = 0; i < bch_bits; ++i)
			CODE::set_le_bit(out, i, code[i] < 0);
		int ecnt = 0;
		for (int i = 0; i < bch_bits; ++i) {
			if (!code[i]) {
				if (ecnt < 24) {
					erasures[ecnt++] = i;
				} else {
					std::cerr << "payload LDPC produced more than 24 erasures." << std::endl;
					return;
				}
			}
		}
		if (ecnt)
			std::cerr << "payload LDPC produced " << ecnt << " erasures." << std::endl;
		int ret = bchdec1(out, out+data_bits/8, erasures, ecnt, data_bits);
		if (ret < 0) {
			std::cerr << "payload BCH error." << std::endl;
			return;
		}
		if (ret)
			std::cerr << "payload BCH corrected " << ret << " errors." << std::endl;
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
	const char *input_name = argv[2];

	DSP::ReadWAV<value> input_file(input_name);

	if (input_file.channels() < 1 || input_file.channels() > 2) {
		std::cerr << "Only real or analytic signal (one or two channels) supported." << std::endl;
		return 1;
	}

	int skip_count = 1;
	if (argc > 3)
		skip_count = std::atoi(argv[3]);

	const int code_len = 64800 / 8;
	uint8_t *output_data = new uint8_t[code_len];

	switch (input_file.rate()) {
	case 8000:
		delete new Decoder<value, cmplx, 8000>(output_data, &input_file, skip_count);
		break;
	case 16000:
		delete new Decoder<value, cmplx, 16000>(output_data, &input_file, skip_count);
		break;
	case 44100:
		delete new Decoder<value, cmplx, 44100>(output_data, &input_file, skip_count);
		break;
	case 48000:
		delete new Decoder<value, cmplx, 48000>(output_data, &input_file, skip_count);
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
	const int data_len = code_len - (12 * 16 + 21600) / 8;
	CODE::Xorshift32 scrambler;
	for (int i = 0; i < data_len; ++i)
		output_data[i] ^= scrambler();
	for (int i = 0; i < data_len; ++i)
		output_file.put(output_data[i]);
	delete []output_data;
	return 0;
}


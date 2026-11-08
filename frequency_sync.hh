/*
Frequency domain correlator for the two identical sync symbols

Copyright 2026 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

#include "fft.hh"

// Takes an FFT every symbol_len samples. As the two sync symbols are back to back,
// one of these windows always falls between them and only sees a cyclic shift of the
// sync symbol: the ratio of adjacent tones is then the known one, turned by the same
// phase, which comes from the shift. As in SchmidlCox, adjacent tones are divided,
// not multiplied, so a strong tone counts no more than any other. Only the tones of
// the kernel take part, so unlike SchmidlCox, which divides all bins of the symbol,
// there are no empty bins to erase by their power, only pairs whose ratio is far off.
// Their sum, against the number of pairs, peaks there. It knows of no integer CFO:
// with FM or AM there is little CFO besides that from the SFO, a small fraction of a
// tone.
template<typename value, typename cmplx, int search_pos, int symbol_len, int guard_len>
class FrequencySync {
	typedef DSP::Const<value> Const;
	DSP::FastFourierTransform<symbol_len, cmplx, -1> fwd;
	cmplx tmp0[symbol_len];
	// the conjugate of the known ratio of adjacent tones, only where it isn't zero,
	// with the bins of the tone and the one before it
	cmplx kern[symbol_len];
	int curr_bin[symbol_len], prev_bin[symbol_len];
	int pairs = 0;
	value threshold;
	value metric_max = 0;
	cmplx sum_max = 0;
	int hop = 0;

	// keeps a single tone next to a faded one from dominating the sum
	static cmplx demod_or_erase(cmplx curr, cmplx prev) {
		if (norm(prev) > value(0)) {
			cmplx cons = curr / prev;
			if (norm(cons) < value(4))
				return cons;
		}
		return 0;
	}

	// correlation of the symbol at pos with the one after it
	static cmplx repetition(const cmplx *samples, int pos) {
		cmplx P = 0;
		for (int i = 0; i < symbol_len; ++i)
			P += samples[pos + i] * conj(samples[pos + symbol_len + i]);
		return P;
	}

public:
	int symbol_pos = 0;
	value cfo_rad = 0;
	value metric = 0;

	// takes the known ratio of adjacent tones, as SchmidlCox, which could also be complex, as with a CAZAC sequence
	FrequencySync(const cmplx *sequence, value threshold = value(0.2)) : threshold(threshold) {
		// the leftmost tone has none before it, so the ratio is zero there
		for (int i = 0; i < symbol_len; ++i) {
			if (norm(sequence[i]) > value(0)) {
				kern[pairs] = conj(sequence[i]);
				curr_bin[pairs] = i;
				prev_bin[pairs++] = (i + symbol_len - 1) % symbol_len;
			}
		}
	}

	bool operator()(const cmplx *samples) {
		if (++hop < symbol_len)
			return false;
		hop = 0;

		fwd(tmp0, samples + search_pos);
		cmplx sum = 0;
		for (int i = 0; i < pairs; ++i)
			sum += kern[i] * demod_or_erase(tmp0[curr_bin[i]], tmp0[prev_bin[i]]);
		metric = abs(sum) / value(pairs);

		// keep the best window until the next one is worse
		if (metric > threshold && metric > metric_max) {
			metric_max = metric;
			sum_max = sum;
			return false;
		}
		if (metric_max == value(0))
			return false;
		metric_max = 0;

		// that window started shift samples after the start of a sync symbol, modulo symbol_len
		int shift = std::nearbyint(arg(sum_max) * symbol_len / Const::TwoPi());
		shift = (shift + symbol_len) % symbol_len;
		// and a hop ago. With the ratio of adjacent tones, a window that reaches a few guard
		// intervals beyond the sync symbols still looks clean, so it started within half a
		// symbol of them: only the first sync symbol is followed by the second
		int first = search_pos - 5 * symbol_len / 2, last = search_pos - symbol_len / 2;
		int pos = 0;
		cmplx P = 0;
		for (int cand = first + (search_pos - shift - first) % symbol_len; cand <= last; cand += symbol_len) {
			cmplx Q = repetition(samples, cand);
			if (norm(Q) > norm(P)) {
				pos = cand;
				P = Q;
			}
		}
		if (norm(P) <= value(0))
			return false;
		symbol_pos = pos;
		cfo_rad = - arg(P) / value(symbol_len);
		return true;
	}
};

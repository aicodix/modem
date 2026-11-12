/*
Sync detection on the raw input in the frequency domain

Copyright 2026 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

#include "fft.hh"

// Like FrequencySync, but on the real input before mixing and filtering, so nothing
// else has to run while waiting for a sync. Every symbol_len samples it takes a Hann
// window over the last two symbols and a real FFT twice as long as a symbol. Within the
// two identical sync symbols the content repeats every symbol, so the tones only fall on
// even bins, and the window, which only spreads a bin to its neighbours, leaves them
// apart while keeping strong out-of-band tones, hum and DC from leaking in. Mixing is a
// shift of the bins, and searching the shifts by whole bins, half a tone, finds the CFO
// within a quarter of a tone, up to range bins away. As in FrequencySync, adjacent tones
// are divided, not multiplied, and the phase of the sum gives the timing modulo a symbol.
// refine() then picks the start among the candidates with the correlation of consecutive
// symbols on the baseband, and the fraction of a tone of the CFO.
template<typename value, typename cmplx, int symbol_len, int carrier, int range>
class RawSync {
	typedef DSP::Const<value> Const;
	static const int fft_len = 2 * symbol_len;
	DSP::RealToHalfComplexTransform<fft_len, cmplx> fwd;
	value win[fft_len], tmp[fft_len];
	cmplx spec[symbol_len + 1];
	cmplx ratio[symbol_len + 1];
	// the conjugate of the known ratio of adjacent tones and the bin of the tone without CFO
	cmplx kern[symbol_len];
	int tone_bin[symbol_len];
	int pairs = 0, low = symbol_len, high = 0;
	value threshold;
	value metric_max = 0;
	cmplx sum_max = 0;
	int cfo_max = 0;
	int hop = 0;

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
	// the best window started shift samples after the start of a sync symbol, modulo symbol_len
	int shift = 0;
	// the CFO in bins of half a tone
	int cfo_bins = 0;

	// takes the known ratio of adjacent tones on the baseband, as FrequencySync, with carrier the tone at the center
	RawSync(const cmplx *sequence, value threshold = value(0.2)) : threshold(threshold) {
		for (int i = 0; i < fft_len; ++i)
			win[i] = value(0.5) - value(0.5) * DSP::UnitCircle<value>::cos(i, fft_len);
		for (int i = 0; i < symbol_len; ++i) {
			if (norm(sequence[i]) > value(0)) {
				int tone = i < symbol_len / 2 ? i : i - symbol_len;
				kern[pairs] = conj(sequence[i]);
				tone_bin[pairs] = 2 * (tone + carrier);
				low = std::min(low, tone_bin[pairs]);
				high = std::max(high, tone_bin[pairs]);
				++pairs;
			}
		}
		low -= range;
		high += range;
	}

	// takes the last fft_len samples of the raw input
	bool operator()(const value *samples) {
		if (++hop < symbol_len)
			return false;
		hop = 0;

		for (int i = 0; i < fft_len; ++i)
			tmp[i] = win[i] * samples[i];
		fwd(spec, tmp);
		// the tone before is a whole tone, two bins, lower
		for (int i = low; i <= high; ++i)
			ratio[i] = demod_or_erase(spec[i], spec[i - 2]);
		cmplx best = 0;
		int cfo = 0;
		for (int c = -range; c <= range; ++c) {
			cmplx sum = 0;
			for (int i = 0; i < pairs; ++i)
				sum += kern[i] * ratio[tone_bin[i] + c];
			if (norm(sum) > norm(best)) {
				best = sum;
				cfo = c;
			}
		}
		metric = abs(best) / value(pairs);

		// keep the best window until the next one is worse
		if (metric > threshold && metric > metric_max) {
			metric_max = metric;
			sum_max = best;
			cfo_max = cfo;
			return false;
		}
		if (metric_max == value(0))
			return false;
		metric_max = 0;
		shift = std::nearbyint(arg(sum_max) * symbol_len / Const::TwoPi());
		shift = (shift + symbol_len) % symbol_len;
		cfo_bins = cfo_max;
		return true;
	}

	// takes the baseband, in which the best window started at start, and picks the start of the first sync
	// symbol within a symbol before or after it: silence before the sync looks as clean as the guard interval
	// after it, so two windows a hop apart can both look clean, the best one early or late by half a symbol
	bool refine(const cmplx *samples, int start) {
		int first = start - symbol_len, last = start + symbol_len;
		int pos = 0;
		cmplx P = 0;
		for (int cand = first + ((start - shift - first) % symbol_len + symbol_len) % symbol_len; cand <= last; cand += symbol_len) {
			cmplx Q = repetition(samples, cand);
			if (norm(Q) > norm(P)) {
				pos = cand;
				P = Q;
			}
		}
		if (norm(P) <= value(0))
			return false;
		symbol_pos = pos;
		// the correlation only knows the CFO modulo a tone, the search knows it within a quarter of a tone
		value tone = Const::TwoPi() / value(symbol_len);
		value fine = - arg(P) / value(symbol_len);
		value coarse = cfo_bins * tone / value(2);
		cfo_rad = fine + tone * std::nearbyint((coarse - fine) / tone);
		return true;
	}
};

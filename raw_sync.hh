/*
Sync detection on the raw input in the frequency domain

Copyright 2026 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

#include "fft.hh"
#include "mls.hh"

// Sync detection on the real input before mixing and filtering, so nothing else has to
// run while waiting for a sync. Every symbol_len samples it takes a Hann window over the
// last two symbols and a real FFT twice as long as a symbol. Within the two identical
// sync symbols the content repeats every symbol, so the tones only fall on even bins, and
// the window, which only spreads a bin to its neighbours, leaves them apart. Mixing is a
// shift of the bins, and searching the shifts by whole bins, half a tone, finds the CFO
// within a quarter of a tone, up to range bins away. Adjacent tones are divided, not
// multiplied, so a strong tone counts no more than any other, which keeps hum, DC and
// other tones out, and the phase of the sum gives the timing modulo a symbol.
// refine() then picks the start among the candidates with the correlation of consecutive
// symbols on the baseband, and the fraction of a tone of the CFO.
// The tone_count tones of the sync symbols are BPSK from the MLS of poly and seed, centered
// on carrier, the tone at the center, so the known ratio of adjacent tones is +1 or -1.
template<typename value, typename cmplx, int symbol_len, int carrier, int range, int tone_count, int poly, int seed>
class RawSync {
	typedef DSP::Const<value> Const;
	static const int fft_len = 2 * symbol_len;
	// each tone and the one before it, a whole tone or two bins lower
	static const int pairs = tone_count - 1;
	// the bin of the first tone that has one before it, without CFO
	static const int first = 2 * (carrier - tone_count / 2 + 1);
	// the bins whose ratio the search over the shifts reads
	static const int low = first - range, high = first + 2 * (pairs - 1) + range;
	static_assert(low >= 2 && high <= symbol_len, "the shifts reach beyond the spectrum");
	DSP::RealToHalfComplexTransform<fft_len, cmplx> fwd;
	value win[fft_len], tmp[fft_len];
	cmplx spec[symbol_len + 1];
	// the real and imaginary parts of the ratios of the even and of the odd bins, so a shift reads every second
	// bin as contiguous arrays of reals and the search is two plain dot products with the real kernel
	value ratio_re[2][symbol_len / 2 + 1], ratio_im[2][symbol_len / 2 + 1];
	// the known ratio of adjacent tones, real as the tones are BPSK
	value kern[pairs];
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

	RawSync(value threshold = value(0.2)) : threshold(threshold) {
		// The FFT treats the window as one period of a periodic signal. Only what repeats
		// exactly within it, as the two sync symbols in a window that lines up with them,
		// stays on its bins. A window off by up to half a symbol, a CFO between the bins and
		// SFO, which stretches the second symbol against the first and moves every tone off
		// its bin by a different amount, all leave a step where the end wraps around to the
		// start. Without a taper, that step leaks every tone into its neighbours, the very
		// tones the ratio compares. Hann goes to zero at both ends, so the step is gone and
		// what leaks stays on the empty odd bins next to each tone. Measured on a clean frame,
		// the metric stays between 0.97 and 1.05 for any alignment and fraction of a bin, and
		// above 0.96 up to 2000 ppm, where a rectangular window drops to 0.62 and to 0.71.
		// Decoded frames didn't change, but Hann syncs at about 1 dB more noise and gives
		// preamble() a third less timing error to start from.
		for (int i = 0; i < fft_len; ++i)
			win[i] = value(0.5) - value(0.5) * DSP::UnitCircle<value>::cos(i, fft_len);
		CODE::MLS seq(poly, seed);
		value prev = 1 - 2 * seq();
		for (int i = 0; i < pairs; ++i) {
			value curr = 1 - 2 * seq();
			kern[i] = prev * curr;
			prev = curr;
		}
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
		for (int i = low; i <= high; ++i) {
			cmplx r = demod_or_erase(spec[i], spec[i - 2]);
			ratio_re[i & 1][i >> 1] = r.real();
			ratio_im[i & 1][i >> 1] = r.imag();
		}
		cmplx best = 0;
		int cfo = 0;
		for (int c = -range; c <= range; ++c) {
			const value *re = ratio_re[(first + c) & 1] + ((first + c) >> 1);
			const value *im = ratio_im[(first + c) & 1] + ((first + c) >> 1);
			value sum_re = 0, sum_im = 0;
			for (int i = 0; i < pairs; ++i) {
				sum_re += kern[i] * re[i];
				sum_im += kern[i] * im[i];
			}
			cmplx sum(sum_re, sum_im);
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

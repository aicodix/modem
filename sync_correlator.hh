/*
Frequency domain correlator for repeated sync symbols

Copyright 2026 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

#include "fft.hh"

// the sync symbol is sent several times back to back after one guard interval,
// so a window of one symbol length inside them sees it cyclically shifted
template<typename value, typename cmplx, int symbol_len, int guard_len, int buffer_len, int tone_count, int tone_off>
class SyncCorrelator {
	static const int extended_len = symbol_len + guard_len;
	// leaves room behind the window for one more copy, a guard interval and the metadata symbol
	static const int search_pos = buffer_len - symbol_len - extended_len;
	DSP::FastFourierTransform<symbol_len, cmplx, -1> fwd;
	DSP::FastFourierTransform<symbol_len, cmplx, 1> bwd;
	cmplx tmp0[symbol_len], tmp1[symbol_len];
	value seq[tone_count];
	int hop = 0, aligned = -1;
	static int bin(int carrier) {
		return (carrier + symbol_len) % symbol_len;
	}
	// correlation with the sync symbol at the shift with the largest value,
	// relative to the energies of both, so a clean copy gives one and noise about 1 / tone_count
	value search(const cmplx *samples, int *shift) {
		fwd(tmp0, samples + search_pos);
		for (int i = 0; i < symbol_len; ++i)
			tmp1[i] = 0;
		value pwr = 0;
		for (int i = 0; i < tone_count; ++i) {
			cmplx v = tmp0[bin(i+tone_off)];
			tmp1[bin(i+tone_off)] = seq[i] * v;
			pwr += norm(v);
		}
		bwd(tmp0, tmp1);
		value peak = 0;
		for (int i = 0; i < symbol_len; ++i) {
			value power = norm(tmp0[i]);
			if (power > peak) {
				peak = power;
				*shift = i;
			}
		}
		return peak / std::max(value(tone_count) * pwr, value(1e-30));
	}
	// the same at a known position, without shift
	value match(const cmplx *samples, int pos) {
		fwd(tmp0, samples + pos);
		cmplx sum = 0;
		value pwr = 0;
		for (int i = 0; i < tone_count; ++i) {
			cmplx v = tmp0[bin(i+tone_off)];
			sum += seq[i] * v;
			pwr += norm(v);
		}
		return norm(sum) / std::max(value(tone_count) * pwr, value(1e-30));
	}
public:
	static constexpr value threshold = 0.05;
	// start of the second to last copy, the last copy follows and then the metadata symbol after a guard interval
	int symbol_pos = 0;
	value metric = 0;
	SyncCorrelator(const value *sequence) {
		for (int i = 0; i < tone_count; ++i)
			seq[i] = sequence[i];
	}
	// call once for every sample, returns true when the metadata symbol follows the copies
	bool operator()(const cmplx *samples) {
		if (++hop < symbol_len)
			return false;
		hop = 0;
		if (aligned < 0) {
			int shift = 0;
			value m = search(samples, &shift);
			if (m < threshold)
				return false;
			// the window started this many samples after the start of a copy
			int offset = (symbol_len - shift) % symbol_len;
			aligned = search_pos - offset;
			metric = m;
		}
		// one hop later, the same position holds the next copy
		if (match(samples, aligned) < threshold) {
			aligned = -1;
			return false;
		}
		if (match(samples, aligned + symbol_len) >= threshold)
			return false;
		symbol_pos = aligned - symbol_len;
		aligned = -1;
		return true;
	}
};

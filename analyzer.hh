/*
Spectrum analyzer and X-Y mode oscilloscope of the baseband signal

Copyright 2026 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

#include "common.hh"
#include "window.hh"
#include "coeffs.hh"

template <typename value, typename cmplx, int rate>
struct Analyzer
{
	// as long as the symbols of the decoder, so the bins are as wide as the tones
	static const int symbol_len = (rate / 300) * 40;
	static const int hop = symbol_len / 2;
	// the middle bins, from 0 to 3000 Hz
	static const int bins = 400;
	static const int width = bins;
	static const int scope_height = bins;
	static const int spectrum_height = bins / 4;
	static const int spectrogram_height = bins / 2;
	static const int subchannels_height = bins / 2;
	DSP::Hann<value> hann;
	DSP::Coeffs<symbol_len, value, true> win;
	DSP::FastFourierTransform<symbol_len, cmplx, -1> fwd;
	// the newest hop samples also go to the scope
	cmplx hist[symbol_len], tmp[symbol_len], fdom[symbol_len];
	value phosphor[width * scope_height];
	value dB[bins], dB_avg[bins], avg_max, power;
	int scope[width * scope_height];
	int spectrum[width * spectrum_height];
	int spectrogram[width * spectrogram_height];
	int subchannels[width * subchannels_height];
	int pos;
	bool fresh, spectrum_on, scope_on, subchannels_on, show_rainbow, spectrum_phosphor, subchannels_rainbow;
	value intensity, decay;

	Analyzer() : win(&hann), avg_max(0), power(0), pos(0), fresh(true),
		spectrum_on(false), scope_on(false), subchannels_on(false),
		show_rainbow(true), spectrum_phosphor(true), subchannels_rainbow(true), intensity(0.5), decay(0.5)
	{
		for (int i = 0; i < symbol_len; ++i)
			hist[i] = 0;
		for (int i = 0; i < width * scope_height; ++i)
			phosphor[i] = scope[i] = 0;
		for (int i = 0; i < width * spectrum_height; ++i)
			spectrum[i] = 0;
		for (int i = 0; i < width * spectrogram_height; ++i)
			spectrogram[i] = 0;
		for (int i = 0; i < width * subchannels_height; ++i)
			subchannels[i] = 0;
	}
	static int rgba(value r, value g, value b, value a)
	{
		auto byte = [](value v) { return int(nearbyint(255 * std::sqrt(std::min(std::max(v, value(0)), value(1))))); };
		int A = nearbyint(255 * std::min(std::max(a, value(0)), value(1)));
		return (A << 24) | (byte(b) << 16) | (byte(g) << 8) | byte(r);
	}
	static int alpha(value a)
	{
		return rgba(1, 1, 1, a);
	}
	static int rainbow(value v)
	{
		v = std::min(std::max(v, value(0)), value(1));
		return rgba(4 * v - 2, 1 - 4 * std::abs(v - value(0.5)), 2 - 4 * v, 4 * v);
	}
	void line(int x0, int y0, int x1, int y1, value a)
	{
		int dx = std::abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
		int dy = -std::abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
		int err = dx + dy;
		a /= std::max(dx, -dy) + 1;
		while (true) {
			if (0 <= x0 && x0 < width && 0 <= y0 && y0 < scope_height)
				phosphor[width * y0 + x0] += a;
			if (x0 == x1 && y0 == y1)
				break;
			int e2 = 2 * err;
			if (e2 >= dy) {
				err += dy;
				x0 += sx;
			}
			if (e2 <= dx) {
				err += dx;
				y0 += sy;
			}
		}
	}
	// vertical span from y0 to y1, the first half in column x0 and the second in x0+1
	void span(int x0, int y0, int y1, int color)
	{
		int mid = (y0 + y1) / 2;
		for (int y = std::min(y0, mid); y <= std::max(y0, mid); ++y)
			spectrum[width * y + x0] = color;
		for (int y = std::min(mid, y1); y <= std::max(mid, y1); ++y)
			spectrum[width * y + x0 + 1] = color;
	}
	void draw_scope()
	{
		const cmplx *trace = hist + symbol_len - hop;
		value sum = 0;
		for (int i = 0; i < hop; ++i)
			sum += norm(trace[i]);
		power = DSP::lerp(power, sum / hop, value(0.1));
		value half = (scope_height - 1) / value(2);
		// three times the RMS reaches the edge
		value scale = half / std::max(3 * std::sqrt(power), value(0.000001));
		// a decay of 100% clears the scope with each update
		for (int i = 0; i < width * scope_height; ++i)
			phosphor[i] *= decay;
		// the brightness of each update, spread over a quarter of the image
		value a = intensity * (width * scope_height / 4) / hop;
		for (int i = 0, x0 = 0, y0 = 0; i < hop; ++i) {
			int x1 = nearbyint(half + scale * trace[i].real());
			int y1 = nearbyint(half - scale * trace[i].imag());
			if (i)
				line(x0, y0, x1, y1, a);
			x0 = x1;
			y0 = y1;
		}
		for (int i = 0; i < width * scope_height; ++i)
			scope[i] = alpha(phosphor[i]);
	}
	void draw_spectrum()
	{
		for (int i = 0; i < symbol_len; ++i)
			tmp[i] = win[i] * hist[i];
		fwd(fdom, tmp);
		// center the spectrum at 1500 Hz
		for (int i = 0; i < bins; ++i)
			dB[i] = DSP::decibel(norm(fdom[(i - bins / 2 + symbol_len) % symbol_len]));
		value tmp_max = dB[0];
		for (int i = 1; i < bins; ++i)
			tmp_max = std::max(tmp_max, dB[i]);
		if (fresh) {
			for (int i = 0; i < bins; ++i)
				dB_avg[i] = dB[i];
			avg_max = tmp_max;
			fresh = false;
		}
		for (int i = 0; i < bins; ++i)
			dB_avg[i] = DSP::lerp(dB_avg[i], dB[i], value(0.05));
		avg_max = DSP::lerp(avg_max, tmp_max, avg_max < tmp_max ? value(0.5) : value(0.05));
		value tmp_min = dB_avg[0];
		for (int i = 1; i < bins; ++i)
			tmp_min = std::min(tmp_min, dB_avg[i]);
		value dB_max = avg_max;
		value dB_min = std::min(std::max(tmp_min, value(-120)), dB_max - 10);
		value range = 1 / (dB_max - dB_min);
		if (spectrum_phosphor) {
			for (int i = 0; i < width * spectrum_height; ++i)
				spectrum[i] = 0xffffff | ((((unsigned)spectrum[i] >> 24) * 7) >> 3) << 24;
		} else {
			for (int i = 0; i < width * spectrum_height; ++i)
				spectrum[i] = 0xffffff;
		}
		for (int b = 0, y0 = 0; b < bins; ++b) {
			int y1 = nearbyint((spectrum_height - 1) * range * (dB_max - dB[b]));
			y1 = std::min(std::max(y1, 0), spectrum_height - 1);
			if (b)
				span(b - 1, y0, y1, -1);
			y0 = y1;
		}
		for (int i = width * spectrogram_height - 1; i >= width; --i)
			spectrogram[i] = spectrogram[i - width];
		for (int b = 0; b < bins; ++b) {
			value v = 1 - range * (dB_max - dB[b]);
			spectrogram[b] = show_rainbow ? rainbow(v) : alpha(v);
		}
	}
	// takes the FFT of a symbol from the decoder after synchronization and draws one row
	void symbol(const cmplx *symbol_fdom)
	{
		if (!subchannels_on)
			return;
		value row[bins];
		for (int i = 0; i < bins; ++i)
			row[i] = DSP::decibel(norm(symbol_fdom[(i - bins / 2 + symbol_len) % symbol_len]));
		value row_max = row[0], row_min = row[0];
		for (int i = 1; i < bins; ++i) {
			row_max = std::max(row_max, row[i]);
			row_min = std::min(row_min, row[i]);
		}
		row_min = std::min(std::max(row_min, row_max - 120), row_max - 10);
		value range = 1 / (row_max - row_min);
		for (int i = width * subchannels_height - 1; i >= width; --i)
			subchannels[i] = subchannels[i - width];
		for (int b = 0; b < bins; ++b) {
			value v = 1 - range * (row_max - row[b]);
			subchannels[b] = subchannels_rainbow ? rainbow(v) : alpha(v);
		}
	}
	// takes one sample of the baseband, does nothing while turned off
	void operator()(cmplx sample)
	{
		if (!spectrum_on && !scope_on)
			return;
		hist[symbol_len - hop + pos++] = sample;
		if (pos < hop)
			return;
		pos = 0;
		if (spectrum_on)
			draw_spectrum();
		if (scope_on)
			draw_scope();
		for (int i = 0; i < symbol_len - hop; ++i)
			hist[i] = hist[i + hop];
	}
};


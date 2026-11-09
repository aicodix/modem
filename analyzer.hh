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
	typedef DSP::Const<value> Const;
	// as long as the symbols of the decoder, so the bins are as wide as the tones
	static const int symbol_len = (rate / 300) * 40;
	static const int hop = symbol_len / 2;
	// the middle bins, from 0 to 3000 Hz
	static const int bins = 400;
	static const int width = bins;
	static const int scope_height = bins;
	static const int spectrum_height = bins / 4;
	static const int spectrogram_height = bins / 2;
	DSP::Hann<value> hann;
	DSP::Coeffs<symbol_len, value, true> win;
	DSP::FastFourierTransform<symbol_len, cmplx, -1> fwd;
	// the newest hop samples also go to the scope
	cmplx hist[symbol_len], tmp[symbol_len], fdom[symbol_len];
	value phosphor[width * scope_height];
	value dB[bins];
	// the hold of each bin, jumps to new extremes and slowly pulls back to the current value,
	// a fraction of the difference per update, at 15 updates per second about a second
	static constexpr value hold_pull = value(1) / 16;
	value dB_low[bins], dB_high[bins];
	bool hold_valid;
	// fixed scales without gain control: the scope gain, and the spectrum in dB of a full scale sine
	value scope_gain, spectrum_max, spectrum_min;
	int scope[width * scope_height];
	int spectrum[width * spectrum_height];
	int spectrogram[width * spectrogram_height];
	// the channel estimate of the decoder, only for the tones,
	// its power in dB of a full scale sine and its phase from -180 to 180 degrees
	static const int channel_width = Common::tone_count;
	static constexpr value channel_top = 0, channel_range = 100;
	int channel_power[channel_width * spectrum_height];
	int channel_phase[channel_width * spectrum_height];
	// the real Es/N0 of each tone after a decoded frame, from 0 to 50 dB
	static constexpr value snr_top = 50, snr_range = 50;
	int channel_snr[channel_width * spectrum_height];
	// the precision of each tone after each digital symbol, on the same scale as the real Es/N0
	int channel_precision[channel_width * spectrum_height];
	int pos;
	// the equalized tones of the latest symbol, on a fixed scale
	// few but large pixels
	static const int constellation_size = 128;
	static constexpr value constellation_scale = 2;
	int constellation[constellation_size * constellation_size];
	bool spectrum_on, scope_on, channel_on, constellation_on, show_rainbow;
	value intensity, decay;

	Analyzer() : win(&hann), hold_valid(false), scope_gain(1), spectrum_max(0), spectrum_min(-100), pos(0),
		spectrum_on(false), scope_on(false), channel_on(false), constellation_on(false),
		show_rainbow(true), intensity(0.5), decay(0.5)
	{
		for (int i = 0; i < symbol_len; ++i)
			hist[i] = 0;
		for (int i = 0; i < width * scope_height; ++i)
			phosphor[i] = scope[i] = 0;
		for (int i = 0; i < width * spectrum_height; ++i)
			spectrum[i] = 0;
		for (int i = 0; i < width * spectrogram_height; ++i)
			spectrogram[i] = 0;
		for (int i = 0; i < channel_width * spectrum_height; ++i)
			channel_power[i] = channel_phase[i] = channel_snr[i] = channel_precision[i] = 0;
		for (int i = 0; i < constellation_size * constellation_size; ++i)
			constellation[i] = 0;
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
	template <int W>
	static void span(int *image, int x0, int y0, int y1, int color)
	{
		int mid = (y0 + y1) / 2;
		for (int y = std::min(y0, mid); y <= std::max(y0, mid); ++y)
			image[W * y + x0] = color;
		for (int y = std::min(mid, y1); y <= std::max(mid, y1); ++y)
			image[W * y + x0 + 1] = color;
	}
	template <int W>
	static void clear(int *image)
	{
		for (int i = 0; i < W * spectrum_height; ++i)
			image[i] = 0xffffff;
	}
	// draws a line through the W dB values
	template <int W>
	static void curve(int *image, const value *dB, value dB_max, value range, int color)
	{
		for (int b = 0, y0 = 0; b < W; ++b) {
			int y1 = nearbyint((spectrum_height - 1) * range * (dB_max - dB[b]));
			y1 = std::min(std::max(y1, 0), spectrum_height - 1);
			if (b)
				span<W>(image, b - 1, y0, y1, color);
			y0 = y1;
		}
	}
	// a new picture with only the line through the W dB values
	template <int W>
	static void plot(int *image, const value *dB, value dB_max, value range)
	{
		clear<W>(image);
		curve<W>(image, dB, dB_max, range, -1);
	}
	void draw_scope()
	{
		const cmplx *trace = hist + symbol_len - hop;
		value half = (scope_height - 1) / value(2);
		// without gain, an amplitude of one reaches the edge, about three times the RMS of a full scale signal
		value scale = half * scope_gain;
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
		// center the spectrum at 1500 Hz, the window is normalized and
		// a full scale sine gives a tone of half its amplitude after mixing
		for (int i = 0; i < bins; ++i)
			dB[i] = DSP::decibel(std::max(4 * norm(fdom[(i - bins / 2 + symbol_len) % symbol_len]), value(1e-30)));
		for (int i = 0; i < bins; ++i) {
			if (!hold_valid)
				dB_low[i] = dB_high[i] = dB[i];
			dB_low[i] = std::min(dB[i], dB_low[i] + hold_pull * (dB[i] - dB_low[i]));
			dB_high[i] = std::max(dB[i], dB_high[i] + hold_pull * (dB[i] - dB_high[i]));
		}
		hold_valid = true;
		value dB_max = spectrum_max;
		value dB_min = std::min(spectrum_min, dB_max - 10);
		value range = 1 / (dB_max - dB_min);
		// min and max dimmed, the current spectrum on top
		clear<width>(spectrum);
		curve<width>(spectrum, dB_low, dB_max, range, alpha(value(0.5)));
		curve<width>(spectrum, dB_high, dB_max, range, alpha(value(0.5)));
		curve<width>(spectrum, dB, dB_max, range, -1);
		for (int i = width * spectrogram_height - 1; i >= width; --i)
			spectrogram[i] = spectrogram[i - width];
		for (int b = 0; b < bins; ++b) {
			value v = 1 - range * (dB_max - dB[b]);
			spectrogram[b] = show_rainbow ? rainbow(v) : alpha(v);
		}
	}
	// takes the channel estimate of the decoder after each symbol, on a fixed scale without any gain control
	void channel_estimate(const cmplx *chan)
	{
		if (!channel_on)
			return;
		value est[channel_width];
		// a full scale sine gives a tone of half its amplitude after mixing, times symbol_len after the FFT
		const value full_scale = value(symbol_len * symbol_len) / 4;
		for (int i = 0; i < channel_width; ++i)
			est[i] = DSP::decibel(norm(chan[i]) / full_scale);
		plot<channel_width>(channel_power, est, channel_top, 1 / channel_range);
		// one dot per tone, as a line would jump where the phase wraps around
		clear<channel_width>(channel_phase);
		for (int i = 0; i < channel_width; ++i) {
			int y = nearbyint((spectrum_height - 1) * (value(0.5) - arg(chan[i]) / Const::TwoPi()));
			y = std::min(std::max(y, 0), spectrum_height - 1);
			channel_phase[channel_width * y + i] = -1;
		}
	}
	// takes the real Es/N0 in dB of each tone from the decoder after a decoded frame
	void tone_snr(const value *snr)
	{
		if (channel_on)
			plot<channel_width>(channel_snr, snr, snr_top, 1 / snr_range);
	}
	// takes the precision of each tone from the decoder after symbol j, nothing for erased or analog symbols
	template <typename DECODER>
	void tone_precision(DECODER &dec, int j)
	{
		if (!channel_on)
			return;
		if (dec.lost[j] || dec.analog_mode) {
			clear<channel_width>(channel_precision);
			return;
		}
		value est[channel_width];
		for (int i = 0; i < channel_width; ++i)
			est[i] = DSP::decibel(std::max(dec.tone_precision[i], value(1e-30)));
		plot<channel_width>(channel_precision, est, snr_top, 1 / snr_range);
	}
	void point(cmplx c)
	{
		value half = (constellation_size - 1) / value(2);
		int x = nearbyint(half + half / constellation_scale * c.real());
		int y = nearbyint(half - half / constellation_scale * c.imag());
		if (x < 0 || x >= constellation_size || y < 0 || y >= constellation_size)
			return;
		constellation[constellation_size * y + x] = -1;
	}
	// takes the equalized symbol j from the decoder and draws a new picture of it,
	// the BPSK tones of the preamble symbol, of the payload symbols only the tones with all bits
	template <typename DECODER>
	void constellation_symbol(DECODER &dec, int j, int start)
	{
		if (!constellation_on)
			return;
		for (int i = 0; i < constellation_size * constellation_size; ++i)
			constellation[i] = 0;
		if (dec.lost[j])
			return;
		for (int i = 0, l = start; i < Common::tone_count; ++i) {
			if (i % Common::block_length == dec.seed_off)
				continue;
			if (!j || dec.analog_mode) {
				point(dec.demod[i]);
			} else {
				int bits = dec.bits_at(l);
				l += bits;
				if (bits == dec.mod_bits)
					point(dec.demod[i]);
			}
		}
	}
	// takes one sample of the baseband, does nothing while turned off
	void operator()(cmplx sample)
	{
		// the hold starts over when the spectrum is turned on again
		if (!spectrum_on)
			hold_valid = false;
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


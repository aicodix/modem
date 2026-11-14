/*
OFDM modem decoder WebAssembly module

Copyright 2026 Ahmet Inan <inan@aicodix.de>
*/

#include "decoder.hh"
#include "analyzer.hh"

typedef float value;
typedef DSP::Complex<value> cmplx;
typedef Decoder<value, cmplx, RATE> decoder_type;
typedef DSP::Const<value> Const;
static decoder_type decoder;
typedef Analyzer<value, cmplx, RATE> analyzer_type;
static_assert(analyzer_type::symbol_len == decoder_type::symbol_len, "analyzer and decoder must agree on the symbol length");
static analyzer_type analyzer;
static int last_events;

extern "C" {

__attribute__((visibility("default")))
float samples[decoder_type::extended_len];
__attribute__((visibility("default")))
unsigned char payload[Common::data_max];
__attribute__((visibility("default")))
char call_sign[10];
/* cfo Hz, sfo ppm, residual cfo mHz, Es/N0 dB min, median, max, real Es/N0 dB min, median, max */
__attribute__((visibility("default")))
float stats[9];
/* the CFO in Hz and SFO in ppm the decoder tracks, after each symbol */
__attribute__((visibility("default")))
float tracking[2];

static void update(int events)
{
	if (events & decoder_type::found_sync) {
		stats[0] = decoder.cfo_rad * (RATE / Const::TwoPi());
		stats[1] = -1000000 * decoder.sfo_rad / Const::TwoPi();
		stats[2] = 1000 * decoder.residual_rad * RATE / (Const::TwoPi() * decoder_type::symbol_len);
	}
	if (events & decoder_type::read_meta) {
		int64_t val = decoder.call_sign;
		for (int i = 8; i >= 0; --i, val /= 40)
			call_sign[i] = "   /0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"[val%40];
		call_sign[9] = 0;
	}
	if (events & decoder_type::frame_end)
		for (int i = 0; i < 3; ++i)
			stats[3+i] = decoder.es_n0[i];
	if (events & decoder_type::real_snr)
		for (int i = 0; i < 3; ++i)
			stats[6+i] = decoder.real_es_n0[i];
	if (events & decoder_type::decoded_frame)
		for (int i = 0; i < decoder.data_bytes; ++i)
			payload[i] = decoder.data[i];
}

/* feeds samples from offset on, stops after the first sample with events and returns the number of samples used */
__attribute__((visibility("default")))
int feed(int offset, int count)
{
	if (offset < 0 || count < 0 || offset + count > decoder_type::extended_len)
		return -1;
	last_events = 0;
	// the views of the baseband need it also while waiting for a sync
	decoder.keep_baseband = analyzer.spectrum_on || analyzer.scope_on;
	for (int i = 0; i < count; ++i) {
		/* the index of the next symbol, unless it follows a sync */
		int symbol = decoder.symbol_index;
		int events = decoder.feed(samples[offset+i]);
		analyzer(decoder.input_hist()[decoder_type::buffer_len-1]);
		/* only for the analyzer, the page doesn't need to stop for it */
		if (events & decoder_type::new_symbol) {
			tracking[0] = decoder.cfo_rad * (RATE / Const::TwoPi());
			tracking[1] = 1000000 * (decoder.step - 1);
			int j = events & decoder_type::found_sync ? 0 : symbol;
			analyzer.channel_estimate(decoder.chan);
			analyzer.tone_precision(decoder, j);
			analyzer.constellation_symbol(decoder, j);
		}
		if (events & decoder_type::real_snr)
			analyzer.tone_snr(decoder.real_tone_snr);
		events &= ~decoder_type::new_symbol;
		if (events) {
			last_events = events;
			update(events);
			return i + 1;
		}
	}
	return count;
}

__attribute__((visibility("default")))
int events()
{
	return last_events;
}

__attribute__((visibility("default")))
int busy()
{
	return decoder.busy();
}

__attribute__((visibility("default")))
int oper_mode()
{
	return decoder.oper_mode;
}

__attribute__((visibility("default")))
int meta_mode()
{
	return decoder.meta_mode;
}

__attribute__((visibility("default")))
int erased_symbol()
{
	return decoder.erased_index;
}

__attribute__((visibility("default")))
int bit_flips()
{
	return decoder.bit_flips;
}

__attribute__((visibility("default")))
int code_bits()
{
	return decoder.frame_bits;
}

/* the payload symbols, the bits per tone and k of the code rate k / 16 of the latest known mode */
__attribute__((visibility("default")))
int symbol_count()
{
	return decoder.symbol_count;
}

__attribute__((visibility("default")))
int mod_bits()
{
	return decoder.mod_bits;
}

__attribute__((visibility("default")))
int code_rate()
{
	if (!decoder.codewords)
		return 0;
	return decoder.data_bits / decoder.codewords >> (decoder.code_order - 4);
}

__attribute__((visibility("default")))
int payload_length()
{
	return decoder.data_bytes;
}

/* RGBA pixels of the analyzer images, all analyzer_width() wide */
__attribute__((visibility("default")))
int *scope()
{
	return analyzer.scope;
}

__attribute__((visibility("default")))
int *spectrum()
{
	return analyzer.spectrum;
}

__attribute__((visibility("default")))
int *spectrogram()
{
	return analyzer.spectrogram;
}

/* the decoder's channel estimate after each symbol, channel_width() wide and spectrum_height() high */
__attribute__((visibility("default")))
int *channel_power()
{
	return analyzer.channel_power;
}

__attribute__((visibility("default")))
int *channel_phase()
{
	return analyzer.channel_phase;
}

/* the real Es/N0 of each tone after a decoded frame, channel_width() wide and spectrum_height() high */
__attribute__((visibility("default")))
int *channel_snr()
{
	return analyzer.channel_snr;
}

/* the precision of each tone after each digital symbol, channel_width() wide and spectrum_height() high */
__attribute__((visibility("default")))
int *channel_precision()
{
	return analyzer.channel_precision;
}

__attribute__((visibility("default")))
int channel_width()
{
	return analyzer_type::channel_width;
}

__attribute__((visibility("default")))
int enable_channel(int on)
{
	if (on < 0 || on > 1)
		return -1;
	analyzer.channel_on = on;
	return 0;
}

/* the equalized symbols of the latest frame, constellation_size() wide and high */
__attribute__((visibility("default")))
int *constellation()
{
	return analyzer.constellation;
}

__attribute__((visibility("default")))
int constellation_size()
{
	return analyzer_type::constellation_size;
}

__attribute__((visibility("default")))
int enable_constellation(int on)
{
	if (on < 0 || on > 1)
		return -1;
	analyzer.constellation_on = on;
	return 0;
}

__attribute__((visibility("default")))
int analyzer_width()
{
	return analyzer_type::width;
}

__attribute__((visibility("default")))
int scope_height()
{
	return analyzer_type::scope_height;
}

__attribute__((visibility("default")))
int spectrum_height()
{
	return analyzer_type::spectrum_height;
}

__attribute__((visibility("default")))
int spectrogram_height()
{
	return analyzer_type::spectrogram_height;
}

/* spectrum and spectrogram */
__attribute__((visibility("default")))
int enable_spectrum(int on)
{
	if (on < 0 || on > 1)
		return -1;
	analyzer.spectrum_on = on;
	return 0;
}

__attribute__((visibility("default")))
int enable_scope(int on)
{
	if (on < 0 || on > 1)
		return -1;
	analyzer.scope_on = on;
	return 0;
}

__attribute__((visibility("default")))
int show_rainbow(int on)
{
	if (on < 0 || on > 1)
		return -1;
	analyzer.show_rainbow = on;
	return 0;
}

/* in percent */
__attribute__((visibility("default")))
int scope_intensity(int percent)
{
	if (percent < 0 || percent > 100)
		return -1;
	analyzer.intensity = value(0.01) * percent;
	return 0;
}

/* in dB, at 0 dB an amplitude of one reaches the edge */
__attribute__((visibility("default")))
int scope_gain(int dB)
{
	if (dB < -20 || dB > 60)
		return -1;
	analyzer.scope_gain = std::sqrt(DSP::idecibel(value(dB)));
	return 0;
}

/* top and bottom of spectrum and spectrogram in dB of a full scale sine, at least 10 dB apart */
__attribute__((visibility("default")))
int spectrum_range(int max, int min)
{
	if (max < -100 || max > 0 || min < -160 || min > max - 10)
		return -1;
	analyzer.spectrum_max = max;
	analyzer.spectrum_min = min;
	return 0;
}

/* in percent per update */
__attribute__((visibility("default")))
int phosphor_decay(int percent)
{
	if (percent < 0 || percent > 100)
		return -1;
	analyzer.decay = 1 - value(0.01) * percent;
	return 0;
}

}

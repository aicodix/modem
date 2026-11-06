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
static analyzer_type analyzer;
static int last_events;

extern "C" {

__attribute__((visibility("default")))
float samples[decoder_type::extended_len];
__attribute__((visibility("default")))
unsigned char payload[Common::data_max];
__attribute__((visibility("default")))
char call_sign[10];
__attribute__((visibility("default")))
float values[2*Common::analog_max];
/* cfo Hz, sfo ppm, residual cfo mHz, Es/N0 dB min, median, max, real Es/N0 dB min, median, max */
__attribute__((visibility("default")))
float stats[9];

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
	if ((events & decoder_type::decoded_frame) && !decoder.analog_mode)
		for (int i = 0; i < decoder.data_bytes; ++i)
			payload[i] = decoder.data[i];
	if ((events & decoder_type::decoded_frame) && decoder.analog_mode)
		for (int i = 0; i < decoder.symbol_count * decoder.data_tones; ++i) {
			values[2*i] = decoder.values[i].real();
			values[2*i+1] = decoder.values[i].imag();
		}
}

/* feeds samples from offset on, stops after the first sample with events and returns the number of samples used */
__attribute__((visibility("default")))
int feed(int offset, int count)
{
	if (offset < 0 || count < 0 || offset + count > decoder_type::extended_len)
		return -1;
	last_events = 0;
	for (int i = 0; i < count; ++i) {
		int events = decoder.feed(samples[offset+i]);
		analyzer(decoder.input_hist()[decoder_type::buffer_len-1]);
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
	return 1 << decoder.code_order;
}

__attribute__((visibility("default")))
int values_length()
{
	if (!decoder.analog_mode)
		return -1;
	return decoder.symbol_count * decoder.data_tones;
}

__attribute__((visibility("default")))
int analog_rate()
{
	return Common::analog_rate;
}

__attribute__((visibility("default")))
int payload_length()
{
	if (decoder.analog_mode)
		return -1;
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

__attribute__((visibility("default")))
int show_phosphor(int on)
{
	if (on < 0 || on > 1)
		return -1;
	analyzer.show_phosphor = on;
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

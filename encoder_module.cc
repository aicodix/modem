/*
OFDM modem encoder WebAssembly module

Copyright 2026 Ahmet Inan <inan@aicodix.de>
*/

#include "encoder.hh"

typedef float value;
typedef DSP::Complex<value> cmplx;
typedef Encoder<value, cmplx, RATE> encoder_type;
static encoder_type encoder;
static bool analog;

extern "C" {

__attribute__((visibility("default")))
unsigned char payload[Common::data_max];
__attribute__((visibility("default")))
float audio[encoder_type::audio_max];
__attribute__((visibility("default")))
float values[2*Common::analog_max];

static int copy_audio()
{
	for (int i = 0; i < encoder.audio_len; ++i)
		audio[i] = encoder.audio[i];
	return encoder.audio_len;
}

__attribute__((visibility("default")))
int start(int oper_mode)
{
	if (oper_mode < 0 || oper_mode > 128)
		return -1;
	analog = oper_mode == 128;
	if (!encoder.start(oper_mode))
		return -1;
	return copy_audio();
}

__attribute__((visibility("default")))
int payload_length()
{
	return encoder.data_bytes;
}

__attribute__((visibility("default")))
int symbol_count()
{
	if (analog)
		return -1;
	return encoder.symbol_count;
}

__attribute__((visibility("default")))
int begin(int length, long long call_sign)
{
	if (analog)
		return -1;
	if (length < 0 || length > encoder.data_bytes)
		return -1;
	if (call_sign <= 0 || call_sign >= 262144000000000LL)
		return -1;
	if (!encoder.begin(payload, length, call_sign))
		return -1;
	return 0;
}

__attribute__((visibility("default")))
int values_length()
{
	return Common::analog_max;
}

__attribute__((visibility("default")))
int begin_analog(int count, long long call_sign)
{
	if (!analog)
		return -1;
	if (count < 1 || count > Common::analog_max)
		return -1;
	if (call_sign <= 0 || call_sign >= 262144000000000LL)
		return -1;
	if (!encoder.begin(reinterpret_cast<cmplx *>(values), count, call_sign))
		return -1;
	return 0;
}

__attribute__((visibility("default")))
int next()
{
	encoder.next();
	return copy_audio();
}

__attribute__((visibility("default")))
int stop()
{
	encoder.stop();
	return copy_audio();
}

}


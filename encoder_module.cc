/*
OFDM modem encoder WebAssembly module

Copyright 2026 Ahmet Inan <inan@aicodix.de>
*/

#include "encoder.hh"

typedef float value;
typedef DSP::Complex<value> cmplx;
typedef Encoder<value, cmplx, RATE> encoder_type;
static encoder_type encoder;

extern "C" {

__attribute__((visibility("default")))
unsigned char payload[Common::data_max];
__attribute__((visibility("default")))
float audio[encoder_type::audio_max];

static int copy_audio()
{
	for (int i = 0; i < encoder.audio_len; ++i)
		audio[i] = encoder.audio[i];
	return encoder.audio_len;
}

__attribute__((visibility("default")))
int start(int oper_mode)
{
	if (oper_mode < 0 || oper_mode > 127)
		return -1;
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
int begin(int length, long long call_sign)
{
	if (length < 0 || length > encoder.data_bytes)
		return -1;
	if (call_sign <= 0 || call_sign >= 262144000000000LL)
		return -1;
	if (!encoder.begin(payload, length, call_sign))
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


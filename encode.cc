/*
OFDM modem encoder

Copyright 2021 Ahmet Inan <inan@aicodix.de>
*/

#include "encoder.hh"
#include "native.hh"

int read_digital(uint8_t *payload, int length, const char *input_name)
{
	std::ifstream input_file(input_name, std::ios::binary);
	if (input_file.bad()) {
		std::cerr << "Couldn't open file \"" << input_name << "\" for reading." << std::endl;
		return -1;
	}
	int count = 0;
	for (int c; count < length && (c = input_file.get()) >= 0; ++count)
		payload[count] = c;
	return count;
}

template <typename value, typename cmplx, int rate>
void encode(DSP::WritePCM<value> *pcm, const char *const *input_names, int input_count, int64_t call_sign, int seconds)
{
	auto *encoder = new Encoder<value, cmplx, rate>();
	uint8_t payload[Common::data_max];
	bool started = false;
	for (int input_index = 0; input_index < input_count; ++input_index) {
		const char *input_name = input_names[input_index];
		if (input_count == 1 && input_name[0] == '-' && input_name[1] == 0)
			input_name = "/dev/stdin";
		int length = read_digital(payload, Common::data_max, input_name);
		if (length < 0)
			continue;
		// the most robust mode that carries all of the input in about the given seconds
		int oper_mode = encoder->pick_mode(seconds, length);
		if (oper_mode < 0) {
			std::cerr << "Input \"" << input_name << "\" too long for about " << seconds << " seconds." << std::endl;
			continue;
		}
		if (!started) {
			started = encoder->start(oper_mode);
			pcm->write(encoder->audio, encoder->audio_len);
		}
		print_setup(*encoder, oper_mode);
		if (!encoder->begin(payload, length, call_sign))
			continue;
		while (encoder->next())
			pcm->write(encoder->audio, encoder->audio_len);
		int last = encoder->symbol_count;
		DSP::quick_sort(encoder->papr, last + 1);
		std::cerr << "PAPR (dB): " << DSP::decibel(encoder->papr[0]) << " .. " << DSP::decibel(encoder->papr[last/2]) << " .. " << DSP::decibel(encoder->papr[last]) << std::endl;
	}
	if (started) {
		encoder->stop();
		pcm->write(encoder->audio, encoder->audio_len);
	}
	delete encoder;
}

template <typename value, typename cmplx, int rate>
void ping(DSP::WritePCM<value> *pcm, int64_t call_sign, int ping_value)
{
	auto *encoder = new Encoder<value, cmplx, rate>();
	// a frame without payload symbols carries the value in its mode
	if (encoder->start(ping_value)) {
		pcm->write(encoder->audio, encoder->audio_len);
		print_setup(*encoder, ping_value);
		encoder->begin(nullptr, 0, call_sign);
		while (encoder->next())
			pcm->write(encoder->audio, encoder->audio_len);
		encoder->stop();
		pcm->write(encoder->audio, encoder->audio_len);
	}
	delete encoder;
}

int64_t base40_encoder(const char *str)
{
	int64_t acc = 0;
	for (char c = *str++; c; c = *str++) {
		acc *= 40;
		if (c == '/')
			acc += 3;
		else if (c >= '0' && c <= '9')
			acc += c - '0' + 4;
		else if (c >= 'a' && c <= 'z')
			acc += c - 'a' + 14;
		else if (c >= 'A' && c <= 'Z')
			acc += c - 'A' + 14;
		else if (c != ' ')
			return -1;
	}
	return acc;
}

int main(int argc, char **argv)
{
	if (argc < 7) {
		std::cerr << "usage: " << argv[0] << " OUTPUT RATE BITS CALLSIGN SECONDS INPUT.." << std::endl;
		std::cerr << "       " << argv[0] << " OUTPUT RATE BITS CALLSIGN 0 PING" << std::endl;
		return 1;
	}

	const char *output_name = argv[1];
	if (output_name[0] == '-' && output_name[1] == 0)
		output_name = "/dev/stdout";
	int output_rate = std::atoi(argv[2]);
	int output_bits = std::atoi(argv[3]);
	int64_t call_sign = base40_encoder(argv[4]);
	if (call_sign <= 0 || call_sign >= 262144000000000L) {
		std::cerr << "Unsupported call sign." << std::endl;
		return 1;
	}
	// how long a frame may take, the encoder picks the closest number of symbols, or a ping without any
	int seconds = std::atoi(argv[5]);
	if (seconds < 0 || seconds > 60) {
		std::cerr << "Unsupported number of seconds." << std::endl;
		return 1;
	}
	int ping_value = std::atoi(argv[6]);
	if (!seconds && (argc != 7 || ping_value < 0 || ping_value > 31)) {
		std::cerr << "A ping takes a single value from 0 to 31." << std::endl;
		return 1;
	}
	int input_index = 6;
	std::cerr << std::fixed << std::setprecision(1);
	typedef float value;
	typedef DSP::Complex<value> cmplx;
	DSP::WriteWAV<value> output_file(output_name, output_rate, output_bits, 1);
	output_file.silence(output_rate);
	int input_count = argc - input_index;
	switch (output_rate) {
	case 44100:
		if (seconds)
			encode<value, cmplx, 44100>(&output_file, argv + input_index, input_count, call_sign, seconds);
		else
			ping<value, cmplx, 44100>(&output_file, call_sign, ping_value);
		break;
	case 48000:
		if (seconds)
			encode<value, cmplx, 48000>(&output_file, argv + input_index, input_count, call_sign, seconds);
		else
			ping<value, cmplx, 48000>(&output_file, call_sign, ping_value);
		break;
	default:
		std::cerr << "Unsupported sample rate." << std::endl;
		return 1;
	}
	output_file.silence(output_rate);

	return 0;
}


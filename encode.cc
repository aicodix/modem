/*
OFDM modem encoder

Copyright 2021 Ahmet Inan <inan@aicodix.de>
*/

#include "encoder.hh"
#include "native.hh"

template <typename value, typename cmplx>
int read_analog(cmplx *values, const char *input_name)
{
	const int analog_max = Common::analog_max;
	DSP::ReadWAV<value> input_file(input_name);
	if (!input_file.good()) {
		std::cerr << "Couldn't open file \"" << input_name << "\" for reading." << std::endl;
		return 0;
	}
	if (input_file.channels() != 2 || input_file.bits() != 32) {
		std::cerr << "Only complex values (two channels, 32 bit float) supported." << std::endl;
		return 0;
	}
	int count = 0;
	while (count < analog_max) {
		input_file.read(reinterpret_cast<value *>(values + count), 1);
		if (!input_file.good())
			break;
		++count;
	}
	if (count == analog_max) {
		cmplx tmp;
		input_file.read(reinterpret_cast<value *>(&tmp), 1);
		if (input_file.good())
			std::cerr << "Ignoring complex values beyond " << analog_max << " in file \"" << input_name << "\"." << std::endl;
	}
	if (!count)
		std::cerr << "No complex values in file \"" << input_name << "\"." << std::endl;
	return count;
}

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
void encode(DSP::WritePCM<value> *pcm, const char *const *input_names, int input_count, int64_t call_sign, int oper_mode)
{
	auto *encoder = new Encoder<value, cmplx, rate>();
	auto *values = new cmplx[Common::analog_max];
	uint8_t payload[Common::data_max];
	bool started = encoder->start(oper_mode);
	if (!(oper_mode & 128))
		print_setup(*encoder, oper_mode);
	if (started) {
		pcm->write(encoder->audio, encoder->audio_len);
		for (int input_index = 0; input_index < input_count; ++input_index) {
			const char *input_name = input_names[input_index];
			if (input_count == 1 && input_name[0] == '-' && input_name[1] == 0)
				input_name = "/dev/stdin";
			if (oper_mode & 128) {
				int count = read_analog<value>(values, input_name);
				if (!count || !encoder->begin(values, count, call_sign))
					continue;
				print_setup(*encoder, encoder->oper_mode);
			} else {
				int length = read_digital(payload, encoder->data_bytes, input_name);
				if (length < 0 || !encoder->begin(payload, length, call_sign))
					continue;
			}
			while (encoder->next())
				pcm->write(encoder->audio, encoder->audio_len);
			int last = encoder->symbol_count;
			DSP::quick_sort(encoder->papr, last + 1);
			std::cerr << "PAPR (dB): " << DSP::decibel(encoder->papr[0]) << " .. " << DSP::decibel(encoder->papr[last/2]) << " .. " << DSP::decibel(encoder->papr[last]) << std::endl;
		}
		encoder->stop();
		pcm->write(encoder->audio, encoder->audio_len);
	}
	delete[] values;
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
	bool analog_mode = argc > 6 && !strcmp(argv[5], "ANALOG");
	if (!analog_mode && argc < 9) {
		std::cerr << "usage: " << argv[0] << " OUTPUT RATE BITS CALLSIGN MODULATION CODERATE FRAMESIZE INPUT.." << std::endl;
		std::cerr << "       " << argv[0] << " OUTPUT RATE BITS CALLSIGN ANALOG INPUT.." << std::endl;
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
	int oper_mode = 0;
	int input_index = 8;
	if (analog_mode) {
		oper_mode = 128;
		input_index = 6;
	} else {
		char *modulation = argv[5];
		if (!strcmp(modulation, "BPSK")) {
			oper_mode |= 0 << 4;
		} else if (!strcmp(modulation, "QPSK")) {
			oper_mode |= 1 << 4;
		} else if (!strcmp(modulation, "8PSK")) {
			oper_mode |= 2 << 4;
		} else if (!strcmp(modulation, "QAM16")) {
			oper_mode |= 3 << 4;
		} else if (!strcmp(modulation, "QAM64")) {
			oper_mode |= 4 << 4;
		} else if (!strcmp(modulation, "QAM256")) {
			oper_mode |= 5 << 4;
		} else if (!strcmp(modulation, "QAM1024")) {
			oper_mode |= 6 << 4;
		} else if (!strcmp(modulation, "QAM4096")) {
			oper_mode |= 7 << 4;
		} else {
			std::cerr << "Unsupported modulation." << std::endl;
			return 1;
		}
		char *code_rate = argv[6];
		if (!strcmp(code_rate, "1/2")) {
			oper_mode |= 0 << 1;
		} else if (!strcmp(code_rate, "2/3")) {
			oper_mode |= 1 << 1;
		} else if (!strcmp(code_rate, "3/4")) {
			oper_mode |= 2 << 1;
		} else if (!strcmp(code_rate, "5/6")) {
			oper_mode |= 3 << 1;
		} else {
			std::cerr << "Unsupported code rate." << std::endl;
			return 1;
		}
		char *frame_size = argv[7];
		if (!strcmp(frame_size, "short")) {
			oper_mode |= 0;
		} else if (!strcmp(frame_size, "normal")) {
			oper_mode |= 1;
		} else {
			std::cerr << "Unsupported frame size." << std::endl;
			return 1;
		}
	}
	std::cerr << std::fixed << std::setprecision(1);
	typedef float value;
	typedef DSP::Complex<value> cmplx;
	DSP::WriteWAV<value> output_file(output_name, output_rate, output_bits, 1);
	output_file.silence(output_rate);
	int input_count = argc - input_index;
	switch (output_rate) {
	case 44100:
		encode<value, cmplx, 44100>(&output_file, argv + input_index, input_count, call_sign, oper_mode);
		break;
	case 48000:
		encode<value, cmplx, 48000>(&output_file, argv + input_index, input_count, call_sign, oper_mode);
		break;
	default:
		std::cerr << "Unsupported sample rate." << std::endl;
		return 1;
	}
	output_file.silence(output_rate);

	return 0;
}


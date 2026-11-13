/*
OFDM modem decoder

Copyright 2021 Ahmet Inan <inan@aicodix.de>
*/

#include "decoder.hh"
#include "native.hh"

void base40_decoder(char *str, int64_t val, int len)
{
	for (int i = len-1; i >= 0; --i, val /= 40)
		str[i] = "   /0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"[val%40];
}

template <typename value, typename cmplx, int rate>
void report(Decoder<value, cmplx, rate> *decoder, int events)
{
	typedef Decoder<value, cmplx, rate> decoder_type;
	typedef DSP::Const<value> Const;
	if (events & decoder_type::found_sync) {
		std::cerr << "symbol pos: " << decoder->sample_count - decoder_type::buffer_len + decoder->symbol_pos << std::endl;
		std::cerr << "coarse cfo: " << decoder->cfo_rad * (rate / Const::TwoPi()) << " Hz " << std::endl;
		std::cerr << "coarse sfo: " << -1000000 * decoder->sfo_rad / Const::TwoPi() << " ppm" << std::endl;
		std::cerr << "residual cfo: " << 1000 * decoder->residual_rad * rate / (Const::TwoPi() * decoder_type::symbol_len) << " mHz" << std::endl;
	}
	if (events & decoder_type::damaged_seed)
		std::cerr << "seed value damaged" << std::endl;
	if (events & decoder_type::preamble_error)
		std::cerr << "preamble decoding error." << std::endl;
	if (events & decoder_type::unsupported_call)
		std::cerr << "call sign unsupported." << std::endl;
	if (events & decoder_type::read_meta) {
		char call_sign[10];
		base40_decoder(call_sign, decoder->call_sign, 9);
		call_sign[9] = 0;
		std::cerr << "call sign: " << call_sign << std::endl;
		print_setup(*decoder, decoder->meta_mode);
	}
	if (events & decoder_type::known_mode)
		std::cerr << "oper mode: " << decoder->oper_mode << std::endl;
	if (events & decoder_type::erased_symbol)
		std::cerr << "seed value damaged, erasing symbol " << decoder->erased_index << std::endl;
	if (events & decoder_type::frame_end)
		std::cerr << "Es/N0 (dB): " << decoder->es_n0[0] << " .. " << decoder->es_n0[1] << " .. " << decoder->es_n0[2] << std::endl;
	if (events & decoder_type::payload_error)
		std::cerr << "payload decoding error." << std::endl;
	if (events & decoder_type::reencoded)
		std::cerr << "bit flips: " << decoder->bit_flips << " of " << (1 << decoder->code_order) << std::endl;
	if (events & decoder_type::real_snr)
		std::cerr << "real Es/N0 (dB): " << decoder->real_es_n0[0] << " .. " << decoder->real_es_n0[1] << " .. " << decoder->real_es_n0[2] << std::endl;
}

template <typename value, typename cmplx, int rate>
void write_output(Decoder<value, cmplx, rate> *decoder, const char *output_name, int output_count)
{
	if (output_count == 1 && output_name[0] == '-' && output_name[1] == 0)
		output_name = "/dev/stdout";
	std::ofstream output_file(output_name, std::ios::binary | std::ios::trunc);
	if (output_file.bad()) {
		std::cerr << "Couldn't open file \"" << output_name << "\" for writing." << std::endl;
		return;
	}
	for (int i = 0; i < decoder->data_bytes; ++i)
		output_file.put(decoder->data[i]);
}

template <typename value, typename cmplx, int rate>
void decode(DSP::ReadPCM<value> *pcm, const char *const *output_names, int output_count)
{
	auto *decoder = new Decoder<value, cmplx, rate>();
	int output_index = 0;
	bool end = false;
	while (output_index < output_count) {
		value sample = 0;
		if (!end) {
			pcm->read(&sample, 1);
			if (!pcm->good()) {
				end = true;
				sample = 0;
			}
		}
		if (end && !decoder->busy())
			break;
		int events = decoder->feed(sample);
		if (events)
			report(decoder, events);
		if (events & decoder->decoded_frame)
			write_output(decoder, output_names[output_index++], output_count);
	}
	delete decoder;
}

int main(int argc, char **argv)
{
	if (argc < 3) {
		std::cerr << "usage: " << argv[0] << " INPUT OUTPUT.." << std::endl;
		return 1;
	}

	typedef float value;
	typedef DSP::Complex<value> cmplx;

	const char *input_name = argv[1];
	if (input_name[0] == '-' && input_name[1] == 0)
		input_name = "/dev/stdin";

	DSP::ReadWAV<value> input_file(input_name);

	if (input_file.channels() != 1) {
		std::cerr << "Only real signal (one channel) supported." << std::endl;
		return 1;
	}

	std::cerr << std::fixed << std::setprecision(1);
	int output_count = argc - 2;
	switch (input_file.rate()) {
	case 44100:
		decode<value, cmplx, 44100>(&input_file, argv+2, output_count);
		break;
	case 48000:
		decode<value, cmplx, 48000>(&input_file, argv+2, output_count);
		break;
	default:
		std::cerr << "Unsupported sample rate." << std::endl;
		return 1;
	}
	return 0;
}


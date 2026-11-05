/*
OFDM modem encoder

Copyright 2021 Ahmet Inan <inan@aicodix.de>
*/

#include "common.hh"
#include "utils.hh"
#include "polar_encoder.hh"

template <typename value, typename cmplx, int rate>
struct Encoder : public Common
{
	typedef int8_t code_type;
	static const int guard_len = rate / 300;
	static const int symbol_len = guard_len * 40;
	static const int tone_off = (1500 * symbol_len) / rate - tone_count / 2;
	static const int audio_max = symbol_len + guard_len;
	DSP::FastFourierTransform<symbol_len, cmplx, -1> fwd;
	DSP::FastFourierTransform<symbol_len, cmplx, 1> bwd;
	CODE::PolarEncoder<code_type> polar_encoder;
	CODE::MLS seq1;
	code_type code[bits_max], perm[bits_max], mesg[bits_max], meta[data_tones];
	cmplx fdom[symbol_len];
	cmplx tdom[symbol_len];
	cmplx test[symbol_len];
	cmplx kern[symbol_len];
	cmplx guard[guard_len];
	cmplx tone[tone_count];
	cmplx temp[tone_count];
	cmplx analog[analog_max];
	value weight[guard_len];
	value papr[symbols_max];
	value audio[audio_max];
	int audio_len;
	int symbol_index, tone_index, meta_index;

	static int bin(int carrier)
	{
		return (carrier + symbol_len) % symbol_len;
	}
	static int nrz(bool bit)
	{
		return 1 - 2 * bit;
	}
	void emit(const cmplx *buf, int len)
	{
		for (int i = 0; i < len; ++i)
			audio[audio_len++] = buf[i].real();
	}
	void clipping_and_filtering(value scale)
	{
		for (int i = 0; i < symbol_len; ++i) {
			value pwr = norm(tdom[i]);
			if (pwr > value(1))
				tdom[i] /= sqrt(pwr);
		}
		fwd(fdom, tdom);
		for (int i = 0; i < symbol_len; ++i) {
			int j = bin(i + tone_off);
			if (i >= tone_count)
				fdom[j] = 0;
			else
				fdom[j] *= 1 / (scale * symbol_len);
		}
		bwd(tdom, fdom);
		for (int i = 0; i < symbol_len; ++i)
			tdom[i] *= scale;
		auto clamp = [](value v){ return v < value(-1) ? value(-1) : v > value(1) ? value(1) : v; };
		for (int i = 0; i < symbol_len; ++i)
			tdom[i] = cmplx(clamp(tdom[i].real()), clamp(tdom[i].imag()));
	}
	void symbol(int symbol_number)
	{
		value scale = value(0.5) / std::sqrt(value(tone_count));
		if (symbol_number < 0) {
			for (int i = 0; i < symbol_len; ++i)
				fdom[i] = 0;
			for (int i = 0; i < tone_count; ++i)
				fdom[bin(i+tone_off)] = tone[i];
			bwd(tdom, fdom);
			for (int i = 0; i < symbol_len; ++i)
				tdom[i] *= scale;
		} else {
			value best_papr = 1000;
			for (int seed_value = 0; seed_value < 128; ++seed_value) {
				for (int i = 0; i < tone_count; ++i)
					temp[i] = tone[i];
				hadamard_encoder(seed, seed_value);
				for (int i = 0; i < seed_tones; ++i)
					temp[i*block_length+seed_off] *= seed[i];
				if (seed_value) {
					CODE::MLS seq(mls2_poly, seed_value);
					for (int i = 0; i < tone_count; ++i)
						if (i % block_length != seed_off)
							temp[i] *= nrz(seq());
				}
				for (int i = 0; i < symbol_len; ++i)
					fdom[i] = 0;
				for (int i = 0; i < tone_count; ++i)
					fdom[bin(i+tone_off)] = temp[i];
				bwd(test, fdom);
				for (int i = 0; i < symbol_len; ++i)
					test[i] *= scale;
				value peak = 0, mean = 0;
				for (int i = 0; i < symbol_len; ++i) {
					value power(norm(test[i]));
					peak = std::max(peak, power);
					mean += power;
				}
				mean /= symbol_len;
				value test_papr(peak / mean);
				if (test_papr < best_papr) {
					best_papr = test_papr;
					papr[symbol_number] = test_papr;
					for (int i = 0; i < symbol_len; ++i)
						tdom[i] = test[i];
					if (test_papr < 5)
						break;
				}
			}
		}
		clipping_and_filtering(scale);
		if (symbol_number != -1) {
			for (int i = 0; i < guard_len; ++i)
				guard[i] = DSP::lerp(guard[i], tdom[i+symbol_len-guard_len], weight[i]);
			emit(guard, guard_len);
		}
		for (int i = 0; i < guard_len; ++i)
			guard[i] = tdom[i];
		emit(tdom, symbol_len);
	}
	void finish()
	{
		for (int i = 0; i < guard_len; ++i)
			guard[i] *= 1 - weight[i];
		emit(guard, guard_len);
		for (int i = 0; i < guard_len; ++i)
			guard[i] = 0;
	}
	void leading_noise(int num = 1)
	{
		CODE::MLS noise(mls2_poly);
		for (int j = 0; j < num; ++j) {
			for (int i = 0; i < tone_count; ++i)
				tone[i] = nrz(noise());
			symbol(-3);
		}
	}
	void schmidl_cox()
	{
		CODE::MLS seq0(mls0_poly, mls0_seed);
		for (int i = 0; i < tone_count; ++i)
			tone[i] = nrz(seq0());
	}
	void meta_data(uint64_t md)
	{
		for (int i = 0; i < 56; ++i)
			mesg[i] = nrz((md >> i) & 1);
		crc0.reset();
		crc0(md << 8);
		for (int i = 0; i < 16; ++i)
			mesg[i+56] = nrz((crc0() >> i) & 1);
		polar_encoder(code, mesg, frozen_256_72, 8);
		shuffle(meta, code, 8, true);
	}
	cmplx map_bits(code_type *b, int bits)
	{
		switch (bits) {
		case 1:
			return PhaseShiftKeying<1, cmplx, code_type>::map(b);
		case 2:
			return PhaseShiftKeying<2, cmplx, code_type>::map(b);
		case 3:
			return PhaseShiftKeying<3, cmplx, code_type>::map(b);
		case 4:
			return QuadratureAmplitudeModulation<4, cmplx, code_type>::map(b);
		case 6:
			return QuadratureAmplitudeModulation<6, cmplx, code_type>::map(b);
		case 8:
			return QuadratureAmplitudeModulation<8, cmplx, code_type>::map(b);
		case 10:
			return QuadratureAmplitudeModulation<10, cmplx, code_type>::map(b);
		case 12:
			return QuadratureAmplitudeModulation<12, cmplx, code_type>::map(b);
		}
		return 0;
	}
	void guard_interval_weights()
	{
		for (int i = 0; i < guard_len / 4; ++i)
			weight[i] = 0;
		for (int i = guard_len / 4; i < guard_len / 4 + guard_len / 2; ++i) {
			value x = value(i - guard_len / 4) / value(guard_len / 2 - 1);
			weight[i] = value(0.5) * (value(1) - std::cos(DSP::Const<value>::Pi() * x));
		}
		for (int i = guard_len / 4 + guard_len / 2; i < guard_len; ++i)
			weight[i] = 1;
	}
	Encoder() : seq1(mls1_poly), audio_len(0), symbol_index(-3)
	{
		guard_interval_weights();
	}
	bool start(int oper_mode)
	{
		audio_len = 0;
		bool analog_mode = oper_mode & 128;
		if (!analog_mode && !setup(oper_mode))
			return false;
		leading_noise();
		return true;
	}
	void prepare(int64_t call_sign)
	{
		meta_data((call_sign << 8) | oper_mode);
		seq1.reset();
		symbol_index = -2;
		tone_index = 0;
		meta_index = 0;
	}
	bool begin(const uint8_t *payload, int length, int64_t call_sign)
	{
		symbol_index = -3;
		if (length < 0)
			return false;
		for (int i = 0; i < data_bytes; ++i)
			data[i] = i < length ? payload[i] : 0;
		CODE::Xorshift32 scrambler;
		for (int i = 0; i < data_bytes; ++i)
			data[i] ^= scrambler();
		for (int i = 0; i < data_bits; ++i)
			mesg[i] = nrz(CODE::get_le_bit(data, i));
		crc1.reset();
		for (int i = 0; i < data_bytes; ++i)
			crc1(data[i]);
		for (int i = 0; i < 32; ++i)
			mesg[i+data_bits] = nrz((crc1()>>i)&1);
		polar_encoder(code, mesg, frozen_bits, code_order);
		shuffle(perm, code, code_order, true);
		prepare(call_sign);
		return true;
	}
	bool begin(const cmplx *values, int count, int64_t call_sign)
	{
		symbol_index = -3;
		if (count < 1 || count > analog_max)
			return false;
		if (!setup(128 | (count - 1) / data_tones))
			return false;
		for (int i = 0; i < count; ++i)
			analog[i] = values[i];
		for (int i = count; i % data_tones; ++i)
			analog[i] = 0;
		prepare(call_sign);
		return true;
	}
	int next()
	{
		audio_len = 0;
		if (symbol_index < -2 || symbol_index > symbol_count)
			return 0;
		if (symbol_index == -2)
			schmidl_cox();
		if (symbol_index < 0) {
			symbol(symbol_index++);
			return audio_len;
		}
		bool analog_mode = oper_mode & 128;
		int j = symbol_index++;
		seed_off = (block_skew * j + first_seed) % block_length;
		for (int i = 0; i < tone_count; ++i) {
			if (i % block_length == seed_off) {
				tone[i] = nrz(seq1());
			} else if (j && analog_mode) {
				tone[i] = analog[tone_index++];
			} else if (j) {
				int bits = bits_at(tone_index);
				tone[i] = map_bits(perm+tone_index, bits);
				tone_index += bits;
			} else {
				tone[i] = map_bits(meta+meta_index++, 1);
			}
		}
		symbol(j);
		if (j == symbol_count) {
			DSP::quick_sort(papr, symbol_count + 1);
			std::cerr << "PAPR (dB): " << DSP::decibel(papr[0]) << " .. " << DSP::decibel(papr[symbol_count/2]) << " .. " << DSP::decibel(papr[symbol_count]) << std::endl;
		}
		return audio_len;
	}
	void stop()
	{
		audio_len = 0;
		finish();
	}
};

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
	if (encoder->start(oper_mode)) {
		pcm->write(encoder->audio, encoder->audio_len);
		for (int input_index = 0; input_index < input_count; ++input_index) {
			const char *input_name = input_names[input_index];
			if (input_count == 1 && input_name[0] == '-' && input_name[1] == 0)
				input_name = "/dev/stdin";
			if (oper_mode & 128) {
				int count = read_analog<value>(values, input_name);
				if (!count || !encoder->begin(values, count, call_sign))
					continue;
			} else {
				int length = read_digital(payload, encoder->data_bytes, input_name);
				if (length < 0 || !encoder->begin(payload, length, call_sign))
					continue;
			}
			while (encoder->next())
				pcm->write(encoder->audio, encoder->audio_len);
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


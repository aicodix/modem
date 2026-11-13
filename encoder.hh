/*
OFDM modem encoder

Copyright 2021 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

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
		}
		return 0;
	}
	void guard_interval_weights()
	{
		for (int i = 0; i < guard_len / 4; ++i)
			weight[i] = 0;
		for (int i = guard_len / 4; i < guard_len / 4 + guard_len / 2; ++i)
			weight[i] = value(0.5) * (value(1) - DSP::UnitCircle<value>::cos(i - guard_len / 4, 2 * (guard_len / 2 - 1)));
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
		if (!setup(oper_mode))
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
		int j = symbol_index++;
		seed_off = (block_skew * j + first_seed) % block_length;
		for (int i = 0; i < tone_count; ++i) {
			if (i % block_length == seed_off) {
				tone[i] = nrz(seq1());
			} else if (j) {
				int bits = bits_at(tone_index);
				tone[i] = map_bits(perm+tone_index, bits);
				tone_index += bits;
			} else {
				tone[i] = map_bits(meta+meta_index++, 1);
			}
		}
		symbol(j);
		return audio_len;
	}
	void stop()
	{
		audio_len = 0;
		finish();
	}
};

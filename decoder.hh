/*
OFDM modem decoder

Copyright 2021 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

#include "common.hh"
#include "schmidl_cox.hh"
#include "frequency_sync.hh"
#include "bip_buffer.hh"
#include "theil_sen.hh"
#include "filter.hh"
#include "window.hh"
#include "fir.hh"
#include "phasor.hh"
#include "spline.hh"
#include "delay.hh"
#include "polar_list_decoder.hh"
#include "polar_encoder.hh"
#include "hadamard_decoder.hh"

template <typename value, typename cmplx, int rate>
struct Decoder : Common
{
	typedef int16_t code_type;
	typedef SIMD<code_type, 32> mesg_type;
	typedef DSP::Const<value> Const;
	static const int guard_len = rate / 300;
	static const int symbol_len = guard_len * 40;
	static const int filter_len = rate / 100 | 1;
	static const int extended_len = symbol_len + guard_len;
	static const int buffer_len = 5 * extended_len;
	static const int search_pos = extended_len;
	// the frequency correlator looks for the sync up to two and a half symbols back
	static const int freq_pos = 5 * symbol_len / 2 + guard_len;
	static const int tone_off = - tone_count / 2;
	// what happened during a call of feed(), as reported by its return value
	static const int found_sync = 1 << 0;
	static const int damaged_seed = 1 << 1;
	static const int preamble_error = 1 << 2;
	static const int unsupported_call = 1 << 3;
	static const int read_meta = 1 << 4;
	static const int known_mode = 1 << 5;
	static const int erased_symbol = 1 << 6;
	static const int frame_end = 1 << 7;
	static const int payload_error = 1 << 8;
	static const int reencoded = 1 << 9;
	static const int real_snr = 1 << 10;
	static const int decoded_frame = 1 << 11;
	static const int new_symbol = 1 << 12;
	DSP::FastFourierTransform<symbol_len, cmplx, -1> fwd;
	DSP::Phasor<cmplx> mix;
	DSP::FIR<cmplx, value, filter_len> lowpass;
	DSP::BipBuffer<cmplx, buffer_len> input_hist;
	DSP::TheilSenEstimator<value, tone_count> tse;
	SchmidlCox<value, cmplx, search_pos, symbol_len, guard_len> correlator;
	FrequencySync<value, cmplx, freq_pos, symbol_len, guard_len> freq_correlator;
	CODE::HadamardDecoder<7> hadamard_decoder;
	CODE::PolarListDecoder<mesg_type, code_max> polar_decoder;
	CODE::PolarEncoder<code_type> polar_encoder;
	mesg_type mesg[bits_max];
	code_type code[bits_max], perm[bits_max], meta[data_tones];
	cmplx demod[tone_count], chan[tone_count], tone[tone_count];
	cmplx knots[seed_tones];
	value innov[seed_tones], sorted[seed_tones];
	value gain[tone_count];
	// equalized values of the data tones
	cmplx values[analog_max];
	code_type reenc[bits_max];
	bool lost[symbols_max];
	cmplx fdom[symbol_len], tdom[symbol_len];
	value index[tone_count], phase[tone_count];
	value snr[symbols_max];
	value cfo_rad, sfo_rad, residual_rad;
	value es_n0[3], real_es_n0[3];
	// real Es/N0 in dB of each tone, averaged over the payload symbols of the frame
	value real_tone_snr[tone_count];
	value tone_sp[tone_count], tone_np[tone_count];
	// precision of each tone of the latest digital symbol, its estimated Es/N0, zero when erased
	value tone_precision[tone_count];
	// noise of each tone before equalization, from the preamble and then the seed tones
	value tone_noise[tone_count], local_noise[tone_count], sorted_noise[tone_count];
	value error_power[tone_count];
	int64_t call_sign;
	int meta_mode, erased_index, bit_flips;
	int events;
	int symbol_pos;
	int crc_bits;
	DSP::Phasor<cmplx> osc;
	CODE::MLS seq1;
	int sample_count, countdown;
	int symbol_index, bit_index, snr_count;
	bool analog_mode;
	// sync detection in the frequency instead of the time domain, for FM or AM without CFO
	bool frequency_sync;

	static int bin(int carrier)
	{
		return (carrier + symbol_len) % symbol_len;
	}
	static value nrz(bool bit)
	{
		return 1 - 2 * bit;
	}
	static cmplx demod_or_erase(cmplx curr, cmplx prev, bool limit = true)
	{
		if (norm(prev) > 0) {
			cmplx demod = curr / prev;
			if (!limit || norm(demod) < 4)
				return demod;
		}
		return 0;
	}
	const cmplx *mls0_seq()
	{
		CODE::MLS seq0(mls0_poly, mls0_seed);
		for (int i = 0; i < symbol_len; ++i)
			fdom[i] = 0;
		value cur = 0, prv = 0;
		for (int i = 0; i < tone_count; ++i, prv = cur)
			fdom[bin(i+tone_off)] = prv * (cur = nrz(seq0()));
		return fdom;
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
	void demap_soft(code_type *b, cmplx c, value precision, int bits)
	{
		switch (bits) {
		case 1:
			return PhaseShiftKeying<1, cmplx, code_type>::soft(b, c, precision);
		case 2:
			return PhaseShiftKeying<2, cmplx, code_type>::soft(b, c, precision);
		case 3:
			return PhaseShiftKeying<3, cmplx, code_type>::soft(b, c, precision);
		case 4:
			return QuadratureAmplitudeModulation<4, cmplx, code_type>::soft(b, c, precision);
		case 6:
			return QuadratureAmplitudeModulation<6, cmplx, code_type>::soft(b, c, precision);
		case 8:
			return QuadratureAmplitudeModulation<8, cmplx, code_type>::soft(b, c, precision);
		case 10:
			return QuadratureAmplitudeModulation<10, cmplx, code_type>::soft(b, c, precision);
		case 12:
			return QuadratureAmplitudeModulation<12, cmplx, code_type>::soft(b, c, precision);
		}
	}
	void demap_hard(code_type *b, cmplx c, int bits)
	{
		switch (bits) {
		case 1:
			return PhaseShiftKeying<1, cmplx, code_type>::hard(b, c);
		case 2:
			return PhaseShiftKeying<2, cmplx, code_type>::hard(b, c);
		case 3:
			return PhaseShiftKeying<3, cmplx, code_type>::hard(b, c);
		case 4:
			return QuadratureAmplitudeModulation<4, cmplx, code_type>::hard(b, c);
		case 6:
			return QuadratureAmplitudeModulation<6, cmplx, code_type>::hard(b, c);
		case 8:
			return QuadratureAmplitudeModulation<8, cmplx, code_type>::hard(b, c);
		case 10:
			return QuadratureAmplitudeModulation<10, cmplx, code_type>::hard(b, c);
		case 12:
			return QuadratureAmplitudeModulation<12, cmplx, code_type>::hard(b, c);
		}
	}
	void reencode(int best)
	{
		for (int i = 0; i < crc_bits; ++i)
			perm[i] = mesg[i].v[best] < 0 ? -1 : 1;
		polar_encoder(reenc, perm, frozen_bits, code_order);
		int code_bits = 1 << code_order;
		bit_flips = 0;
		for (int i = 0; i < code_bits; ++i)
			bit_flips += code[i] * reenc[i] <= 0;
		events |= reencoded;
		shuffle(perm, reenc, code_order, true);
		int count = 0;
		for (int i = 0; i < tone_count; ++i)
			tone_sp[i] = tone_np[i] = 0;
		for (int j = 1, k = 0, n = 0; j <= symbol_count; ++j) {
			seed_off = (block_skew * j + first_seed) % block_length;
			value sp = 0, np = 0;
			for (int i = 0; i < tone_count; ++i) {
				if (i % block_length != seed_off) {
					int bits = bits_at(k);
					cmplx hard = map_bits(perm+k, bits);
					k += bits;
					value s = norm(hard), e = norm(values[n++] - hard);
					sp += s;
					np += e;
					if (!lost[j]) {
						tone_sp[i] += s;
						tone_np[i] += e;
					}
				}
			}
			if (!lost[j])
				snr[count++] = sp / std::max(np, value(1e-30) * sp);
		}
		if (!count)
			return;
		// a tone that carried no payload, as a seed tone in every symbol of a short frame, gets 0 dB
		for (int i = 0; i < tone_count; ++i)
			real_tone_snr[i] = tone_sp[i] > 0 ? DSP::decibel(tone_sp[i] / std::max(tone_np[i], value(1e-30) * tone_sp[i])) : 0;
		DSP::quick_sort(snr, count);
		real_es_n0[0] = DSP::decibel(snr[0]);
		real_es_n0[1] = DSP::decibel(snr[(count-1)/2]);
		real_es_n0[2] = DSP::decibel(snr[count-1]);
		events |= real_snr;
	}
	int64_t meta_data()
	{
		shuffle(code, perm, 8, false);
		polar_decoder(nullptr, mesg, code, frozen_256_72, 8);
		int best = -1;
		for (int k = 0; k < mesg_type::SIZE; ++k) {
			crc0.reset();
			for (int i = 0; i < 72; ++i)
				crc0(mesg[i].v[k] < 0);
			if (crc0() == 0) {
				best = k;
				break;
			}
		}
		if (best < 0)
			return -1;
		for (int i = 0; i < 72; ++i)
			perm[i] = mesg[i].v[best] < 0 ? -1 : 1;
		polar_encoder(code, perm, frozen_256_72, 8);
		shuffle(meta, code, 8, true);
		uint64_t md = 0;
		for (int i = 0; i < 56; ++i)
			md |= uint64_t(mesg[i].v[best] < 0) << i;
		return md;
	}
	void process(const cmplx *buf)
	{
		auto clamp = [](int v){ return v < -127 ? -127 : v > 127 ? 127 : v; };
		int j = symbol_index;
		events |= new_symbol;
		seed_off = (block_skew * j + first_seed) % block_length;
		if (j) {
			for (int i = 0; i < symbol_len; ++i)
				tdom[i] = buf[i] * osc();
			for (int i = 0; i < guard_len; ++i)
				osc();
			fwd(fdom, tdom);
		}
		for (int i = 0; i < tone_count; ++i)
			tone[i] = fdom[bin(i+tone_off)];
		for (int i = seed_off; i < tone_count; i += block_length)
			tone[i] *= nrz(seq1());
		for (int i = 0; i < tone_count; ++i)
			demod[i] = demod_or_erase(tone[i], chan[i], !analog_mode);
		for (int i = 0; i < seed_tones; ++i)
			seed[i] = clamp(std::nearbyint(127 * demod[i*block_length+seed_off].real()));
		int seed_value = hadamard_decoder(seed);
		bool erased = seed_value < 0;
		lost[j] = erased;
		if (erased && !j) {
			events |= damaged_seed;
			oper_mode = -1;
			return;
		}
		if (erased) {
			events |= erased_symbol;
			erased_index = j;
			for (int i = 0; i < tone_count; ++i)
				demod[i] = 0;
		} else {
			hadamard_encoder(seed, seed_value);
			for (int i = 0; i < seed_tones; ++i) {
				tone[block_length*i+seed_off] *= seed[i];
				demod[block_length*i+seed_off] *= seed[i];
			}
			for (int i = 0; i < seed_tones; ++i) {
				index[i] = tone_off + block_length * i + seed_off;
				phase[i] = arg(demod[block_length*i+seed_off]);
			}
			tse.compute(index, phase, seed_tones);
			//std::cerr << "Theil-Sen slope = " << tse.slope() << std::endl;
			//std::cerr << "Theil-Sen yint = " << tse.yint() << std::endl;
			for (int i = 0; i < tone_count; ++i)
				demod[i] *= DSP::polar<value>(1, -tse(i+tone_off));
			for (int i = 0; i < tone_count; ++i)
				chan[i] *= DSP::polar<value>(1, tse(i+tone_off));
			if (j) {
				// update all tones with the seed tones of this symbol, interpolated values with a lower weight
				// keep the estimate where a seed tone changed far more than the median, as with a carrier
				for (int i = 0; i < seed_tones; ++i)
					sorted[i] = innov[i] = norm(tone[block_length*i+seed_off] - chan[block_length*i+seed_off]);
				value limit = 10 * DSP::quick_select(sorted, seed_tones / 2, seed_tones);
				for (int i = 0; i < seed_tones; ++i)
					knots[i] = innov[i] <= limit ? tone[block_length*i+seed_off] : chan[block_length*i+seed_off];
				DSP::UniformNaturalCubicSpline<seed_tones, cmplx, value> spline(knots, seed_off, block_length);
				for (int i = 0; i < tone_count; ++i) {
					if (i % block_length == seed_off) {
						chan[i] = DSP::lerp(chan[i], knots[i/block_length], value(0.5));
					} else {
						chan[i] = DSP::lerp(chan[i], spline(i), value(1) / value(3));
						demod[i] = demod_or_erase(tone[i], chan[i], !analog_mode);
					}
				}
			}
			if (seed_value) {
				CODE::MLS seq(mls2_poly, seed_value);
				for (int i = 0; i < tone_count; ++i)
					if (i % block_length != seed_off)
						demod[i] *= nrz(seq());
			}
		}
		if (analog_mode) {
			value sp = 0, np = 0;
			for (int i = seed_off; !erased && i < tone_count; i += block_length) {
				sp += 1;
				np += norm(demod[i] - cmplx(1, 0));
			}
			if (!erased)
				snr[snr_count++] = sp / std::max(np, value(1e-9) * sp);
			for (int i = 0; i < tone_count; ++i)
				if (i % block_length != seed_off)
					values[bit_index++] = demod[i];
		} else {
			value sp = 0, np = 0;
			for (int i = 0, l = bit_index; !erased && i < tone_count; ++i) {
				cmplx hard(1, 0);
				if (i % block_length != seed_off) {
					int bits = bits_at(l);
					demap_hard(perm+l, demod[i], bits);
					hard = map_bits(perm+l, bits);
					l += bits;
				}
				cmplx error = demod[i] - hard;
				sp += norm(hard);
				np += norm(error);
				// the error before equalization, flat for white noise
				gain[i] = norm(chan[i]);
				error_power[i] = norm(error) * gain[i];
			}
			if (erased) {
				// erased symbols get zero LLRs
				for (int i = 0; i < tone_count; ++i)
					tone_precision[i] = 0;
			} else {
				snr[snr_count++] = sp / np;
				// the noise of each tone from known symbols, all tones of the preamble, then the seed tones
				if (!j) {
					for (int i = 0; i < tone_count; ++i)
						tone_noise[i] = error_power[i];
				} else {
					for (int i = seed_off; i < tone_count; i += block_length)
						tone_noise[i] = DSP::lerp(tone_noise[i], error_power[i], value(0.5));
				}
				for (int i = 0; i < tone_count; ++i) {
					value sum = 0;
					int count = 0;
					for (int k = std::max(i - 8, 0); k <= std::min(i + 8, tone_count - 1); ++k, ++count)
						sum += tone_noise[k];
					sorted_noise[i] = local_noise[i] = sum / value(count);
				}
				// tones with more than twice the median noise, as from interference, get less precision
				value median = DSP::quick_select(sorted_noise, tone_count / 2, tone_count);
				value limit = 2 * median;
				// the white noise of this symbol, without those tones
				value sum = 0;
				int count = 0;
				for (int i = 0; i < tone_count; ++i) {
					if (local_noise[i] <= limit) {
						sum += error_power[i];
						++count;
					}
				}
				value es = sp / value(tone_count);
				value noise = std::max(sum / value(count), value(1e-30) * es);
				for (int i = 0; i < tone_count; ++i) {
					tone_precision[i] = es * gain[i] / noise;
					if (local_noise[i] > limit)
						tone_precision[i] *= median / local_noise[i];
				}
			}
			for (int i = 0, n = (j - 1) * data_tones; i < tone_count; ++i) {
				if (i % block_length != seed_off) {
					if (j)
						values[n++] = demod[i];
					int bits = bits_at(bit_index);
					// more resolution after rounding, quantize() caps the slope
					demap_soft(perm+bit_index, demod[i], value(300) * tone_precision[i], bits);
					bit_index += bits;
				}
			}
		}
		if (!j) {
			int64_t meta_info = meta_data();
			if (meta_info < 0) {
				events |= preamble_error;
				return;
			}
			call_sign = meta_info >> 8;
			if (call_sign == 0 || call_sign >= 262144000000000L) {
				events |= unsupported_call;
				return;
			}
			int mode = meta_info & 255;
			meta_mode = mode;
			events |= read_meta;
			if (!setup(mode))
				return;
			events |= known_mode;
			analog_mode = mode & 128;
			bit_index = 0;
			CODE::MLS seq(mls2_poly, seed_value);
			for (int i = 0, m = 0; i < tone_count; ++i) {
				if (i % block_length != seed_off) {
					tone[i] *= value(meta[m++]);
					if (seed_value)
						tone[i] *= nrz(seq());
				}
			}
			for (int i = 0; i < tone_count; ++i)
				chan[i] = DSP::lerp(chan[i], tone[i], value(0.5));
			symbol_index = 1;
			countdown = symbol_pos + symbol_len + extended_len + extended_len;
			return;
		}
		if (j < symbol_count) {
			++symbol_index;
			countdown = extended_len;
			return;
		}
		finish();
	}
	void finish()
	{
		DSP::quick_sort(snr, snr_count);
		es_n0[0] = DSP::decibel(snr[0]);
		es_n0[1] = DSP::decibel(snr[(snr_count-1)/2]);
		es_n0[2] = DSP::decibel(snr[snr_count-1]);
		events |= frame_end;
		if (analog_mode) {
			events |= decoded_frame;
			return;
		}
		crc_bits = data_bits + 32;
		shuffle(code, perm, code_order, false);
		polar_decoder(nullptr, mesg, code, frozen_bits, code_order);
		int best = -1;
		for (int k = 0; k < mesg_type::SIZE; ++k) {
			crc1.reset();
			for (int i = 0; i < crc_bits; ++i)
				crc1(mesg[i].v[k] < 0);
			if (crc1() == 0) {
				best = k;
				break;
			}
		}
		if (best < 0) {
			events |= payload_error;
			return;
		}
		for (int i = 0; i < data_bits; ++i)
			CODE::set_le_bit(data, i, mesg[i].v[best] < 0);
		reencode(best);
		CODE::Xorshift32 scrambler;
		for (int i = 0; i < data_bytes; ++i)
			data[i] ^= scrambler();
		events |= decoded_frame;
	}
	void preamble(const cmplx *buf)
	{
		symbol_pos = frequency_sync ? freq_correlator.symbol_pos : correlator.symbol_pos;
		cfo_rad = frequency_sync ? freq_correlator.cfo_rad : correlator.cfo_rad;
		events |= found_sync;

		osc.omega(-cfo_rad);
		for (int i = 0; i < symbol_len; ++i)
			tdom[i] = buf[i+symbol_pos] * osc();
		fwd(fdom, tdom);
		for (int i = 0; i < tone_count; ++i)
			tone[i] = fdom[bin(i+tone_off)];
		for (int i = 0; i < symbol_len; ++i)
			tdom[i] = buf[i+symbol_pos+symbol_len] * osc();
		for (int i = 0; i < guard_len; ++i)
			osc();
		fwd(fdom, tdom);
		for (int i = 0; i < tone_count; ++i)
			chan[i] = fdom[bin(i+tone_off)];
		for (int i = 0; i < tone_count; ++i) {
			index[i] = tone_off + i;
			phase[i] = arg(demod_or_erase(chan[i], tone[i]));
		}
		tse.compute(index, phase, tone_count);
		sfo_rad = tse.slope();
		residual_rad = tse.yint();
		for (int i = 0; i < tone_count; ++i)
			tone[i] *= DSP::polar<value>(1, tse(i+tone_off));
		for (int i = 0; i < tone_count; ++i)
			chan[i] = DSP::lerp(chan[i], tone[i], value(0.5));
		CODE::MLS seq0(mls0_poly, mls0_seed);
		for (int i = 0; i < tone_count; ++i)
			chan[i] *= nrz(seq0());
		for (int i = 0; i < symbol_len; ++i)
			tdom[i] = buf[i+symbol_pos+symbol_len+extended_len] * osc();
		for (int i = 0; i < guard_len; ++i)
			osc();
		fwd(fdom, tdom);
		seq1.reset();
		mod_bits = 1;
		oper_mode = -1;
		symbol_count = 0;
		snr_count = 0;
		analog_mode = false;
		symbol_index = 0;
		bit_index = 0;
		process(buf);
	}
	Decoder() : lowpass(DSP::LowPass<value>(value(1275) / value(rate)), DSP::Kaiser<value>(2)), correlator(mls0_seq()), freq_correlator(mls0_seq()), seq1(mls1_poly), sample_count(0), countdown(0), frequency_sync(false)
	{
		mix.omega(-1500, rate);
	}
	int feed(value sample)
	{
		events = 0;
		const cmplx *buf = input_hist(lowpass(sample * mix()));
		bool found = frequency_sync ? freq_correlator(buf) : correlator(buf);
		if (countdown) {
			if (!--countdown)
				process(buf);
		} else {
			++sample_count;
			if (found)
				preamble(buf);
		}
		return events;
	}
	bool busy()
	{
		return countdown > 0;
	}
};

/*
OFDM modem decoder

Copyright 2021 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

#include "common.hh"
#include "raw_sync.hh"
#include "bip_buffer.hh"
#include "theil_sen.hh"
#include "filter.hh"
#include "window.hh"
#include "down_conversion.hh"
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
	// the windows start this many samples into the buffer, so resampling them can also reach before their start
	static const int window_margin = 32;
	// the raw input, enough to filter the whole buffer again after a sync
	static const int raw_len = buffer_len + filter_len;
	// the tone at 1500 Hz, where the signal is centered
	static const int carrier = 1500 * symbol_len / rate;
	static_assert(carrier * rate == 1500 * symbol_len, "1500 Hz must be a whole tone");
	// the sync detection searches a CFO of up to 250 Hz, in bins of half a tone
	static const int cfo_range = 250 * 2 * symbol_len / rate;
	// where the best window of the sync detection starts in the buffer when preamble() takes over. The sync can
	// start up to a symbol before or after it, as silence before the sync looks as clean as the guard interval after it
	static const int sync_start = symbol_len + guard_len;
	// the window ended a hop before the detection fired, and the filter delays the baseband by half its length,
	// so wait this long after it for the sync and the symbol after it, with some room for preamble() to move the windows
	static const int sync_delay = buffer_len - 3 * symbol_len - sync_start + (filter_len - 1) / 2;
	static_assert(sync_start - symbol_len >= window_margin, "sync too early in the buffer");
	static_assert(sync_start + symbol_len + 3 * symbol_len + 3 * guard_len <= buffer_len, "sync too late in the buffer");
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
	DSP::DownConversion<cmplx, filter_len> ddc;
	DSP::BipBuffer<cmplx, buffer_len> input_hist;
	DSP::BipBuffer<value, raw_len> raw_hist;
	DSP::TheilSenEstimator<value, tone_count> tse;
	RawSync<value, cmplx, symbol_len, carrier, cfo_range> raw_sync;
	CODE::HadamardDecoder<7> hadamard_decoder;
	CODE::PolarListDecoder<mesg_type, code_max> polar_decoder;
	CODE::PolarEncoder<code_type> polar_encoder;
	mesg_type mesg[bits_max];
	code_type code[bits_max], perm[bits_max], meta[data_tones];
	cmplx demod[tone_count], chan[tone_count], tone[tone_count];
	cmplx knots[seed_tones];
	value innov[seed_tones], sorted[seed_tones];
	// equalized values of the data tones
	cmplx values[analog_max];
	code_type reenc[bits_max];
	bool lost[symbols_max];
	cmplx fdom[symbol_len], tdom[symbol_len];
	value index[tone_count], phase[tone_count];
	value snr[symbols_max];
	value cfo_rad, sfo_rad, residual_rad;
	// the change of the timing since the last window, as the slope of the phase over the tones, in radians per tone
	value drift_rad;
	// the change of the phase since the last window, which a residual CFO leaves
	value phase_rad;
	// the received samples per sent sample, and where the next window starts after the whole samples of the countdown
	value step, frac;
	value es_n0[3], real_es_n0[3];
	// real Es/N0 in dB of each tone, averaged over the payload symbols of the frame
	value real_tone_snr[tone_count];
	value tone_sp[tone_count], tone_np[tone_count];
	// precision of each tone of the latest digital symbol, its estimated Es/N0, zero when erased
	value tone_precision[tone_count];
	// noise of each tone before equalization, from the sync symbols, then the seed tones
	value tone_noise[tone_count], local_noise[tone_count], sorted_noise[tone_count];
	int64_t call_sign;
	int meta_mode, erased_index, bit_flips;
	int events;
	int symbol_pos;
	int crc_bits;
	DSP::Phasor<cmplx> osc;
	CODE::MLS seq1;
	int sample_count, countdown, pending;
	int symbol_index, bit_index, snr_count;
	bool analog_mode;
	// no longer used, the sync detection works on the raw input for everything
	bool frequency_sync;
	// keep mixing and filtering while waiting for a sync, for the baseband views of the page
	bool keep_baseband;

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
			window(buf, window_margin + frac + symbol_len / 2);
			for (int i = 0; i < guard_len; ++i)
				osc();
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
			drift_rad += tse.slope();
			phase_rad += tse.yint();
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
			value sp = 0, np = 0, noise = 0;
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
				value power = norm(error) * norm(chan[i]);
				noise += power;
				// the noise of each tone from the known seed tones, which take turns over all tones
				if (i % block_length == seed_off)
					tone_noise[i] = DSP::lerp(tone_noise[i], power, value(0.5));
			}
			if (erased) {
				// erased symbols get zero LLRs
				for (int i = 0; i < tone_count; ++i)
					tone_precision[i] = 0;
			} else {
				snr[snr_count++] = sp / np;
				for (int i = 0; i < tone_count; ++i) {
					value sum = 0;
					int count = 0;
					for (int k = std::max(i - 8, 0); k <= std::min(i + 8, tone_count - 1); ++k, ++count)
						sum += tone_noise[k];
					sorted_noise[i] = local_noise[i] = sum / value(count);
				}
				// tones with more than twice the median noise, as from interference, get less precision
				value median = DSP::quick_select(sorted_noise, tone_count / 2, tone_count);
				// the white noise of this symbol, constellations have unit power
				noise = std::max(noise / value(tone_count), value(1e-30));
				for (int i = 0; i < tone_count; ++i) {
					tone_precision[i] = norm(chan[i]) / noise;
					if (local_noise[i] > 2 * median)
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
			countdown = symbol_pos + symbol_len + extended_len + extended_len - window_margin + resync();
			return;
		}
		if (j < symbol_count) {
			++symbol_index;
			countdown = extended_len + resync();
			return;
		}
		finish();
	}
	// cubic Lagrange interpolation, plenty for the oversampled baseband
	static cmplx interpolate(const cmplx *buf, value x)
	{
		int k = x;
		value t = x - value(k);
		const cmplx *y = buf + k - 1;
		value a = t + 1, b = t, c = t - 1, d = t - 2;
		return (- b * c * d / value(6)) * y[0] + (a * c * d / value(2)) * y[1] + (- a * b * d / value(2)) * y[2] + (a * b * c / value(6)) * y[3];
	}
	// the FFT of the window around middle, resampled by step, where the sync measures the timing
	void window(const cmplx *buf, value middle)
	{
		for (int i = 0; i < symbol_len; ++i)
			tdom[i] = interpolate(buf, middle + (i - symbol_len / 2) * step) * osc();
		fwd(fdom, tdom);
	}
	// predict where the next window starts from the SFO, correct it by how far the timing drifted from the last
	// prediction, as the seed tones measure it by the slope of their phase, and correct the SFO by a part of it
	int resync()
	{
		value drift = - drift_rad * symbol_len / Const::TwoPi();
		drift_rad = 0;
		step += value(0.25) * drift / value(extended_len);
		step = std::min(std::max(step, value(0.995)), value(1.005));
		// the same for the CFO, the phase changes by the residual CFO from one window to the next
		cfo_rad += value(0.25) * phase_rad / value(extended_len);
		phase_rad = 0;
		osc.omega(-cfo_rad);
		value advance = extended_len * (step - 1) + drift;
		// the window moved against the tones by the drift, and against the oscillator, which only counts samples, by all of it
		for (int i = 0; i < tone_count; ++i)
			chan[i] *= DSP::polar<value>(1, (i+tone_off) * drift * Const::TwoPi() / value(symbol_len) + advance * cfo_rad);
		value next = frac + advance;
		int shift = std::nearbyint(next);
		frac = next - value(shift);
		return shift;
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
		symbol_pos = raw_sync.symbol_pos;
		cfo_rad = raw_sync.cfo_rad;
		// resampling reaches before the window
		if (symbol_pos < window_margin)
			return;
		events |= found_sync;

		// the two sync symbols give the SFO, and the CFO the correlators got wrong. The CFO of the correlators
		// comes from all tones weighted by their power, and the phase of the tones also turns by the SFO, so
		// a channel stronger on one side gets it wrong. Take them again resampled with the SFO and with the
		// corrected CFO, so their interference between the tones doesn't stay in the channel estimate
		// The second time also measures where the first sync symbol starts, by the phase of the ratio of adjacent
		// tones against their known ratio, as the frequency correlator estimates it with the CFO still there, and
		// when the windows are off by more than an eighth of the guard interval, the third time moves them there
		step = 1;
		value offset = 0;
		for (int pass = 0; pass < 3; ++pass) {
			value middle = symbol_pos + offset + symbol_len / 2;
			osc.omega(-cfo_rad);
			window(buf, middle);
			for (int i = 0; i < tone_count; ++i)
				tone[i] = fdom[bin(i+tone_off)];
			window(buf, middle + symbol_len * step);
			for (int i = 0; i < tone_count; ++i)
				chan[i] = fdom[bin(i+tone_off)];
			for (int i = 0; i < tone_count; ++i) {
				index[i] = tone_off + i;
				phase[i] = arg(demod_or_erase(chan[i], tone[i]));
			}
			tse.compute(index, phase, tone_count);
			step -= tse.slope() / Const::TwoPi();
			step = std::min(std::max(step, value(0.995)), value(1.005));
			residual_rad = tse.yint();
			cfo_rad += residual_rad / value(symbol_len);
			if (pass == 1) {
				CODE::MLS seq0(mls0_poly, mls0_seed);
				// both sync symbols, which halves the error of the estimate
				cmplx sum = 0;
				value mag = 0;
				for (int i = 0, prev = nrz(seq0()), curr; i < tone_count - 1; ++i, prev = curr) {
					value known = prev * (curr = nrz(seq0()));
					cmplx ratio1 = known * demod_or_erase(tone[i+1], tone[i]);
					cmplx ratio2 = known * demod_or_erase(chan[i+1], chan[i]);
					sum += ratio1 + ratio2;
					mag += abs(ratio1) + abs(ratio2);
				}
				value late = arg(sum) * symbol_len / Const::TwoPi();
				// a radian is a thousand samples, so the estimate scatters by tens of samples until most ratios agree,
				// and a few samples are well within the guard interval, so only move when it is reliable and needed.
				// Keep the windows within the buffer
				if (abs(sum) < value(0.88) * mag || abs(late) <= value(guard_len / 8) || abs(late) > 2 * guard_len || symbol_pos + offset - late < window_margin)
					break;
				offset -= late;
			}
		}
		// the whole samples of the offset go into the position of the sync, as the correlators report it
		int shift = std::nearbyint(offset);
		symbol_pos += shift;
		offset -= shift;
		value middle = symbol_pos + offset + symbol_len / 2;
		sfo_rad = Const::TwoPi() * (1 - step);
		osc.omega(-cfo_rad);
		for (int i = 0; i < guard_len; ++i)
			osc();
		for (int i = 0; i < tone_count; ++i)
			tone[i] *= DSP::polar<value>(1, tse(i+tone_off));
		// the noise of each tone from the difference of the two sync symbols, which carries it twice
		for (int i = 0; i < tone_count; ++i)
			tone_noise[i] = norm(chan[i] - tone[i]) / 2;
		for (int i = 0; i < tone_count; ++i)
			chan[i] = DSP::lerp(chan[i], tone[i], value(0.5));
		CODE::MLS seq0(mls0_poly, mls0_seed);
		for (int i = 0; i < tone_count; ++i)
			chan[i] *= nrz(seq0());
		window(buf, middle + (symbol_len + extended_len) * step);
		for (int i = 0; i < guard_len; ++i)
			osc();
		seq1.reset();
		mod_bits = 1;
		oper_mode = -1;
		symbol_count = 0;
		snr_count = 0;
		analog_mode = false;
		symbol_index = 0;
		bit_index = 0;
		drift_rad = 0;
		phase_rad = 0;
		// how far the first symbol is from where the countdown takes the next one
		frac = offset + (symbol_len + extended_len) * (step - 1);
		process(buf);
	}
	Decoder() : ddc(-1500, rate, DSP::LowPass<value>(value(1275) / value(rate)), DSP::Kaiser<value>(2)), raw_sync(mls0_seq()), seq1(mls1_poly), sample_count(0), countdown(0), pending(0), frequency_sync(false), keep_baseband(false)
	{
	}
	// mixes and filters the whole raw input again, as nothing did while waiting for a sync
	const cmplx *catch_up()
	{
		const value *raw = raw_hist();
		const cmplx *buf = input_hist();
		for (int i = raw_len - buffer_len - filter_len + 1; i < raw_len; ++i)
			buf = input_hist(ddc(raw[i]));
		return buf;
	}
	int feed(value sample)
	{
		events = 0;
		const value *raw = raw_hist(sample);
		if (countdown) {
			const cmplx *buf = input_hist(ddc(sample));
			if (!--countdown)
				process(buf);
			return events;
		}
		if (keep_baseband)
			input_hist(ddc(sample));
		++sample_count;
		// the detection rests while waiting for the preamble and during a frame, so it keeps what it found
		if (pending) {
			if (!--pending) {
				const cmplx *buf = catch_up();
				if (raw_sync.refine(buf, sync_start))
					preamble(buf);
			}
		} else if (raw_sync(raw + raw_len - 2 * symbol_len)) {
			pending = sync_delay;
		}
		return events;
	}
	// also while waiting for the preamble after a sync
	bool busy()
	{
		return countdown > 0 || pending > 0;
	}
};

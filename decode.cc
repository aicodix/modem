/*
OFDM modem decoder

Copyright 2021 Ahmet Inan <inan@aicodix.de>
*/

#include "common.hh"
#include "schmidl_cox.hh"
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
	static const int tone_off = - tone_count / 2;
	DSP::ReadPCM<value> *pcm;
	DSP::FastFourierTransform<symbol_len, cmplx, -1> fwd;
	DSP::Phasor<cmplx> mix;
	DSP::FIR<cmplx, value, filter_len> lowpass;
	DSP::BipBuffer<cmplx, buffer_len> input_hist;
	DSP::TheilSenEstimator<value, tone_count> tse;
	SchmidlCox<value, cmplx, search_pos, symbol_len, guard_len> correlator;
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
	value cfo_rad, sfo_rad;
	int symbol_pos;
	int crc_bits;

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
	static void base40_decoder(char *str, int64_t val, int len)
	{
		for (int i = len-1; i >= 0; --i, val /= 40)
			str[i] = "   /0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"[val%40];
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
		int code_bits = 1 << code_order, flips = 0;
		for (int i = 0; i < code_bits; ++i)
			flips += code[i] * reenc[i] <= 0;
		std::cerr << "bit flips: " << flips << " of " << code_bits << std::endl;
		shuffle(perm, reenc, code_order, true);
		int count = 0;
		for (int j = 1, k = 0, n = 0; j <= symbol_count; ++j) {
			seed_off = (block_skew * j + first_seed) % block_length;
			value sp = 0, np = 0;
			for (int i = 0; i < tone_count; ++i) {
				if (i % block_length != seed_off) {
					int bits = bits_at(k);
					cmplx hard = map_bits(perm+k, bits);
					k += bits;
					sp += norm(hard);
					np += norm(values[n++] - hard);
				}
			}
			if (!lost[j])
				snr[count++] = sp / std::max(np, value(1e-30) * sp);
		}
		if (!count)
			return;
		DSP::quick_sort(snr, count);
		std::cerr << "real Es/N0 (dB): " << DSP::decibel(snr[0]) << " .. " << DSP::decibel(snr[(count-1)/2]) << " .. " << DSP::decibel(snr[count-1]) << std::endl;
	}
	const cmplx *next_sample()
	{
		value tmp;
		pcm->read(&tmp, 1);
		if (!pcm->good())
			tmp = 0;
		return input_hist(lowpass(tmp * mix()));
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
	Decoder(DSP::ReadPCM<value> *pcm, const char *const *output_names, int output_count) : pcm(pcm), lowpass(DSP::LowPass<value>(value(1275) / value(rate)), DSP::Kaiser<value>(2)), correlator(mls0_seq())
	{
		mix.omega(-1500, rate);
		DSP::Phasor<cmplx> osc;
		const cmplx *buf;
		int output_index = 0;
		int sample_count = 0;
		while (output_index < output_count) {
			do {
				buf = next_sample();
				if (!pcm->good())
					return;
				++sample_count;
			} while (!correlator(buf));

			symbol_pos = correlator.symbol_pos;
			cfo_rad = correlator.cfo_rad;
			std::cerr << "symbol pos: " << sample_count - buffer_len + symbol_pos << std::endl;
			std::cerr << "coarse cfo: " << cfo_rad * (rate / Const::TwoPi()) << " Hz " << std::endl;

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
			std::cerr << "coarse sfo: " << -1000000 * tse.slope() / Const::TwoPi() << " ppm" << std::endl;
			std::cerr << "residual cfo: " << 1000 * tse.yint() * rate / (Const::TwoPi() * symbol_len) << " mHz" << std::endl;
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
			CODE::MLS seq1(mls1_poly);
			auto clamp = [](int v){ return v < -127 ? -127 : v > 127 ? 127 : v; };
			mod_bits = 1;
			oper_mode = -1;
			symbol_count = 0;
			int snr_count = 0;
			bool analog_mode = false;
			for (int j = 0, k = 0; j < symbol_count + 1; ++j) {
				seed_off = (block_skew * j + first_seed) % block_length;
				if (j) {
					for (int i = 0; i < extended_len; ++i)
						correlator(buf = next_sample());
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
					std::cerr << "seed value damaged" << std::endl;
					oper_mode = -1;
					break;
				}
				if (erased) {
					std::cerr << "seed value damaged, erasing symbol " << j << std::endl;
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
							values[k++] = demod[i];
				} else {
					value sp = 0, np = 0, mg = 0;
					for (int i = 0, l = k; !erased && i < tone_count; ++i) {
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
						// the precision of each tone follows its channel power
						gain[i] = norm(chan[i]);
						mg += gain[i];
					}
					// erased symbols get zero LLRs
					value precision = 0;
					if (!erased) {
						snr[snr_count++] = precision = sp / np;
						// relative to the average channel power of this symbol
						precision *= value(tone_count) / std::max(mg, value(1e-30));
					}
					// more resolution after rounding, quantize() caps the slope
					precision *= value(300);
					for (int i = 0, n = (j - 1) * data_tones; i < tone_count; ++i) {
						if (i % block_length != seed_off) {
							if (j)
								values[n++] = demod[i];
							int bits = bits_at(k);
							demap_soft(perm+k, demod[i], precision * gain[i], bits);
							k += bits;
						}
					}
				}
				if (!j) {
					int64_t meta_info = meta_data();
					if (meta_info < 0) {
						std::cerr << "preamble decoding error." << std::endl;
						break;
					}
					int64_t call = meta_info >> 8;
					if (call == 0 || call >= 262144000000000L) {
						std::cerr << "call sign unsupported." << std::endl;
						break;
					}
					char call_sign[10];
					base40_decoder(call_sign, call, 9);
					call_sign[9] = 0;
					std::cerr << "call sign: " << call_sign << std::endl;
					int mode = meta_info & 255;
					bool known = setup(mode);
					print_setup(mode);
					if (!known)
						break;
					analog_mode = mode & 128;
					k = 0;
					for (int i = 0; i < symbol_pos+symbol_len+extended_len; ++i)
						correlator(buf = next_sample());
					std::cerr << "oper mode: " << oper_mode << std::endl;
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
				}
			}
			if (oper_mode < 0)
				continue;
			DSP::quick_sort(snr, snr_count);
			std::cerr << "Es/N0 (dB): " << DSP::decibel(snr[0]) << " .. " << DSP::decibel(snr[(snr_count-1)/2]) << " .. " << DSP::decibel(snr[snr_count-1]) << std::endl;
			if (analog_mode) {
				const char *output_name = output_names[output_index++];
				if (output_count == 1 && output_name[0] == '-' && output_name[1] == 0)
					output_name = "/dev/stdout";
				DSP::WriteWAV<value> output_file(output_name, analog_rate, 32, 2);
				if (!output_file.good()) {
					std::cerr << "Couldn't open file \"" << output_name << "\" for writing." << std::endl;
					continue;
				}
				output_file.write(reinterpret_cast<value *>(values), symbol_count * data_tones, 2);
				continue;
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
				std::cerr << "payload decoding error." << std::endl;
				continue;
			}
			for (int i = 0; i < data_bits; ++i)
				CODE::set_le_bit(data, i, mesg[i].v[best] < 0);
			reencode(best);

			const char *output_name = output_names[output_index++];
			if (output_count == 1 && output_name[0] == '-' && output_name[1] == 0)
				output_name = "/dev/stdout";
			std::ofstream output_file(output_name, std::ios::binary | std::ios::trunc);
			if (output_file.bad()) {
				std::cerr << "Couldn't open file \"" << output_name << "\" for writing." << std::endl;
				continue;
			}
			CODE::Xorshift32 scrambler;
			for (int i = 0; i < data_bytes; ++i)
				data[i] ^= scrambler();
			for (int i = 0; i < data_bytes; ++i)
				output_file.put(data[i]);
		}
	}
};

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
		delete new Decoder<value, cmplx, 44100>(&input_file, argv+2, output_count);
		break;
	case 48000:
		delete new Decoder<value, cmplx, 48000>(&input_file, argv+2, output_count);
		break;
	default:
		std::cerr << "Unsupported sample rate." << std::endl;
		return 1;
	}
	return 0;
}


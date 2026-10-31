/*
Quadrature amplitude modulation

Copyright 2018 Ahmet Inan <xdsopl@gmail.com>
*/

#pragma once

#include "quantize.hh"

template <int BITS, typename TYPE, typename CODE>
struct QuadratureAmplitudeModulation
{
	static_assert(BITS >= 4 && BITS <= 12 && BITS % 2 == 0, "unsupported number of bits");
	static const int NUM = 1 << BITS;
	typedef TYPE complex_type;
	typedef typename TYPE::value_type value_type;
	typedef CODE code_type;

	// sqrt(2 * (NUM - 1) / 3) / (sqrt(NUM) - 1)
	static constexpr value_type FAC =
		BITS == 4 ? 1.0540925533894596 :
		BITS == 6 ? 0.9258200997725516 :
		BITS == 8 ? 0.8692269873603529 :
		BITS == 10 ? 0.8424235391742344 :
		0.8293555858802010;
	static constexpr value_type RCP = ((1 << BITS / 2) - 1) * FAC;
	static constexpr value_type AMP = 1 / RCP;
	static constexpr value_type DIST = 2 * AMP;

	static constexpr value_type amp(int i)
	{
		return AMP * i;
	}

	static void hard(code_type *b, complex_type c)
	{
		value_type re = c.real(), im = c.imag();
		b[0] = re < amp(0) ? code_type(-1) : code_type(1);
		b[1] = im < amp(0) ? code_type(-1) : code_type(1);
		for (int i = 2; i < BITS; i += 2) {
			re = std::abs(re) - amp(1 << ((BITS - i) / 2));
			im = std::abs(im) - amp(1 << ((BITS - i) / 2));
			b[i] = re < amp(0) ? code_type(-1) : code_type(1);
			b[i+1] = im < amp(0) ? code_type(-1) : code_type(1);
		}
	}

	static void soft(code_type *b, complex_type c, value_type precision)
	{
		value_type slope = DIST * precision;
		value_type re = c.real(), im = c.imag();
		b[0] = quantize<code_type>(slope, re);
		b[1] = quantize<code_type>(slope, im);
		for (int i = 2; i < BITS; i += 2) {
			re = std::abs(re) - amp(1 << ((BITS - i) / 2));
			im = std::abs(im) - amp(1 << ((BITS - i) / 2));
			b[i] = quantize<code_type>(slope, re);
			b[i+1] = quantize<code_type>(slope, im);
		}
	}

	static complex_type map(code_type *b)
	{
		value_type re = b[BITS-2], im = b[BITS-1];
		for (int i = BITS - 4; i >= 0; i -= 2) {
			re = b[i] * (re + value_type(1 << ((BITS - i) / 2 - 1)));
			im = b[i+1] * (im + value_type(1 << ((BITS - i) / 2 - 1)));
		}
		return AMP * complex_type(re, im);
	}
};

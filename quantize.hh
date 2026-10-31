/*
Soft decision quantization

Copyright 2026 Ahmet Inan <inan@aicodix.de>
*/

#pragma once

#include <cmath>
#include <cstdint>
#include <algorithm>
#include <type_traits>

// the slope is capped, leaving two bits of headroom for |value| < 2
template <typename CODE, typename VALUE>
CODE quantize(VALUE slope, VALUE value)
{
	if constexpr (std::is_same<CODE, int8_t>::value) {
		value = std::nearbyint(value * std::min<VALUE>(slope, 16));
		return std::min<VALUE>(std::max<VALUE>(value, -127), 127);
	} else if constexpr (std::is_same<CODE, int16_t>::value) {
		value = std::nearbyint(value * std::min<VALUE>(slope, 4096));
		return std::min<VALUE>(std::max<VALUE>(value, -32767), 32767);
	} else if constexpr (std::is_integral<CODE>::value) {
		return std::nearbyint(value * slope);
	} else {
		return value * slope;
	}
}

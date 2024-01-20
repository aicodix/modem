/*
Table generator for frozen bits

Copyright 2021 Ahmet Inan <inan@aicodix.de>
*/

#include <cmath>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <functional>
#include "polar_freezer.hh"

template <int M>
void code(int N, int K)
{
	long double erasure_probability = (long double)(N - K) / N;
	double design_SNR = 10 * std::log10(-std::log(erasure_probability));
	std::cerr << "design SNR: " << design_SNR << std::endl;
	auto freeze = new CODE::PolarCodeConst0<M>;
	double better_SNR = design_SNR + 1.59175;
	std::cerr << "better SNR: " << better_SNR << std::endl;
	long double better_probability = std::exp(-pow(10.0, better_SNR / 10));
	auto frozen = new uint32_t[1<<(M-5)];
	(*freeze)(frozen, M, K+(1<<M)-N, better_probability);
	delete freeze;

	std::cerr << "Polar(" << N << ", " << K << ")" << std::endl;
	std::cout << "static const uint32_t frozen_" << std::dec << N << "_" << K << "[" << (1<<(M-5)) << "] = { " << std::hex;
	for (int i = 0; i < 1<<(M-5); ++i)
		std::cout << "0x" << frozen[i] << ", ";
	std::cout << "};" << std::endl;
}

int main()
{
	code<12>(4096, 2048+32+67);
	code<13>(8192, 4096+32+133);
	code<14>(16384, 8192+32+265);
	return 0;
}

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
void code(int N, int K, double P)
{
	auto freeze = new CODE::PolarCodeConst0<M>;
	auto frozen = new uint32_t[1<<(M-5)];
	(*freeze)(frozen, M, K+(1<<M)-N, P);
	delete freeze;

	std::cout << "static const uint32_t frozen_" << std::dec << N << "_" << K << "[" << (1<<(M-5)) << "] = { " << std::hex;
	for (int i = 0; i < 1<<(M-5); ++i)
		std::cout << "0x" << frozen[i] << ", ";
	std::cout << "};" << std::endl;
}

int main()
{
	code<16>(65536, 32768+32, 0.42);
	code<15>(32768, 16384+32, 0.40);
	code<14>(16384, 8192+32, 0.39);
	code<13>(8192, 4096+32, 0.38);
	code<12>(4096, 2048+32, 0.37);
	return 0;
}

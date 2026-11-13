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
	auto freeze = new CODE::BhattacharyyaFreezer<M>;
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
	// call sign and mode with 16 bit CRC
	code<8>(256, 48+8+16, 0.5);
	return 0;
}

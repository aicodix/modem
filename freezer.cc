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
	// call sign and mode with 16 bit CRC
	code<8>(256, 47+8+16, 0.5);
	// 1/2-rate payload with 32 bit CRC
	code<11>(2048, 1024+32, 0.33);
	code<12>(4096, 2048+32, 0.34);
	code<13>(8192, 4096+32, 0.36);
	code<14>(16384, 8192+32, 0.39);
	code<15>(32768, 16384+32, 0.40);
	code<16>(65536, 32768+32, 0.42);
	// 2/3-rate payload with 32 bit CRC
	code<11>(2048, 1368+32, 0.19);
	code<12>(4096, 2736+32, 0.20);
	code<13>(8192, 5472+32, 0.21);
	code<14>(16384, 10944+32, 0.23);
	code<15>(32768, 21888+32, 0.24);
	code<16>(65536, 43776+32, 0.26);
	// 3/4-rate payload with 32 bit CRC
	code<11>(2048, 1536+32, 0.11);
	code<12>(4096, 3072+32, 0.12);
	code<13>(8192, 6144+32, 0.15);
	code<14>(16384, 12288+32, 0.15);
	code<15>(32768, 24576+32, 0.17);
	code<16>(65536, 49152+32, 0.19);
	// 5/6-rate payload with 32 bit CRC
	code<11>(2048, 1704+32, 0.07);
	code<12>(4096, 3408+32, 0.08);
	code<13>(8192, 6816+32, 0.09);
	code<14>(16384, 13632+32, 0.11);
	code<15>(32768, 27264+32, 0.11);
	code<16>(65536, 54528+32, 0.11);
	return 0;
}

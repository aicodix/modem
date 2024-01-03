
CXXFLAGS = -std=c++17 -W -Wall -Ofast -fno-exceptions -fno-rtti -I../dsp -I../code
CXX = clang++ -stdlib=libc++ -march=native
#CXX = g++ -march=native

#CXX = armv7a-hardfloat-linux-gnueabi-g++ -static -mfpu=neon -march=armv7-a
#QEMU = qemu-arm

#CXX = aarch64-unknown-linux-gnu-g++ -static -march=armv8-a+crc+simd -mtune=cortex-a72
#QEMU = qemu-aarch64

.PHONY: all

all: encode decode

test: encode decode
	$(QEMU) ./encode encoded.wav 8000 8 1 /dev/urandom 25
	$(QEMU) ./decode /dev/null encoded.wav

encode: encode.cc
	$(CXX) $(CXXFLAGS) $< -o $@

decode: decode.cc
	$(CXX) $(CXXFLAGS) $< -o $@

freezer: freezer.cc
	$(CXX) $(CXXFLAGS) $< -o $@

.PHONY: clean

clean:
	rm -f encode decode freezer


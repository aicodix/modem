
CXXFLAGS = -std=c++17 -W -Wall -O3 -ffast-math -fno-exceptions -fno-rtti -I../dsp -I../code
#CXXFLAGS += -g -fsanitize=address
CXX = clang++ -stdlib=libc++ -march=native
#CXX = g++ -march=native

#CXX = armv7a-hardfloat-linux-gnueabi-g++ -static -mfpu=neon -march=armv7-a
#QEMU = qemu-arm

#CXX = aarch64-unknown-linux-gnu-g++ -static -march=armv8-a+crc+simd -mtune=cortex-a72
#QEMU = qemu-aarch64

.PHONY: all

all: encode decode

test: encode decode
	$(QEMU) ./encode audio.wav 48000 16 1 1500 ANONYMOUS 48 /dev/urandom
	$(QEMU) ./decode audio.wav /dev/null

encode: encode.cc common.hh
	$(CXX) $(CXXFLAGS) $< -o $@

decode: decode.cc common.hh schmidl_cox.hh
	$(CXX) $(CXXFLAGS) $< -o $@

freezer: freezer.cc
	$(CXX) $(CXXFLAGS) $< -o $@

.PHONY: clean

clean:
	rm -f encode decode freezer


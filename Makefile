
CXXFLAGS = -std=c++17 -W -Wall -Ofast -fno-exceptions -fno-rtti -I../dsp -I../code
CXX = clang++ -stdlib=libc++ -march=native
#CXX = g++ -march=native

#CXX = armv7a-hardfloat-linux-gnueabi-g++ -static -mfpu=neon -march=armv7-a
#QEMU = qemu-arm

.PHONY: all

all: encode decode

test: encode decode
	$(QEMU) ./encode encoded.wav 8000 8 1 /dev/urandom
	$(QEMU) ./decode /dev/null encoded.wav

encode: encode.cc
	$(CXX) $(CXXFLAGS) $< -o $@

decode: decode.cc
	$(CXX) $(CXXFLAGS) $< -o $@

.PHONY: clean

clean:
	rm -f encode decode


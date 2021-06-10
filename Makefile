
CXXFLAGS = -std=c++17 -W -Wall -Ofast -fno-exceptions -fno-rtti -march=native -I../dsp -I../code
CXX = clang++ -stdlib=libc++
#CXX = g++

.PHONY: all

all: encode decode

encode: encode.cc
	$(CXX) $(CXXFLAGS) $< -o $@

decode: decode.cc
	$(CXX) $(CXXFLAGS) $< -o $@

.PHONY: clean

clean:
	rm -f encode decode


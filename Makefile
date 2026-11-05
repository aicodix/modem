
CXXFLAGS = -std=c++17 -W -Wall -O3 -ffast-math -fno-exceptions -fno-rtti -I../dsp -I../code
#CXXFLAGS += -g -fsanitize=address
CXX = clang++ -stdlib=libc++ -march=native
#CXX = g++ -march=native

#CXX = armv7a-hardfloat-linux-gnueabi-g++ -static -mfpu=neon -march=armv7-a
#QEMU = qemu-arm

#CXX = aarch64-unknown-linux-gnu-g++ -static -march=armv8-a+crc+simd -mtune=cortex-a72
#QEMU = qemu-aarch64

WASMCXX = clang++
WASMFLAGS = -std=c++17 -W -Wall -O3 -ffast-math -fno-exceptions -fno-rtti -ffreestanding -fvisibility=hidden --target=wasm32 -Iquirks -I../dsp -I../code
WASMLDFLAGS = -nostdlib -Wl,--export-dynamic,--no-entry,--export=__wasm_call_ctors
QUIRKS = quirks/algorithm quirks/cassert quirks/cmath quirks/cstdint quirks/cstdlib quirks/cstring quirks/functional quirks/type_traits

.PHONY: all

all: encode decode

test: encode decode
	$(QEMU) ./encode audio.wav 48000 32 ANONYMOUS QAM16 1/2 short /dev/urandom
	$(QEMU) ./decode audio.wav /dev/null

encode: encode.cc encoder.hh native.hh common.hh psk.hh qam.hh quantize.hh
	$(CXX) $(CXXFLAGS) $< -o $@

decode: decode.cc decoder.hh native.hh common.hh schmidl_cox.hh psk.hh qam.hh quantize.hh
	$(CXX) $(CXXFLAGS) $< -o $@

.PHONY: wasm

wasm: temp/encoder_44100.wasm temp/encoder_48000.wasm temp/decoder_44100.wasm temp/decoder_48000.wasm

temp/encoder_%.wasm: encoder_module.cc encoder.hh common.hh psk.hh qam.hh quantize.hh $(QUIRKS)
	$(WASMCXX) $(WASMFLAGS) $(WASMLDFLAGS) -DRATE=$* -o $@ $<

temp/decoder_%.wasm: decoder_module.cc decoder.hh common.hh schmidl_cox.hh psk.hh qam.hh quantize.hh $(QUIRKS)
	$(WASMCXX) $(WASMFLAGS) $(WASMLDFLAGS) -DRATE=$* -o $@ $<

.PHONY: check-wasm

check-wasm: temp/encoder_44100.wasm temp/encoder_48000.wasm temp/decoder_44100.wasm temp/decoder_48000.wasm encode decode
	node check_encoder.mjs
	node check_decoder.mjs

freezer: freezer.cc
	$(CXX) $(CXXFLAGS) $< -o $@

polar_tables.hh: freezer
	./$< > $@

.PHONY: clean

clean:
	rm -f encode decode freezer temp/*.wasm


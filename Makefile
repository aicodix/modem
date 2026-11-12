
CXXFLAGS = -std=c++17 -W -Wall -O3 -ffast-math -fno-exceptions -fno-rtti -I../dsp -I../code
#CXXFLAGS += -g -fsanitize=address
CXX = clang++ -stdlib=libc++ -march=native
#CXX = g++ -march=native

#CXX = armv7a-hardfloat-linux-gnueabi-g++ -static -mfpu=neon -march=armv7-a
#QEMU = qemu-arm

#CXX = aarch64-unknown-linux-gnu-g++ -static -march=armv8-a+crc+simd -mtune=cortex-a72
#QEMU = qemu-aarch64

WASMCXX = clang++
WASMFLAGS = -std=c++17 -W -Wall -O3 -ffast-math -fno-exceptions -fno-rtti -ffreestanding -fvisibility=hidden --target=wasm32 -msimd128 -Iquirks -I../dsp -I../code
WASMLDFLAGS = -nostdlib -Wl,--export-dynamic,--no-entry,--export=__wasm_call_ctors
QUIRKS = quirks/algorithm quirks/cassert quirks/cmath quirks/cstdint quirks/cstdlib quirks/cstring quirks/functional quirks/limits quirks/type_traits

.PHONY: all

all: encode decode wasm

test: encode decode
	$(QEMU) ./encode audio.wav 48000 32 ANONYMOUS QAM16 1/2 short /dev/urandom
	$(QEMU) ./decode audio.wav /dev/null

encode: encode.cc encoder.hh native.hh common.hh psk.hh qam.hh quantize.hh
	$(CXX) $(CXXFLAGS) $< -o $@

decode: decode.cc decoder.hh native.hh common.hh schmidl_cox.hh frequency_sync.hh psk.hh qam.hh quantize.hh
	$(CXX) $(CXXFLAGS) $< -o $@

.PHONY: play-modes

play-modes: encode
	for m in BPSK QPSK 8PSK QAM16 QAM64 QAM256 QAM1024 QAM4096; do for r in 1/2 2/3 3/4 5/6; do for f in short normal; do echo "$$m $$r $$f"; echo "$$m $$r $$f" | $(QEMU) ./encode - 48000 16 ANONYMOUS $$m $$r $$f - 2>/dev/null | play -q -; done; done; done

.PHONY: wasm

wasm: temp/encoder_44100.wasm temp/encoder_48000.wasm temp/decoder_44100.wasm temp/decoder_48000.wasm

temp/encoder_%.wasm: encoder_module.cc encoder.hh common.hh psk.hh qam.hh quantize.hh $(QUIRKS)
	$(WASMCXX) $(WASMFLAGS) $(WASMLDFLAGS) -DRATE=$* -o $@ $<

temp/decoder_%.wasm: decoder_module.cc decoder.hh analyzer.hh common.hh schmidl_cox.hh frequency_sync.hh psk.hh qam.hh quantize.hh $(QUIRKS)
	$(WASMCXX) $(WASMFLAGS) $(WASMLDFLAGS) -DRATE=$* -o $@ $<

.PHONY: check-wasm

check-wasm: temp/encoder_44100.wasm temp/encoder_48000.wasm temp/decoder_44100.wasm temp/decoder_48000.wasm encode decode
	node check_encoder.mjs
	node check_decoder.mjs

.PHONY: pages

# a snapshot of the page and its modules for GitHub Pages, served from the docs directory
# the page shows the commit it was built from, "-dirty" if files outside of docs differ from it
pages: index.html manifest.webmanifest sw.js icon.svg icon-192.png icon-512.png temp/encoder_44100.wasm temp/encoder_48000.wasm temp/decoder_44100.wasm temp/decoder_48000.wasm
	mkdir -p docs/temp
	hash=$$(git rev-parse --short HEAD) && \
	version=$$hash$$(git diff --quiet HEAD -- . ':!docs' || echo -dirty) && \
	grep -q 'tree/next">next<' index.html && \
	sed 's#tree/next">next<#tree/'$$hash'">'$$version'<#' index.html > docs/index.html
	cp manifest.webmanifest sw.js icon.svg icon-192.png icon-512.png docs/
	cp temp/encoder_44100.wasm temp/encoder_48000.wasm temp/decoder_44100.wasm temp/decoder_48000.wasm docs/temp/
	touch docs/.nojekyll

.PHONY: webserver

webserver: wasm
	python3 -m http.server 8080 --bind 127.0.0.1

freezer: freezer.cc
	$(CXX) $(CXXFLAGS) $< -o $@

polar_tables.hh: freezer
	./$< > $@

.PHONY: clean

clean:
	rm -f encode decode freezer temp/*.wasm


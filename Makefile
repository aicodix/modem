
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
	head -c 1024 /dev/urandom | $(QEMU) ./encode audio.wav 48000 32 ANONYMOUS 5 -
	$(QEMU) ./decode audio.wav /dev/null

encode: encode.cc encoder.hh native.hh common.hh psk.hh qam.hh quantize.hh
	$(CXX) $(CXXFLAGS) $< -o $@

decode: decode.cc decoder.hh native.hh common.hh raw_sync.hh psk.hh qam.hh quantize.hh
	$(CXX) $(CXXFLAGS) $< -o $@

.PHONY: play-modes

# a frame of each number of payload symbols, at the fewest seconds that allow it, with the payload of the most and of the
# least robust mode at that many symbols, 12 and 240 bytes per symbol, so the encoder can't pick fewer symbols
play-modes: encode
	for p in "1 4" "2 8" "3 16" "5 32" "10 64" "18 128" "36 256"; do set -- $$p; for b in 12 240; do n=$$(($$2 * $$b)); echo "$$1 s $$n B"; (echo "$$1 s $$n B"; head -c $$n /dev/zero) | head -c $$n | $(QEMU) ./encode - 48000 16 ANONYMOUS $$1 - 2>/dev/null | play -q -; done; done

.PHONY: wasm

wasm: temp/encoder_44100.wasm temp/encoder_48000.wasm temp/decoder_44100.wasm temp/decoder_48000.wasm

temp/encoder_%.wasm: encoder_module.cc encoder.hh common.hh psk.hh qam.hh quantize.hh $(QUIRKS)
	$(WASMCXX) $(WASMFLAGS) $(WASMLDFLAGS) -DRATE=$* -o $@ $<

temp/decoder_%.wasm: decoder_module.cc decoder.hh analyzer.hh common.hh raw_sync.hh psk.hh qam.hh quantize.hh $(QUIRKS)
	$(WASMCXX) $(WASMFLAGS) $(WASMLDFLAGS) -DRATE=$* -o $@ $<

.PHONY: check-wasm

check-wasm: temp/encoder_44100.wasm temp/encoder_48000.wasm temp/decoder_44100.wasm temp/decoder_48000.wasm encode decode
	node check_encoder.mjs
	node check_decoder.mjs
	node check_frozen.mjs

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

nested: ../code/tools/nested.cc
	$(CXX) $(CXXFLAGS) $< -o $@

polar_tables.hh: freezer nested
	./freezer > $@
	for n in 10 11 12 13 14 15 16; do ./nested $$n 4 32 >> $@; done

.PHONY: clean

clean:
	rm -f encode decode freezer nested temp/*.wasm


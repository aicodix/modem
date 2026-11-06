# CLAUDE.md

C++17 OFDM modem: a command-line encoder that turns data into a WAV audio signal
and a decoder that recovers it, with polar codes and PSK/QAM modulations in a
2400 Hz bandwidth. See README.md for usage, modes and simulation.

## Layout

- `encode.cc`, `decode.cc`: the two programs, native drivers around the
  `Encoder` struct in `encoder.hh` and the `Decoder` struct in `decoder.hh`,
  both templated on sample rate.
- `native.hh`: what only the native programs need, such as iostream, WAV
  file I/O and `print_setup()`.
- `encoder_module.cc`, `decoder_module.cc`, `quirks/`: the WebAssembly
  modules and the stand-ins for standard headers they are built with (see
  "WebAssembly module").
- `check_encoder.mjs`, `check_decoder.mjs`, `index.html`: the Node.js checks
  of the modules and the browser page that sends text and receives it from
  the microphone or a WAV file. For analog frames, the page takes and gives
  their complex values as stereo WAV files, like the native programs.
- `common.hh`: what both sides share, such as tone layout, sequences, the
  interleaver (`shuffle()`), bits per tone (`bits_at()`) and the operation
  mode table in `Common::setup()`.
- `psk.hh`, `qam.hh`: constellation mapping, hard and soft demapping, templated
  on bits per symbol (`PhaseShiftKeying<1..3>`, `QuadratureAmplitudeModulation<4..12>`).
- `quantize.hh`: turns soft values into the code type of the LLRs.
- `schmidl_cox.hh`: synchronization (timing, CFO, SFO).
- `analyzer.hh`: spectrum, spectrogram and X-Y oscilloscope of the decoder's
  baseband for the page, fed by the decoder module and drawn as RGBA pixels.
  It takes FFTs as long as the decoder's symbols, so the bins are as wide as
  the tones, and draws the middle 400 bins, from 0 to 3000 Hz. Spectrum and
  scope are off until the page turns them on, and take no time while off.
- `freezer.cc`: generates `polar_tables.hh`. Don't edit the tables by hand;
  regenerate them with `make polar_tables.hh`.

Most building blocks come from two sibling checkouts found by `-I../dsp -I../code`:
[dsp](https://github.com/aicodix/dsp) (namespace `DSP`) and
[code](https://github.com/aicodix/code) (namespace `CODE`). The simulation tools
come from [disorders](https://github.com/aicodix/disorders), also expected next
to this repository. When a change belongs in one of those libraries, say so
instead of working around it here.

## Building and testing

```
make
make test
```

`make test` encodes random data to `audio.wav` and decodes it. There is no test
suite; check a change by encoding, decoding and comparing, for example with the
simulation pipeline in the README ("Simulating"), which adds multipath, CFO,
SFO and AWGN with the disorders tools and ends with `diff -s`.
The decoder reports the operation mode and Es/N0 statistics on stderr. After a
frame decodes, it re-encodes the message and also reports the number of code
bits whose soft value had the wrong sign (erased bits count too) and the real
Es/N0 against the re-encoded symbols, which is exact also near threshold.

The default compiler is `clang++ -stdlib=libc++ -march=native`. If libc++ is
not installed, build with `make CXX="g++ -march=native"`. The Makefile has
commented-out lines for `g++` and for ARM cross-compilation under QEMU.
Run `make clean` to remove the binaries.
When a header joins the build, add it to the dependencies of `encode` and
`decode` in the Makefile, but not `polar_tables.hh`, which has its own rule.

`make wasm` builds the encoder and decoder modules `temp/encoder_44100.wasm`,
`temp/encoder_48000.wasm`, `temp/decoder_44100.wasm` and
`temp/decoder_48000.wasm` with clang for `wasm32`. `make check-wasm` runs
`check_encoder.mjs` and `check_decoder.mjs` with Node.js: they encode random
data in all modes and both rates with the encoder modules and decode it with
the native `decode`, and the other way around with the native `encode` and
the decoder modules. To try `index.html`, serve the repository with any
static web server that sends `.wasm` as `application/wasm`, for example
`python3 -m http.server`, and open `index.html`.

Everything is compiled with `-std=c++17 -O3 -ffast-math -fno-exceptions -fno-rtti`,
so code must:

- not use exceptions or RTTI;
- stay correct under `-ffast-math`. Don't rely on signed zeros, NaN or infinity
  checks, or a strict order of floating-point operations. Guard divisions that
  can approach zero, such as noise estimates on a clean signal.

## Code style

- Every file starts with a short C comment block: a one-line description,
  a blank line, then `Copyright <year> Ahmet Inan <inan@aicodix.de>`.
- Headers use `#pragma once`. Modem code lives in the global namespace and
  refers to the libraries as `DSP::` and `CODE::`.
- Indent with tabs and use K&R braces. Keep the code compact and avoid
  unnecessary abstraction.
- Sizes are fixed at compile time with fixed-size member arrays. Avoid heap
  allocation in the signal path.

## Conventions worth knowing

- Encoder and decoder must stay in step. A change to anything on the air, such
  as the mode table, preamble, sequences, scrambler, frozen bits or constellation
  mapping, changes both sides, and the README tables when durations, payloads
  or bitrates change.
- Soft demapping produces int16 LLRs for a min-sum polar list decoder, which
  ignores their overall scale. What matters is resolution after rounding and
  headroom against saturation. `quantize()` caps the slope `DIST * precision`
  at 4096, so every constellation gets the same headroom, and `decode.cc`
  scales the precision by 300 for resolution.
- Each tone's precision is weighted by its channel power, relative to the
  average of its symbol. The logged Es/N0 stays unweighted.
- The precision comes from hard decisions and overestimates Es/N0 by several dB
  for dense QAM near threshold. Within a frame that is a common scale and does
  not hurt decoding, but the logged Es/N0 is not exact there. The real Es/N0
  after decoding is.
- The signal is always a real (one channel) audio signal centered at 1500 Hz.
  The decoder mixes it down by 1500 Hz and lowpass filters it, flat to ±1050 Hz
  and at least 65 dB down from ±1500 Hz, which removes DC and out-of-band noise
  before synchronization. This limits the CFO to about ±250 Hz.
- Below roughly 7 dB in-band SNR the preamble no longer synchronizes, so BPSK and
  QPSK at low code rates fail at sync before the code is tested. To compare
  decoder changes, use denser modulations or higher code rates near their
  threshold, many noise seeds per point, and report decoded frames out of trials.

## Workflow

- Work on the `next` branch. Keep changes to this modem minimal; larger
  redesigns are for its successor.
- Commit messages are a single short lowercase line, such as
  "moved common header includes to common.hh". There is no body and
  no trailers: no `Co-Authored-By`, no session links.
- Don't push. Deliver changes as a `git format-patch` patch for the
  maintainer to review and apply.
- When changing the encoder, check that its output is unchanged or as
  intended. When changing the decoder, check that decoding still succeeds and
  compare old and new versions over many noise seeds near the threshold of the
  affected modes. A refactoring should be bit-exact: compare the outputs of
  old and new code over many inputs, including decision boundaries.
  Say how you checked it.

## WebAssembly module

The goal is to also build the modem as a WebAssembly module for the browser,
the same way as [HYDRA](https://github.com/xdsopl/hydra) and the old
WebAssembly examples on aicodix.de, while the native `encode` and `decode`
programs keep working and stay the reference:

- There are four modules: an encoder and a decoder, each for 44100 Hz and
  for 48000 Hz, built from the same sources with the rate fixed at compile
  time. The rate is a template parameter, so buffer sizes, FFT and filter
  lengths stay compile-time constants, and a page that only sends never
  loads the decoder. The page loads the module matching
  `AudioContext.sampleRate`, or asks for a rate with
  `new AudioContext({ sampleRate: 48000 })` and lets the browser resample.
- Compile straight to `--target=wasm32` with clang. No Emscripten and no libc:
  `-ffreestanding -fvisibility=hidden -fno-exceptions -fno-rtti`, linked with
  `-nostdlib -Wl,--export-dynamic,--no-entry`. Only what is marked
  `__attribute__((visibility("default")))` is exported, wrapped in `extern "C"`.
- Don't link with `--allow-undefined`. A missing function must fail the link,
  not become a wasm import that only fails when the browser loads the module.
- The encoder and decoder build filters, windows and sequences in their
  constructors, so the module needs its static constructors run once before
  first use: export `__wasm_call_ctors` and call it from JavaScript after
  instantiating.
- No heap. The interface is fixed-size exported buffers, such as the payload,
  the audio samples and a status record, which JavaScript reads and writes
  through typed arrays on `memory.buffer`. Exported functions take and return
  plain integers or floats, and return -1 on bad arguments or failure.
- The modem, `dsp` and `code` include standard headers and call `std::abs`,
  `std::min`, `std::nearbyint` and the like. For the module, a `quirks/`
  directory found first with `-Iquirks` holds minimal stand-ins for these
  headers (`cmath`, `algorithm`, `functional`, `type_traits`, `cstdint`,
  ...). They define only what the code uses, in namespace `std`, and map math
  to the `__builtin_*` functions that are single wasm instructions, such as
  `sqrt`, `abs`, `min`, `max`, `nearbyint` and `copysign`, like `quirks.hh`
  in the old examples. `sin` and `cos` are declared there but not defined,
  so a use in a module fails at link time. Small libc pieces like `memset`
  and `memcpy` go there too when the compiler needs them.
- Everything else that would come from libm comes from `dsp`, which already
  avoids it: `decibel()` needs only `sqrt`, `UnitCircle` has its own
  polynomials, and `common.hh` includes `sincos.hh` and `atan2.hh` before
  `complex.hh`, so `Phasor::omega(value)`, `polar()` and `arg()` use them,
  natively as well as in the modules. Use `UnitCircle` for angles that are
  fractions of a circle, as the encoder's guard interval weights do. When
  a module needs another such function, add it to `dsp`, not here.
- The modules are built with `-msimd128`, so `SIMD<>` from `code` uses its
  `simd128*.hh` backend instead of the generic loops. That roughly halves the
  decoding time, most of which goes to the polar list decoder.
- The signal path must not depend on the host: keep iostream, `std::cerr`
  logging, `wav.hh`/`pcm.hh` file I/O and argument parsing in the native
  programs. The module gets samples from buffers and reports through return
  values and the status buffer. The encoder already works this way: its
  native driver in `encode.cc` reads the input, calls `start()`, `begin()`,
  `next()` once per symbol and `stop()`, writes each symbol's samples from
  `audio` and does the logging, including `print_setup()` and the PAPR line.
- The decoder is a streaming state machine: `Decoder::feed()` takes one
  sample and returns a bit mask of events, such as a sync, the metadata, a
  damaged seed or a decoded frame, with their values in members. The native
  driver in `decode.cc` prints them in `report()` and writes the output
  files. The module's `feed()` stops after each sample with events, so
  JavaScript reads them with their values before the next sample can change
  them. On the page, audio comes from an `AudioWorklet`, not the deprecated
  `ScriptProcessorNode`.
- Every step must stay testable natively: the same sources build for the host
  with sanitizers. A change on the way to the module doesn't have to be bit
  exact, since `-ffast-math` and `UnitCircle` instead of `std::cos` move
  samples by rounding errors. It is enough that the digital data decodes
  back unchanged, over all modes and both rates. Node.js can run the modules
  without a browser: a script writes their samples to a WAV file for the
  native `decode`, and the other way around.
- A plain HTML page with JavaScript, without a framework, loads the module with
  `WebAssembly.instantiateStreaming()` and runs the tests and examples in the
  browser.
- Work in small, reviewable steps, one logical change per patch. Before
  editing, say in one or two sentences what will change, and ask first if a
  step would touch more than about three files or the scope is unclear.

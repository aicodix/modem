# CLAUDE.md

C++17 OFDM modem: a command-line encoder that turns data into a WAV audio signal
and a decoder that recovers it, with polar codes and PSK/QAM modulations in a
2400 Hz bandwidth. See README.md for usage, modes and simulation.

## Layout

- `encode.cc`, `decode.cc`: the two programs, each one `Encoder` or `Decoder`
  struct templated on sample rate.
- `common.hh`: what both sides share, such as tone layout, sequences and the
  operation mode table in `Common::setup()`.
- `psk.hh`, `qam.hh`: constellation mapping, hard and soft demapping.
- `schmidl_cox.hh`: synchronization (timing, CFO, SFO).
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
The decoder reports the operation mode and Es/N0 statistics on stderr.

The default compiler is `clang++ -stdlib=libc++ -march=native`. If libc++ is
not installed, build with `make CXX="g++ -march=native"`. The Makefile has
commented-out lines for `g++` and for ARM cross-compilation under QEMU.
Run `make clean` to remove the binaries.

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
  headroom against saturation.
- Below roughly 7 dB in-band SNR the preamble no longer synchronizes, so BPSK and
  QPSK at low code rates fail at sync before the code is tested. To compare
  decoder changes, use denser modulations or higher code rates near their
  threshold, many noise seeds per point, and report decoded frames out of trials.

## Workflow

- Work on the `next` branch.
- Commit messages are a single short lowercase line, such as
  "moved common header includes to common.hh". There is no body and
  no trailers: no `Co-Authored-By`, no session links.
- Don't push. Deliver changes as a `git format-patch` patch for the
  maintainer to review and apply.
- When changing the encoder, check that its output is unchanged or as
  intended. When changing the decoder, check that decoding still succeeds and
  compare old and new versions over many noise seeds near the threshold of the
  affected modes. Say how you checked it.

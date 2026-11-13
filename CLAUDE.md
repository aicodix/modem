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
  their complex values as stereo WAV files, like the native programs. Its
  meters show the CFO and SFO the decoder tracks after each symbol, from the
  module's `tracking` buffer, on symmetric logarithmic scales, dimmed between
  frames, and the RMS and peak level of the microphone, which the page
  measures itself. Like the views of the analyzer, the meters are off until
  the page turns them on and take no time while off, and the page only
  redraws while at least one of them is on.
- `common.hh`: what both sides share, such as tone layout, sequences, the
  interleaver (`shuffle()`), bits per tone (`bits_at()`) and the operation
  mode table in `Common::setup()`.
- `psk.hh`, `qam.hh`: constellation mapping, hard and soft demapping, templated
  on bits per symbol (`PhaseShiftKeying<1..3>`, `QuadratureAmplitudeModulation<4..12>`).
- `quantize.hh`: turns soft values into the code type of the LLRs.
- `raw_sync.hh`: sync detection in the frequency domain on the raw input,
  before mixing and filtering, so while waiting for a sync the decoder only
  keeps the raw samples and nothing else runs. Every `symbol_len` samples it
  takes a Hann window over the last two symbols and a real FFT twice as long
  as a symbol. Within the two identical sync symbols the tones only fall on
  even bins, and the window only spreads a bin to its neighbours, so it keeps
  the tones apart. It also hides the step where the end of the window wraps
  around to its start, which a misaligned window, a CFO between the bins and
  SFO leave, and which would otherwise leak every tone into its neighbours
  (see the comment in `RawSync`). Mixing
  is a shift of the bins: it sums the ratio of adjacent tones times the
  conjugate of their known ratio, against the number of tone pairs, for every
  shift by half a tone up to 250 Hz, which finds the CFO within a quarter of
  a tone. Dividing instead of multiplying keeps strong tones, as from a
  carrier or speech, from dominating. Its phase gives the timing modulo a
  symbol. After it fires, the decoder waits until the sync and the symbol
  after it are in the buffer, mixes and filters the whole buffer of raw input
  it kept (`catch_up()`), and `refine()` picks the start among the candidates
  within a symbol of the best window by the correlation between consecutive
  symbols, which also gives the fraction of a tone of the CFO. `preamble()`
  then estimates SFO and residual CFO. The detection rests while waiting, so
  it keeps what it found, and `busy()` counts the waiting. It keeps running
  during a frame: the decoder lags a few symbols behind the input, and the
  sync of a frame right after another would be gone by its end.
- `analyzer.hh`: spectrum, spectrogram and X-Y oscilloscope of the decoder's
  baseband for the page, fed by the decoder module and drawn as RGBA pixels.
  It takes FFTs as long as the decoder's symbols, so the bins are as wide as
  the tones, and draws the middle 400 bins, from 0 to 3000 Hz. None of its
  views has gain control: the page sets the scope gain and the top and bottom
  of spectrum and spectrogram in dB of a full scale sine. The spectrum shows
  three curves: the current one and the dimmed min and max of each bin over
  the last 32 updates, about two seconds, with the moving extrema of
  `movext.hh` from `dsp`. It also plots the
  power and phase of the decoder's channel estimate of the 320 tones after
  each symbol after synchronization, when the decoder reports the `new_symbol`
  event, which the module keeps from the page. Both have fixed scales: the
  power from 0 to -100 dB of a full scale sine, the phase from -180 to 180
  degrees, one dot per tone. After each decoded frame it plots the real Es/N0
  of each tone, `real_tone_snr` of the decoder, from 0 to 50 dB, and after
  each digital symbol the precision of each tone the soft demapper uses,
  `tone_precision` of the decoder, its estimated Es/N0, on the same scale,
  blank for erased or analog symbols. Its
  constellation shows the equalized tones of the latest symbol on a fixed
  scale: the BPSK tones of the preamble symbol, then only the payload tones
  that carry all bits of the mode, without seed tones, the tones with fewer
  bits that pad the code to whole symbols, or erased symbols. All views are
  off until the page turns them on, and take no time while off.
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

`make test` encodes random data to `audio.wav` and decodes it. `make play-modes`
plays a frame in every digital mode with sox `play`, each carrying the name of its
mode, for testing the decoder over the air. There is no test
suite; check a change by encoding, decoding and comparing, for example with the
simulation pipeline in the README ("Simulating"), which adds multipath, CFO,
SFO and AWGN with the disorders tools and ends with `diff -s`.
The decoder reports the operation mode and Es/N0 statistics on stderr. After a
frame decodes, it re-encodes the message and also reports the number of code
bits whose soft value had the wrong sign (erased bits count too) and the real
Es/N0 against the re-encoded symbols, which is exact also near threshold.
It also keeps the real Es/N0 of each tone, averaged over the payload symbols
of the frame, in `real_tone_snr`, which the page plots.

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
the decoder modules. To try `index.html`, run `make webserver`, which builds
the modules and serves the repository on http://127.0.0.1:8080/ only, and
open `index.html` there. Any other static web server works too, as long as it
sends `.wasm` as `application/wasm`.
`make pages` builds the modules and copies the page and the modules to `docs/`,
the snapshot GitHub Pages serves, so work on `next` doesn't change the public
page until `docs/` is rebuilt and committed.

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

## Design decisions

The signal is made for audio channels such as room acoustics, audio equipment
and FM or AM radio, not for the channels textbook OFDM is tuned to. Measured
with the disorders tools: a Watterson channel with 2 ms delay and 1 Hz spread
has a notch every 500 Hz, about 20 dB deep, its phase turns about five times
across the band, and the tones next to the notches change a lot from one
symbol to the next as the notches move. The LTE channel models EPA, EVA and
ETU, by contrast, are smooth over the bandwidth of LTE's sync signal and
barely change within a symbol. Results that rely on such channels don't carry
over. Before proposing to change any of the following, test it on such a
channel, also near threshold:

- The sync symbol is BPSK from an MLS, not a CAZAC sequence such as
  Zadoff-Chu. With 320 tones inside a larger FFT, the best Zadoff-Chu root
  has only about 2 dB less PAPR than the MLS, many roots have more, and the
  real audio signal adds 3 dB to all of them. The ratio of adjacent
  Zadoff-Chu tones turns by the same phase from tone to tone, so a shift by
  whole tones only turns it as a whole, and the CFO search of the sync
  detection can't tell the shifts apart for any root. The sync symbols get
  no SLM: the receiver has to know them exactly.
- There are two identical sync symbols back to back. They let the sync
  detection hop once per symbol with a window over two symbols, which then
  sees mostly the repeating sync, and their correlation picks the start and
  the fraction of a tone of the CFO. A third symbol isn't needed for
  detection.
- Sync detection uses the ratio of adjacent tones, a differential across
  frequency, not across time and not a correlation with the known tones.
  - On the Watterson channel the ratio stays flat apart from the tones next
    to a notch. A ratio across time fails exactly at the notches, as they
    move between symbols.
  - Correlating with the known tones, with an inverse FFT to search all
    delays, detects about 3 to 4 dB lower, as a ratio multiplies noise by
    noise. But the preamble already fails before the differential sync does,
    so this gains no decoded frames.
  - The ratio treats every channel alike. Its noise penalty is won back by
    many tone pairs, which longer symbols provide.
- Symbols are long: 7.5 Hz tones, a 133 ms symbol and a 3.3 ms guard
  interval. Neighbouring tones are close enough that even a 2 ms echo only
  turns the phase by about 5 degrees from one tone to the next, and there are
  320 tones to average. The price is a search over 133 shifts by half a
  tone for a CFO of 250 Hz, which is cheap once per symbol. A
  room's reverberation lasts far longer than the guard interval, which only
  covers the early reflections.
- Most of the CFO is SFO in disguise. FM and AM demodulation don't pass an
  offset of the RF carrier to the audio, unlike SSB. What remains are the
  sound card clocks, whose SFO scales all frequencies, 147 ppm looks like a
  CFO of 0.22 Hz at 1500 Hz, and acoustic Doppler. Where sound travels
  through the air, as from a speaker to the microphone of a walkie-talkie,
  moving either changes the length of the path: with sound at 343 m/s, 10 cm/s
  is already about 290 ppm, more than typical sound card crystals, and a
  casual 0.5 m/s about 1500 ppm, 2 Hz at 1500 Hz. It changes over time, a
  time-varying SFO, which is why people have to hold the device still.
  - The fractional CFO from the sync is still needed. The seed tones track
    slow phase drift from symbol to symbol, but an offset of a fraction of a
    tone causes interference between the tones: 0.44 Hz, 6 percent of a tone,
    puts it at about -19 dB, which would cap the dense QAM modes. The
    correlation of the two sync symbols weights the tones by their power,
    and with SFO their phase also turns along the band, so a channel
    stronger on one side biases it. `preamble()` corrects it by the Theil-Sen
    intercept of the phases between the two sync symbols, at the middle tone,
    and gets the SFO from their slope. Then it takes the two sync symbols
    again, resampled with that SFO and with the corrected CFO, and the first
    symbol too, as their interference between the tones would otherwise stay
    in the channel estimate: it took the first symbols of a frame to track it
    out, and the first symbol had about 12 dB Es/N0 at 1000 ppm. The second
    time also measures the timing, by the phase of the ratio of adjacent
    tones of both sync symbols against their known ratio, now without the
    CFO, which leaves the timing of the sync detection late by up to
    half a guard interval at 2000 ppm. It moves the windows and takes them a
    third time only when the ratios agree to 88 percent and the timing is off
    by more than an eighth of the guard interval: a radian is a thousand
    samples, so near threshold the estimate scatters by tens of samples, and
    moving there cost frames. Both sync symbols halve its error against the
    first alone, which let the agreement drop from 95 percent and gained
    frames near threshold with sync detection in the frequency domain and
    SFO. The sync detection stays cheap, all of this runs once per sync.
  - The decoder resamples every symbol after the preamble, without a separate
    estimator. `step`, the received samples per sent sample, starts from the
    SFO of the sync. Each symbol, the Theil-Sen slope and intercept of the
    seed tones give how far the timing drifted and how far the phase turned
    since the last window. `resync()` corrects the start of the next window by
    the whole drift, `step` and the CFO by a quarter of it, a second order
    loop each, and turns `chan` by what the window moved. `process()` reads
    the window with cubic Lagrange interpolation around its middle, where
    the unresampled sync windows measure the timing, starting
    `window_margin` samples into the buffer to reach before it. Without
    it, SFO caps the Es/N0 by the interference between the tones, about
    27 dB at 300 ppm and 19 dB at 1000 ppm, and the timing drifts beyond
    the guard interval on long frames.
- Adjacent tones are divided, not multiplied, so a strong tone or carrier
  can't dominate, and ratios that are far off, as next to a notch, are erased.
- The decoder doesn't trust the channel estimate of the sync for the whole
  frame. It tracks the channel with the seed tones of every symbol and
  weights the soft bits of each tone by its own channel power and noise.
- The noise of each tone is measured on known symbols only: the difference
  of the two sync symbols, then the seed tones, which take turns over all
  tones. Hard decisions on the payload get the noise of faded
  tones wrong, as their errors wrap around at the decision boundaries, which
  cost frames near threshold on multipath.
  Few symbols give a noisy estimate per tone, so it only counts where it is
  more than twice the median, as with interference in part of the band,
  and white noise keeps the weights of the channel power alone.

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
- Each tone's precision is its channel power over the white noise of its
  symbol before equalization, the mean over all tones of their errors times
  their channel power, and less where its own noise is more than twice the
  median. Interference raises that mean and lowers all precisions alike, by
  about 3 dB for noise across a sixth of the band, which the decoder ignores
  as a common scale, so the page shows the shape of the real Es/N0 but not
  always its level. The logged Es/N0 stays unweighted.
- The level of the precision and the logged Es/N0 come from hard decisions
  and overestimate Es/N0 by several dB for dense QAM near threshold. Within a
  frame that is a common scale and does not hurt decoding, but they are not
  exact there. The real Es/N0 after decoding is.
- The signal is always a real (one channel) audio signal centered at 1500 Hz.
  The sync detection works on the raw input, where its window keeps DC, hum
  and out-of-band tones out of the tones. Only after a sync, during a frame,
  and while the page shows its spectrum or scope (`keep_baseband`) does the
  decoder mix the signal down by 1500 Hz and lowpass filter it with
  `DSP::DownConversion`, flat to ±1050 Hz and at least 65 dB down from
  ±1500 Hz. With the search, this limits the CFO to about ±250 Hz, where the
  outer tones are already attenuated.
- Near threshold, the coarse timing of the sync detection scatters by tens
  of samples, which `preamble()` corrects. BPSK 1/2 short then decodes all frames
  from about 3 dB in-band SNR, also with multipath, CFO and SFO, and fails at
  about the same point as its code does. To compare decoder changes, use many
  noise seeds per point near the threshold of the affected modes, and report
  decoded frames out of trials. Run old and new on the same seeds and count
  the trials only one of them decodes: two totals that differ by a few
  frames are usually noise.
- The guard interval of the symbol after the two sync symbols starts as their
  cyclic continuation, as the encoder's guard interval weights are zero for
  its first quarter. A window that reaches past the second sync symbol can
  therefore look clean, and the sync detection, whose ratio of adjacent
  tones hardly notices a few guard intervals of other symbols, may pick such
  a window. Silence before the sync looks just as clean, so two windows a hop
  apart can both look clean, and the best one can be early or late by half a
  symbol: `refine()` looks a whole symbol either way. Don't assume the best
  window lies within the sync symbols.
- When testing anything that works on a grid of hops or blocks, such as the
  sync detection, sweep the offset of the signal over a whole hop, for
  example by prepending 0 to `symbol_len` zero samples to a clean WAV. Native
  runs of the same file always see the same phase; the Node checks feed
  their files back to back and see others, which is how a failure at a few
  phases showed up.
- Test false syncs with structured signals, such as tones, chirps or speech,
  not only with white noise, which a detector weighted by power never
  mistakes for a sync but a strong tone can.
- Don't take the QAM1024 and QAM4096 modes too seriously. They are only there
  for fun and might be replaced with something else in the future, so don't
  tune or design anything around them.

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

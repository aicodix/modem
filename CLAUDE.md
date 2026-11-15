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
- `check_encoder.mjs`, `check_decoder.mjs`, `check_frozen.mjs`, `frozen/`,
  `index.html`: the Node.js checks of the modules, the recordings of the
  frozen air interface (see "Frozen air interface") and the browser page
  that sends text and receives it from the microphone or a WAV file. Its
  meters show the CFO and SFO the decoder tracks after each symbol, from the
  module's `tracking` buffer, on symmetric logarithmic scales, dimmed
  between frames, and the RMS and peak level of the microphone, which the
  page measures itself. Like the views of the analyzer, the meters are off
  until the page turns them on and take no time while off, and the page only
  redraws while at least one of them is on.
- `common.hh`: what both sides share, such as tone layout, sequences, the
  interleaver (`shuffle()`) and the operation mode in `Common::setup()`.
  The mode byte carries the number of payload symbols, 2 << (mode >> 5)
  from 4 to 256, and an index into the table `efficiency` of 32 pairs of
  bits per tone and code rate k / 16, each set from the most robust to the
  least. Indices 0 to 25 hold the pairs from QPSK 3/16 to QAM256 15/16
  that no other pair beats in both bits per tone and threshold on AWGN,
  without those that add less than about 0.6 dB, listed under "Design
  decisions". From `bpsk_first`, 26 to 31 hold BPSK 1/16, 2/16, 3/16,
  4/16, 6/16 and 8/16 for bad channels, which `pick_mode()` only picks
  when asked for, for a negative SECONDS of `encode` or "BPSK only" on the
  page, and never otherwise.
  The code bits of a frame are a power of two, or three times one for
  QAM64, split into codewords of up to order 16 with a 32 bit CRC each,
  and the interleaver spreads them across
  the whole frame, skipping the indices beyond its size below a power of
  two. A frame only decodes if all its codewords do. `pick_mode()` picks
  the lowest index of the asked for set that carries the payload in a frame
  of at most the given seconds, and of those the fewest symbols, as fewer
  symbols give about the same thresholds. A mode with
  mode >> 5 of zero is a ping, a frame without payload symbols, that only
  carries the value mode & 31.
  The xorshift
  triples of the interleaver were chosen from all with full period by the
  rank of the low and of the top bits of an index together with those of
  the indices up to eight steps later, which tells how evenly neighbours
  spread, with `tools/rank.cc` from `code`. Ties went to the triple whose
  pairs of indices look most like those of a random permutation.
  `seed_shuffle()` spreads the 64 bits of the Hadamard codeword of the seed
  across the seed tones with the 6 bit xorshift triple that, of all six with
  full period, keeps the most Hamming distance between the codewords on the
  seed tones left over when a burst of adjacent seed tones is lost: up to 57
  of 64, instead of 31 without it.
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
  sync of a frame right after another would be gone by its end. As the
  detection hops by a symbol, but frames by a symbol and its guard interval,
  it fires up to about half a symbol early or late, depending on the length
  of the frame before. If that frame hasn't ended when the waiting is over,
  the preamble waits for it up to a symbol and `refine()` searches that much
  earlier in the buffer, but not before the windows.
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
  each symbol the precision of each tone the soft demapper uses,
  `tone_precision` of the decoder, its estimated Es/N0, on the same scale,
  blank for erased symbols. Its constellation shows the equalized tones of
  the latest symbol on a fixed scale: the BPSK tones of the preamble symbol,
  then the payload tones, without seed tones or erased symbols. All views are
  off until the page turns them on, and take no time while off.
- `freezer.cc`: generates the frozen bits of the metadata code in
  `polar_tables.hh`, and `tools/nested.cc` from `code` the nested tables of
  the payload codes of length 1024 to 65536 with a 32 bit CRC, where bit i
  is frozen for a payload of k / 16 of the code if its 4 bit entry is at
  least k. Don't edit the tables by hand; regenerate them with
  `make polar_tables.hh`.

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
plays frames of every number of payload symbols, each with a small payload
and with one close to the largest, with sox `play`, each saying its seconds
and size, for testing the decoder over the air. There is no test
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
the decoder modules. `check_frozen.mjs` decodes the recordings in `frozen/`
with the decoder modules and the native `decode` and compares them with
`frozen/frames.json` and the texts next to them, padded with zeros. To try
`index.html`, run `make webserver`, which builds
the modules and serves the repository on http://127.0.0.1:8080/ only, and
open `index.html` there. Any other static web server works too, as long as it
sends `.wasm` as `application/wasm`.
`make pages` builds the modules and copies the page and the modules to `docs/`,
the snapshot GitHub Pages serves, so work on `next` doesn't change the public
page until `docs/` is rebuilt and committed. Robusto isn't published for now:
don't run `make pages` here, and leave `docs/` as the snapshot of `next`, as
deleting it would conflict with every rebase onto `next`.

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

## Frozen air interface

The scheme is in use by a BBS and is frozen: frames sent today have to decode
with every later decoder. Don't change anything on the air. A change that
can't be avoided needs a different sync sequence, as `next` and `robusto`
have, so neither finds the frames of the other: the sync is the version.
`frozen/` holds recordings that the encoder of the frozen scheme made, 16
bit and with a quarter second of silence before and after, of memes from the
call signs RICK, CATS and DOGE, that every decoder has to get: a ping, BPSK
1/16 and 8/16, QPSK 3/16 and 8/16, QAM64 11/16 and QAM256 15/16, at 44100
and 48000 Hz, with codes from 1024 to 8192 bits. `make check-wasm` decodes
them with every decoder, and every later decoder has to decode them. Never
replace them with new ones.

On the air, and frozen:

- The signal: real audio at 44100 or 48000 Hz, 320 tones 7.5 Hz apart
  centered at 1500 Hz, symbols of `guard_len * 40` samples with a guard
  interval of `guard_len = rate / 300`, so a symbol and its guard last
  41 / 300 seconds, and the shape of the guard interval from
  `guard_interval_weights()`.
- The frame: one symbol of leading noise from the MLS of `mls2_poly`, two
  identical sync symbols of BPSK from the MLS of `mls0_poly` 0x2d1 from
  `mls0_seed` 214, the preamble symbol, 2 << (mode >> 5) payload symbols,
  none for a ping, and the closing guard interval.
- The preamble: 48 bits of the base 40 call sign and 8 bits of the mode,
  with the CRC-16 of poly 0xA8F4, polar coded with `frozen_256_72`,
  interleaved by `shuffle()` over 256 bits and sent as BPSK on its data
  tones.
- The mode: the meaning of the mode byte and the table `efficiency`, its
  32 pairs of modulation and rate in their order, and `setup()`, which
  derives the code length, the number of codewords and the frozen bits.
- The seed tones: every fifth tone from `seed_off`, which moves by
  `block_skew` from `first_seed` each symbol, carrying the codeword of
  `HadamardEncoder<7>` of the seed value, spread by `seed_shuffle()`, and
  the scrambling of the other tones of a symbol by the MLS of `mls2_poly`
  0x163 from the seed value, none for seed value 0.
- The payload: zero padding to the payload of the mode, the `Xorshift32`
  scrambler from its default state, the split into codewords of
  `data_bits / codewords` bits, each with its CRC-32 of poly 0x8F6E37A0,
  the bits of each byte from the lowest, the nested polar tables and the
  polar encoder of `code`, `shuffle()` with its xorshift triples over the
  whole frame, and the order of the code bits over the data tones and
  symbols.
- The constellations of `psk.hh` and `qam.hh`, with their bit mapping.
- The tables in `polar_tables.hh`, even if `freezer.cc` or `tools/nested.cc`
  in `code` change: regenerating them must not change them.

Not on the air, and free to change:

- Which seed value SLM picks and how it searches, and the clipping and
  filtering of the encoder, as long as the tones stay those of the frame.
- Which mode `pick_mode()` picks, and with it the options of `encode` and
  the page.
- Everything in the decoder: sync detection, timing, tracking, demapping,
  the polar list decoder.

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
- The sync sequence is the first 320 bits of the MLS of poly 0x2d1 from seed
  214. Nothing correlates with the tones themselves: the sync detection and
  the timing in `preamble()` only use the ratio of adjacent tones, which is
  again the same MLS, just further along. So a sequence has to be chosen for
  two things at once. Its PAPR follows the correlation of the tones
  themselves, as the power over time is the Fourier transform of it, while
  the ratios have to correlate well against the shifts of the CFO search.
  The two don't depend on each other, so a search for low PAPR alone, which
  gave the poly 0x331 of next, can land on poor ratios: up to 29 of 319
  within 250 Hz. Of all windows of all degree 9 MLS, 0x2d1 from 214 has
  7.2 dB PAPR as the real audio signal, 0.3 dB less, and 19 of 319. Longer
  MLS have about 1 dB more PAPR. Golay sequences of length 320 have at most
  3 dB PAPR at baseband, but their ratios repeat, with sidelobes from 77 to
  209 of 319, and Legendre sequences get about 6 dB of audio PAPR, but no
  better ratios than an MLS. A search over all 320 bits gets 6.5 dB at 17 of
  319, but needs a table. The sync is also a different sequence than that of
  next, so neither finds the frames of the other, though the metric of one
  on the sync of the other can reach about 0.2, the threshold.
- There are two identical sync symbols back to back. They let the sync
  detection hop once per symbol with a window over two symbols, which then
  sees mostly the repeating sync, and their correlation picks the start and
  the fraction of a tone of the CFO. A third symbol isn't needed for
  detection.
  - A train of more sync symbols in place of the leading noise was tried.
    Averaging them into the channel estimate decoded the meta data of 193
    instead of 139 of 200 frames on static multipath near threshold, but
    nothing more on a Watterson channel with 1 Hz spread. Finding the end of
    a fading train needed the coherence of consecutive symbols in
    `refine()`, a bigger buffer and more, and a variable number of sync
    symbols works against concatenating frames. It was dropped.
  - A single sync symbol was modeled. A Hann window over one period leaks
    every tone into its neighbours, so the ratio of adjacent tones mixes
    them, and shorter windows only smear them more: about 3 dB less
    detection, and nothing repeats for the fine CFO and the SFO. Tones on
    every other subcarrier, a symbol of two identical halves, keep the
    scheme at half the scale, but still detect about 2 dB worse. Either way
    the channel estimate comes from one symbol, about 1.2 dB worse for the
    meta data, which already fails before the sync detection does.
  - The window is exactly two symbols long, so its leakage lands on the odd
    bins, between the tones. Shorter Hann windows leak into the neighbouring
    tones and detect worse. Hops shorter than a symbol, or the worst
    alignment of the windows against the sync, change nothing measurable
    near threshold.
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
  - Each symbol moves the estimate halfway to the seed tones and a third of
    the way to the spline between them, an average over a few symbols. On a
    Watterson channel with 1 Hz spread that lags behind the fading, which
    left a real Es/N0 of about 0 dB even without noise. Only BPSK is meant
    for such skywave channels, so only BPSK follows each symbol more: all
    the way to the seed tones and two thirds to the spline. At 32 symbols,
    64 seeds, CFO and SFO, at 1 Hz without noise BPSK 1/16, 2/16, 3/16,
    4/16, 6/16 and 8/16 decoded 42, 37, 33, 32, 27 and 27 frames instead of
    25, 21, 18, 14, 12 and 13, at 1 Hz and an `awgn` level of -10 27 to 17
    instead of 20 to 6, and at 0.5 Hz and -3, where noise decides, about as
    many as before. On AWGN it cost at most 0.2 dB for 1/16 to 4/16, which
    the preamble limits anyway, and 0.5 to 0.6 dB for 6/16 and 8/16.
    Following only the current symbol gained more without noise, but cost
    up to 1.3 dB on AWGN and lost frames where noise decides. The other
    modes keep the averaging: no averaging would have cost them 1.3 to
    1.9 dB on AWGN, and they aren't meant for such channels.
- The noise of each tone is measured on known symbols only: the difference
  of the two sync symbols, then the seed tones, which take turns over all
  tones. Hard decisions on the payload get the noise of faded
  tones wrong, as their errors wrap around at the decision boundaries, which
  cost frames near threshold on multipath.
  Few symbols give a noisy estimate per tone, so it only counts where it is
  more than twice the median, as with interference in part of the band,
  and white noise keeps the weights of the channel power alone.

- The seed tones carry the codeword of the seed without a scrambler. As the
  xorshift of `seed_shuffle()` is linear, every permuted codeword on the
  seed tones after the first is a shift of the MLS of poly 0x6d or its
  complement, apart from the all equal ones of seed values 0 and 64, which
  SLM avoids. Without the permutation, the Walsh codes on the evenly spaced
  seed tones are peaky and needed scrambling by the MLS of poly 0x43,
  running across the symbols of a frame: without it, the median PAPR after
  SLM over 256 symbols of QPSK to QAM256 rose from 7.06 to 7.88 dB, and all
  128 seed values were tried. With the permutation, that scrambler raised it
  from 7.01 to 7.08 dB, 68 instead of 55 percent of the symbols stayed above
  the 7 dB that ends the search early, and each symbol tried about 106
  instead of 96 seed values. Other MLS of degree 6 to 11 did the same, apart
  from 0x6d, which only turns codewords into codewords. Worst cases stayed
  at about 7.5 to 7.7 dB either way.
- The search of SLM starts after the seed value of the symbol before, from 0
  at each frame, and ends at the first seed value below 7 dB. Starting from
  0 at each symbol, BPSK, which ends early, used only about 84 of the 126
  seed values that SLM picks, one of them for 7 percent of the symbols, and
  gave a symbol the seed value of the one before 3.4 percent of the time.
  Now it uses all 126 about evenly and never repeats one, with the same
  PAPR. Denser modulations mostly search all 128 seed values, so where they
  start doesn't matter. The decoder doesn't change.
- The interleaver triples ranked by `tools/rank.cc` replaced hand-picked
  ones whose pairs of indices clustered far beyond those of random
  permutations. Near threshold on static multipath and the Watterson
  channel, 336 paired trials decoded 23 frames only with the new and 27 only
  with the old interleavers: no measurable difference.
- Thresholds of all 60 pairs of QPSK, QAM16, QAM64 and QAM256 with rates
  from 1/16 to 15/16, measured on the current robusto with the mode forced,
  at 32 symbols, as the `awgn` level at which half of the frames decode,
  bisected to 0.5 dB steps of 32 seeds each and interpolated, on AWGN with a
  CFO of 234.567 Hz and an SFO of 147 ppm. A higher level is more noise.
  Those that no other pair beats in both bits per tone and threshold:
  QPSK 1/16 -0.2, 2/16 -0.5, 3/16 -1.0, 4/16 -1.9, 5/16 -2.8, 6/16 -3.6,
  7/16 -4.2, 8/16 -4.8, 9/16 -5.5, 10/16 -6.2, 11/16 -6.8, 12/16 -7.6,
  13/16 -8.3;
  QAM16 7/16 -9.1, 8/16 -9.8, 9/16 -10.8, 10/16 -11.7, 11/16 -12.6,
  12/16 -13.4, 13/16 -14.6, 14/16 -15.8;
  QAM64 9/16 -15.3, 10/16 -16.3, 11/16 -17.7, 12/16 -18.8, 13/16 -20.1,
  14/16 -21.5, 15/16 -23.4;
  QAM256 10/16 -20.7, 11/16 -22.2, 12/16 -23.7, 13/16 -25.2, 14/16 -26.8,
  15/16 -28.9.
  - Unlike the earlier measurement, the highest rates of each modulation now
    beat the next modulation at the same bits per tone: QPSK 12/16 beats
    QAM16 6/16 by 0.6 dB, QAM16 13/16 and 14/16 beat QAM64 at the same bits
    per tone by about 0.5 dB, and QAM64 12/16 beats QAM256 9/16 by 0.5 dB.
    Every curve falls from all to none of the frames without an error floor,
    and the 26 pairs of the table before came out within 0.1 dB of the run
    before this one.
  - The table leaves out QPSK 1/16 and 2/16, only 0.6 to 0.8 dB beyond QPSK
    3/16 where the preamble limits, as the BPSK modes cover bad channels,
    and QPSK 7/16 and 11/16, QAM64 9/16, 10/16 and 15/16 and QAM256 10/16,
    as each adds less than about 0.6 dB to its neighbours. Its steps are 0.7
    to 1.5 dB, 1.9 dB from QAM16 14/16 to QAM64 11/16 and 2.1 dB to QAM256
    15/16. Where two modulations came close, it keeps the lower one, which
    suffers less from phase noise. Code lengths differ between the
    modulations, so the order may shift a little at other numbers of
    symbols.
  - Multipath is left out, as that of a room depends too much on the room.
    BPSK is measured on a Watterson channel, where fading, not noise,
    decides.
- BPSK and QPSK from 1/16 and 3/16 to 15/16 and 8/16 on a Watterson channel,
  before BPSK tracked the channel faster (see above), with 2 ms delay, 64
  seeds of fading and noise, at 32 symbols with CFO and SFO: at 1 Hz spread
  without noise, BPSK 1/16, 2/16, 3/16 and 4/16 decode 31, 23, 16 and 12
  frames, so lower rates pay off where the fading limits, and only BPSK goes
  that low. At 0.5 Hz they decode 55 to 49, at 0.1 Hz all of them, where
  noise again decides. At the same bits per tone QPSK decodes as many frames
  as BPSK or more, also on fading, as QPSK 3/16 against BPSK 6/16 and QPSK
  4/16 against BPSK 8/16. So BPSK above 4/16 only pays off for phase noise,
  which the disorders tools don't simulate, and from 10/16 on it falls off
  quickly even at 0.1 Hz. The BPSK set keeps 1/16 to 4/16 for fading and
  6/16 and 8/16 for phase noise.
  - Unlike on AWGN, the payload, not the preamble, limits there: without
    noise, sync never failed, and the preamble lost the same 5 of 64 frames
    at 1 Hz in every mode and none at 0.5 Hz, so it would let about 92
    percent through at 1 Hz, where BPSK 1/16 decodes 48. All other frames
    failed in the payload, which is why each lower rate still gains.
  - The gain of 1/16 and 2/16 holds with some noise: at 1 Hz and an `awgn`
    level of -15, BPSK 1/16 to 4/16 decode 30, 20, 16 and 14 frames, at -10
    17, 17, 13 and 13, and at 0.5 Hz and -3, near where noise decides, 22,
    14, 14 and 10. Their flatter FER curves keep them from the error rates
    of the higher rates on good channels, but where fading limits they get
    the most frames through.
- Earlier thresholds of 54 pairs of modulation and rate k / 16, from BPSK
  1/16 to QAM256 15/16, with the mode forced, at 32 symbols, as the `awgn`
  level of the disorders tools at which half of the frames decode,
  interpolated from 16 seeds per 1 dB step on static multipath with CFO and
  SFO, before the changes since. A higher level is more noise. Those that
  no other pair beats in both bits per tone and threshold:
  BPSK 1/16 -0.7, 2/16 -0.8, 4/16 -1.0;
  QPSK 3/16 -1.5, 4/16 -2.5, 5/16 -3.5, 6/16 -4.5, 7/16 -5.3, 8/16 -5.6,
  9/16 -6.5, 10/16 -7.5, 11/16 -8.5;
  QAM16 6/16 -8.7, 7/16 -9.8, 8/16 -10.6, 9/16 -11.6, 10/16 -12.5,
  11/16 -13.5, 12/16 -14.9;
  QAM64 9/16 -16.4, 10/16 -17.5, 11/16 -18.7, 12/16 -20.5;
  QAM256 9/16 -20.5, 10/16 -22.5, 11/16 -24.4, 12/16 -26.3, 13/16 -28.5,
  14/16 -31.5, 15/16 -38.2.
  - Below about 0.4 bits per tone the preamble, not the payload, fails
    first: 70 to 100 percent of the failures of BPSK 1/16 to 5/16 are sync
    or preamble errors, and their thresholds lie within 1 dB. Lower rates
    only help where the channel itself limits, see below.
  - QPSK beats BPSK at the same bits per tone, by 1 dB at 0.375 and 0.5 and
    by 2 to 4 dB from 0.625 to 0.875, so BPSK above 4/16 is never worth it.
    The higher rates of each modulation lose to the next one at a lower
    rate: QAM16 6/16 beats QPSK 12/16 by 0.8 dB, QAM64 5/16 QPSK 15/16 by
    1.9 dB, QAM16 9/16 QAM64 6/16 by 0.9 dB. At 3 bits per tone QAM16 12/16
    and QAM64 8/16 are within 0.4 dB, and QAM64 9/16 carries more at the
    threshold of QAM16 13/16. QAM64 12/16 and QAM256 9/16 tie.
  - 4 symbols, codes of a quarter or an eighth the length, give the same
    thresholds within 0.5 dB, apart from QAM256 15/16, 3.4 dB better.
  - On a Watterson channel with 2 ms delay and 0.1 Hz spread, 8 seeds per
    2 dB step, the order stays, QAM256 9/16 beats QAM64 12/16 by 1 dB, and
    QAM256 12/16 and above never decode half the frames even without
    noise: the fading alone leaves a real Es/N0 of about 14 to 22 dB for
    32 symbols. It drops to about 8 dB at 0.5 Hz, and at 1 Hz to about
    0 dB, where without noise BPSK 1/16 decodes 21 percent, BPSK 2/16 and
    QPSK 3/16 6 percent and rates from 1/2 none of 32 frames. An older build
    of robusto with the previous modes gave the same Es/N0, so this isn't a
    regression of the new modes.

## Conventions worth knowing

- Encoder and decoder must stay in step, and the air interface is frozen,
  see "Frozen air interface". A change to anything on the air, such as the
  mode table, preamble, sequences, scrambler, frozen bits or constellation
  mapping, would break every frame already sent.
- Soft demapping produces int16 LLRs for a min-sum polar list decoder, which
  ignores their overall scale. What matters is resolution after rounding and
  headroom against saturation. `quantize()` caps the slope `DIST * precision`
  at 4096, so every constellation gets the same headroom, and `process()`
  scales the precision by 300 for resolution.
  - int8 LLRs were measured and lose to int16 where dense QAM meets high
    rates. Translating the int16 scaling (cap 16, two bits of headroom)
    costs about 0.5 dB on QAM256 5/6 and 0.25 dB on QAM64 and QAM256 1/2.
    Larger caps saturate the sums inside the polar decoder: a cap of 32 had
    an error floor on QAM64 1/2 2 dB above its threshold, and clipping the
    channel LLRs lower didn't help. The best found, scale 1.75 and cap 24,
    still lost 0.1 dB at rate 1/2 and 0.3 to 0.5 dB on QAM256 13/16 and 5/6.
    PSK and QAM16 lost nothing. The lowest and highest rates of the nested
    codes need the range even more.
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
  of samples, which `preamble()` corrects. BPSK at rate 1/2 over 8 symbols then decodes all frames
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
  phases showed up. Frames right after another need the same sweep, as half
  of the phases once lost the second frame.
- Test false syncs with structured signals, such as tones, chirps or speech,
  not only with white noise, which a detector weighted by power never
  mistakes for a sync but a strong tone can.

## Workflow

- Work on the `robusto` branch, the successor of `next`, whose air
  interface is now frozen (see "Frozen air interface"); `next` only gets
  minimal changes. A redesign of the air interface needs a new branch with
  its own sync sequence.
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

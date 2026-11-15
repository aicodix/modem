
### OFDM MODEM

Quick start:

Create file ```uncoded.dat``` with ```2048``` bits of random data:

```
dd if=/dev/urandom of=uncoded.dat bs=1 count=256
```

Encode file ```uncoded.dat``` to ```encoded.wav``` [WAV](https://en.wikipedia.org/wiki/WAV) audio file with 48000 Hz sample rate and 32 bit float samples, in a frame of about 2 seconds:

```
./encode encoded.wav 48000 32 ANONYMOUS 2 uncoded.dat
```

Start recording to ```recorded.wav``` audio file and stop after 5 seconds:

```
sox -d -c 1 -e float recorded.wav
```

Start playing ```encoded.wav``` audio file:

```
play encoded.wav
```

Decode ```recorded.wav``` audio file to ```decoded.dat``` file:

```
./decode recorded.wav decoded.dat
```

Compare original ```uncoded.dat``` with ```decoded.dat```:

```
diff -s uncoded.dat decoded.dat
```

### Stream Audio

This currently works on macOS and Linux only.

Stream live audio into the decoder:

```
sox -q -d -c 1 -e float -t wav - | ./decode - decoded.dat
```

Stream audio out of the encoder:

```
./encode - 48000 32 ANONYMOUS 2 uncoded.dat | play -q -
```

### Supported Modes

All modes need a bandwidth of 2400 Hz.
The encoder picks the mode from SECONDS and the length of the input: first the number of payload symbols whose frame lasts closest to SECONDS, then the lowest efficiency index, the most robust mode, that carries all of the input, from the BPSK modes only for a negative SECONDS.
An input that no mode carries in about SECONDS is refused.
The decoder writes the whole payload of the mode, so a shorter input comes back padded with zeros.
The code bits of a frame are split into polar codewords of up to 65536 bits, each with its own 32 bit CRC, and interleaved across the whole frame.
QAM64 frames always have three codewords, or six in the longest frames, the others a power of two of them.

These are the frames for each number of payload symbols:

| Symbols | Duration |  SECONDS |
| ------- | -------- | -------- |
|       4 |     1.0s |        1 |
|       8 |     1.5s |        2 |
|      16 |     2.6s |        3 |
|      32 |     4.8s |   4 .. 6 |
|      64 |     9.2s |  7 .. 13 |
|     128 |    17.9s | 14 .. 26 |
|     256 |    35.4s | 27 .. 60 |

These are the efficiencies, the payload grows with the number of symbols.
The threshold is the noise level of `awgn` from [disorders](https://github.com/aicodix/disorders) at which half of the frames of 32 symbols decode on AWGN with CFO and SFO, a higher level is more noise.

| Index | Modulation | Code rate | Bits per tone | Payload per symbol | Threshold |
| ----- | ---------- | --------- | ------------- | ------------------ | --------- |
|     0 | QPSK       |      3/16 |         0.375 |                12B |      -1.0 |
|     1 | QPSK       |      4/16 |           0.5 |                16B |      -1.9 |
|     2 | QPSK       |      5/16 |         0.625 |                20B |      -2.8 |
|     3 | QPSK       |      6/16 |          0.75 |                24B |      -3.6 |
|     4 | QPSK       |      8/16 |             1 |                32B |      -4.8 |
|     5 | QPSK       |      9/16 |         1.125 |                36B |      -5.5 |
|     6 | QPSK       |     10/16 |          1.25 |                40B |      -6.2 |
|     7 | QPSK       |     12/16 |           1.5 |                48B |      -7.6 |
|     8 | QPSK       |     13/16 |         1.625 |                52B |      -8.3 |
|     9 | QAM16      |      7/16 |          1.75 |                56B |      -9.1 |
|    10 | QAM16      |      8/16 |             2 |                64B |      -9.8 |
|    11 | QAM16      |      9/16 |          2.25 |                72B |     -10.8 |
|    12 | QAM16      |     10/16 |           2.5 |                80B |     -11.7 |
|    13 | QAM16      |     11/16 |          2.75 |                88B |     -12.6 |
|    14 | QAM16      |     12/16 |             3 |                96B |     -13.4 |
|    15 | QAM16      |     13/16 |          3.25 |               104B |     -14.6 |
|    16 | QAM16      |     14/16 |           3.5 |               112B |     -15.8 |
|    17 | QAM64      |     11/16 |         4.125 |               132B |     -17.7 |
|    18 | QAM64      |     12/16 |           4.5 |               144B |     -18.8 |
|    19 | QAM64      |     13/16 |         4.875 |               156B |     -20.1 |
|    20 | QAM64      |     14/16 |          5.25 |               168B |     -21.5 |
|    21 | QAM256     |     11/16 |           5.5 |               176B |     -22.2 |
|    22 | QAM256     |     12/16 |             6 |               192B |     -23.7 |
|    23 | QAM256     |     13/16 |           6.5 |               208B |     -25.2 |
|    24 | QAM256     |     14/16 |             7 |               224B |     -26.8 |
|    25 | QAM256     |     15/16 |           7.5 |               240B |     -28.9 |

The BPSK modes are for bad channels and only picked on request, with a negative SECONDS for the encoder, or with "BPSK only" on the page:

```
./encode encoded.wav 48000 32 ANONYMOUS -5 uncoded.dat
```

Their low rates survive fast fading, 6/16 and 8/16 are for phase noise.
Decoded are the frames of 32 symbols out of 64 on a [Watterson](https://github.com/aicodix/disorders) channel with a delay of 2 ms and a spread of 0.5 or 1 Hz, without noise:

| Index | Code rate | Bits per tone | Payload per symbol | Decoded at 0.5 Hz | Decoded at 1 Hz |
| ----- | --------- | ------------- | ------------------ | ----------------- | --------------- |
|    26 |      1/16 |        0.0625 |                 2B |               86% |             48% |
|    27 |      2/16 |         0.125 |                 4B |               81% |             36% |
|    28 |      3/16 |        0.1875 |                 6B |               78% |             25% |
|    29 |      4/16 |          0.25 |                 8B |               77% |             19% |
|    30 |      6/16 |         0.375 |                12B |               73% |             17% |
|    31 |      8/16 |           0.5 |                16B |               72% |             19% |

### Pings

A ping is a frame without payload, only the call sign and a value from 0 to 31, as robust as a frame can be and about half a second long:

```
./encode ping.wav 48000 32 ANONYMOUS 0 7
```

The decoder reports it as ```ping 7 received``` and writes no output file for it.

### Simulating

Prerequisite: [disorders](https://github.com/aicodix/disorders)

Encode ```uncoded.dat``` to audio signal, add [multipath](https://en.wikipedia.org/wiki/Multipath_propagation), [CFO, SFO](https://en.wikipedia.org/wiki/Carrier_frequency_offset), [AWGN](https://en.wikipedia.org/wiki/Additive_white_Gaussian_noise), decode and compare:

```
./encode - 48000 32 ANONYMOUS 2 uncoded.dat | ../disorders/multipath - - ../disorders/multipath.txt 10 | ../disorders/cfo - - 234.567 | ../disorders/sfo - - 147 | ../disorders/awgn - - -20 | ./decode - - | diff -q -s uncoded.dat -
```

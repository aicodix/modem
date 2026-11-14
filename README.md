
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
The encoder picks the mode from SECONDS and the length of the input: first the number of payload symbols whose frame lasts closest to SECONDS, then the lowest efficiency index, the most robust mode, that carries all of the input.
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

These are the efficiencies, the payload grows with the number of symbols:

| Index | Modulation | Code rate | Bits per tone | Payload per symbol |
| ----- | ---------- | --------- | ------------- | ------------------ |
|     0 | BPSK       |      1/16 |        0.0625 |                 2B |
|     1 | BPSK       |      2/16 |         0.125 |                 4B |
|     2 | BPSK       |      3/16 |        0.1875 |                 6B |
|     3 | BPSK       |      4/16 |          0.25 |                 8B |
|     4 | BPSK       |      5/16 |        0.3125 |                10B |
|     5 | BPSK       |      6/16 |         0.375 |                12B |
|     6 | BPSK       |      7/16 |        0.4375 |                14B |
|     7 | BPSK       |      8/16 |           0.5 |                16B |
|     8 | BPSK       |      9/16 |        0.5625 |                18B |
|     9 | BPSK       |     10/16 |         0.625 |                20B |
|    10 | BPSK       |     11/16 |        0.6875 |                22B |
|    11 | BPSK       |     12/16 |          0.75 |                24B |
|    12 | BPSK       |     14/16 |         0.875 |                28B |
|    13 | QPSK       |      8/16 |             1 |                32B |
|    14 | QPSK       |      9/16 |         1.125 |                36B |
|    15 | QPSK       |     10/16 |          1.25 |                40B |
|    16 | QPSK       |     11/16 |         1.375 |                44B |
|    17 | QPSK       |     12/16 |           1.5 |                48B |
|    18 | QPSK       |     14/16 |          1.75 |                56B |
|    19 | QAM16      |      8/16 |             2 |                64B |
|    20 | QAM16      |      9/16 |          2.25 |                72B |
|    21 | QAM16      |     10/16 |           2.5 |                80B |
|    22 | QAM16      |     11/16 |          2.75 |                88B |
|    23 | QAM64      |      8/16 |             3 |                96B |
|    24 | QAM64      |      9/16 |         3.375 |               108B |
|    25 | QAM64      |     10/16 |          3.75 |               120B |
|    26 | QAM256     |      8/16 |             4 |               128B |
|    27 | QAM256     |      9/16 |           4.5 |               144B |
|    28 | QAM256     |     10/16 |             5 |               160B |
|    29 | QAM256     |     11/16 |           5.5 |               176B |
|    30 | QAM256     |     13/16 |           6.5 |               208B |
|    31 | QAM256     |     15/16 |           7.5 |               240B |

### Simulating

Prerequisite: [disorders](https://github.com/aicodix/disorders)

Encode ```uncoded.dat``` to audio signal, add [multipath](https://en.wikipedia.org/wiki/Multipath_propagation), [CFO, SFO](https://en.wikipedia.org/wiki/Carrier_frequency_offset), [AWGN](https://en.wikipedia.org/wiki/Additive_white_Gaussian_noise), decode and compare:

```
./encode - 48000 32 ANONYMOUS 2 uncoded.dat | ../disorders/multipath - - ../disorders/multipath.txt 10 | ../disorders/cfo - - 234.567 | ../disorders/sfo - - 147 | ../disorders/awgn - - -20 | ./decode - - | diff -q -s uncoded.dat -
```

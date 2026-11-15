
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

These are the efficiencies, the payload grows with the number of symbols.
The threshold is the noise level of `awgn` from [disorders](https://github.com/aicodix/disorders) at which half of the frames of 32 symbols decode on AWGN with CFO and SFO, a higher level is more noise. BPSK is meant for fading channels and has no AWGN threshold.
The indices from 27 to 31 are reserved.

| Index | Modulation | Code rate | Bits per tone | Payload per symbol | Threshold |
| ----- | ---------- | --------- | ------------- | ------------------ | --------- |
|     0 | BPSK       |      4/16 |          0.25 |                 8B |         - |
|     1 | QPSK       |      3/16 |         0.375 |                12B |      -1.0 |
|     2 | QPSK       |      4/16 |           0.5 |                16B |      -2.0 |
|     3 | QPSK       |      5/16 |         0.625 |                20B |      -2.8 |
|     4 | QPSK       |      6/16 |          0.75 |                24B |      -3.6 |
|     5 | QPSK       |      7/16 |         0.875 |                28B |      -4.2 |
|     6 | QPSK       |      8/16 |             1 |                32B |      -4.8 |
|     7 | QPSK       |      9/16 |         1.125 |                36B |      -5.5 |
|     8 | QPSK       |     10/16 |          1.25 |                40B |      -6.2 |
|     9 | QPSK       |     11/16 |         1.375 |                44B |      -6.8 |
|    10 | QAM16      |      6/16 |           1.5 |                48B |      -8.0 |
|    11 | QAM16      |      7/16 |          1.75 |                56B |      -9.0 |
|    12 | QAM16      |      8/16 |             2 |                64B |      -9.9 |
|    13 | QAM16      |      9/16 |          2.25 |                72B |     -10.8 |
|    14 | QAM16      |     10/16 |           2.5 |                80B |     -11.7 |
|    15 | QAM16      |     11/16 |          2.75 |                88B |     -12.5 |
|    16 | QAM16      |     12/16 |             3 |                96B |     -13.6 |
|    17 | QAM64      |      9/16 |         3.375 |               108B |     -15.2 |
|    18 | QAM64      |     10/16 |          3.75 |               120B |     -16.3 |
|    19 | QAM64      |     11/16 |         4.125 |               132B |     -17.6 |
|    20 | QAM256     |      9/16 |           4.5 |               144B |     -19.3 |
|    21 | QAM256     |     10/16 |             5 |               160B |     -20.8 |
|    22 | QAM256     |     11/16 |           5.5 |               176B |     -22.2 |
|    23 | QAM256     |     12/16 |             6 |               192B |     -23.6 |
|    24 | QAM256     |     13/16 |           6.5 |               208B |     -25.2 |
|    25 | QAM256     |     14/16 |             7 |               224B |     -26.8 |
|    26 | QAM256     |     15/16 |           7.5 |               240B |     -29.0 |

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

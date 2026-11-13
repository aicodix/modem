
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
The decoder writes the whole payload of the mode, so a shorter input comes back padded with zeros.
A frame has a single polar code of at most 65536 bits, so the long frames only have the lower efficiencies.
An input too large for them gets fewer symbols, and an input that no mode carries is refused.

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

These are the efficiencies, the payload grows with the number of symbols up to the longest frame:

| Index | Modulation | Code rate | Bits per tone | Payload per symbol | Longest frame |
| ----- | ---------- | --------- | ------------- | ------------------ | ------------- |
|     0 | BPSK       |      1/16 |        0.0625 |                 2B |   256 symbols |
|     1 | BPSK       |      2/16 |         0.125 |                 4B |   256 symbols |
|     2 | BPSK       |      3/16 |        0.1875 |                 6B |   256 symbols |
|     3 | BPSK       |      4/16 |          0.25 |                 8B |   256 symbols |
|     4 | BPSK       |      5/16 |        0.3125 |                10B |   256 symbols |
|     5 | BPSK       |      6/16 |         0.375 |                12B |   256 symbols |
|     6 | BPSK       |      7/16 |        0.4375 |                14B |   256 symbols |
|     7 | BPSK       |      8/16 |           0.5 |                16B |   256 symbols |
|     8 | BPSK       |      9/16 |        0.5625 |                18B |   256 symbols |
|     9 | BPSK       |     10/16 |         0.625 |                20B |   256 symbols |
|    10 | BPSK       |     11/16 |        0.6875 |                22B |   256 symbols |
|    11 | BPSK       |     12/16 |          0.75 |                24B |   256 symbols |
|    12 | BPSK       |     14/16 |         0.875 |                28B |   256 symbols |
|    13 | QPSK       |      8/16 |             1 |                32B |   128 symbols |
|    14 | QPSK       |      9/16 |         1.125 |                36B |   128 symbols |
|    15 | QPSK       |     10/16 |          1.25 |                40B |   128 symbols |
|    16 | QPSK       |     11/16 |         1.375 |                44B |   128 symbols |
|    17 | QPSK       |     12/16 |           1.5 |                48B |   128 symbols |
|    18 | QPSK       |     14/16 |          1.75 |                56B |   128 symbols |
|    19 | QAM16      |      8/16 |             2 |                64B |    64 symbols |
|    20 | QAM16      |      9/16 |          2.25 |                72B |    64 symbols |
|    21 | QAM16      |     10/16 |           2.5 |                80B |    64 symbols |
|    22 | QAM16      |     11/16 |          2.75 |                88B |    64 symbols |
|    23 | QAM16      |     12/16 |             3 |                96B |    64 symbols |
|    24 | QAM16      |     14/16 |           3.5 |               112B |    64 symbols |
|    25 | QAM256     |      8/16 |             4 |               128B |    32 symbols |
|    26 | QAM256     |      9/16 |           4.5 |               144B |    32 symbols |
|    27 | QAM256     |     10/16 |             5 |               160B |    32 symbols |
|    28 | QAM256     |     11/16 |           5.5 |               176B |    32 symbols |
|    29 | QAM256     |     12/16 |             6 |               192B |    32 symbols |
|    30 | QAM256     |     13/16 |           6.5 |               208B |    32 symbols |
|    31 | QAM256     |     15/16 |           7.5 |               240B |    32 symbols |

### Simulating

Prerequisite: [disorders](https://github.com/aicodix/disorders)

Encode ```uncoded.dat``` to audio signal, add [multipath](https://en.wikipedia.org/wiki/Multipath_propagation), [CFO, SFO](https://en.wikipedia.org/wiki/Carrier_frequency_offset), [AWGN](https://en.wikipedia.org/wiki/Additive_white_Gaussian_noise), decode and compare:

```
./encode - 48000 32 ANONYMOUS 2 uncoded.dat | ../disorders/multipath - - ../disorders/multipath.txt 10 | ../disorders/cfo - - 234.567 | ../disorders/sfo - - 147 | ../disorders/awgn - - -20 | ./decode - - | diff -q -s uncoded.dat -
```

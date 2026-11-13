
### OFDM MODEM

Quick start:

Create file ```uncoded.dat``` with ```2048``` bits of random data:

```
dd if=/dev/urandom of=uncoded.dat bs=1 count=256
```

Encode file ```uncoded.dat``` to ```encoded.wav``` [WAV](https://en.wikipedia.org/wiki/WAV) audio file with 48000 Hz sample rate and 32 bit float samples:

```
./encode encoded.wav 48000 32 ANONYMOUS QAM16 1/2 short uncoded.dat
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
./encode - 48000 32 ANONYMOUS QAM16 1/2 short uncoded.dat | play -q -
```

### Supported Modes

All modes need a bandwidth of 2400 Hz and there are two frame sizes.

These are the durations for each modulation scheme:

| Modulation | Short | Normal |
| ---------- | ----- | ------ |
| BPSK       |  1.5s |   2.6s |
| QPSK       |  1.0s |   2.6s |
| QAM16      |  1.0s |   2.6s |
| QAM64      |  1.9s |   3.4s |
| QAM256     |  1.5s |   2.6s |

These are the numbers for the short modes with a 1/2 code rate:

| Modulation | Payload | Bitrate |
| ---------- | ------- | ------- |
| BPSK       |    128B | 0.7kb/s |
| QPSK       |    128B | 1.1kb/s |
| QAM16      |    256B | 2.1kb/s |
| QAM64      |   1024B | 4.3kb/s |
| QAM256     |   1024B | 5.4kb/s |

These are the numbers for the normal modes with a 1/2 code rate:

| Modulation | Payload | Bitrate |
| ---------- | ------- | ------- |
| BPSK       |    256B | 0.8kb/s |
| QPSK       |    512B | 1.6kb/s |
| QAM16      |   1024B | 3.2kb/s |
| QAM64      |   2048B | 4.8kb/s |
| QAM256     |   2048B | 6.3kb/s |

These are the numbers for the short modes with an 11/16 code rate:

| Modulation | Payload | Bitrate |
| ---------- | ------- | ------- |
| BPSK       |    176B | 0.9kb/s |
| QPSK       |    176B | 1.5kb/s |
| QAM16      |    352B | 2.9kb/s |
| QAM64      |   1408B | 5.9kb/s |
| QAM256     |   1408B | 7.5kb/s |

These are the numbers for the normal modes with an 11/16 code rate:

| Modulation | Payload | Bitrate |
| ---------- | ------- | ------- |
| BPSK       |    352B | 1.1kb/s |
| QPSK       |    704B | 2.2kb/s |
| QAM16      |   1408B | 4.3kb/s |
| QAM64      |   2816B | 6.6kb/s |
| QAM256     |   2816B | 8.7kb/s |

These are the numbers for the short modes with a 3/4 code rate:

| Modulation | Payload | Bitrate |
| ---------- | ------- | ------- |
| BPSK       |    192B | 1.0kb/s |
| QPSK       |    192B | 1.6kb/s |
| QAM16      |    384B | 3.2kb/s |
| QAM64      |   1536B | 6.4kb/s |
| QAM256     |   1536B | 8.2kb/s |

These are the numbers for the normal modes with a 3/4 code rate:

| Modulation | Payload | Bitrate |
| ---------- | ------- | ------- |
| BPSK       |    384B | 1.2kb/s |
| QPSK       |    768B | 2.4kb/s |
| QAM16      |   1536B | 4.7kb/s |
| QAM64      |   3072B | 7.2kb/s |
| QAM256     |   3072B | 9.5kb/s |

These are the numbers for the short modes with a 13/16 code rate:

| Modulation | Payload | Bitrate |
| ---------- | ------- | ------- |
| BPSK       |    208B | 1.1kb/s |
| QPSK       |    208B | 1.7kb/s |
| QAM16      |    416B | 3.5kb/s |
| QAM64      |   1664B | 7.0kb/s |
| QAM256     |   1664B | 8.9kb/s |

These are the numbers for the normal modes with a 13/16 code rate:

| Modulation | Payload | Bitrate |
| ---------- | ------- | ------- |
| BPSK       |    416B | 1.3kb/s |
| QPSK       |    832B | 2.6kb/s |
| QAM16      |   1664B | 5.1kb/s |
| QAM64      |   3328B | 7.8kb/s |
| QAM256     |   3328B |10.3kb/s |

### Simulating

Prerequisite: [disorders](https://github.com/aicodix/disorders)

Encode ```uncoded.dat``` to audio signal, add [multipath](https://en.wikipedia.org/wiki/Multipath_propagation), [CFO, SFO](https://en.wikipedia.org/wiki/Carrier_frequency_offset), [AWGN](https://en.wikipedia.org/wiki/Additive_white_Gaussian_noise), decode and compare:

```
./encode - 48000 32 ANONYMOUS QAM16 1/2 short uncoded.dat | ../disorders/multipath - - ../disorders/multipath.txt 10 | ../disorders/cfo - - 234.567 | ../disorders/sfo - - 147 | ../disorders/awgn - - -20 | ./decode - - | diff -q -s uncoded.dat -
```

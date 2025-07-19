
### OFDM MODEM

Quick start:

Create file ```uncoded.dat``` with ```2048``` bits of random data:

```
dd if=/dev/urandom of=uncoded.dat bs=1 count=256
```

Encode file ```uncoded.dat``` to ```encoded.wav``` [WAV](https://en.wikipedia.org/wiki/WAV) audio file with 48000 Hz sample rate, 16 bits and only 1 (real) channel:

```
./encode encoded.wav 48000 16 1 1500 ANONYMOUS QAM16 1/2 short uncoded.dat
```

Start recording to ```recorded.wav``` audio file and stop after 5 seconds:

```
arecord -c 1 -f S16_LE -r 48000 -d 5 recorded.wav
```

Start playing ```encoded.wav``` audio file:

```
aplay encoded.wav
```

Decode ```recorded.wav``` audio file to ```decoded.dat``` file:

```
./decode recorded.wav decoded.dat
```

Compare original ```uncoded.dat``` with ```decoded.dat```:

```
diff -s uncoded.dat decoded.dat
```

### Supported Modes

All modes need a bandwidth of 2400 Hz.

These are the durations for each differential modulation scheme:

| Modulation | Short | Normal |
| ---------- | ----- | ------ |
| DBPSK      |  1.5s |   2.6s |
| DQPSK      |  1.0s |   2.6s |
| D8PSK      |  1.9s |   3.4s |

And these are the durations for each coherent modulation scheme:

| Modulation | Short | Normal |
| ---------- | ----- | ------ |
| QAM16      |  1.0s |   2.6s |
| QAM64      |  1.9s |   3.4s |
| QAM256     |  1.5s |   2.6s |
| QAM1024    |  2.2s |   4.0s |
| QAM4096    |  1.9s |   3.4s |

Currently all digital modes use a 1/2-rate forward error correction code.

Therefore we have the following byte payloads for each modulation:

| Modulation | Short | Normal |
| ---------- | ----- | ------ |
| DBPSK      |  128B |   256B |
| DQPSK      |  128B |   512B |
| D8PSK      |  512B |  1024B |
| QAM16      |  256B |  1024B |
| QAM64      | 1024B |  2048B |
| QAM256     | 1024B |  2048B |
| QAM1024    | 2048B |  4096B |
| QAM4096    | 2048B |  4096B |

Which give us the following bitrates for each modulation:

| Modulation |   Short |  Normal |
| ---------- | ------- | ------- |
| DBPSK      |  681b/s |  789b/s |
| DQPSK      | 1070b/s | 1577b/s |
| D8PSK      | 2141b/s | 2398b/s |
| QAM16      | 2141b/s | 3155b/s |
| QAM64      | 4282b/s | 4795b/s |
| QAM256     | 5449b/s | 6310b/s |
| QAM1024    | 7493b/s | 8268b/s |
| QAM4096    | 8563b/s | 9591b/s |

### Simulating

Prerequisite: [disorders](https://github.com/aicodix/disorders)

Encode ```uncoded.dat``` to [analytic](https://en.wikipedia.org/wiki/Analytic_signal) audio signal, add [multipath](https://en.wikipedia.org/wiki/Multipath_propagation), [CFO, SFO](https://en.wikipedia.org/wiki/Carrier_frequency_offset), [AWGN](https://en.wikipedia.org/wiki/Additive_white_Gaussian_noise), decode and compare:

```
./encode - 48000 16 2 1500 ANONYMOUS QAM16 1/2 short uncoded.dat | ../disorders/multipath - - ../disorders/multipath.txt 10 | ../disorders/cfo - - 234.567 | ../disorders/sfo - - 147 | ../disorders/awgn - - -20 | ./decode - - | diff -q -s uncoded.dat -
```

### Reading

* Robust frequency and timing synchronization for OFDM  
by Timothy M. Schmidl and Donald C. Cox - 1997
* On Timing Offset Estimation for OFDM Systems  
by H. Minn, M. Zeng, and V. K. Bhargava - 2000


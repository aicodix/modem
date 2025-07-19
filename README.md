
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

All modes need a bandwidth of 2400 Hz and there are two frame sizes.

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

These are the numbers for the short modes with a 1/2 code rate:

| Modulation | Payload | Bitrate |
| ---------- | ------- | ------- |
| DBPSK      |    128B | 0.7kb/s |
| DQPSK      |    128B | 1.1kb/s |
| D8PSK      |    512B | 2.1kb/s |
| QAM16      |    256B | 2.1kb/s |
| QAM64      |   1024B | 4.3kb/s |
| QAM256     |   1024B | 5.4kb/s |
| QAM1024    |   2048B | 7.5kb/s |
| QAM4096    |   2048B | 8.6kb/s |

These are the numbers for the normal modes with a 1/2 code rate:

| Modulation | Payload | Bitrate |
| ---------- | ------- | ------- |
| DBPSK      |    256B | 0.8kb/s |
| DQPSK      |    512B | 1.6kb/s |
| D8PSK      |   1024B | 2.4kb/s |
| QAM16      |   1024B | 3.2kb/s |
| QAM64      |   2048B | 4.8kb/s |
| QAM256     |   2048B | 6.3kb/s |
| QAM1024    |   4096B | 8.3kb/s |
| QAM4096    |   4096B | 9.6kb/s |

These are the numbers for the short modes with a 2/3 code rate:

| Modulation | Payload | Bitrate |
| ---------- | ------- | ------- |
| DBPSK      |    171B | 0.9kb/s |
| DQPSK      |    171B | 1.4kb/s |
| D8PSK      |    684B | 2.9kb/s |
| QAM16      |    342B | 2.9kb/s |
| QAM64      |   1368B | 5.7kb/s |
| QAM256     |   1368B | 7.3kb/s |
| QAM1024    |   2736B |10.0kb/s |
| QAM4096    |   2736B |11.4kb/s |

These are the numbers for the normal modes with a 2/3 code rate:

| Modulation | Payload | Bitrate |
| ---------- | ------- | ------- |
| DBPSK      |    342B | 1.1kb/s |
| DQPSK      |    684B | 2.1kb/s |
| D8PSK      |   1368B | 3.2kb/s |
| QAM16      |   1368B | 4.2kb/s |
| QAM64      |   2736B | 6.4kb/s |
| QAM256     |   2736B | 8.4kb/s |
| QAM1024    |   5472B |11.0kb/s |
| QAM4096    |   5472B |12.8kb/s |

These are the numbers for the short modes with a 3/4 code rate:

| Modulation | Payload | Bitrate |
| ---------- | ------- | ------- |
| DBPSK      |    192B | 1.0kb/s |
| DQPSK      |    192B | 1.6kb/s |
| D8PSK      |    768B | 3.2kb/s |
| QAM16      |    384B | 3.2kb/s |
| QAM64      |   1536B | 6.4kb/s |
| QAM256     |   1536B | 8.2kb/s |
| QAM1024    |   3072B |11.2kb/s |
| QAM4096    |   3072B |12.8kb/s |

These are the numbers for the normal modes with a 3/4 code rate:

| Modulation | Payload | Bitrate |
| ---------- | ------- | ------- |
| DBPSK      |    384B | 1.2kb/s |
| DQPSK      |    768B | 2.4kb/s |
| D8PSK      |   1536B | 3.6kb/s |
| QAM16      |   1536B | 4.7kb/s |
| QAM64      |   3072B | 7.2kb/s |
| QAM256     |   3072B | 9.5kb/s |
| QAM1024    |   6144B |12.4kb/s |
| QAM4096    |   6144B |14.4kb/s |

These are the numbers for the short modes with a 5/6 code rate:

| Modulation | Payload | Bitrate |
| ---------- | ------- | ------- |
| DBPSK      |    213B | 1.1kb/s |
| DQPSK      |    213B | 1.8kb/s |
| D8PSK      |    852B | 3.6kb/s |
| QAM16      |    426B | 3.6kb/s |
| QAM64      |   1704B | 7.1kb/s |
| QAM256     |   1704B | 9.1kb/s |
| QAM1024    |   3408B |12.5kb/s |
| QAM4096    |   3408B |14.2kb/s |

These are the numbers for the normal modes with a 5/6 code rate:

| Modulation | Payload | Bitrate |
| ---------- | ------- | ------- |
| DBPSK      |    426B | 1.3kb/s |
| DQPSK      |    852B | 2.6kb/s |
| D8PSK      |   1704B | 4.0kb/s |
| QAM16      |   1704B | 5.2kb/s |
| QAM64      |   3408B | 8.0kb/s |
| QAM256     |   3408B |10.5kb/s |
| QAM1024    |   6816B |13.8kb/s |
| QAM4096    |   6816B |16.0kb/s |

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


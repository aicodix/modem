
### OFDM MODEM

Quick start:

Create file ```uncoded.dat``` with ```2048``` bits of random data:

```
dd if=/dev/urandom of=uncoded.dat bs=1 count=256
```

Encode file ```uncoded.dat``` to ```encoded.wav``` [WAV](https://en.wikipedia.org/wiki/WAV) audio file with 8000 Hz sample rate, 16 bits and only 1 (real) channel:

```
./encode encoded.wav 8000 16 1 1500 5 uncoded.dat
```

Start recording to ```recorded.wav``` audio file and stop after 5 seconds:

```
arecord -c 1 -f S16_LE -r 8000 -d 5 recorded.wav
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

All modes need a bandwidth of 2000 Hz and use a 1/2-rate forward error correction code

These use a differential modulation scheme:
* Mode 0: DBPSK, 1.70 seconds and 128 bytes
* Mode 1: DBPSK, 3.06 seconds and 256 bytes
* Mode 2: DBPSK, 5.78 seconds and 512 bytes
* Mode 3: DQPSK, 1.02 seconds and 128 bytes
* Mode 4: DQPSK, 1.70 seconds and 256 bytes
* Mode 5: DQPSK, 3.06 seconds and 512 bytes
* Mode 6: DQPSK, 5.78 seconds and 1024 bytes
* Mode 7: D8PSK, 2.21 seconds and 512 bytes
* Mode 8: D8PSK, 4.08 seconds and 1024 bytes
* Mode 9: D8PSK, 7.82 seconds and 2048 bytes

And these a coherent modulation scheme:
* Mode 10: BPSK, 1.70 seconds and 128 bytes
* Mode 11: BPSK, 3.06 seconds and 256 bytes
* Mode 12: BPSK, 5.78 seconds and 512 bytes
* Mode 13: QPSK, 1.02 seconds and 128 bytes
* Mode 14: QPSK, 1.70 seconds and 256 bytes
* Mode 15: QPSK, 3.06 seconds and 512 bytes
* Mode 16: QPSK, 5.78 seconds and 1024 bytes
* Mode 17: QAM16, 1.02 seconds and 256 bytes
* Mode 18: QAM16, 1.70 seconds and 512 bytes
* Mode 19: QAM16, 3.06 seconds and 1024 bytes
* Mode 20: QAM16, 5.78 seconds and 2048 bytes
* Mode 21: QAM64, 2.21 seconds and 1024 bytes
* Mode 22: QAM64, 4.08 seconds and 2048 bytes
* Mode 23: QAM64, 7.82 seconds and 4096 bytes
* Mode 24: QAM256, 1.02 seconds and 512 bytes
* Mode 25: QAM256, 1.70 seconds and 1024 bytes
* Mode 26: QAM256, 3.06 seconds and 2048 bytes
* Mode 27: QAM256, 5.78 seconds and 4096 bytes

### Simulating

Prerequisite: [disorders](https://github.com/aicodix/disorders)

Encode ```uncoded.dat``` to [analytic](https://en.wikipedia.org/wiki/Analytic_signal) audio signal, add [multipath](https://en.wikipedia.org/wiki/Multipath_propagation), [CFO, SFO](https://en.wikipedia.org/wiki/Carrier_frequency_offset), [AWGN](https://en.wikipedia.org/wiki/Additive_white_Gaussian_noise), decode and compare:

```
./encode - 8000 16 2 1500 5 uncoded.dat | ../disorders/multipath - - ../disorders/multipath.txt 10 | ../disorders/cfo - - 234.567 | ../disorders/sfo - - 147 | ../disorders/awgn - - -30 | ./decode - - | diff -q -s uncoded.dat -
```

### Reading

* Robust frequency and timing synchronization for OFDM  
by Timothy M. Schmidl and Donald C. Cox - 1997
* On Timing Offset Estimation for OFDM Systems  
by H. Minn, M. Zeng, and V. K. Bhargava - 2000



### OFDM MODEM

Quick start:

Create file ```uncoded.dat``` with ```5440``` bits of random data:

```
dd if=/dev/urandom of=uncoded.dat bs=1 count=680
```

Encode file ```uncoded.dat``` to ```encoded.wav``` [WAV](https://en.wikipedia.org/wiki/WAV) audio file with 8000 Hz sample rate, 16 bits and only 1 (real) channel:

```
./encode encoded.wav 8000 16 1 uncoded.dat
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
./decode decoded.dat recorded.wav
```

Compare original ```uncoded.dat``` with ```decoded.dat```:

```
diff -s uncoded.dat decoded.dat
```

### Simulating

Prerequisite: [disorders](https://github.com/aicodix/disorders)

Encode ```uncoded.dat``` to [analytic](https://en.wikipedia.org/wiki/Analytic_signal) audio signal, add [multipath](https://en.wikipedia.org/wiki/Multipath_propagation), [CFO, SFO](https://en.wikipedia.org/wiki/Carrier_frequency_offset), [AWGN](https://en.wikipedia.org/wiki/Additive_white_Gaussian_noise), decode and compare:

```
./encode - 8000 16 2 uncoded.dat 2000 | ../disorders/multipath - - ../disorders/multipath.txt 10 | ../disorders/cfo - - 234.567 | ../disorders/sfo - - 147 | ../disorders/awgn - - -33 | ./decode - - | diff -q -s uncoded.dat -
```

### Reading

* Robust frequency and timing synchronization for OFDM  
by Timothy M. Schmidl and Donald C. Cox - 1997
* On Timing Offset Estimation for OFDM Systems  
by H. Minn, M. Zeng, and V. K. Bhargava - 2000


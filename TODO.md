### TODO

* Find the error floor of about 40 dB in the decoder: on a clean loopback the real Es/N0 of the symbols of a frame spreads from about 37 to 47 dB, with or without the clipping of the encoder and the averaging of the channel estimate, so the precision jumps from symbol to symbol. It only shows above about 37 dB SNR and doesn't hurt decoding

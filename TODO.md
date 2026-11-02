### TODO
* implement analog mode (most significant mode bit is 1), where only complex values are the payload
* compute SNR per subchannel using sync and preamble and update with each symbol
* re-encode entire message and output number of bitflips and real SNR
* bandpass 300..2700 Hz (or lowpass 1200 Hz after mixing down by 1500 Hz) only the correlator input, with relaxed edges leaving room for the CFO, a FIR filter keeps the timing aligned, or add band limited noise to the awgn tool of disorders


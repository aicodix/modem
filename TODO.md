### TODO
* implement analog mode (most significant mode bit is 1), where only complex values are the payload
* re-encode preamble to use entire symbol to update the channel estimation
* compute SNR per subchannel using sync and preamble and update with each symbol
* re-encode entire message and output number of bitflips and real SNR


### TODO
* refine Bhattacharyya parameters for each payload code before freezing the specification
* erase symbol or abort if either head or tail data is damaged, as we can't recover anything when that happens
* implement analog mode (most significant mode bit is 1), where only complex values are the payload
* re-encode preamble to use entire symbol to update the channel estimation
* re-encode entire message and output number of bitflips and real SNR


### ROOFTOP

Ideas for the successor of this modem: a decentralized BBS for local
communities, over line of sight FM radio on VHF and UHF, with repeaters on
mountain tops. HYDRA is one of its components. Not HF and not SSB anymore.

### The channel

FM hands over audio, not a carrier, so there is practically no frequency
offset left. What matters instead:

* Key up time: PTT, VOX and squelch cut off the start of a transmission.
* Pre-emphasis and de-emphasis tilt the band, or don't, depending on whether
  the radio is fed through its microphone or a flat data port. Without
  de-emphasis, discriminator noise rises with frequency, about 20 dB across
  the band.
* The deviation limiter clips peaks, so PAPR matters even more.
* Repeaters pass the signal through two FM links and their own audio chain:
  CTCSS filters, limiters, emphasis, hang time and delay.
* Sound card clocks still need SFO correction.

The first step is an FM channel model for the simulations: emphasis, noise
that rises with frequency, a limiter, a cut off start, a squelch tail and a
repeater chain.

### Remove

* Analog mode.
* QAM1024 and QAM4096.
* The Schmidl & Cox correlator.
* The large CFO search and its ±250 Hz limit, which only SSB needs.

### Sync: the lead in is the sync

Send k copies of one known pseudo noise symbol, as the leading noise symbol
already is, then the metadata symbol, then the payload. No separate sync pair.

* Detect with a frequency domain correlator that hops by one symbol length:
  one FFT, a multiply with the known spectrum and one inverse FFT per hop.
  Two surviving copies are enough, as a window of one symbol length then lies
  fully inside them, and the peak gives the timing within a symbol.
* Find the end, not the beginning: correlate each following symbol aligned
  window with the known symbol. The copies give a high value, the metadata
  symbol about 1/sqrt(320) of it. The first low value is the metadata symbol.
* It does not matter how many copies the radio cuts off, as long as two
  survive. Choose k from real radios, likely 3 to 5.
* The copies also give SFO and residual CFO from their phase changes, an
  averaged first channel estimate and a per tone noise estimate from their
  differences, before the metadata symbol.

To measure first: end detection and timing at low SNR, with 0, 1 or 2 copies
cut off.

### Payload sizes

Use the nested frozen bit tables of code/tools/nested.cc instead of four code
rates and two frame sizes. One byte per bit position holds all rates: bit i is
frozen for a payload of 256 * j bits, if threshold[i] >= j. Each rate keeps its
own design erasure probability, chosen for the same union bound, and the sets
are still nested: `./nested 8 32` reports no violations for N = 65536.

* Steps of 256 bits are coarse for short frames. Either one table per code
  order, or accept the padding. It does not need to be finer than that.
* Air time should follow the payload, so short frames like acks are short.
* The metadata needs the payload size instead of code rate and frame size.

### Per tone SNR

Today each tone's LLR precision is the SNR of its symbol, weighted by its
channel power. That is only right for white noise: a weak tone can be almost
perfect and a strong one can sit under a whistle.

Measured on next, decoding with a running per tone estimate: squared errors of
the preamble, the seed tones and hard decisions on the payload, summed per tone
from the preamble on and averaged over 16 neighbouring tones. Today's estimate
counts as 64 observations, so it dominates while a tone has few of its own:

    noise = (64 * predicted + n * measured) / (64 + n)

* Noise rising with frequency, as on a flat data port: the point where half the
  frames decode moves down by about 1.2 dB in BPSK, QPSK, QAM16 and QAM64.
* With multipath on top: about 0.5 dB. With a narrow interferer: 0.3 to 0.5 dB.
* White noise and room multipath: about the same as today.
* Each tone's own estimate alone, without today's estimate mixed in, loses up
  to 1 dB on white noise and multipath, mostly with short QAM frames.

Shorter symbols are not needed for this: averaging over neighbouring tones
gives the same number of observations, and the measurements prefer averaging
over 8 to 16 tones, while 30 Hz tones would fix it at about 4.

### Sync limit of next

With the position check after the correlator widened to two guard intervals,
BPSK 1/2 short decodes all frames from about 3 dB in band SNR, also with
multipath, CFO and SFO. The sync is no longer what fails first.

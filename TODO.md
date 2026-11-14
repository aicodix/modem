### TODO

* Measure thresholds of the candidate pairs of modulation and code rate on static multipath and on a Watterson channel, and replace the draft table `efficiency` with them:
  * QAM16 12/16 against QAM64 8/16 at 3 bits per tone
  * BPSK against QPSK at the low rates, where phase noise matters
  * whether the steps of BPSK at the lowest rates, 3 dB from 1/16 to 2/16, are fine enough
* Time the decoding of long frames near threshold in the browser, where up to eight codewords decode one after another at the end of a frame
* Remove the unused `Common::mod_distance()`
* `polar_freezer.hh` in `code` needs `#include <cstdint>` to build `freezer.cc` with g++ and libstdc++

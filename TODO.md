### TODO

* Replace the draft table `efficiency` with the measured pairs in CLAUDE.md, 29 of them, and decide what the other three indices carry
* Find out why a Watterson channel with 1 Hz spread leaves a real Es/N0 of about 0 dB even without noise, and whether faster tracking of the channel helps
* Time the decoding of long frames near threshold in the browser, where up to eight codewords decode one after another at the end of a frame
* Remove the unused `Common::mod_distance()`
* `polar_freezer.hh` in `code` needs `#include <cstdint>` to build `freezer.cc` with g++ and libstdc++

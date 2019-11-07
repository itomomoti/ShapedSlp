#include "Common.hpp"

void padVLine
(
 const uint64_t pad
 ) {
  for (uint64_t i = 0; i < pad; ++i) {
    std::cout << "|";
  }
}


uint32_t ceilLog2(uint64_t x) {
  if (x == 0) {
    return 1;
  }
  return 64 - __builtin_clzll(x);
}

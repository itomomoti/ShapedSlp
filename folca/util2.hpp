#ifndef INCLUDE_GUARD_UTIL_HPP_
#define INCLUDE_GUARD_UTIL_HPP_

#include <stdint.h>
#include <limits.h>

namespace comp {

  typedef uint64_t var_t;
  typedef uint64_t len_t;
  const uint64_t kDummy = (18446744073709551615ULL);
  const uint64_t kAlphabetSize = 256;
  const uint64_t kOne = 0x0000000000000001ULL;
  const bool kOP = 1;
  const bool kCP = 0;

  static const uint64_t PRIMES[] = {
    /* 0*/  8 + 3,
    /* 1*/  16 + 3,
    /* 2*/  32 + 5,
    /* 3*/  64 + 3,
    /* 4*/  128 + 3,
    /* 5*/  256 + 27,
    /* 6*/  512 + 9,
    /* 7*/  1024 + 9,
    /* 8*/  2048 + 5,
    /* 9*/  4096 + 3,
    /*10*/  8192 + 27,
    /*11*/  16384 + 43,
    /*12*/  32768 + 3,
    /*13*/  65536 + 45,
    /*14*/  131072 + 29,
    /*15*/  262144 + 3,
    /*16*/  524288 + 21,
    /*17*/  1048576 + 7,
    /*18*/  2097152 + 17,
    /*19*/  4194304 + 15,
    /*20*/  8388608 + 9,
    /*21*/  16777216 + 43,
    /*22*/  33554432 + 35,
    /*23*/  67108864 + 15,
    /*24*/  134217728 + 29,
    /*25*/  268435456 + 3,
    /*26*/  536870912 + 11,
    /*27*/  1073741824 + 85,
  };

} // namespace comp

#endif // INCLUDE_GUARD_UTIL_HPP_

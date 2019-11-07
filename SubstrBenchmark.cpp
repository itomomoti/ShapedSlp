#include <iostream>
#include <string>
#include <queue>
#include "cmdline.h"
#include "Common.hpp"
#include "PoSlp.hpp"
#include "ShapedSlp.hpp"
#include "DirectAccessibleGammaCode.hpp"
#include "DirectAccessibleGammaCodeWithRrr.hpp"
#include "VlcVec.hpp"

using namespace std;

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;


template<class SlpT>
void measureExpand
(
 SlpT slp,
 std::string in,
 const uint64_t numItr,
 const uint64_t lenExpand,
 const uint64_t firstPos,
 const uint64_t jump,
 const bool dummy_flag
) {
  auto start = timer::now();
  ifstream fs(in);
  slp.load(fs);
  auto stop = timer::now();
  cout << "time to load (ms): " << duration_cast<milliseconds>(stop-start).count() << endl;
  // slp.printStatus();

  cout << "numItr = " << numItr << ", lenExpand = " << lenExpand << ", firstPos = " << firstPos << ", jump = " << jump << endl;
  const uint64_t textLen = slp.getLen();
  if (firstPos >= textLen or lenExpand >= textLen or lenExpand >= textLen) {
    cout << "firstPos, lenExpand and jump should be smaller than text length, which is " << textLen << endl;
    exit(1);
  }

  string substr;
  substr.resize(lenExpand);
  start = timer::now();
  uint64_t beg = firstPos;
  for (uint64_t i = 0; i < numItr; ++i) {
    // substr[0] = slp.charAt(beg);
    slp.expandSubstr(beg, lenExpand, substr.data());
    beg += jump;
    if (beg > textLen - lenExpand) {
      beg -= (textLen - lenExpand);
    }
  }
  stop = timer::now();
  cout << "time to expand (micro sec per query): " << (double)duration_cast<microseconds>(stop-start).count() / numItr << endl;

  if (dummy_flag) {
    cout << substr << endl;
  }
}


int main(int argc, char* argv[])
{
  cmdline::parser parser;
  parser.add<std::string>("input", 'i', "input file name in which ShapedSlp data structure is written.", true);
  parser.add<uint32_t>("encoding", 'e', "encoding: [0] PoSlp. [1] ShapedSlp. [2] ShapedSlpDagc", true);
  parser.add<uint64_t>("numItr", 'n', "number of iterations", true);
  parser.add<uint64_t>("lenExpand", 'l', "length to expand", true);
  parser.add<uint64_t>("firstPos", 'f', "first position to access", false, 0);
  parser.add<uint64_t>("jump", 'j', "amount of jump when determining the next position to access", false, 38201); // default 38201 is a prime number
  parser.add<bool>("dummy_flag", 0, "this is dummy flag to prevent that optimization deletes codes", false, false);
  parser.parse_check(argc, argv);
  const std::string in = parser.get<std::string>("input");
  const uint32_t encoding = parser.get<uint32_t>("encoding");
  const uint64_t numItr = parser.get<uint64_t>("numItr");
  const uint64_t lenExpand = parser.get<uint64_t>("lenExpand");
  const uint64_t firstPos = parser.get<uint64_t>("firstPos");
  const uint64_t jump = parser.get<uint64_t>("jump");
  const bool dummy_flag = parser.get<bool>("dummy_flag");

  if (numItr == 0) {
    cout << "numItr should be > 0." << endl;
    exit(1);
  }

  using var_t = uint32_t;

  // { // correctness check
  //   PoSlp<var_t> poslp;
  //   ShapedSlp<var_t> sslp;
  //   ShapedSlpDagc<var_t> sslpdagc;
  //   {
  //     ifstream fs("temp0");
  //     poslp.load(fs);
  //   }
  //   {
  //     ifstream fs("temp1");
  //     sslp.load(fs);
  //   }
  //   {
  //     ifstream fs("temp2");
  //     sslpdagc.load(fs);
  //   }

  //   cout << "numItr = " << numItr << ", lenExpand = " << lenExpand << ", firstPos = " << firstPos << ", jump = " << jump << endl;
  //   const uint64_t textLen = poslp.getLen();
  //   if (firstPos >= textLen or lenExpand >= textLen or lenExpand >= textLen) {
  //     cout << "firstPos, lenExpand and jump should be smaller than text length, which is " << textLen << endl;
  //     exit(1);
  //   }

  //   string substr0, substr1, substr2;
  //   substr0.resize(lenExpand);
  //   substr1.resize(lenExpand);
  //   substr2.resize(lenExpand);
  //   uint64_t beg = firstPos;
  //   for (uint64_t i = 0; i < numItr; ++i) {
  //     cout << beg << endl;
  //     // substr0[0] = poslp.charAt(beg);
  //     // substr1[0] = sslp.charAt(beg);
  //     // substr2[0] = sslpdagc.charAt(beg);
  //     poslp.expandSubstr(beg, lenExpand, substr0.data());
  //     sslp.expandSubstr(beg, lenExpand, substr1.data());
  //     sslpdagc.expandSubstr(beg, lenExpand, substr2.data());
  //     if ((substr0 != substr1) or (substr1 != substr2) or (substr0 != substr2)) {
  //       cout << "wrong: beg = " << beg << endl;
  //       cout << substr0 << endl;
  //       cout << substr1 << endl;
  //       cout << substr2 << endl;
  //       exit(1);
  //     }
  //     beg += jump;
  //     if (beg > textLen - lenExpand) {
  //       beg -= (textLen - lenExpand);
  //     }
  //   }
  // }

  if (encoding == 0) {
    PoSlp<var_t> poslp;
    measureExpand(poslp, in, numItr, lenExpand, firstPos, jump, dummy_flag);
  } else if (encoding == 1) {
    ShapedSlp<var_t, DirectAccessibleGammaCode> sslp;
    measureExpand(sslp, in, numItr, lenExpand, firstPos, jump, dummy_flag);
  } else if (encoding == 2) {
    ShapedSlp<var_t, DirectAccessibleGammaCodeWithRrr> sslp;
    measureExpand(sslp, in, numItr, lenExpand, firstPos, jump, dummy_flag);
  } else if (encoding == 3) {
    ShapedSlp<var_t, VlcVec<sdsl::coder::elias_delta, 64> > sslp;
    measureExpand(sslp, in, numItr, lenExpand, firstPos, jump, dummy_flag);
  } else if (encoding == 4) {
    ShapedSlp<var_t, VlcVec<sdsl::coder::elias_delta, 128> > sslp;
    measureExpand(sslp, in, numItr, lenExpand, firstPos, jump, dummy_flag);
  }

  return 0;
}

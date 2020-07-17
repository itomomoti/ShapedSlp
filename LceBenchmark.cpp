#include <iostream>
#include <string>
#include <queue>
#include <thread>
#include <random>
#include "cmdline.h"
#include "Common.hpp"
#include "PlainSlp.hpp"
#include "PoSlp.hpp"
#include "ShapedSlp_Status.hpp"
#include "ShapedSlp.hpp"
#include "ShapedSlpV2.hpp"
#include "SelfShapedSlp.hpp"
#include "SelfShapedSlpV2.hpp"
#include "DirectAccessibleGammaCode.hpp"
#include "IncBitLenCode.hpp"
#include "FixedBitLenCode.hpp"
#include "SelectType.hpp"
#include "VlcVec.hpp"

using namespace std;

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

using var_t = uint32_t;


template<class SlpT>
uint64_t naiveLceToR
(
 const SlpT & slp,
 uint64_t p1,
 uint64_t p2
) {
  if (p1 > p2) {
    std::swap(p1, p2);
  }
  // Now p1 \le p2
  const uint64_t rem = p2;
  for ( ; p2 < slp.getLen(); ++p1, ++p2) {
    // std::cout << "naive: p1 = " << p1 << ": " << slp.charAt(p1) << ", p2 = " << p2 << ": " << slp.charAt(p2) << std::endl;
    if (slp.charAt(p1) != slp.charAt(p2)) {
      break;
    }
  }

  return p2 - rem;
}


template<class SlpT>
void measure
(
 std::string in,
 const uint64_t numItr,
 const uint64_t givenSeed
) {
  SlpT slp;

  auto start = timer::now();
  ifstream fs(in);
  slp.load(fs);
  auto stop = timer::now();
  cout << "time to load (ms): " << duration_cast<milliseconds>(stop-start).count() << endl;
  // slp.printStatus();

  random_device rnd;
  const uint64_t seed = (givenSeed) ? givenSeed : rnd(); // if givenSeed == 0, choose it randomly by rnd()
  const uint64_t textLen = slp.getLen();
  cout << "numItr = " << numItr << ", text len = " << textLen << ", seed = " << seed << endl;
  if (numItr == 0) {
    cout << "numItr should be > 0." << endl;
    exit(1);
  }

  // {//debug
  //   printDerivationTree(slp);
  // }
  const uint64_t numLoop = 11;
  uint64_t checksum0 = 0;
  uint64_t checksum1 = 0;

  {
    std::vector<double> times(numLoop);
    for (uint64_t loop = 0; loop < numLoop; ++loop) {
      uniform_int_distribution<uint64_t> rndUniform(0, textLen - 1);
      mt19937_64 mt(seed);
      start = timer::now();
      for (uint64_t i = 0; i < numItr; ++i) {
        const uint64_t p1 = rndUniform(mt);
        const uint64_t p2 = rndUniform(mt);
        checksum0 += lceToR(slp, p1, p2);
      }
      stop = timer::now();
      times[loop] = (double)duration_cast<microseconds>(stop-start).count() / numItr;
    }
    std::sort(times.begin(), times.end());
    cout << "time to lce queries (micro sec per query): " << times[numLoop / 2] << endl;
  }

  {
    std::vector<double> times(numLoop);
    for (uint64_t loop = 0; loop < numLoop; ++loop) {
      uniform_int_distribution<uint64_t> rndUniform(0, textLen - 1);
      mt19937_64 mt(seed);
      start = timer::now();
      for (uint64_t i = 0; i < numItr; ++i) {
        const uint64_t p1 = rndUniform(mt);
        const uint64_t p2 = rndUniform(mt);
        checksum1 += lceToR(slp, p1, p2);
      }
      stop = timer::now();
      times[loop] = (double)duration_cast<microseconds>(stop-start).count() / numItr;
    }
    std::sort(times.begin(), times.end());
    cout << "time to naive lce queries (micro sec per query): " << times[numLoop / 2] << endl;
  }

  cout << "averagec LCE length = " << checksum0 / (numItr * numLoop) << endl;
  if (checksum0 == checksum1) {
    cout << "checksum match: checksum0 = " << checksum0 << ", checksum1 = " << checksum1 << endl;
  } else {
    cout << "checksum ERROR: checksum0 = " << checksum0 << ", checksum1 = " << checksum1 << endl;
  }
}


int main(int argc, char* argv[])
{
  using Fblc = FixedBitLenCode<>;
  using SelSd = SelectSdvec<>;
  using SelMcl = SelectMcl<>;
  using DagcSd = DirectAccessibleGammaCode<SelSd>;
  using DagcMcl = DirectAccessibleGammaCode<SelMcl>;
  using Vlc64 = VlcVec<sdsl::coder::elias_delta, 64>;
  using Vlc128 = VlcVec<sdsl::coder::elias_delta, 128>;
  using funcs_type = map<string,
                         void(*)
                         (
                          std::string in,
                          const uint64_t numItr,
                          const uint64_t givenSeed
                          )>;
  funcs_type funcs;

  //// PlainSlp
  funcs.insert(make_pair("PlainSlp_FblcFblc", measure<PlainSlp<var_t, Fblc, Fblc>>));
  funcs.insert(make_pair("PlainSlp_IblcFblc", measure<PlainSlp<var_t, IncBitLenCode, Fblc>>));
  funcs.insert(make_pair("PlainSlp_32Fblc", measure<PlainSlp<var_t, FixedBitLenCode<32>, Fblc>>));

  //// PoSlp: Post-order SLP
  //// Sometimes PoSlp_Sd is better than PoSlp_Iblc
  funcs.insert(make_pair("PoSlp_Iblc", measure<PoSlp<var_t, IncBitLenCode>>));
  funcs.insert(make_pair("PoSlp_Sd", measure<PoSlp<var_t, DagcSd>>));

  //// ShapedSlp: plain implementation of slp encoding that utilizes shape-tree grammar
  //// Since bit length to represent slp element is small, SelMcl is good for them.
  //// For stg and bal element, SelSd is better
  // funcs.insert(make_pair("ShapedSlp_SdMclSd_SdMcl", measure<ShapedSlp<var_t, DagcSd, DagcMcl, DagcSd, SelSd, SelMcl>>));
  // funcs.insert(make_pair("ShapedSlp_SdSdSd_SdMcl", measure<ShapedSlp<var_t, DagcSd, DagcSd, DagcSd, SelSd, SelMcl>>));

  //// ShapedSlpV2: all vlc vectors are merged into one.
  //// Generally encoding size gets worse than ShapedSlp_SdMclSd_SdMcl because
  //// - Since bit length to represnet stg and bal element is large, DagcSd is a good choice.
  //// - On the other hand, bit size to represent slp element is significantly small, and so SelMcl should be used
  // funcs.insert(make_pair("ShapedSlpV2_Sd_SdMcl", measure<ShapedSlpV2<var_t, DagcSd, SelSd, SelMcl>>));

  //// SelfShapedSlp: ShapedSlp that does not use shape-tree grammar
  funcs.insert(make_pair("SelfShapedSlp_SdSd_Sd", measure<SelfShapedSlp<var_t, DagcSd, DagcSd, SelSd>>));
  funcs.insert(make_pair("SelfShapedSlp_SdSd_Mcl", measure<SelfShapedSlp<var_t, DagcSd, DagcSd, SelMcl>>));

  //// SelfShapedSlpV2:
  //// attempted to asign smaller offsets to frequent variables by giving special seats for hi-frequent ones
  // funcs.insert(make_pair("SelfShapedSlpV2_SdSd_Sd", measure<SelfShapedSlpV2<var_t, DagcSd, DagcSd, SelSd>>));

  string methodList;
  for (auto itr = funcs.begin(); itr != funcs.end(); ++itr) {
    methodList += itr->first + ". ";
  }


  cmdline::parser parser;
  parser.add<string>("input", 'i', "input file name in which ShapedSlp data structure is written.", true);
  parser.add<string>("encoding", 'e', "encoding: " + methodList, true);
  parser.add<uint64_t>("numItr", 'n', "number of iterations", true);
  parser.add<uint64_t>("seed", 's', "seed for random function", false, 0);
  parser.add<uint64_t>("numThreads", 't', "number of threads", false, 0);
  parser.add<bool>("dummy_flag", 0, "this is dummy flag to prevent that optimization deletes codes", false, false);
  parser.parse_check(argc, argv);
  const string in = parser.get<string>("input");
  const string encoding = parser.get<string>("encoding");
  const uint64_t numItr = parser.get<uint64_t>("numItr");
  const uint64_t seed = parser.get<uint64_t>("seed");
  const uint64_t nt = parser.get<uint64_t>("numThreads");
  const bool dummy_flag = parser.get<bool>("dummy_flag");

  if (numItr == 0) {
    cout << "numItr should be > 0." << endl;
    exit(1);
  }


  if (encoding.compare("All") == 0) {
    for (auto itr = funcs.begin(); itr != funcs.end(); ++itr) {
      cout << itr->first << ": BEGIN" << std::endl;
      itr->second(in + itr->first, numItr, seed);
      cout << itr->first << ": END" << std::endl;
    }
  } else {
    auto itr = funcs.find(encoding);
    if (itr != funcs.end()) {
      cout << itr->first << ": BEGIN" << std::endl;
      itr->second(in, numItr, seed);
      cout << itr->first << ": END" << std::endl;
    } else {
      cerr << "error: specify a valid encoding name in " + methodList << endl;
      exit(1);
    }
  }

  return 0;
}

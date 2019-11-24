#include <iostream>
#include <string>
#include <queue>
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
#include "SelectType.hpp"
#include "VlcVec.hpp"


using namespace std;

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

using var_t = uint32_t;


template
<
  class SlpT
  >
void measure
(
 const NaiveSlp<var_t> & slp,
 std::string out
) {
  auto start = timer::now();
  SlpT sslp(slp);
  auto stop = timer::now();
  cout << "time to encode (ms): " << duration_cast<milliseconds>(stop-start).count() << endl;
  sslp.printStatus();

  start = timer::now();
  ofstream fs(out);
  sslp.serialize(fs);
  stop = timer::now();
  cout << "time to serialize (ms): " << duration_cast<milliseconds>(stop-start).count() << endl;
}


template
<
  class SlpT
  >
void measure_PoSlp
(
 const NaiveSlp<var_t> & slp,
 std::string out
) {
  NaiveSlp<var_t> temp(slp);
  auto start = timer::now();
  temp.makeBinaryTree();
  SlpT poslp;
  poslp.makePoSlpWithLen(temp);
  auto stop = timer::now();
  cout << "time to encode (ms): " << duration_cast<milliseconds>(stop-start).count() << endl;
  poslp.printStatus();

  start = timer::now();
  ofstream fs(out);
  poslp.serialize(fs);
  stop = timer::now();
  cout << "time to serialize (ms): " << duration_cast<milliseconds>(stop-start).count() << endl;
}

template
<
  class SlpT
  >
void measure_PlainSlp
(
 const NaiveSlp<var_t> & slp,
 std::string out
) {
  NaiveSlp<var_t> temp(slp);
  auto start = timer::now();
  temp.makeBinaryTree();
  SlpT pslp;
  pslp.init(temp);
  auto stop = timer::now();
  cout << "time to encode (ms): " << duration_cast<milliseconds>(stop-start).count() << endl;
  pslp.printStatus();

  start = timer::now();
  ofstream fs(out);
  pslp.serialize(fs);
  stop = timer::now();
  cout << "time to serialize (ms): " << duration_cast<milliseconds>(stop-start).count() << endl;
}


int main(int argc, char* argv[])
{
  using Fiv = IntVec<>;
  using SelSd = SelectSdvec<>;
  using SelMcl = SelectMcl<>;
  using DagcSd = DirectAccessibleGammaCode<SelSd>;
  using DagcMcl = DirectAccessibleGammaCode<SelMcl>;
  using Vlc64 = VlcVec<sdsl::coder::elias_delta, 64>;
  using Vlc128 = VlcVec<sdsl::coder::elias_delta, 128>;
  using funcs_type = map<string,
                         void(*)
                         (
                          const NaiveSlp<var_t> & slp,
                          std::string out
                          )>;
  funcs_type funcs;

  //// PlainSlp
  funcs.insert(make_pair("PlainSlp_FivFiv", measure_PlainSlp<PlainSlp<var_t, Fiv, Fiv>>));
  funcs.insert(make_pair("PlainSlp_IblcFiv", measure_PlainSlp<PlainSlp<var_t, IncBitLenCode, Fiv>>));
  funcs.insert(make_pair("PlainSlp_32Fiv", measure_PlainSlp<PlainSlp<var_t, IntVec<32>, Fiv>>));

  //// PoSlp: Post-order SLP
  //// Sometimes PoSlp_Sd is better than PoSlp_Iblc
  funcs.insert(make_pair("PoSlp_Iblc", measure_PoSlp<PoSlp<var_t, IncBitLenCode>>));
  funcs.insert(make_pair("PoSlp_Sd", measure_PoSlp<PoSlp<var_t, DagcSd>>));
  // funcs.insert(make_pair("PoSlp_Mcl", measure_PoSlp<PoSlp<var_t, DagcMcl>>));

  //// ShapedSlp: plain implementation of slp encoding that utilizes shape-tree grammar
  //// Since bit length to represent slp element is small, SelMcl is good for them.
  //// For stg and bal element, SelSd is better
  funcs.insert(make_pair("ShapedSlp_SdMclSd_SdMcl", measure<ShapedSlp<var_t, DagcSd, DagcMcl, DagcSd, SelSd, SelMcl>>));
  funcs.insert(make_pair("ShapedSlp_SdSdSd_SdMcl", measure<ShapedSlp<var_t, DagcSd, DagcSd, DagcSd, SelSd, SelMcl>>));

  //// ShapedSlpV2: all vlc vectors are merged into one.
  //// Generally encoding size gets worse than ShapedSlp_SdMclSd_SdMcl because
  //// - bit length to represnet stg and bal element is large and DagcSd is a good choice, whereas
  //// - bit size to represent slp element is significantly small and DagcMcl should be used
  funcs.insert(make_pair("ShapedSlpV2_Sd_SdMcl", measure<ShapedSlpV2<var_t, DagcSd, SelSd, SelMcl>>));
  // funcs.insert(make_pair("ShapedSlpV2_SdSdSd", measure<ShapedSlp<var_t, DagcSd, SelSd, SelSd>>));
  // funcs.insert(make_pair("ShapedSlpV2_SdMclMcl", measure<ShapedSlp<var_t, DagcSd, SelMcl, SelMcl>>));
  // funcs.insert(make_pair("ShapedSlpV2_Vlc128SdSd", measure<ShapedSlp<var_t, Vlc128, SelSd, SelSd>>));

  //// SelfShapedSlp: ShapedSlp that does not use shape-tree grammar
  funcs.insert(make_pair("SelfShapedSlp_SdSd_Sd", measure<SelfShapedSlp<var_t, DagcSd, DagcSd, SelSd>>));
  funcs.insert(make_pair("SelfShapedSlp_SdSd_Mcl", measure<SelfShapedSlp<var_t, DagcSd, DagcSd, SelMcl>>));
  // funcs.insert(make_pair("SelfShapedSlp_MclMcl_Sd", measure<SelfShapedSlp<var_t, DagcMcl, DagcMcl, SelSd>>));
  // funcs.insert(make_pair("SelfShapedSlp_SdMcl_Sd", measure<SelfShapedSlp<var_t, DagcSd, DagcMcl, SelSd>>));

  //// SelfShapedSlpV2:
  //// attempted to asign smaller offsets to frequent variables by giving special seats for hi-frequent ones
  funcs.insert(make_pair("SelfShapedSlpV2_SdSd_Sd", measure<SelfShapedSlpV2<var_t, DagcSd, DagcSd, SelSd>>));
  // funcs.insert(make_pair("SelfShapedSlpV2_SdSd_Mcl", measure<SelfShapedSlpV2<var_t, DagcSd, DagcSd, SelMcl>>));

  string methodList;
  for (auto itr = funcs.begin(); itr != funcs.end(); ++itr) {
    methodList += itr->first + ". ";
  }


  cmdline::parser parser;
  parser.add<string>("input", 'i', "input file name. <input>.C and <input>.R in Navarro's RePair format", true);
  parser.add<string>("output", 'o', "output file to which data structure is written", true);
  parser.add<string>("format", 'f', "format of input: NavarroRepair. Bigrepair.", true);
  parser.add<string>("encoding", 'e', "encoding: " + methodList, true);
  parser.parse_check(argc, argv);
  const string in = parser.get<string>("input");
  const string out = parser.get<string>("output");
  const string format = parser.get<string>("format");
  const string encoding = parser.get<string>("encoding");

  NaiveSlp<var_t> slp;
  if (format.compare("NavarroRepair") == 0) {
    slp.load_NavarroRepair(in.data());
  } else if (format.compare("Bigrepair") == 0) {
    slp.load_Bigrepair(in.data());
  } else {
    cerr << "error: specify a valid format name: NavarroRepair. Bigrepair" << endl;
    exit(1);
  }

  if (encoding.compare("All") == 0) {
    for (auto itr = funcs.begin(); itr != funcs.end(); ++itr) {
      if ((itr->first).compare("ShapedSlp_Status_SdSdSd_SdMcl") != 0) {
        cout << itr->first << ": BEGIN" << std::endl;
        itr->second(slp, out + "_" + itr->first);
        cout << itr->first << ": END" << std::endl;
      }
    }
  } else {
    auto itr = funcs.find(encoding);
    if (itr != funcs.end()) {
      itr->second(slp, out);
    } else {
      cerr << "error: specify a valid encoding name: " + methodList << endl;
      exit(1);
    }
  }

  return 0;
}

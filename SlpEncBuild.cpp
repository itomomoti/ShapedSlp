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


int main(int argc, char* argv[])
{
  cmdline::parser parser;
  parser.add<std::string>("input", 'i', "input file name. <input>.C and <input>.R in Navarro's RePair format", true);
  parser.add<std::string>("output", 'o', "output file to which data structure is written", true);
  parser.add<uint32_t>("format", 'f', "format of input: [0] Navarro's Repair. [1] Bigrepair.", true);
  parser.add<uint32_t>("encoding", 'e', "encoding: [0] PoSlp. [1] ShapedSlpDagc. [2] ShapedSlpDagcWithRrr [3] ShapedSlp_64. [4] ShapedSlp_128. ", true);
  parser.parse_check(argc, argv);
  const std::string in = parser.get<std::string>("input");
  const std::string out = parser.get<std::string>("output");
  const uint32_t format = parser.get<uint32_t>("format");
  const uint32_t encoding = parser.get<uint32_t>("encoding");

  using var_t = uint32_t;

  NaiveSlp<var_t> slp;
  if (format == 0) {
    slp.load_NavarroRepair(in.data());
  } else if (format == 1) {
    slp.load_Bigrepair(in.data());
  }
  cout << "height = " << slp.calcHeight() << std::endl;

  if (encoding == 0) {
    cout << "PoSlp:" << endl;
    auto start = timer::now();
    slp.makeBinaryTree();
    PoSlp<var_t> poslp;
    poslp.makePoSlpWithLen(slp);
    auto stop = timer::now();
    cout << "time to encode (ms): " << duration_cast<milliseconds>(stop-start).count() << endl;
    poslp.printStatus();

    start = timer::now();
    ofstream fs(out);
    poslp.serialize(fs);
    stop = timer::now();
    cout << "time to serialize (ms): " << duration_cast<milliseconds>(stop-start).count() << endl;
  } else if (encoding == 1) {
    cout << "ShapedSlp_Dac:" << endl;
    auto start = timer::now();
    ShapedSlp<var_t, DirectAccessibleGammaCode> sslp(slp);
    auto stop = timer::now();
    cout << "time to encode (ms): " << duration_cast<milliseconds>(stop-start).count() << endl;
    sslp.printStatus();

    start = timer::now();
    ofstream fs(out);
    sslp.serialize(fs);
    stop = timer::now();
    cout << "time to serialize (ms): " << duration_cast<milliseconds>(stop-start).count() << endl;
  } else if (encoding == 2) {
    cout << "ShapedSlp_Rrr:" << endl;
    auto start = timer::now();
    ShapedSlp<var_t, DirectAccessibleGammaCodeWithRrr> sslp(slp);
    auto stop = timer::now();
    cout << "time to encode (ms): " << duration_cast<milliseconds>(stop-start).count() << endl;
    sslp.printStatus();

    start = timer::now();
    ofstream fs(out);
    sslp.serialize(fs);
    stop = timer::now();
    cout << "time to serialize (ms): " << duration_cast<milliseconds>(stop-start).count() << endl;
  } else if (encoding == 3) {
    using VlcT = VlcVec<sdsl::coder::elias_delta, 64>;
    cout << "ShapedSlp_64:" << endl;
    auto start = timer::now();
    ShapedSlp<var_t, VlcT> sslp(slp);
    auto stop = timer::now();
    cout << "time to encode (ms): " << duration_cast<milliseconds>(stop-start).count() << endl;
    sslp.printStatus();

    start = timer::now();
    ofstream fs(out);
    sslp.serialize(fs);
    stop = timer::now();
    cout << "time to serialize (ms): " << duration_cast<milliseconds>(stop-start).count() << endl;
  } else if (encoding == 4) {
    using VlcT = VlcVec<sdsl::coder::elias_delta, 128>;
    cout << "ShapedSlp_128:" << endl;
    auto start = timer::now();
    ShapedSlp<var_t, VlcT> sslp(slp);
    auto stop = timer::now();
    cout << "time to encode (ms): " << duration_cast<milliseconds>(stop-start).count() << endl;
    sslp.printStatus();

    start = timer::now();
    ofstream fs(out);
    sslp.serialize(fs);
    stop = timer::now();
    cout << "time to serialize (ms): " << duration_cast<milliseconds>(stop-start).count() << endl;
  }

  return 0;
}

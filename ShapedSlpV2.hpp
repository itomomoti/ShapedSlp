#ifndef INCLUDE_GUARD_ShapedSlpV2
#define INCLUDE_GUARD_ShapedSlpV2

#include <sys/stat.h>
#include <iostream>
#include <string>
#include <stdint.h> // include uint64_t etc.
#include <map>
#include <set>
#include "Common.hpp"
#include "NaiveSlp.hpp"
#include "RecSplit.hpp"
#include <sdsl/bit_vectors.hpp>
#include <sdsl/vlc_vector.hpp>
#include <sdsl/coder.hpp>


/*!
 * @file ShapedSlpV2.hpp
 * @brief An SLP encoding that utilizes its shape-tree grammar
 * @author Tomohiro I
 * @date 2019-11-07
 */
template
<
  typename tparam_var_t,
  class DacT,
  class StgDivSelT,
  class SlpDivSelT
  >
class ShapedSlpV2
{
public:
  //// Public constant, alias etc.
  using var_t = tparam_var_t;


private:
  //// parameter of RecSplit
  static constexpr size_t kBucketSize = 100;
  static constexpr size_t kLeaf = 8;

  std::vector<char> alph_;
  sdsl::sd_vector<> seqSBV_;
  sdsl::sd_vector<>::rank_1_type seqRank_;
  sdsl::sd_vector<>::select_1_type seqSel_;
  StgDivSelT stgDivSel_;
  SlpDivSelT slpDivSel_;
  DacT vlcSeq_;
  DacT vlc_;
  sux::function::RecSplit<kLeaf> * rs_; // minimal perfect hash: map from "expansion lengths" to IDs for them


public:
  ShapedSlpV2
  () : rs_(nullptr)
  {}


  ShapedSlpV2
  (
   const NaiveSlp<var_t> & slp,
   const bool freqSort = true
   ) : rs_(nullptr)
  {
    makeShapedSlp(slp);
  }


  ~ShapedSlpV2() {
    delete(rs_);
  }


  size_t getAlphSize() const {
    return alph_.size();
  }


  size_t getLen() const {
    return seqSBV_.size() - 1;
  }


  size_t getLenSeq() const {
    return vlcSeq_.size() / 2;
  }


  size_t getNumRulesOfStg() const {
    return stgDivSel_.size();
  }


  size_t getNumRulesOfSlp() const {
    return slpDivSel_.size();
  }


  char charAt
  (
   const uint64_t pos //!< 0-based position
   ) const {
    assert(pos < getLen());

    const uint64_t seqPos = seqRank_(pos + 1);
    const uint64_t varLen = lenOfSeqAt(seqPos);
    const uint64_t prevSum = (seqPos > 0) ? seqSel_(seqPos) : 0;
    return charAt(pos - prevSum, varLen, vlcSeq_[2 * seqPos], vlcSeq_[2 * seqPos + 1]);
  }


  char charAt
  (
   const uint64_t pos, //!< 0-based relative position in a variable
   const uint64_t varLen, //!< expansion length of the variable
   const var_t stgOffset, //!< stg offset for the variable
   const var_t slpOffset //!< slp offset for the variable 
   ) const {
    assert(pos < varLen);
    // std::cout << "pos = " << pos << ", varLen = " << varLen << ", stgOffset = " << stgOffset << ", slpOffset = " << slpOffset << std::endl;

    if (varLen == 1) {
      return alph_[slpOffset];
    }
    const uint64_t h = hashLen(varLen);
    const uint64_t stgNum = stgOffset + ((h == 0) ? 0 : stgDivSel_(h) + 1);
    const uint64_t slpNum = (stgNum == 0) ? 0 : slpDivSel_(stgNum) + 1;
    const uint64_t vlcPos = 3 * stgNum + 2 * slpNum;
    const uint64_t leftLen = decLeftVarLen(varLen, vlc_[vlcPos]);
    if (pos < leftLen) {
      return charAt(pos, leftLen, vlc_[vlcPos + 1], vlc_[vlcPos + 3 + 2 * slpOffset]);
    } else {
      return charAt(pos - leftLen, varLen - leftLen, vlc_[vlcPos + 2], vlc_[vlcPos + 4 + 2 * slpOffset]);
    }
  }


  void expandSubstr
  (
   const uint64_t pos, //!< 0-based beginning position
   uint64_t len, //!< length to expand
   char * str //!< [out] must have length at least 'len'
   ) const {
    assert(pos < getLen());
    assert(len > 0);
    assert(len <= getLen() - pos);

    uint64_t seqPos = seqRank_(pos + 1);
    const uint64_t varLen = lenOfSeqAt(seqPos);
    const uint64_t prevSum = (seqPos > 0) ? seqSel_(seqPos) : 0;
    expandSubstr(pos - prevSum, len, str, varLen, vlcSeq_[2 * seqPos], vlcSeq_[2 * seqPos + 1]);
    for (uint64_t maxExLen = prevSum + varLen - pos; maxExLen < len; ) {
      len -= maxExLen;
      str += maxExLen;
      maxExLen = lenOfSeqAt(++seqPos);
      expandPref(len, str, maxExLen, vlcSeq_[2 * seqPos], vlcSeq_[2 * seqPos + 1]);
    }
  }


  void expandSubstr
  (
   const uint64_t pos, //!< 0-based beginning position (relative in variable)
   const uint64_t len, //!< length to expand
   char * str, //!< [out] must have length at least 'len'
   const uint64_t varLen, //!< expansion length of the variable
   const var_t stgOffset, //!< stg offset for the variable
   const var_t slpOffset //!< slp offset for the variable 
   ) const {
    // std::cout << "pos = " << pos << ", len = " << len << ", varLen = " << varLen << ", stgOffset = " << stgOffset << ", slpOffset = " << slpOffset << std::endl;
    assert(pos < varLen);

    if (varLen == 1) {
      *str = alph_[slpOffset];
      return;
    }
    const uint64_t h = hashLen(varLen);
    const uint64_t stgNum = stgOffset + ((h == 0) ? 0 : stgDivSel_(h) + 1);
    const uint64_t slpNum = (stgNum == 0) ? 0 : slpDivSel_(stgNum) + 1;
    const uint64_t vlcPos = 3 * stgNum + 2 * slpNum;
    const uint64_t leftLen = decLeftVarLen(varLen, vlc_[vlcPos]);
    if (pos < leftLen) {
      expandSubstr(pos, len, str, leftLen, vlc_[vlcPos + 1], vlc_[vlcPos + 3 + 2 * slpOffset]);
      if (leftLen - pos < len) {
        expandPref(len - (leftLen - pos), str + (leftLen - pos), varLen - leftLen, vlc_[vlcPos + 2], vlc_[vlcPos + 4 + 2 * slpOffset]);
      }
    } else {
      expandSubstr(pos - leftLen, len, str, varLen - leftLen, vlc_[vlcPos + 2], vlc_[vlcPos + 4 + 2 * slpOffset]);
    }
  }


  void expandPref
  (
   const uint64_t len, //!< length to expand
   char * str, //!< [out] must have length at least 'len'
   const uint64_t varLen, //!< expansion length of the variable
   const var_t stgOffset, //!< stg offset for the variable
   const var_t slpOffset //!< slp offset for the variable 
   ) const {
    // std::cout << "len = " << len << ", varLen = " << varLen << ", stgOffset = " << stgOffset << ", slpOffset = " << slpOffset << std::endl;
    assert(len > 0);

    if (varLen == 1) {
      *str = alph_[slpOffset];
      return;
    }
    const uint64_t h = hashLen(varLen);
    const uint64_t stgNum = stgOffset + ((h == 0) ? 0 : stgDivSel_(h) + 1);
    const uint64_t slpNum = (stgNum == 0) ? 0 : slpDivSel_(stgNum) + 1;
    const uint64_t vlcPos = 3 * stgNum + 2 * slpNum;
    const uint64_t leftLen = decLeftVarLen(varLen, vlc_[vlcPos]);
    expandPref(len, str, leftLen, vlc_[vlcPos + 1], vlc_[vlcPos + 3 + 2 * slpOffset]);
    if (len > leftLen) {
      expandPref(len - leftLen, str + leftLen, varLen - leftLen, vlc_[vlcPos + 2], vlc_[vlcPos + 4 + 2 * slpOffset]);
    }
  }


  void printStatus
  (
   const bool verbose = false
   ) const {
    std::cout << "ShapedSlpV2 object (" << this << ") " << __func__ << "(" << verbose << ") BEGIN" << std::endl;
    const size_t len = getLen();
    const size_t alphSize = getAlphSize();
    const size_t lenSeq = getLenSeq();
    const size_t numRulesOfStg = getNumRulesOfStg();
    const size_t numRulesOfSlp = getNumRulesOfSlp();
    const size_t numDistLen = rs_->size();
    std::cout << "alphSize = " << alphSize << ", len = " << len << ", lenSeq = " << lenSeq
              << ", numRulesOfSlp = " << numRulesOfSlp
              << ", numRulesOfStg = " << numRulesOfStg
              << ", numDistLen = " << numDistLen
              << std::endl;

    const size_t bytesAlphSize = sizeof(std::vector<char>) + (sizeof(char) * alph_.size());
    const size_t bytesSeqSBV = sdsl::size_in_bytes(seqSBV_);
    const size_t bytesVlcSeq = vlcSeq_.calcMemBytes();
    const size_t bytesVlc = vlc_.calcMemBytes();
    const size_t bytesStgDivSel = stgDivSel_.calcMemBytes();
    const size_t bytesSlpDivSel = slpDivSel_.calcMemBytes();

    const size_t bytesMph = calcMemBytesOfMph();
    const size_t bytesTotal = bytesMph + bytesSeqSBV + bytesVlcSeq + bytesVlc + bytesStgDivSel + bytesSlpDivSel + bytesAlphSize;
    const size_t bytesEstStgWithLen = estimateEncSizeOfStgWithLen();
    const size_t bytesEstSlpWithLen = estimateEncSizeOfSlpWithLen();
    const size_t bytesEstSlp = estimateEncSizeOfSlp();
    std::cout << "Sizes (bytes) for various approach (small o() term is ignored for the ones with est.)" << std::endl;
    std::cout << "New encoding = " << bytesTotal << std::endl;
    std::cout << "| MPH = " << bytesMph << std::endl;
    std::cout << "| | MPH / numDistLen = " << (double) bytesMph / numDistLen << std::endl;
    std::cout << "| seqSBV = " << bytesSeqSBV << std::endl;
    std::cout << "| vlcSeq = " << bytesVlcSeq << std::endl;
    std::cout << "| | vlcSeq per entry = " << (double)bytesVlcSeq / (2 * lenSeq) << std::endl;
    std::cout << "| vlcRules = " << bytesVlc << std::endl;
    std::cout << "| | vlcRules per entry = " << (double)bytesVlc / (3 * numRulesOfStg + 2 * numRulesOfSlp) << std::endl;
    std::cout << "| stgDiv = " << bytesStgDivSel << std::endl;
    std::cout << "| slpDiv = " << bytesSlpDivSel << std::endl;
    std::cout << "| alph = " << bytesAlphSize << std::endl;
    std::cout << "MaruyamaEnc of Stg + POSLP (est.) = " << bytesEstStgWithLen + bytesEstSlp << std::endl;
    std::cout << "| MaruyamaEnc of Stg (est.) = " << bytesEstStgWithLen << std::endl;
    std::cout << "| POSLP (est.) = " << bytesEstSlp << std::endl;
    std::cout << "MaruyamaEnc of Slp (est.) = " << bytesEstSlpWithLen << std::endl;
    if (verbose) {
      std::cout << "alph_" << std::endl;
      printVec(alph_);
      std::cout << std::endl;
      std::cout << "vlc_" << std::endl;
      printVec(vlc_);
      std::cout << "vlcSeq_" << std::endl;
      printVec(vlcSeq_);
      std::cout << "hash_" << std::endl;
      for (uint64_t i = 0; i < getLen(); ++i) {
        std::cout << "(" << i << ":" << hashLen(i) << ") ";
      }
      std::cout << std::endl;
      std::cout << "stgDiv_" << std::endl;
      printVec(stgDivSel_);
      std::cout << "stgDivSel_" << std::endl;
      for (uint64_t i = 1; i <= rs_->size(); ++i) {
        std::cout << "(" << i << ":" << stgDivSel_(i) << ") ";
      }
      std::cout << std::endl;
      std::cout << "slpDiv_" << std::endl;
      printVec(slpDivSel_);
      std::cout << "slpDivSel_" << std::endl;
      for (uint64_t i = 1; i <= getNumRulesOfStg(); ++i) {
        std::cout << "(" << i << ":" << slpDivSel_(i) << ") ";
      }
      std::cout << std::endl;
    }
    std::cout << "ShapedSlpV2 object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
  }


  size_t calcMemBytesOfMph() const {
    char fname[] = "rs_temp_output"; // temp
    fstream fs;
    fs.exceptions(fstream::failbit | fstream::badbit);
    fs.open(fname, fstream::out | fstream::binary | fstream::trunc);
    fs << (*rs_);
    struct stat s;
    stat(fname, &s);
    return s.st_size;
  }


  size_t estimateEncSizeOfStgWithLen() const {
    //// |G| lg |S| + |G| lg(1 + σ/|G|) + 5|G| + o(|G|) bits
    //// in Theorem 3 of the paper "Fully-Online Grammar Compression”, SPIRE 2013
    uint64_t ls = getLenSeq();
    uint64_t nr = getNumRulesOfStg();
    uint64_t as = 1;
    return ((ls + nr) * ceilLog2(ls + nr + as) + 2 * (ls + nr)) / 8 +
      ((ls + nr) * log(getLen() / (ls + nr)) / log(2.0) + 3 * (ls + nr)) / 8;
  }


  size_t estimateEncSizeOfSlpWithLen() const {
    //// |G| lg |S| + |G| lg(1 + σ/|G|) + 5|G| + o(|G|) bits
    //// in Theorem 3 of the paper "Fully-Online Grammar Compression”, SPIRE 2013
    uint64_t ls = getLenSeq();
    uint64_t nr = getNumRulesOfSlp();
    uint64_t as = getAlphSize();
    // return ((ls + nr) * ceilLog2(getLen()) + (ls + nr) * log(1 + (double)as/(ls + nr)) / log(2.0) + 5 * (ls + nr)) / 8;
    return estimateEncSizeOfSlp() +
      ((ls + nr) * log(getLen() / (ls + nr)) / log(2.0) + 3 * (ls + nr)) / 8;
  }


  size_t estimateEncSizeOfSlp() const {
    //// |G| lg(|G| + σ) + 2|G| + o(|G|) bits
    //// in Theorem 3 of the paper "Fully-Online Grammar Compression”, SPIRE 2013
    uint64_t ls = getLenSeq();
    uint64_t nr = getNumRulesOfSlp();
    uint64_t as = getAlphSize();
    return ((ls + nr) * ceilLog2(ls + nr + as) + 2 * (ls + nr)) / 8; // little cheat on representation of leaves
  }


  void load
  (
   std::istream & in
   ) {
    uint64_t alphSize = 0;
    in.read((char*) & alphSize, sizeof(alphSize));
    alph_.resize(alphSize);
    in.read((char*) alph_.data(), alphSize * sizeof(alph_[0]));
    rs_ = new sux::function::RecSplit<kLeaf>();
    in >> (*rs_);
    seqSBV_.load(in);
    seqRank_.load(in);
    seqSel_.load(in);
    stgDivSel_.load(in);
    slpDivSel_.load(in);
    vlcSeq_.load(in);
    vlc_.load(in);

    seqRank_.set_vector(&seqSBV_);
    seqSel_.set_vector(&seqSBV_);
  }


  void serialize
  (
   std::ostream & out
   ) const {
    assert(rs_ != nullptr);

    uint64_t alphSize = getAlphSize();
    out.write((char*) & alphSize, sizeof(alphSize));
    out.write((char*) alph_.data(), alphSize * sizeof(alph_[0]));
    out << (*rs_);
    seqSBV_.serialize(out);
    seqRank_.serialize(out);
    seqSel_.serialize(out);
    stgDivSel_.serialize(out);
    slpDivSel_.serialize(out);
    vlcSeq_.serialize(out);
    vlc_.serialize(out);
  }

  
private:
  //// used to create input for RecSplit
  string uint2Str(const uint64_t n) const {
    string ret;
    ret.resize(8);
    for (uint64_t i = 0; i < 8; ++i) {
      ret[i] = (n >> (8 * i)) & 0xFF;
    }
    return ret;
  }


  uint64_t hashLen(uint64_t len) const {
    return (*rs_)(uint2Str(len));
  }


  uint64_t lenOfSeqAt(uint64_t i) const {
    assert(i < getLenSeq());
    return (i > 0) ? seqSel_(i+1) - seqSel_(i) : seqSel_(i+1);
  }


  uint64_t encBal
  (
   const uint64_t varlen,
   const uint64_t leftvarlen
   ) const {
    if (varlen/2 <= leftvarlen) {
      return ((leftvarlen - varlen/2) << 1); // lsb is 0
    } else {
      return ((varlen/2 - leftvarlen) << 1) + 1; // lsb is 1
    }
  }


  uint64_t decLeftVarLen
  (
   const uint64_t varlen,
   const uint64_t balenc
   ) const {
    if ((balenc & 1) == 0) { // lsb is 0
      return varlen/2 + (balenc >> 1);
    } else { // lsb is 1
      return varlen/2 - (balenc >> 1);
    }
  }


  void makeShapedSlp
  (
   const NaiveSlp<var_t> & slp,
   const bool freqSort = true
   ) {
    alph_.resize(slp.getAlphSize());
    for (uint64_t i = 0; i < slp.getAlphSize(); ++i) {
      alph_[i] = slp.getChar(i);
    }
    NaiveSlp<var_t> stg;
    std::vector<uint64_t> slp2stg(slp.getNumRules()); // map from slp variable to stg variable

    { // make stg
      stg.setLenSeq(slp.getLenSeq());
      stg.setAlphSize(1); // in shape grammar all leaves are labeled with 0
      stg.setChar(0, '!'); // for printing
      uint64_t numRules = 0;

      std::map<PairT<var_t>, var_t> p2stg;
      for (uint64_t i = 0; i < slp.getNumRules(); ++i) { // assumption: each variable refers to smaller variables
        PairT<var_t> p;
        p.left = (slp.getLeft(i) < slp.getAlphSize()) ? 0 : slp2stg[slp.getLeft(i) - slp.getAlphSize()];
        p.right = (slp.getRight(i) < slp.getAlphSize()) ? 0 : slp2stg[slp.getRight(i) - slp.getAlphSize()];
        uint64_t val;
        auto itr = p2stg.find(p);
        if (itr != p2stg.end()) { // pair found
          val = (*itr).second;
        } else { // new shape grammar variable
          val = ++numRules;
          stg.setNumRules(numRules);
          stg.setRule(val - 1, p); // alphabet size of stg is canceled by -1
          p2stg.insert(std::make_pair(p, static_cast<var_t>(val)));
        }
        slp2stg[i] = val;
      }
      for (uint64_t i = 0; i < slp.getLenSeq(); ++i) {
        const auto v = slp.getSeq(i);
        stg.setSeq(i, (v < slp.getAlphSize()) ? 0 : slp2stg[v - slp.getAlphSize()]);
      }
    }

    std::vector<uint64_t> stglen(stg.getNumRules());
    stg.makeLenVec(stglen); // stglen is now: expansion lengths

    { // construct prefix sum data structure
      sdsl::int_vector<64> psum(slp.getLenSeq());
      uint64_t s = 0;
      for (uint64_t i = 0; i < slp.getLenSeq(); ++i) {
        s += stg.getLenOfVar(stg.getSeq(i), stglen);
        psum[i] = s;
      }
      seqSBV_ = std::move(sdsl::sd_vector<>(psum.begin(), psum.end()));
      seqRank_.set_vector(&seqSBV_);
      seqSel_.set_vector(&seqSBV_);
    }

    { // build minimal perfect hash
      std::set<uint64_t> distLen;
      for (uint64_t i = 0; i < stglen.size(); ++i) {
        distLen.insert(stglen[i]);
      }
      std::vector<std::string> keys;
      for (auto itr = distLen.begin(); itr != distLen.end(); ++itr) {
        keys.push_back(uint2Str(*itr));
      }
      rs_ = new sux::function::RecSplit<kLeaf>(keys, kBucketSize);
    }

    std::vector<var_t> stgOrder(stg.getNumRules());
    std::vector<var_t> stgOffset(stg.getNumRules());
    {
      for (uint64_t i = 0; i < stg.getNumRules(); ++i) {
        stgOrder[i] = i;
      }
      if (freqSort) {
        std::vector<uint64_t> ruleFreq(stg.getNumRules());
        std::vector<uint64_t> alphFreq(stg.getAlphSize());
        stg.makeFreqInRulesVec(ruleFreq, alphFreq);
        std::stable_sort
          (
           stgOrder.begin(),
           stgOrder.end(),
           [&](uint64_t x, uint64_t y) { return ruleFreq[x] > ruleFreq[y]; }
           );
      }
      std::vector<var_t> hashVal(stg.getNumRules());
      for (uint64_t i = 0; i < stg.getNumRules(); ++i) {
        hashVal[i] = hashLen(stglen[i]);
      }
      std::stable_sort
        (
         stgOrder.begin(),
         stgOrder.end(),
         [&](uint64_t x, uint64_t y) { return hashVal[x] < hashVal[y]; }
         );
      sdsl::bit_vector stgDiv(stg.getNumRules(), 0);
      uint64_t offset = 0;
      for (uint64_t i = 0; i < stg.getNumRules() - 1; ++i) {
        stgOffset[stgOrder[i]] = offset++;
        if (hashVal[stgOrder[i]] == hashVal[stgOrder[i+1]]) {
          stgDiv[i] = 0;
        } else {
          stgDiv[i] = 1;
          offset = 0;
        }
      }
      stgDiv[stg.getNumRules() - 1] = 1;
      stgOffset[stgOrder[stg.getNumRules() - 1]] = offset;
      stgDivSel_.init(std::move(stgDiv));
    }

    std::vector<var_t> slpOrder(slp.getNumRules());
    std::vector<var_t> slpOffset(slp.getNumRules());
    {
      for (uint64_t i = 0; i < slpOrder.size(); ++i) {
        slpOrder[i] = i;
      }
      if (freqSort) {
        std::vector<uint64_t> ruleFreq(slp.getNumRules());
        std::vector<uint64_t> alphFreq(slp.getAlphSize()); // TODO: sorting alphabet based on frequency
        slp.makeFreqInRulesVec(ruleFreq, alphFreq);
        std::stable_sort
          (
           slpOrder.begin(),
           slpOrder.end(),
           [&](uint64_t x, uint64_t y) { return ruleFreq[x] > ruleFreq[y]; }
           );
      }
      std::vector<var_t> stgRank(stg.getNumRules());
      for (uint64_t i = 0; i < stgOrder.size(); ++i) {
        stgRank[stgOrder[i]] = i;
      }
      for (uint64_t i = 0; i < slpOrder.size(); ++i) {
        slp2stg[i] = stgRank[slp2stg[i] - 1]; // replace stg with rank of stg (-1 to cancel a letter)
      }
      std::stable_sort
        (
         slpOrder.begin(),
         slpOrder.end(),
         [&](uint64_t x, uint64_t y) { return slp2stg[x] < slp2stg[y]; }
         );
      sdsl::bit_vector slpDiv(slp.getNumRules(), 0);
      uint64_t offset = 0;
      for (uint64_t i = 0; i < slp.getNumRules() - 1; ++i) {
        slpOffset[slpOrder[i]] = offset++;
        if (slp2stg[slpOrder[i]] == slp2stg[slpOrder[i+1]]) {
          slpDiv[i] = 0;
        } else {
          slpDiv[i] = 1;
          offset = 0;
        }
      }
      slpDiv[slp.getNumRules() - 1] = 1;
      slpOffset[slpOrder[slp.getNumRules() - 1]] = offset;
      slpDivSel_.init(std::move(slpDiv));
    }

    {
      const uint64_t dfsize = 3 * stg.getNumRules() + 2 * slp.getNumRules();
      std::vector<uint64_t> df(dfsize);
      uint64_t dfpos = 0;
      uint64_t slppos = 0;
      for (uint64_t stgpos = 0; stgpos < stg.getNumRules(); ++stgpos) {
        const uint64_t stgRuleId = stgOrder[stgpos]; // in [0..stg.getNumRules())
        const uint64_t stgLeftVar = stg.getLeft(stgRuleId); // in [0..stg.getNumRules() + stg.getAlphSize())
        const uint64_t stgRightVar = stg.getRight(stgRuleId); // in [0..stg.getNumRules() + stg.getAlphSize())
        const uint64_t varlen = stglen[stgRuleId];
        const uint64_t leftvarlen = stg.getLenOfVar(stgLeftVar, stglen);
        df[dfpos++] = encBal(varlen, leftvarlen);
        df[dfpos++] = (stgLeftVar < stg.getAlphSize()) ? stgLeftVar : stgOffset[stgLeftVar - stg.getAlphSize()];
        df[dfpos++] = (stgRightVar < stg.getAlphSize()) ? stgRightVar : stgOffset[stgRightVar - stg.getAlphSize()];        
        for (const uint64_t nextSlpPos = slpDivSel_(stgpos + 1) + 1; slppos < nextSlpPos; ++slppos) {
          const uint64_t slpRuleId = slpOrder[slppos]; // in [0..slp.getNumRules())
          const uint64_t slpLeftVar = slp.getLeft(slpRuleId); // in [0..slp.getNumRules() + slp.getAlphSize())
          const uint64_t slpRightVar = slp.getRight(slpRuleId); // in [0..slp.getNumRules() + slp.getAlphSize())
          df[dfpos++] = (slpLeftVar < slp.getAlphSize()) ? slpLeftVar : slpOffset[slpLeftVar - slp.getAlphSize()];
          df[dfpos++] = (slpRightVar < slp.getAlphSize()) ? slpRightVar : slpOffset[slpRightVar - slp.getAlphSize()];
        }
      }
      vlc_.init(df);
    }

    {
      const uint64_t dfsize = 2 * slp.getLenSeq();
      std::vector<uint64_t> df(dfsize);
      for (uint64_t pos = 0; pos < dfsize / 2; ++pos) {
        const uint64_t stgVar = stg.getSeq(pos);
        df[2 * pos] = (stgVar < stg.getAlphSize()) ? stgVar : stgOffset[stgVar - stg.getAlphSize()];
        const uint64_t slpVar = slp.getSeq(pos);
        df[2 * pos + 1] = (slpVar < slp.getAlphSize()) ? slpVar : slpOffset[slpVar - slp.getAlphSize()];
      }
      vlcSeq_.init(df);
    }
  }
};

#endif

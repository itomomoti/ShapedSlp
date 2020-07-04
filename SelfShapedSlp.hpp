#ifndef INCLUDE_GUARD_SelfShapedSlp
#define INCLUDE_GUARD_SelfShapedSlp

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


// #define PRINT_STATUS_SelfShapedSlp


template
<
  typename tparam_var_t,
  class DacT,
  class BalDacT,
  class DivSelT
  >
class SelfShapedSlp
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
  DivSelT slpDivSel_;
  sdsl::bit_vector balBv_;
  sdsl::rank_support_v5<> balBvRank_;
  DacT vlcSeq_;
  DacT vlc_;
  BalDacT bal_;
  sux::function::RecSplit<kLeaf> * rs_; // minimal perfect hash: from "expansion lengths" to IDs for them


public:
  SelfShapedSlp
  () : rs_(nullptr)
  {}


  SelfShapedSlp
  (
   const NaiveSlp<var_t> & slp,
   const bool freqSort = true
   ) : rs_(nullptr)
  {
    makeShapedSlp(slp);
  }


  ~SelfShapedSlp() {
    delete(rs_);
  }


  size_t getAlphSize() const {
    return alph_.size();
  }


  size_t getLen() const {
    return seqSBV_.size() - 1;
  }


  size_t getLenSeq() const {
    return vlcSeq_.size();
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
    return charAt(pos - prevSum, varLen, vlcSeq_[seqPos]);
  }


  char charAt
  (
   const uint64_t pos, //!< relative position in a variable
   const uint64_t varLen, //!< expansion length of the variable
   const var_t slpOffset //!< slp offset for the variable 
   ) const {
    assert(pos < varLen);
    // std::cout << "pos = " << pos << ", varLen = " << varLen << ", slpOffset = " << slpOffset << std::endl;

    if (varLen == 1) {
      return alph_[slpOffset];
    }
    const uint64_t h = hashLen(varLen);
    const uint64_t slpId = slpDivSel_(h + 1) + slpOffset;
    const uint64_t balPos = h + balBvRank_(slpId - h);
    const uint64_t leftLen = decLeftVarLen(varLen, bal_[balPos]);
    if (pos < leftLen) {
      return charAt(pos, leftLen, vlc_[2 * slpId]);
    } else {
      return charAt(pos - leftLen, varLen - leftLen, vlc_[2 * slpId + 1]);
    }
  }


  void expandSubstr
  (
   const uint64_t pos, //!< beginning position
   uint64_t len, //!< length to expand
   char * str //!< [out] must have length at least 'len'
   ) const {
    assert(pos < getLen());
    assert(len > 0);
    assert(len <= getLen() - pos);

    uint64_t seqPos = seqRank_(pos + 1);
    const uint64_t varLen = lenOfSeqAt(seqPos);
    const uint64_t prevSum = (seqPos > 0) ? seqSel_(seqPos) : 0;
    expandSubstr(pos - prevSum, len, str, varLen, vlcSeq_[seqPos]);
    for (uint64_t maxExLen = prevSum + varLen - pos; maxExLen < len; ) {
      len -= maxExLen;
      str += maxExLen;
      maxExLen = lenOfSeqAt(++seqPos);
      expandPref(len, str, maxExLen, vlcSeq_[seqPos]);
    }
  }


  void expandSubstr
  (
   const uint64_t pos, //!< beginning position (relative in variable)
   const uint64_t len, //!< length to expand
   char * str, //!< [out] must have length at least 'len'
   const uint64_t varLen, //!< expansion length of the variable
   const var_t slpOffset //!< slp offset for the variable 
   ) const {
    // std::cout << "pos = " << pos << ", len = " << len << ", varLen = " << varLen << ", slpOffset = " << slpOffset << std::endl;
    assert(pos < varLen);

    if (varLen == 1) {
      *str = alph_[slpOffset];
      return;
    }
    const uint64_t h = hashLen(varLen);
    const uint64_t slpId = slpDivSel_(h + 1) + slpOffset;
    const uint64_t balPos = h + balBvRank_(slpId - h);
    const uint64_t leftLen = decLeftVarLen(varLen, bal_[balPos]);
    if (pos < leftLen) {
      expandSubstr(pos, len, str, leftLen, vlc_[2 * slpId]);
      if (leftLen - pos < len) {
        expandPref(len - (leftLen - pos), str + (leftLen - pos), varLen - leftLen, vlc_[2 * slpId + 1]);
      }
    } else {
      expandSubstr(pos - leftLen, len, str, varLen - leftLen, vlc_[2 * slpId + 1]);
    }
  }


  void expandPref
  (
   const uint64_t len, //!< length to expand
   char * str, //!< [out] must have length at least 'len'
   const uint64_t varLen, //!< expansion length of the variable
   const var_t slpOffset //!< slp offset for the variable 
   ) const {
    // std::cout << "len = " << len << ", varLen = " << varLen << ", stgOffset = " << stgOffset << ", slpOffset = " << slpOffset << std::endl;
    assert(len > 0);

    if (varLen == 1) {
      *str = alph_[slpOffset];
      return;
    }
    const uint64_t h = hashLen(varLen);
    const uint64_t slpId = slpDivSel_(h + 1) + slpOffset;
    const uint64_t balPos = h + balBvRank_(slpId - h);
    const uint64_t leftLen = decLeftVarLen(varLen, bal_[balPos]);
    expandPref(len, str, leftLen, vlc_[2 * slpId]);
    if (len > leftLen) {
      expandPref(len - leftLen, str + leftLen, varLen - leftLen, vlc_[2 * slpId + 1]);
    }
  }


  void printStatus
  (
   const bool verbose = false
   ) const {
    std::cout << "SelfShapedSlp object (" << this << ") " << __func__ << "(" << verbose << ") BEGIN" << std::endl;
    const size_t len = getLen();
    const size_t alphSize = getAlphSize();
    const size_t lenSeq = getLenSeq();
    const size_t numRulesOfSlp = getNumRulesOfSlp();
    const size_t numDistLen = rs_->size();
    const size_t numDistLenPairs = numDistLen + balBvRank_(balBv_.size());
    std::cout << "alphSize = " << alphSize << ", len = " << len << ", lenSeq = " << lenSeq
              << ", numRulesOfSlp = " << numRulesOfSlp
              << ", numDistLenPairs = " << numDistLenPairs
              << ", numDistLen = " << numDistLen
              << std::endl;

    const size_t bytesAlphSize = sizeof(std::vector<char>) + (sizeof(char) * alph_.size());
    const size_t bytesSeqSBV = sdsl::size_in_bytes(seqSBV_);
    const size_t bytesVlcSeq = vlcSeq_.calcMemBytes();
    const size_t bytesVlc = vlc_.calcMemBytes();
    const size_t bytesDivSel = slpDivSel_.calcMemBytes();
    const size_t bytesBal = bal_.calcMemBytes();
    const size_t bytesBalBv = sdsl::size_in_bytes(balBv_);
    const size_t bytesBalBvRank = sdsl::size_in_bytes(balBvRank_);

    const size_t bytesMph = calcMemBytesOfMph();
    const size_t bytesTotal = bytesAlphSize + bytesMph + bytesSeqSBV + bytesVlcSeq + bytesVlc + bytesDivSel + bytesBal + bytesBalBv + bytesBalBvRank;
    const size_t bytesEstSlpWithLen = estimateEncSizeOfSlpWithLen();
    const size_t bytesEstSlp = estimateEncSizeOfSlp();
    std::cout << "Sizes (bytes) for various approach (small o() term is ignored for the ones with est.)" << std::endl;
    std::cout << "New encoding = " << bytesTotal << std::endl;
    std::cout << "| MPH = " << bytesMph << std::endl;
    std::cout << "| * MPH / numDistLen = " << (double) bytesMph / numDistLen << std::endl;
    std::cout << "| seqSBV = " << bytesSeqSBV << std::endl;
    std::cout << "| vlcSeq = " << bytesVlcSeq << std::endl;
    std::cout << "| * vlcSeq per entry = " << (double)bytesVlcSeq / lenSeq << std::endl;
    std::cout << "| vlcRules = " << bytesVlc << std::endl;
    std::cout << "| * vlcRules per entry = " << (double)bytesVlc / (2 * numRulesOfSlp) << std::endl;
    std::cout << "| bal = " << bytesBal + bytesBalBv + bytesBalBvRank << std::endl;
    std::cout << "| | balvlc = " << bytesBal << std::endl;
    std::cout << "| | bitv = " << bytesBalBv << std::endl;
    std::cout << "| | rank = " << bytesBalBvRank << std::endl;
    std::cout << "| * per entry = " << (double)(bytesBal + bytesBalBv + bytesBalBvRank) / numDistLenPairs << std::endl;
    std::cout << "| slpDiv = " << bytesDivSel << std::endl;
    std::cout << "| alph = " << bytesAlphSize << std::endl;
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
      std::cout << "slpDiv_" << std::endl;
      printVec(slpDivSel_);
    }
    std::cout << "SelfShapedSlp object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
  }


  size_t calcMemBytesOfMph() const {
    char fname[] = "rs_temp_output"; // temp
    std::fstream fs;
    fs.exceptions(std::fstream::failbit | std::fstream::badbit);
    fs.open(fname, std::fstream::out | std::fstream::binary | std::fstream::trunc);
    fs << (*rs_);
    struct stat s;
    stat(fname, &s);
    return s.st_size;
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
    seqRank_.load(in, &seqSBV_);
    seqSel_.load(in, &seqSBV_);
    slpDivSel_.load(in);
    balBv_.load(in);
    balBvRank_.load(in, &balBv_);
    vlcSeq_.load(in);
    vlc_.load(in);
    bal_.load(in);
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
    slpDivSel_.serialize(out);
    balBv_.serialize(out);
    balBvRank_.serialize(out);
    vlcSeq_.serialize(out);
    vlc_.serialize(out);
    bal_.serialize(out);
  }

  
private:
  //// used to create input for RecSplit
  std::string uint2Str(const uint64_t n) const {
    std::string ret;
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
    return leftvarlen;
  }


  uint64_t decLeftVarLen
  (
   const uint64_t varlen,
   const uint64_t balenc
   ) const {
    return balenc;
  }


  // uint64_t encBal
  // (
  //  const uint64_t varlen,
  //  const uint64_t leftvarlen
  //  ) const {
  //   if (varlen/2 <= leftvarlen) {
  //     return ((leftvarlen - varlen/2) << 1); // lsb is 0
  //   } else {
  //     return ((varlen/2 - leftvarlen) << 1) + 1; // lsb is 1
  //   }
  // }


  // uint64_t decLeftVarLen
  // (
  //  const uint64_t varlen,
  //  const uint64_t balenc
  //  ) const {
  //   if ((balenc & 1) == 0) { // lsb is 0
  //     return varlen/2 + (balenc >> 1);
  //   } else { // lsb is 1
  //     return varlen/2 - (balenc >> 1);
  //   }
  // }


  // std::vector<uint64_t> calcOrder
  // (
  //  const std::vector<std::pair<uint64_t, uint64_t>> & v
  //  ) {
  //   static constexpr uint64_t limit = 1024;

  //   const uint64_t num = v.size();
  //   const std::vector<uint64_t> beg;
  //   const std::vector<uint64_t> len;
  //   beg.push_back(0);
  //   len.push_back(v[0].first);
  //   for (uint64_t i = 1; i < num; ++i) {
  //     if (v[i-1].first != v[i].first) {
  //       beg.push_back(i);
  //       len.push_back(v[i].first);
  //     }
  //   }
  //   beg.push_back(num);
  //   const uint64_t numDistPair = len.size();

  //   for (uint64_t i = 0; i < )
  //   std::map<uint64_t, std::uint64_t> w;
  //   std::map<uint64_t, uint64_t> n;
  //   for (uint64_t i = 0; i < num; ++i) {
  //     auto itr = w.find(v[i].first);
  //     if (itr == w.end()) {
  //       w.insert(std::pair(1, v[i].second));
  //     } else {
  //       itr->first += 1;
  //       itr->first += (v[i].second);
  //     }
  //   }
  //   for (uint64_t i = 0; i < v.size(); ++i) {
  //     w[v[i].first] += v[i].second;
  //     ++(n[v[i].first]);
  //   }
  //   return order;
  // }


  void makeShapedSlp
  (
   const NaiveSlp<var_t> & slp,
   const bool freqSort = true
   ) {
    alph_.resize(slp.getAlphSize());
    for (uint64_t i = 0; i < slp.getAlphSize(); ++i) {
      alph_[i] = slp.getChar(i);
    }

    std::vector<uint64_t> slplen(slp.getNumRules());
    slp.makeLenVec(slplen); // expansion lengths

    { // construct prefix sum data structure
      sdsl::int_vector<64> psum(slp.getLenSeq());
      uint64_t s = 0;
      for (uint64_t i = 0; i < slp.getLenSeq(); ++i) {
        s += slp.getLenOfVar(slp.getSeq(i), slplen);
        psum[i] = s;
      }
      seqSBV_ = std::move(sdsl::sd_vector<>(psum.begin(), psum.end()));
      seqRank_.set_vector(&seqSBV_);
      seqSel_.set_vector(&seqSBV_);
    }

    { // build minimal perfect hash
      std::set<uint64_t> distLen;
      for (uint64_t i = 0; i < slplen.size(); ++i) {
        distLen.insert(slplen[i]);
      }
      std::vector<std::string> keys;
      for (auto itr = distLen.begin(); itr != distLen.end(); ++itr) {
        keys.push_back(uint2Str(*itr));
      }
      rs_ = new sux::function::RecSplit<kLeaf>(keys, kBucketSize);
    }

    std::vector<var_t> slpOrder(slp.getNumRules());
    std::vector<var_t> slpOffset(slp.getNumRules());
    {
      for (uint64_t i = 0; i < slp.getNumRules(); ++i) {
        slpOrder[i] = i;
      }
      std::vector<uint64_t> ruleFreq(slp.getNumRules());
      std::vector<uint64_t> alphFreq(slp.getAlphSize());
      slp.makeFreqInRulesVec(ruleFreq, alphFreq);
      std::stable_sort
        (
         slpOrder.begin(),
         slpOrder.end(),
         [&](uint64_t x, uint64_t y) { return ruleFreq[x] > ruleFreq[y]; }
         );
      std::stable_sort
        (
         slpOrder.begin(),
         slpOrder.end(),
         [&](uint64_t x, uint64_t y) { return slp.getLenOfVar(slp.getLeft(x), slplen) < slp.getLenOfVar(slp.getLeft(y), slplen); }
         );
      std::vector<var_t> hashVal(slp.getNumRules());
      for (uint64_t i = 0; i < slp.getNumRules(); ++i) {
        hashVal[i] = hashLen(slplen[i]);
      }
      std::stable_sort
        (
         slpOrder.begin(),
         slpOrder.end(),
         [&](uint64_t x, uint64_t y) { return hashVal[x] < hashVal[y]; }
         );

      { // sort expansion-length pairs by the frequencies of the most frequent elements
        std::vector<std::pair<uint64_t, uint64_t>> v;
        uint64_t beg = 0;
        uint64_t end = 0;
        uint64_t offset = 0;
        for (uint64_t i = 0; i < slp.getNumRules(); ++i) {
          const uint64_t id = slpOrder[i];
          if (offset++ == 0) {
            v.push_back(std::make_pair(slp.getLenOfVar(slp.getLeft(id), slplen), ruleFreq[id]));
          }
          ++end;
          if (i < slp.getNumRules() - 1 and hashVal[id] == hashVal[slpOrder[i+1]]) {
            if (slp.getLenOfVar(slp.getLeft(id), slplen) != slp.getLenOfVar(slp.getLeft(slpOrder[i+1]), slplen)) {
              offset = 0;
            }
            continue;
          }

          std::sort
            (
             v.begin(),
             v.end(),
             [](auto x, auto y) { return x.second > y.second; }
             );
          std::map<uint64_t, uint64_t> m;
          for (uint64_t i = 0; i < v.size(); ++i) {
            m.insert(std::make_pair(v[i].first, i));
          }
          std::stable_sort
            (
             slpOrder.begin() + beg,
             slpOrder.begin() + end,
             [&](uint64_t x, uint64_t y) { return m[slp.getLenOfVar(slp.getLeft(x), slplen)] < m[slp.getLenOfVar(slp.getLeft(y), slplen)]; }
             );
          m.clear();
          v.clear();
          offset = 0;
          beg = end = i + 1;
        }
      }

      sdsl::bit_vector slpDiv(slp.getNumRules(), 0);
      uint64_t numZeros = 0;
      uint64_t offset = 0;
      slpDiv[0] = 1;
      slpOffset[slpOrder[0]] = offset++;
      for (uint64_t i = 1; i < slp.getNumRules(); ++i) {
        if (hashVal[slpOrder[i]] == hashVal[slpOrder[i-1]]) {
          slpDiv[i] = 0;
          ++numZeros;
        } else {
          slpDiv[i] = 1;
          offset = 0;
        }
        slpOffset[slpOrder[i]] = offset++;
      }
      slpDivSel_.init(std::move(slpDiv));

      balBv_ = sdsl::bit_vector(numZeros, 0);
      numZeros = 0;
      for (uint64_t i = 1; i < slp.getNumRules(); ++i) {
        const uint64_t prev = slpOrder[i-1];
        const uint64_t cur = slpOrder[i];
        if (hashVal[prev] == hashVal[cur]) {
          balBv_[numZeros++] = (slp.getLenOfVar(slp.getLeft(prev), slplen) != slp.getLenOfVar(slp.getLeft(cur), slplen));
        }
      }
      balBvRank_ = std::move(sdsl::rank_support_v5<>(&balBv_));

#ifdef PRINT_STATUS_SelfShapedSlp
      {
        cout << "height = " << slp.calcHeight() << std::endl;
      }
      {
        std::ofstream ofs("distr_pairs.csv");
        for (uint64_t i = 0; i < slpOrder.size(); ++i) {
          const uint64_t slpId = slpOrder[i];
          ofs << slplen[slpId] << "," << slp.getLenOfVar(slp.getLeft(slpId), slplen) << "," << ruleFreq[slpId] << "," << slpOffset[slpId] << std::endl;
        }
      }      
#endif
    }

    {
      const uint64_t dfsize = 2 * slp.getNumRules();
      std::vector<uint64_t> df(dfsize);
      std::vector<uint64_t> bal;
      uint64_t dfPos = 0;
      uint64_t balPos = 0;
      uint64_t zPos = 0;
      for (uint64_t pos = 0; pos < slp.getNumRules(); ++pos) {
        const uint64_t slpRuleId = slpOrder[pos]; // in [0..slp.getNumRules())
        const uint64_t slpLeftVar = slp.getLeft(slpRuleId); // in [0..slp.getNumRules() + slp.getAlphSize())
        const uint64_t slpRightVar = slp.getRight(slpRuleId); // in [0..slp.getNumRules() + slp.getAlphSize())
        df[dfPos++] = (slpLeftVar < slp.getAlphSize()) ? slpLeftVar : slpOffset[slpLeftVar - slp.getAlphSize()];
        df[dfPos++] = (slpRightVar < slp.getAlphSize()) ? slpRightVar : slpOffset[slpRightVar - slp.getAlphSize()];
        
        if (slpDivSel_[pos] || balBv_[zPos++]) {
          const uint64_t varlen = slplen[slpRuleId];
          const uint64_t leftvarlen = slp.getLenOfVar(slpLeftVar, slplen);
          bal.push_back(encBal(varlen, leftvarlen));
        }
      }
      vlc_.init(df);
      bal_.init(bal);
    }

    {
      const uint64_t dfsize = slp.getLenSeq();
      std::vector<uint64_t> df(dfsize);
      for (uint64_t pos = 0; pos < dfsize; ++pos) {
        const uint64_t slpVar = slp.getSeq(pos);
        df[pos] = (slpVar < slp.getAlphSize()) ? slpVar : slpOffset[slpVar - slp.getAlphSize()];
      }
      vlcSeq_.init(df);
    }
  }
};

#endif

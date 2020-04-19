#ifndef INCLUDE_GUARD_PoSlp
#define INCLUDE_GUARD_PoSlp

#include <stdint.h> // include uint64_t etc.
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <stack>
#include <sdsl/bit_vectors.hpp>
#include "Common.hpp"
#include "NaiveSlp.hpp"
#include "fbtree.h"

// #define DEBUGBUG

template
<
  typename tparam_var_t,
  class DacT
  >
class PoSlp
{
public:
  //// Public constant, alias etc.
  using var_t = tparam_var_t;


private:
  size_t len_;
  size_t numRules_;
  var_t startVar_;
  std::vector<char> alph_;
  comp::FBTREE fbt_;
  DacT leaves_;
  sdsl::sd_vector<> pa_; // position array
  sdsl::sd_vector<>::select_1_type paSel_;


public:
  PoSlp () : len_(0),
             numRules_(0)
  {}


  ~PoSlp()
  {}


  size_t getLen() const {
    return len_;
  }


  size_t getLenSeq() const {
    return 1;
  }


  size_t getNumRules() const {
    return numRules_;
  }


  size_t getAlphSize() const {
    return alph_.size();
  }


  char getChar(uint64_t i) const {
    return alph_[i];
  }


  var_t getLeft
  (
   const uint64_t ruleId
   ) const {
    // std::cout << ruleId << ": insel = " << fbt_.InSelect(ruleId + 1) << ", left_insel = " << fbt_.Left(fbt_.InSelect(ruleId + 1)) << std::endl;
    const int64_t treePos = fbt_.Left(fbt_.InSelect(ruleId + 1));
    if (fbt_.IsLeaf(treePos)) {
      return leaves_[fbt_.LeafRank(treePos) - 1];
    } else {
      return getAlphSize() + fbt_.InRank(treePos) - 1;
    }
  }


  var_t getRight
  (
   const uint64_t ruleId
   ) const {
    // std::cout << ruleId << ": insel = " << fbt_.InSelect(ruleId + 1) << ", right_insel = " << fbt_.Right(fbt_.InSelect(ruleId + 1)) << std::endl;
    const int64_t treePos = fbt_.Right(fbt_.InSelect(ruleId + 1));
    if (fbt_.IsLeaf(treePos)) {
      return leaves_[fbt_.LeafRank(treePos) - 1];
    } else {
      return getAlphSize() + fbt_.InRank(treePos) - 1;
    }
  }


  var_t getSeq(uint64_t i) const {
    return startVar_;
  }


  size_t getNumRulesOfSlp() const {
    return numRules_;
  }


  char charAt
  (
   const uint64_t pos //!< 0-based position
   ) const {
    return charAt(pos, getLen(), 2 * getNumRulesOfSlp());
  }


  char charAt
  (
   const uint64_t pos, //!< relative position in a variable
   const uint64_t len, //!< expansion length of the variable
   uint64_t treePos //!< 0-based position in fbtree of a variable
   ) const {
    assert(pos < len);
    // std::cout << "pos = " << pos << ", len = " << len << ", stgOffset = " << stgOffset << ", slpOffset = " << slpOffset << std::endl;

    if (len == 1) {
      return getChar(leaves_[fbt_.LeafRank(treePos) - 1]);
    }
    if (fbt_.IsLeaf(treePos)) {
      treePos = fbt_.InSelect(leaves_[fbt_.LeafRank(treePos) - 1] - getAlphSize() + 1);
    }
    //// assume that treePos is not a leaf
    const uint64_t leftTreePos = fbt_.Left(treePos);
    const uint64_t leftLen = len - (accPa(treePos) - accPa(leftTreePos));
    if (pos < leftLen) {
      return charAt(pos, leftLen, leftTreePos);
    } else {
      return charAt(pos - leftLen, len - leftLen, treePos - 1);
    }
  }


  uint64_t accPa(const uint64_t treePos) const {
    return paSel_(fbt_.LeafRank(fbt_.RmostLeaf(treePos)));
  }


  void expandSubstr
  (
   const uint64_t pos, //!< beginning position
   const uint64_t lenExpand, //!< length to expand
   char * str //!< [out] must have length at least 'len'
   ) const {
    assert(pos < getLen());
    assert(lenExpand > 0);
    assert(lenExpand <= getLen() - pos);

    expandSubstr(pos, lenExpand, str, getLen(), 2 * getNumRulesOfSlp());
  }


  void expandSubstr
  (
   const uint64_t pos, //!< relative position in a variable
   const uint64_t lenExpand, //!< length to expand
   char * str, //!< [out] must have length at least 'len'
   const uint64_t varLen, //!< expansion length of the variable
   uint64_t treePos //!< 0-based position in fbtree of a variable
   ) const {
    // std::cout << "pos = " << pos << ", len = " << len << ", varLen = " << varLen << ", stgOffset = " << stgOffset << ", slpOffset = " << slpOffset << std::endl;
    assert(lenExpand > 0);
    assert(pos < varLen);

    if (varLen == 1) {
      *str = getChar(leaves_[fbt_.LeafRank(treePos) - 1]);
      return;
    }
    if (fbt_.IsLeaf(treePos)) {
      treePos = fbt_.InSelect(leaves_[fbt_.LeafRank(treePos) - 1] - getAlphSize() + 1);
    }
    //// assume that treePos is not a leaf
    const uint64_t leftTreePos = fbt_.Left(treePos);
    const uint64_t leftLen = varLen - (accPa(treePos) - accPa(leftTreePos));
    if (pos < leftLen) {
      expandSubstr(pos, lenExpand, str, leftLen, leftTreePos);
      if (leftLen - pos < lenExpand) {
        expandPref(lenExpand - (leftLen - pos), str + (leftLen - pos), varLen - leftLen, treePos - 1);
      }
    } else {
      expandSubstr(pos - leftLen, lenExpand, str, varLen - leftLen, treePos - 1);
    }
  }


  void expandPref
  (
   const uint64_t lenExpand, //!< length to expand
   char * str, //!< [out] must have length at least 'len'
   const uint64_t varLen, //!< expansion length of the variable
   uint64_t treePos //!< 0-based position in fbtree of a variable
   ) const {
    // std::cout << "len = " << len << ", varLen = " << varLen << ", stgOffset = " << stgOffset << ", slpOffset = " << slpOffset << std::endl;
    assert(lenExpand > 0);

    if (varLen == 1) {
      *str = getChar(leaves_[fbt_.LeafRank(treePos) - 1]);
      return;
    }
    if (fbt_.IsLeaf(treePos)) {
      treePos = fbt_.InSelect(leaves_[fbt_.LeafRank(treePos) - 1] - getAlphSize() + 1);
    }
    //// assume that treePos is not a leaf
    const uint64_t leftTreePos = fbt_.Left(treePos);
    const uint64_t leftLen = varLen - (accPa(treePos) - accPa(leftTreePos));
    expandPref(lenExpand, str, leftLen, leftTreePos);
    if (lenExpand > leftLen) {
      expandPref(lenExpand - leftLen, str + leftLen, varLen - leftLen, treePos - 1);
    }
  }


  void printStatus
  (
   const bool verbose = false
   ) const {
    std::cout << "PoSlp object (" << this << ") " << __func__ << "(" << verbose << ") BEGIN" << std::endl;
    std::cout << "alphSize = " << getAlphSize() << ", len = " << getLen() << ", lenSeq = " << getLenSeq()
              << ", numRulesOfSlp = " << getNumRulesOfSlp()
              << std::endl;
    const size_t bytesAlph = sizeof(std::vector<char>) + (sizeof(char) * alph_.size());
    const size_t bytesFbt = fbt_.Space();
    const size_t bytesLeaves = leaves_.calcMemBytes();
    const size_t bytesPositionArray = (pa_.size()) ? sdsl::size_in_bytes(pa_) : 0;
    std::cout << "DS sizes (bytes)" << std::endl;
    std::cout << bytesAlph + bytesFbt + bytesLeaves + bytesPositionArray << std::endl;
    std::cout << "| alph = " << bytesAlph << std::endl;
    std::cout << "| fbtree = " << bytesFbt << std::endl;
    std::cout << "| leaves = " << bytesLeaves << std::endl;
    std::cout << "| position = " << bytesPositionArray << std::endl;
    if (verbose) {
      std::cout << "bp: " << std::endl;
      for (uint64_t i = 0; i < fbt_.Length(); ++i) {
        std::cout << fbt_.Get(i);
      }
      std::cout << std::endl;
      std::cout << "leaves: " << std::endl;
      for (uint64_t i = 0; i < getNumRulesOfSlp() + 1; ++i) {
        std::cout << i << "(" << leaves_[i] << ") ";
      }
      std::cout << std::endl;
      std::cout << "rules: " << std::endl;
      for (uint64_t i = 0; i < getNumRulesOfSlp(); ++i) {
        const auto l = getLeft(i);
        const auto r = getRight(i);
        std::cout << i + getAlphSize() << "(";
        if (l < getAlphSize()) {
          std::cout << "'" << getChar(l) << "'";
        } else {
          std::cout << l;
        }
        std::cout << ", ";
        if (r < getAlphSize()) {
          std::cout << "'" << getChar(r) << "'";
        } else {
          std::cout << r;
        }
        std::cout << ") ";
      }
      std::cout << std::endl;
      std::cout << "position array: " << std::endl;
      for (uint64_t i = 0; i < getNumRulesOfSlp() + 1; ++i) {
        std::cout << paSel_(i + 1) << " ";
      }
      std::cout << std::endl;
    }
    std::cout << "PoSlp object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
  }


  size_t calcMemBytes() const {
    size_t ret = 0;
    ret += sizeof(std::vector<char>) + (sizeof(char) * alph_.size());
    ret += fbt_.Space();
    ret += leaves_.calcMemBytes();
    ret += sdsl::size_in_bytes(pa_);
    return ret;
  }


  void load
  (
   std::istream & in
   ) {
    in.read((char*) & len_, sizeof(len_));
    in.read((char*) & numRules_, sizeof(numRules_));
    in.read((char*) & startVar_, sizeof(startVar_));
    uint64_t alphSize = 0;
    in.read((char*) & alphSize, sizeof(alphSize));
    alph_.resize(alphSize);
    in.read((char*) alph_.data(), alphSize * sizeof(alph_[0]));
    fbt_.Load(in);
    leaves_.load(in);
    pa_.load(in);
    paSel_.load(in);

    paSel_.set_vector(&pa_);
  }


  void serialize
  (
   std::ostream & out
   ) const {
    out.write((char*) & len_, sizeof(len_));
    out.write((char*) & numRules_, sizeof(numRules_));
    out.write((char*) & startVar_, sizeof(startVar_));
    uint64_t alphSize = getAlphSize();
    out.write((char*) & alphSize, sizeof(alphSize));
    out.write((char*) alph_.data(), alphSize * sizeof(alph_[0]));
    fbt_.Save(out);
    leaves_.serialize(out);
    pa_.serialize(out);
    paSel_.serialize(out);
  }


  void makePoSlp
  (
   const NaiveSlp<var_t> & slp
   ) {
    assert(slp.getLenSeq() == 1);

    const uint64_t alphSize = slp.getAlphSize();
    alph_.resize(alphSize);
    for (uint64_t i = 0; i < alphSize; ++i) {
      alph_[i] = slp.getChar(i);
    }
    numRules_ = slp.getNumRules();
    std::vector<uint64_t> leaves(slp.getNumRules() + 1);
    fbt_.Build(1024, 2, 2 * slp.getNumRules() + 2);
    std::vector<uint64_t> postRename(slp.getNumRules()); // rename by post-order
    for (uint64_t i = 0; i < postRename.size(); ++i) {
      postRename[i] = UINT64_MAX;
    }
    std::stack<std::pair<var_t, uint64_t> > st;
    st.push(std::make_pair(slp.getSeq(0), 0));
    uint64_t leafPos = 0;
    uint64_t newName = alphSize;
    while (!st.empty()) {
      auto & e = st.top();
      if (e.first < alphSize || postRename[e.first - alphSize] != UINT64_MAX) { // leaf
        fbt_.PushBack(comp::kOP);
        leaves[leafPos++] = (e.first < alphSize) ? e.first : postRename[e.first - alphSize];
        st.pop();
      } else if (e.second == 0) {
        e.second = 1;
        st.push(std::make_pair(slp.getLeft(e.first - alphSize), 0));
      } else if (e.second == 1) {
        e.second = 2;
        st.push(std::make_pair(slp.getRight(e.first - alphSize), 0));
      } else {
        fbt_.PushBack(comp::kCP);
        postRename[e.first - alphSize] = newName++;
        st.pop();
      }
    }
    fbt_.PushBack(comp::kCP);
    leaves_.init(leaves);
    startVar_ = newName - 1;
  }


  void makePoSlpWithLen
  (
   const NaiveSlp<var_t> & slp
   ) {
    assert(slp.getLenSeq() == 1);

    std::vector<uint64_t> lenOfRule(slp.getNumRules());
    slp.makeLenVec(lenOfRule);
    std::vector<uint64_t> psum(slp.getNumRules() + 1);

    const uint64_t alphSize = slp.getAlphSize();
    alph_.resize(alphSize);
    for (uint64_t i = 0; i < alphSize; ++i) {
      alph_[i] = slp.getChar(i);
    }
    numRules_ = slp.getNumRules();
    std::vector<uint64_t> leaves(slp.getNumRules() + 1);
    fbt_.Build(1024, 2, 2 * slp.getNumRules() + 2);
    std::vector<uint64_t> postRename(slp.getNumRules()); // rename by post-order
    for (uint64_t i = 0; i < postRename.size(); ++i) {
      postRename[i] = UINT64_MAX;
    }
    std::stack<std::pair<var_t, uint64_t> > st;
    st.push(std::make_pair(slp.getSeq(0), 0));
    uint64_t leafPos = 0;
    uint64_t newName = alphSize;
    uint64_t sumLen = 0;
    while (!st.empty()) {
      auto & e = st.top();
      if (e.first < alphSize || postRename[e.first - alphSize] != UINT64_MAX) { // leaf
        fbt_.PushBack(comp::kOP);
        sumLen += (e.first < alphSize) ? 1 : lenOfRule[e.first - alphSize];
        psum[leafPos] = sumLen;
        leaves[leafPos++] = (e.first < alphSize) ? e.first : postRename[e.first - alphSize];
        st.pop();
      } else if (e.second == 0) {
        e.second = 1;
        st.push(std::make_pair(slp.getLeft(e.first - alphSize), 0));
      } else if (e.second == 1) {
        e.second = 2;
        st.push(std::make_pair(slp.getRight(e.first - alphSize), 0));
      } else {
        fbt_.PushBack(comp::kCP);
        postRename[e.first - alphSize] = newName++;
        st.pop();
      }
    }
    fbt_.PushBack(comp::kCP);
    leaves_.init(leaves);
    startVar_ = newName - 1;

    // length
    len_ = sumLen;
    pa_ = std::move(sdsl::sd_vector<>(psum.begin(), psum.end()));
    paSel_.set_vector(&pa_);
  }



};

#endif

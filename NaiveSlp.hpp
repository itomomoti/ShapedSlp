/*!
 * Copyright (c) 2017 Tomohiro I
 *
 * This program is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 */
/*!
 * @file NaiveSlp.hpp
 * @brief 
 * @author Tomohiro I
 * @date 2019-10-15
 */
#ifndef INCLUDE_GUARD_NaiveSlp
#define INCLUDE_GUARD_NaiveSlp

#include <sys/stat.h>
#include <stdint.h> // include uint64_t etc.
#include <assert.h>
#include <cmath>
#include <cstring>
#include <iostream>
#include "Common.hpp"


/*!
 * @brief NaiveSlp
 */
template<typename tparam_var_t>
class NaiveSlp
{
public:
  // Public constant, alias etc.
  using var_t = tparam_var_t;


  std::vector<var_t> seq_; // final sequence of text
  std::vector<Tpair<var_t>> rules_; // rules
  std::vector<char> alph_; // map from terminal IDs to actual characters


public:
  NaiveSlp() {};


  ~NaiveSlp() {};


  NaiveSlp
  (
   const NaiveSlp & other
   ) : seq_(other.seq_), rules_(other.rules_), alph_(other.alph_) {
  }


  void load_NavarroRepair
  (
   const char * fname_base
   ) {
    char fname[1024];
    FILE *Rf, *Cf;
    struct stat s;

    strcpy(fname, fname_base);
    strcat(fname, ".R");
    if (stat(fname, &s) != 0) {
      fprintf(stderr, "Error: cannot stat file %s\n", fname);
      exit(1);
    }
    const uint64_t len = s.st_size;
    Rf = fopen(fname, "r");
    if (Rf == NULL) {
      fprintf(stderr, "Error: cannot open file %s for reading\n", fname);
      exit(1);
    }
    unsigned int alphSize;
    if (fread(&alphSize, sizeof(int), 1, Rf) != 1)	{ // read alphabet size
      fprintf(stderr, "Error: cannot read file %s\n", fname);
      exit(1);
    }
    alph_.resize(alphSize);
    if (fread(alph_.data(), sizeof(char), alphSize, Rf) != alphSize) { // read map
      fprintf(stderr, "Error: cannot read file %s\n", fname);
      exit(1);
    }
    const uint64_t numRules = (len - sizeof(int) - alphSize) / sizeof(Tpair<var_t>);
    rules_.resize(numRules);
    if (fread(rules_.data(), sizeof(Tpair<var_t>), numRules, Rf) != numRules) { // read rules
      fprintf(stderr, "Error: cannot read file %s\n", fname);
      exit(1);
    }
    fclose(Rf);

    strcpy(fname, fname_base);
    strcat(fname, ".C");
    if (stat(fname, &s) != 0) {
      fprintf(stderr, "Error: cannot stat file %s\n", fname);
      exit(1);
    }
    const uint64_t lenSeq = s.st_size / sizeof(var_t);
    seq_.resize(lenSeq);

    Cf = fopen(fname, "r");
    if (Cf == NULL) {
      fprintf(stderr, "Error: cannot open file %s for reading\n", fname);
      exit(1);
    }
    if (fread(seq_.data(), sizeof(var_t), lenSeq, Cf) != lenSeq) { // read sequence
      fprintf(stderr, "Error: cannot read file %s\n", fname);
      exit(1);
    }
    fclose(Cf);
  }


  void load_Bigrepair
  (
   const char * fname_base
   ) {
    char fname[1024];
    FILE *Rf, *Cf;
    struct stat s;

    strcpy(fname, fname_base);
    strcat(fname, ".R");
    if (stat(fname, &s) != 0) {
      fprintf(stderr, "Error: cannot stat file %s\n", fname);
      exit(1);
    }
    const uint64_t len = s.st_size;
    Rf = fopen(fname, "r");
    if (Rf == NULL) {
      fprintf(stderr, "Error: cannot open file %s for reading\n", fname);
      exit(1);
    }
    unsigned int alphSize;
    if (fread(&alphSize, sizeof(int), 1, Rf) != 1)	{ // read alphabet size
      fprintf(stderr, "Error: cannot read file %s\n", fname);
      exit(1);
    }
    alphSize = 256; // [0,256) are terminals
    alph_.resize(alphSize);
    for (uint64_t i = 0; i < alph_.size(); ++i) {
      alph_[i] = i;
    }
    const uint64_t numRules = (len - sizeof(int)) / sizeof(Tpair<var_t>);
    rules_.resize(numRules);
    if (fread(rules_.data(), sizeof(Tpair<var_t>), numRules, Rf) != numRules) { // read rules
      fprintf(stderr, "Error: cannot read file %s\n", fname);
      exit(1);
    }
    fclose(Rf);

    strcpy(fname, fname_base);
    strcat(fname, ".C");
    if (stat(fname, &s) != 0) {
      fprintf(stderr, "Error: cannot stat file %s\n", fname);
      exit(1);
    }
    const uint64_t lenSeq = s.st_size / sizeof(var_t);
    seq_.resize(lenSeq);

    Cf = fopen(fname, "r");
    if (Cf == NULL) {
      fprintf(stderr, "Error: cannot open file %s for reading\n", fname);
      exit(1);
    }
    if (fread(seq_.data(), sizeof(var_t), lenSeq, Cf) != lenSeq) { // read sequence
      fprintf(stderr, "Error: cannot read file %s\n", fname);
      exit(1);
    }
    fclose(Cf);
  }


  size_t getLenSeq() const {
    return seq_.size();
  }


  size_t getNumRules() const {
    return rules_.size();
  }


  size_t getAlphSize() const {
    return alph_.size();
  }


  char getChar(uint64_t i) const {
    return alph_[i];
  }


  var_t getLeft(uint64_t i) const {
    return rules_[i].left;
  }


  var_t getRight(uint64_t i) const {
    return rules_[i].right;
  }


  var_t getSeq(uint64_t i) const {
    return seq_[i];
  }


  void setLenSeq(size_t n) {
    return seq_.resize(n);
  }


  void setNumRules(size_t n) {
    return rules_.resize(n);
  }


  void setAlphSize(size_t n) {
    alph_.resize(n);
  }


  void setChar(uint64_t i, char c) {
    alph_[i] = c;
  }


  void pushPair(Tpair<var_t> p) {
    rules_.push_back(p);
  }


  void setLeft(uint64_t i, var_t v) {
    rules_[i].left = v;
  }


  void setRight(uint64_t i, var_t v) {
    rules_[i].right = v;
  }


  void setRule(uint64_t i, Tpair<var_t> p) {
    setLeft(i, p.left);
    setRight(i, p.right);
  }


  void setSeq(uint64_t i, var_t v) {
    seq_[i] = v;
  }


  void printStatus
  (
   const bool verbose = false
   ) const {
    std::cout << "NaiveSlp object (" << this << ") " << __func__ << "(" << verbose << ") BEGIN" << std::endl;

    uint64_t ls = getLenSeq();
    uint64_t nr = getNumRules();
    uint64_t as = getAlphSize();
    std::cout << "lenSeq_ = " << ls << std::endl;
    std::cout << "numRules_ = " << nr << std::endl;
    std::cout << "alphSize_ = " << as << std::endl;
    std::cout << "Estimated space for POSLP encoding (bytes) = "
              << estimateEncSize()
              << std::endl;
    if (verbose) {
      for (uint64_t i = 0; i < getAlphSize(); ++i) {
        std::cout << getChar(i) << " ";
      }
      std::cout << std::endl;
      std::cout << "lenSeq_ = " << getLenSeq() << std::endl;
      for (uint64_t i = 0; i < getLenSeq(); ++i) {
        std::cout << getSeq(i) << " ";
      }
      std::cout << std::endl;
      std::cout << "numRules_ = " << getNumRules() << std::endl;
      for (uint64_t i = 0; i < getNumRules(); ++i) {
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
    }
    std::cout << "NaiveSlp object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
  }


  uint64_t calcHeight() const {
    std::vector<uint64_t> hvec(getNumRules());
    uint64_t ret = 0;
    for (uint64_t i = 0; i < getNumRules(); ++i) {
      const uint64_t lh = (rules_[i].left < getAlphSize()) ? 1 : hvec[rules_[i].left - getAlphSize()] + 1;
      const uint64_t rh = (rules_[i].right < getAlphSize()) ? 1 : hvec[rules_[i].right - getAlphSize()] + 1;
      hvec[i] = std::max(lh, rh);
      ret = std::max(hvec[i], ret);
    }
    return ret;
  }


  void makeLenVec
  (
   std::vector<uint64_t> & lenVec
   ) const {
    for (uint64_t i = 0; i < getNumRules(); ++i) {
      const uint64_t lLen = (rules_[i].left < getAlphSize()) ? 1 : lenVec[rules_[i].left - getAlphSize()];
      const uint64_t rLen = (rules_[i].right < getAlphSize()) ? 1 : lenVec[rules_[i].right - getAlphSize()];
      lenVec[i] = lLen + rLen;
    }
  }


  void makeFreqInRulesVec
  (
   std::vector<uint64_t> & ruleFreqVec,
   std::vector<uint64_t> & alphFreqVec
   ) const {
    for (uint64_t i = 0; i < getNumRules(); ++i) {
      ruleFreqVec[i] = 0;
    }
    for (uint64_t i = 0; i < getAlphSize(); ++i) {
      alphFreqVec[i] = 0;
    }
    for (uint64_t i = 0; i < getLenSeq(); ++i) {
      uint64_t v = getSeq(i);
      if (v < getAlphSize()) {
        alphFreqVec[v]++;
      } else {
        ruleFreqVec[v - getAlphSize()]++;
      }
    }
    for (uint64_t i = 0; i < getNumRules(); ++i) {
      uint64_t v = rules_[i].left;
      if (v < getAlphSize()) {
        alphFreqVec[v]++;
      } else {
        ruleFreqVec[v - getAlphSize()]++;
      }
      v = rules_[i].right;
      if (v < getAlphSize()) {
        alphFreqVec[v]++;
      } else {
        ruleFreqVec[v - getAlphSize()]++;
      }
    }
  }


  void makeBinaryTree() {
    while (seq_.size() > 1) {
      const uint64_t len = seq_.size();
      const uint64_t n = len / 2;
      const bool r = len % 2;
      for (uint64_t i = 0; i < n; ++i) {
        const var_t newId = getNumRules() + getAlphSize();
        Tpair<var_t> p;
        p.left = seq_[2 * i];
        p.right = seq_[2 * i + 1];
        pushPair(p);
        seq_[i] = newId;
      }
      if (r) {
        seq_[n] = seq_[len - 1];
      }
      seq_.resize(n + r);
    }
  }


  uint64_t getLenOfVar
  (
   const uint64_t var,
   const std::vector<uint64_t> & lenVec
   ) const {
    assert(var < getAlphSize() + getNumRules());

    return (var < getAlphSize()) ? 1 : lenVec[var - getAlphSize()];
  }


  void printExpand
  (
   const uint64_t var,
   const std::vector<uint64_t> & lenVec,
   const uint64_t pad
   ) const {
    padVLine(pad);
    std::cout << "var = " << var << "(" << (var < getAlphSize()) << ") len = " << getLenOfVar(var, lenVec) << std::endl;
    if (var < getAlphSize()) {
      padVLine(pad);
      std::cout << getChar(var) << std::endl;
    } else {
      padVLine(pad);
      std::cout << "left: " << std::endl;
      printExpand(getLeft(var - getAlphSize()), lenVec, pad + 1);
      padVLine(pad);
      std::cout << "right: " << std::endl;
      printExpand(getRight(var - getAlphSize()), lenVec, pad + 1);
    }
  }


  size_t estimateEncSize() const {
    // |G| lg(|G| + σ) + 2|G| + o(|G|) bits
    // in Theorem 3 of the paper "Fully-Online Grammar Compression”, SPIRE 2013
    uint64_t ls = getLenSeq();
    uint64_t nr = getNumRules();
    uint64_t as = getAlphSize();
    return ((nr + ls) * ceilLog2(nr + as) + 2 * (nr + ls)) / 8;
  }


  size_t estimateEncSizeWithLen(const size_t totalLen) const {
    // |G| lg |S| + |G| lg(1 + σ/|G|) + 5|G| + o(|G|) bits
    // in Theorem 3 of the paper "Fully-Online Grammar Compression”, SPIRE 2013
    uint64_t ls = getLenSeq();
    uint64_t nr = getNumRules();
    uint64_t as = getAlphSize();
    return ((ls + nr) * ceilLog2(totalLen) + (ls + nr) * ceil(log(1 + (double)as/(ls + nr))) + 5 * (ls + nr)) / 8;
  }
};



#endif

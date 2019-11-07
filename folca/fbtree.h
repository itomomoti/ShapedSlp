#ifndef INCLUDE_GUARD_FBTREE_H_
#define INCLUDE_GUARD_FBTREE_H_

#include <string>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include "util2.hpp"
#include "onlinebp.h"

#pragma pack(1)
namespace comp {

  class FBTREE {
  public:
  FBTREE() : bp_(NULL) {}
    void Build(int L, int alpha, u64 capaBitLen);
    void PushBack(int c);
    i64 Parent(i64 i) const;
    i64 Left(i64 i) const;
    i64 Right(i64 i) const;
    i64 InSelect(i64 i) const;
    i64 LeafSelect(i64 i) const;
    i64 InRank(i64 i) const;
    i64 LeafRank(i64 i) const;    
    i64 LmostLeaf(i64 i) const;
    i64 RmostLeaf(i64 i) const;
    bool IsLeaf(i64 i) const;
    bool IsLeft(i64 i) const;
    bool Get(i64 i) const;
    i64 Space() const;
    i64 Length() const;
    void Clear();
    void Save(std::ostream &ofs) const;
    void Load(std::istream &ifs);
  private:
    onlinebp* bp_;
  }; // class FBTREE

} // namespace comp

#pragma pack()
#endif // INCLUDE_GURAD_FBTREE_H_

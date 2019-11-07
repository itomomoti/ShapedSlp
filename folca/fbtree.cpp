#include <vector>
#include <iostream>
#include <stdint.h>
#include "fbtree.h"
#include "onlinebp.h"

using namespace std;

namespace comp {

  void FBTREE::Build(int L, int alpha, u64 capaBitLen) {
    bp_ = onlinebp_new(L, alpha);
    onlinebp_reserve(bp_, capaBitLen);
  }

  void FBTREE::PushBack(int c) {
    onlinebp_push(bp_, c);
  }

  i64 FBTREE::Parent(i64 i) const {
    return onlinebp_parent(bp_, i);
  }

  i64 FBTREE::Left(i64 i) const {
    i64 s, t;
    onlinebp_leftchild(bp_, i, &s, &t);
    return t;
  }

  i64 FBTREE::Right(i64 i) const {
    i64 s, t;
    onlinebp_rightchild(bp_, i, &s, &t);
    return t;
  }  
  
  i64 FBTREE::InSelect(i64 i) const {
    return onlinebp_select(bp_, i, kCP);
  }

  i64 FBTREE::LeafSelect(i64 i) const {
    return onlinebp_select(bp_, i, kOP);
  }

  i64 FBTREE::InRank(i64 i) const {
    return onlinebp_rank(bp_, i, kCP);
  }

  i64 FBTREE::LeafRank(i64 i) const {
    return onlinebp_rank(bp_, i, kOP);
  }

  i64 FBTREE::LmostLeaf(i64 i) const {
    if (IsLeaf(i)) return i;
    return onlinebp_leftmost_leaf(bp_, i);
  }

  i64 FBTREE::RmostLeaf(i64 i) const {
    if (IsLeaf(i)) return i;
    return onlinebp_rightmost_leaf(bp_, i);
  }

  bool FBTREE::IsLeaf(i64 i) const {
    if (onlinebp_get(bp_, i) == kOP) return true;
    return false;
  }

  bool FBTREE::IsLeft(i64 i) const {
    if(onlinebp_get(bp_,i+1) == kOP) return true;
    return false;
  }

  bool FBTREE::Get(i64 i) const {
    if(onlinebp_get(bp_,i) == kOP) return true;
    return false;
  }

  i64 FBTREE::Space() const {
    return onlinebp_size(bp_);
  }

  i64 FBTREE::Length() const {
    return onlinebp_length(bp_);
  }

  void FBTREE::Clear() {
    onlinebp_free(bp_);
  }

  void FBTREE::Save(ostream &ofs) const {
    
    vector<uint64_t> tmp_bits;
    uint64_t length = Length();
    uint64_t vec_len = (length / 64) + 1;
    
    tmp_bits.resize(vec_len);
    vector<uint64_t> (tmp_bits).swap(tmp_bits);

    for(size_t i = 0; i < vec_len; i++){
      tmp_bits[i] = 0;
    } 
   
    for(size_t i = 0; i < length; i++){
      if (Get(i)){
	tmp_bits[i / 64] |= (1LLU << (i % 64));
      } 
    }
    int L = onlinebp_L(bp_);
    int alpha = onlinebp_alpha(bp_);
    ofs.write((char*)&L,    sizeof(L));
    ofs.write((char*)&alpha, sizeof(alpha));
    ofs.write((char*)&length, sizeof(length));
    ofs.write((char*)&tmp_bits[0], sizeof(tmp_bits[0]) * vec_len);
    vector<uint64_t> ().swap(tmp_bits);
  }
  
  void FBTREE::Load(istream &ifs){
    
    uint64_t length = 0;
    uint64_t vec_len = 0;
    int L;
    int alpha;
    vector<uint64_t> tmp_bits;
    tmp_bits.clear();

    ifs.read((char*)&L, sizeof(L));
    ifs.read((char*)&alpha, sizeof(alpha));
    ifs.read((char*)&length, sizeof(length));

    vec_len = length/64 + 1;

    tmp_bits.resize(vec_len);
    vector<uint64_t> (tmp_bits).swap(tmp_bits);
    ifs.read((char*)&tmp_bits[0], sizeof(tmp_bits[0]) * vec_len); 
  
    Build(L,alpha,length);
    for(size_t i = 0; i < length; i++){
      if((tmp_bits[i / 64] & (1LLU << (i % 64))) > 0){
	PushBack(1);
      }
      else{
	PushBack(0);
      }
    }
    vector<uint64_t> ().swap(tmp_bits);
  }

} // namespace comp

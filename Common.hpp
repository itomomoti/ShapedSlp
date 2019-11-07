#ifndef INCLUDE_GUARD_Common
#define INCLUDE_GUARD_Common

#include <stdint.h> // include uint64_t etc.
#include <iostream>
#include <fstream>
#include <string>
#include <queue>
#include <stack>

template<typename var_t>
struct Tpair
{
  var_t left, right;

  bool operator<(const Tpair & another) const
  {
    return (this->left < another.left) || ((this->left == another.left) && this->right < another.right);
  };
};


void padVLine
(
 const uint64_t pad
 );


template<class type>
void printArray(type * arr, uint64_t n, std::string delimiter = " ")
{
  for (uint64_t i = 0; i < n; ++i) {
    std::cout << arr[i] << delimiter;
  }
  std::cout << std::endl;
}


template<class Container>
void printVec(const Container & vec)
{
  for (uint64_t i = 0; i < vec.size(); ++i) {
    std::cout << "(" << i << ":" << vec[i] << ") ";
  }
  std::cout << std::endl;
}


uint32_t ceilLog2(uint64_t x);


template<class SlpT>
void decompressByCharAt
(
 const SlpT & slp,
 const std::string & ofile
) {
  std::ofstream ofs(ofile);
  for (uint64_t i = 0; i < slp.getLen(); ++i) {
    char c = slp.charAt(i);
    ofs.write(&c, 1);
  }
}


template<class SlpT>
void decompressByRecurse
(
 const SlpT & slp,
 const std::string & ofile
) {
  using var_t = typename SlpT::var_t;
  const uint64_t alphSize = slp.getAlphSize();
  std::ofstream ofs(ofile);
  for (uint64_t i = 0; i < slp.getLenSeq(); ++i) {
    const var_t var = slp.getSeq(i);
    std::stack<std::pair<var_t, uint64_t> > st;
    st.push(std::make_pair(var, 0));
    while (!st.empty()) {
      auto & e = st.top();
      if (e.first < alphSize) { // leaf
        char c = slp.getChar(e.first);
        ofs.write(&c, 1);
        st.pop();
      } else if (e.second == 0) {
        e.second = 1;
        st.push(std::make_pair(slp.getLeft(e.first - alphSize), 0));
      } else if (e.second == 1) {
        e.second = 2;
        st.push(std::make_pair(slp.getRight(e.first - alphSize), 0));
      } else {
        st.pop();
      }
    }
  }
}



template<typename ArrayT>
class PackedArrayTypeValRef
{
  friend ArrayT;


private:
  ArrayT * const obj_;
  const uint64_t idx_;


  PackedArrayTypeValRef
  (
   ArrayT * obj,
   uint64_t idx
   ) :
    obj_(obj),
    idx_(idx)
  {}


public:
  uint64_t operator=
  (
   uint64_t val
   ) {
    obj_->write(val, idx_);
    return val;
  }


  operator uint64_t() const {
    return obj_->read(idx_);
  }
};




/*!
 * @tparam kBucketWidth: Bitwidth of bucket size
 * @tparam elem_t: type of element to be sorted
 * @tparam
 *   Func: Fuction that returns kBucketWidth width integer from an element
 */
template<uint8_t kBucketWidth = 8, class elem_t, class Func>
void my_bucket_sort
(
 elem_t * earray, //!< given array to be sorted by some criterion specified by func
 uint64_t n, //!< length of array
 Func func
 ) {
  const uint64_t kBS = UINT64_C(1) << kBucketWidth; // bucket size
  std::queue<elem_t> bucket[kBS];
  for (uint64_t i = 0; i < n; ++i) {
    auto e = earray[i];
    bucket[func(e)].push(e);
  }
  uint64_t i = 0;
  for (uint64_t k = 0; k < kBS; ++k) {
    while (!bucket[k].empty()) {
      elem_t e = bucket[k].front();
      bucket[k].pop();
      earray[i++] = e;
    }
  }
}


/*!
 * @tparam kBucketWidth: Bitwidth of bucket size
 * @tparam elem_t: type of element to be sorted
 * @tparam
 *   Func: Fuction that returns kBucketWidth width integer from an element
 */
template<uint8_t kBucketWidth = 8, class elem_t, class keys_t>
void my_radix_sort
(
 elem_t * earray, //!< given array to be sorted by some criterion specified by func
 keys_t * keys, //!< i \in [0..n) is sorted based on keys[i]
 uint64_t n, //!< length of array
 uint8_t keyWidth
 ) {
  const uint64_t mask = (UINT64_C(1) << kBucketWidth) - 1; // bucket size
  for (uint64_t k = 0; k < (keyWidth + kBucketWidth - 1) / kBucketWidth; ++k) {
    my_bucket_sort<8>
      (earray, n, 
       [keys, k](uint64_t i){
         return (keys[i] >> (kBucketWidth * k)) & mask;
       }
       );
  }
}



#endif

#ifndef INCLUDE_GUARD_FixedBitLenCode
#define INCLUDE_GUARD_FixedBitLenCode

#include <sdsl/bit_vectors.hpp>
#include "Common.hpp"




template
<
  uint8_t t_width = 0
  >
class FixedBitLenCode
{
private:
  sdsl::int_vector<t_width> vec_;


public:
  FixedBitLenCode
  ()
  {}


  ~FixedBitLenCode()
  {}


  template<class vecT>
  void init
  (
   const vecT & vec
   ) {
    if (t_width == 0) {
      uint8_t w = 1;
      for (uint64_t i = 0; i < vec.size(); ++i) {
        w = std::max(w, static_cast<uint8_t>(sdsl::bits::hi(vec[i]) + 1));
      }
      vec_.width(w);
    }
    vec_.resize(vec.size());
    for (uint64_t i = 0; i < vec.size(); ++i) {
      vec_[i] = vec[i];
    }
  }


  /*!
   * @brief Read only accessor.
   */
  uint64_t operator[]
  (
   const size_t idx //!< in [0, num_).
   ) const {
    assert(idx < num_);

    return this->read(idx);
  }


  /*!
   * @brief Read value at 'idx'.
   */
  uint64_t read
  (
   const size_t idx //!< in [0, num_)
   ) const {
    assert(idx < num_);

    return vec_[idx];
  }


  /*!
   * @brief Get size.
   */
  size_t size() const noexcept {
    return vec_.size();
  }


  /*!
   * @brief Calculate total memory usage in bytes.
   */
  size_t calcMemBytes() const noexcept {
    size_t ret = sizeof(*this);
    ret += sdsl::size_in_bytes(vec_);
    return ret;
  }


  void load
  (
   std::istream & in
   ) {
    vec_.load(in);
  }


  void serialize
  (
   std::ostream & out
   ) const {
    vec_.serialize(out);
  }


  void printStatus
  (
   const bool verbose = false
   ) const noexcept {
    std::cout << "IntVec object (" << this << ") " << __func__ << "(" << verbose << ") BEGIN" << std::endl;
    if (verbose) {
      const auto size = this->size();
      std::cout << "dump stored values" << std::endl;
      for (uint64_t i = 0; i < size; ++i) {
        std::cout << this->read(i) << ", ";
      }
      std::cout << std::endl;
    }
    std::cout << "IntVec object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
  }
};


#endif

#ifndef INCLUDE_GUARD_VlcVec
#define INCLUDE_GUARD_VlcVec

#include <sdsl/bit_vectors.hpp>
#include "Common.hpp"

template
<
  class t_coder = sdsl::coder::elias_delta,
  uint32_t t_dens = 128
  >
class VlcVec
{
private:
  sdsl::vlc_vector<t_coder, t_dens> vlc_;


public:
  VlcVec
  ()
  {}


  ~VlcVec()
  {}


  template<class vecT>
  void init
  (
   const vecT & vec
   ) {
    vlc_ = std::move(sdsl::vlc_vector<t_coder>(vec));
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

    return vlc_[idx];
  }


  /*!
   * @brief Get size.
   */
  size_t size() const noexcept {
    return vlc_.size();
  }


  /*!
   * @brief Calculate total memory usage in bytes.
   */
  size_t calcMemBytes() const noexcept {
    size_t ret = sizeof(*this);
    ret += sdsl::size_in_bytes(vlc_);
    return ret;
  }


  void load
  (
   std::istream & in
   ) {
    vlc_.load(in);
  }


  void serialize
  (
   std::ostream & out
   ) const {
    vlc_.serialize(out);
  }


  void printStatus
  (
   const bool verbose = false
   ) const noexcept {
    std::cout << "VlcVec object (" << this << ") " << __func__ << "(" << verbose << ") BEGIN" << std::endl;
    if (verbose) {
      const auto size = this->size();
      std::cout << "dump stored values" << std::endl;
      for (uint64_t i = 0; i < size; ++i) {
        std::cout << this->read(i) << ", ";
      }
      std::cout << std::endl;
    }
    std::cout << "VlcVec object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
  }
};




template
<
  uint8_t t_width = 0
  >
class IntVec
{
private:
  sdsl::int_vector<t_width> vec_;


public:
  IntVec
  ()
  {}


  ~IntVec()
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

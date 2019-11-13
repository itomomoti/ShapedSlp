#ifndef INCLUDE_GUARD_DirectAccessibleGammaCode
#define INCLUDE_GUARD_DirectAccessibleGammaCode

#include <sdsl/bit_vectors.hpp>
#include "Common.hpp"

template
<
  class SelectT
  >
class DirectAccessibleGammaCode
{
private:
  uint64_t num_;
  uint64_t numWords_;
  uint64_t * array_; //!< Array to store values.
  SelectT sel_;


public:
  DirectAccessibleGammaCode
  () : num_(0),
       numWords_(0),
       array_(nullptr)
  {}


  ~DirectAccessibleGammaCode()
  {
    free(array_);
  }


  template<class vecT>
  void init
  (
   const vecT & vec
   ) {
    num_ = vec.size();
    sdsl::bit_vector bv;
    {
      uint64_t bvSize = 1; // +1 for sentinel
      for (uint64_t i = 0; i < vec.size(); ++i) {
        bvSize += sdsl::bits::hi(vec[i] + 1) + 1;
      }
      bv.resize(bvSize);
      numWords_ = (bvSize - vec.size() + 64) / 64;
      array_ = static_cast<uint64_t *>(malloc(numWords_ * sizeof(uint64_t)));
    }
    {
      uint64_t arrPos = 0;
      uint64_t bvPos = 0;
      for (uint64_t i = 0; i < vec.size(); ++i) {
        const uint64_t hi = sdsl::bits::hi(vec[i] + 1);
        bv[bvPos++] = 1;
        for (uint64_t j = 0; j < hi; ++j) {
          bv[bvPos++] = 0;
        }
        if (hi) {
          const uint64_t val = (vec[i] + 1) ^ (1ULL << hi);
          sdsl::bits::write_int(array_ + (arrPos / 64), val, arrPos % 64, hi);
          arrPos += hi;
        }
      }
      bv[bvPos] = 1;
    }
    sel_.init(std::move(bv));
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

    const uint64_t bvBeg = sel_(idx + 1);
    const uint64_t hi = sel_(idx + 2) - bvBeg - 1;
    // std::cout << "bvBeg = " << bvBeg << ", hi = " << hi << std::endl;
    uint64_t ret = 1ULL << hi;
    if (hi) {
      const uint64_t bitPos = bvBeg - idx;
      // std::cout << "bitPos = " << bitPos << std::endl;
      ret += sdsl::bits::read_int(array_ + (bitPos / 64), bitPos % 64, hi);
    }
    return ret - 1;
  }


  /*!
   * @brief Get size.
   */
  size_t size() const noexcept {
    return num_;
  }


  /*!
   * @brief Calculate total memory usage in bytes.
   */
  size_t calcMemBytes() const noexcept {
    size_t ret = sizeof(*this);
    ret += sel_.calcMemBytes();
    ret += sizeof(uint64_t) * numWords_;
    return ret;
  }


  void load
  (
   std::istream & in
   ) {
    in.read((char*) & num_, sizeof(num_));
    in.read((char*) & numWords_, sizeof(numWords_));
    array_ = static_cast<uint64_t *>(malloc(numWords_ * sizeof(uint64_t)));
    in.read((char*) array_, numWords_ * sizeof(uint64_t));
    sel_.load(in);
  }


  void serialize
  (
   std::ostream & out
   ) const {
    out.write((char*) & num_, sizeof(num_));
    out.write((char*) & numWords_, sizeof(numWords_));
    out.write((char*) array_, numWords_ * sizeof(uint64_t));
    sel_.serialize(out);
  }


  void printStatus
  (
   const bool verbose = false
   ) const noexcept {
    std::cout << "DirectAccessibleGammaCode object (" << this << ") " << __func__ << "(" << verbose << ") BEGIN" << std::endl;
    std::cout << "size = " << this->size() << ", numWords = " << numWords_ << std::endl;
    sel_.printStatus(verbose);
    if (verbose) {
      const auto size = this->size();
      std::cout << "dump stored values" << std::endl;
      for (uint64_t i = 0; i < size; ++i) {
        std::cout << this->read(i) << ", ";
      }
      std::cout << std::endl;
    }
    std::cout << "DirectAccessibleGammaCode object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
  }
};

#endif

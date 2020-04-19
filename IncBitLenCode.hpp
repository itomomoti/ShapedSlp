#ifndef INCLUDE_GUARD_IncBitLenCode
#define INCLUDE_GUARD_IncBitLenCode

#include <sdsl/bits.hpp>
#include "Common.hpp"

class IncBitLenCode
{
private:
  static constexpr uint64_t kOffset = 256;
  static constexpr uint64_t kOffsetW = 8;


  uint64_t num_;
  uint64_t numWords_;
  uint64_t * array_; //!< Array to store values.
  uint64_t * basePos_;


public:
  IncBitLenCode
  (
   size_t num = 0
   ) : num_(0),
       array_(nullptr),
       basePos_(nullptr)
  {
    if (num) {
      initNum(num);
    }
  }


  ~IncBitLenCode()
  {
    free(array_);
    free(basePos_);
  }


  void initNum
  (
   size_t num
   ) {
    assert(num_ == 0);

    num_ = num;
    const uint8_t h = sdsl::bits::hi(num + kOffset);
    basePos_ = static_cast<uint64_t *>(malloc((h - kOffsetW + 1) * sizeof(uint64_t)));
    uint64_t b = 0;
    for (uint8_t i = kOffsetW; i < h; ++i) {
      basePos_[i - kOffsetW] = b;
      b += (i + 1) * (1ULL << i);
    }
    basePos_[h - kOffsetW] = b;
    b += (h + 1) * ((num + kOffset) ^ (1ULL << h));
    numWords_ = (b + 64) / 64;
    array_ = static_cast<uint64_t *>(malloc(numWords_ * sizeof(uint64_t)));
  }


  template<class vecT>
  void init
  (
   const vecT & vec
   ) {
    assert(num_ == 0);

    initNum(vec.size());
    for (uint64_t i = 0; i < vec.size(); ++i) {
      this->write(vec[i], i);
    }
  }


  /*!
   * @brief Array subscript for writing
   */
  PackedArrayTypeValRef<IncBitLenCode> operator[]
  (
   const size_t idx
   ) {
    return PackedArrayTypeValRef<IncBitLenCode>(this, idx);
  }


  /*!
   * @brief Read only accessor.
   */
  uint64_t operator[]
  (
   const size_t idx //!< in [0, capacity_).
   ) const {
    return this->read(idx);
  }


  /*!
   * @brief Read value at 'idx'.
   */
  uint64_t read
  (
   const size_t idx //!< in [0, capacity_).
   ) const {
    const uint8_t h = sdsl::bits::hi(idx + kOffset);
    const uint64_t bitPos = basePos_[h - kOffsetW] + (h + 1) * ((idx + kOffset) ^ (1ULL << h));
    return sdsl::bits::read_int(array_ + (bitPos / 64), bitPos % 64, h + 1);
  }


  /*!
   * @brief Write 'val' at 'idx'.
   */
  void write
  (
   const uint64_t val,
   const size_t idx
   ) {
    const uint8_t h = sdsl::bits::hi(idx + kOffset);
    const uint64_t bitPos = basePos_[h - kOffsetW] + (h + 1) * ((idx + kOffset) ^ (1ULL << h));
    sdsl::bits::write_int(array_ + (bitPos / 64), val, bitPos % 64, h + 1);
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
    const uint8_t h = sdsl::bits::hi(num_ + kOffset);
    return sizeof(*this) + sizeof(uint64_t) * (numWords_ + (h - kOffsetW + 1));
  }


  void load
  (
   std::istream & in
   ) {
    in.read((char*) & num_, sizeof(num_));
    in.read((char*) & numWords_, sizeof(numWords_));
    array_ = static_cast<uint64_t *>(malloc(numWords_ * sizeof(uint64_t)));
    in.read((char*) array_, numWords_ * sizeof(uint64_t));
    const uint8_t h = sdsl::bits::hi(num_ + kOffset);
    basePos_ = static_cast<uint64_t *>(malloc((h - kOffsetW + 1) * sizeof(uint64_t)));
    in.read((char*) basePos_, (h - kOffsetW + 1) * sizeof(uint64_t));
  }


  void serialize
  (
   std::ostream & out
   ) const {
    out.write((char*) & num_, sizeof(num_));
    out.write((char*) & numWords_, sizeof(numWords_));
    out.write((char*) array_, numWords_ * sizeof(uint64_t));
    const uint8_t h = sdsl::bits::hi(num_ + kOffset);
    out.write((char*) basePos_, (h - kOffsetW + 1) * sizeof(uint64_t));
  }


  void printStatus
  (
   const bool verbose = false
   ) const noexcept {
    std::cout << "IncBitLenCode object (" << this << ") " << __func__ << "(" << verbose << ") BEGIN" << std::endl;
    std::cout << "size = " << this->size() << ", numWords = " << numWords_ << std::endl;
    if (verbose) {
      const uint8_t h = sdsl::bits::hi(num_ + kOffset);
      std::cout << "hi = " << h << std::endl;
      std::cout << "basePos_ dump: " << std::endl;
      for (uint64_t i = kOffsetW; i <= h; ++i) {
        std::cout << basePos_[i - kOffsetW] << ", ";
      }
      std::cout << std::endl;
      const auto size = this->size();
      std::cout << "dump stored values" << std::endl;
      for (uint64_t i = 0; i < size; ++i) {
        std::cout << this->read(i) << ", ";
      }
      std::cout << std::endl;
    }
    std::cout << "IncBitLenCode object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
  }
};

#endif

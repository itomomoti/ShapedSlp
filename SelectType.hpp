#ifndef INCLUDE_GUARD_SelectType
#define INCLUDE_GUARD_SelectType

#include <sdsl/bit_vectors.hpp>
#include "Common.hpp"

template<uint8_t t_b=1, uint8_t t_pat_len=1>
class SelectMcl
{
private:
  sdsl::bit_vector bv_;
  sdsl::select_support_mcl<> bvSel_;


public:
  SelectMcl
  ()
  {}


  ~SelectMcl
  ()
  {}


  void init
  (
   sdsl::bit_vector && bv
   ) {
    bv_ = std::move(sdsl::bit_vector(bv));
    bvSel_.init_slow(&bv_);
  }


  //// access to bit
  bool operator[]
  (
   uint64_t idx
   ) const {
    return bv_[idx];
  }


  uint64_t operator()
  (
   const size_t idx
   ) const {
    return bvSel_(idx);
  }


  /*!
   * @brief Get bit size.
   */
  size_t size() const noexcept {
    return bv_.size();
  }


  /*!
   * @brief Calculate total memory usage in bytes.
   */
  size_t calcMemBytes() const noexcept {
    size_t ret = sizeof(*this);
    ret += sdsl::size_in_bytes(bv_);
    ret += sdsl::size_in_bytes(bvSel_);
    return ret;
  }


  void load
  (
   std::istream & in
   ) {
    bv_.load(in);
    bvSel_.load(in);

    bvSel_.set_vector(&bv_);
  }


  void serialize
  (
   std::ostream & out
   ) const {
    bv_.serialize(out);
    bvSel_.serialize(out);
  }


  void printStatus
  (
   const bool verbose = false
   ) const noexcept {
    std::cout << "SelectMcl object (" << this << ") " << __func__ << "(" << verbose << ") BEGIN" << std::endl;
    std::cout << "bit size = " << this->size() << std::endl;
    if (verbose) {
      std::cout << "dump bits" << std::endl;
      printArray(bv_, bv_.size(), "");
    }
    std::cout << "SelectMcl object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
  }
};




template
<
  class t_hi_bit_vector = typename sdsl::bit_vector,
  class t_select_1     = typename t_hi_bit_vector::select_1_type,
  class t_select_0     = typename t_hi_bit_vector::select_0_type
  >
class SelectSdvec
{
private:
  using sdT = sdsl::sd_vector<t_hi_bit_vector, t_select_1, t_select_0>;


  sdT bv_;
  sdT::select_1_type bvSel_;


public:
  SelectSdvec
  ()
  {}


  ~SelectSdvec
  ()
  {}


  void init
  (
   sdsl::bit_vector && bv
   ) {
    bv_ = std::move(sdT(bv));
    bvSel_.set_vector(&bv_);
  }


  //// access to bit
  bool operator[]
  (
   uint64_t idx
   ) const {
    return bv_[idx];
  }


  uint64_t operator()
  (
   const size_t idx
   ) const {
    return bvSel_(idx);
  }


  /*!
   * @brief Get bit size.
   */
  size_t size() const noexcept {
    return bv_.size();
  }


  /*!
   * @brief Calculate total memory usage in bytes.
   */
  size_t calcMemBytes() const noexcept {
    size_t ret = sizeof(*this);
    ret += sdsl::size_in_bytes(bv_);
    ret += sdsl::size_in_bytes(bvSel_);
    return ret;
  }


  void load
  (
   std::istream & in
   ) {
    bv_.load(in);
    bvSel_.load(in);

    bvSel_.set_vector(&bv_);
  }


  void serialize
  (
   std::ostream & out
   ) const {
    bv_.serialize(out);
    bvSel_.serialize(out);
  }


  void printStatus
  (
   const bool verbose = false
   ) const noexcept {
    std::cout << "SelectSdvec object (" << this << ") " << __func__ << "(" << verbose << ") BEGIN" << std::endl;
    std::cout << "bit size = " << this->size() << std::endl;
    if (verbose) {
      std::cout << "dump bits" << std::endl;
      printArray(bv_, bv_.size(), "");
    }
    std::cout << "SelectSdvec object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
  }
};




template
<
  uint16_t t_bs=63,
  class t_rac=sdsl::int_vector<>,
  uint16_t t_k=32
  >
class SelectRrr
{
private:
  using rrrT = sdsl::rrr_vector<t_bs, t_rac, t_k>;
  using rrrSelT = rrrT::select_1_type;


  rrrT bv_;
  rrrSelT bvSel_;


public:
  SelectRrr
  ()
  {}


  ~SelectRrr
  ()
  {}


  void init
  (
   sdsl::bit_vector && bv
   ) {
    bv_ = std::move(rrrT(bv));
    bvSel_ = std::move(rrrSelT(&bv_));
  }


  //// access to bit
  bool operator[]
  (
   uint64_t idx
   ) const {
    return bv_[idx];
  }


  uint64_t operator()
  (
   const size_t idx
   ) const {
    return bvSel_(idx);
  }


  /*!
   * @brief Get bit size.
   */
  size_t size() const noexcept {
    return bv_.size();
  }


  /*!
   * @brief Calculate total memory usage in bytes.
   */
  size_t calcMemBytes() const noexcept {
    size_t ret = sizeof(*this);
    ret += sdsl::size_in_bytes(bv_);
    ret += sdsl::size_in_bytes(bvSel_);
    return ret;
  }


  void load
  (
   std::istream & in
   ) {
    bv_.load(in);
    bvSel_.load(in);

    bvSel_.set_vector(&bv_);
  }


  void serialize
  (
   std::ostream & out
   ) const {
    bv_.serialize(out);
    bvSel_.serialize(out);
  }


  void printStatus
  (
   const bool verbose = false
   ) const noexcept {
    std::cout << "SelectRrr object (" << this << ") " << __func__ << "(" << verbose << ") BEGIN" << std::endl;
    std::cout << "bit size = " << this->size() << std::endl;
    if (verbose) {
      std::cout << "dump bits" << std::endl;
      printArray(bv_, bv_.size(), "");
    }
    std::cout << "SelectRrr object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
  }
};

#endif

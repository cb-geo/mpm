#ifndef MPM_CG_MKL_H_
#define MPM_CG_MKL_H_

#include <cmath>

#include <Eigen/Sparse>
#include <mkl.h>

#include "cg_eigen.h"

namespace mpm {
template <typename Traits>
class CGMKL;
}

//! MPM Eigen CG class
//! \brief Conjugate Gradient solver class using Intel MKL
template <typename Traits>
class mpm::CGMKL : public CGEigen<Traits> {
 public:
  using CGEigen<Traits>::CGEigen;

  bool solve() override;

 protected:
  //! Logger
  using Solver<Traits>::console_;
};

//! Handle raw pointers in MKL
//! Only functionality is:
//!   1) to delete the memory when the object goes out of scope and;
//!   2) implicitly cast to the raw pointer so it can be passed to MKL directly.
//! \tparam T Type
template <typename T>
struct ptr_t : std::unique_ptr<T[]> {
  ptr_t(const size_t& n) : std::unique_ptr<T[]>(new T[n]) {}
  operator T*() { return this->get(); }
};

#include "cg_mkl.tcc"

#endif  // MPM_CG_MKL_H_
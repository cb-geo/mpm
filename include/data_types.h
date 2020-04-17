#ifndef MPM_DATA_TYPES_H_
#define MPM_DATA_TYPES_H_

#ifdef USE_MKL
#define EIGEN_USE_MKL_ALL
#endif

#include <Eigen/Dense>

namespace mpm {

//! Global index type for the node
using Index = unsigned long long;

//! Return zero
template <typename Ttype>
Ttype zero();

//! Zero for vector dimension 2
template <>
inline Eigen::Matrix<double, 2, 1> zero() {
  return Eigen::Matrix<double, 2, 1>::Zero();
}

//! Zero for vector dimension 3
template <>
inline Eigen::Matrix<double, 3, 1> zero() {
  return Eigen::Matrix<double, 3, 1>::Zero();
}

//! Zero for vector dimension 6
template <>
inline Eigen::Matrix<double, 6, 1> zero() {
  return Eigen::Matrix<double, 6, 1>::Zero();
}

//! Zero for a double
template <>
inline double zero() {
  return 0.;
}

}  // namespace mpm

#endif  // MPM_DATA_TYPES_H_

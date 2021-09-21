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

//! Zero
template <>
inline Eigen::Matrix<double, 2, 1> zero() {
  return Eigen::Matrix<double, 2, 1>::Zero();
}

//! Zero
template <>
inline Eigen::Matrix<double, 3, 1> zero() {
  return Eigen::Matrix<double, 3, 1>::Zero();
}

//! Zero
template <>
inline double zero() {
  return 0.;
}

//! Position type
//! None: No position is specified
//! Side: Nodes along model vertical or horizontal sides
//! Corner: Nodes at model corner
enum class Position { None, Side, Corner };

}  // namespace mpm

#endif  // MPM_DATA_TYPES_H_

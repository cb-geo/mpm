#ifndef MPM_DATA_TYPES_H_
#define MPM_DATA_TYPES_H_

#ifdef USE_MKL
#define EIGEN_USE_MKL_ALL
#include <Eigen/Dense>
#endif

namespace mpm {

//! Global index type for the node
using Index = unsigned long long;

}  // namespace mpm

#endif  // MPM_DATA_TYPES_H_

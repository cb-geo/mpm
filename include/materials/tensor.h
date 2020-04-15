#ifndef MPM_TENSOR_H_
#define MPM_TENSOR_H_

#include <cmath>

#include "Eigen/Dense"

namespace mpm {
namespace tensor {

//! Compute mean stress p
//! \param[in] stress Stress in Voigt notation where positive is tension
inline double compute_mean_p(const Eigen::Matrix<double, 6, 1>& stress);

//! Compute deviatoric stress
//! \param[in] stress Stress in Voigt notation where positive is tension
inline Eigen::Matrix<double, 6, 1> compute_deviatoric_stress(
    const Eigen::Matrix<double, 6, 1>& stress);

//! Compute J2 invariant
//! \param[in] stress Stress in Voigt notation where positive is tension
inline double compute_j2(const Eigen::Matrix<double, 6, 1>& stress);

//! Compute J3 invariant
//! \param[in] stress Stress in Voigt notation where positive is tension
inline double compute_j3(const Eigen::Matrix<double, 6, 1>& stress);

//! Compute deviatoric q
//! \param[in] stress Stress in Voigt notation where positive is tension
inline double compute_deviatoric_q(const Eigen::Matrix<double, 6, 1>& stress);

//! Compute Lode angle
//! \param[in] stress Stress in Voigt notation where positive is tension
inline double compute_lode_angle(const Eigen::Matrix<double, 6, 1>& stress);

}  // namespace tensor
}  // namespace mpm

#include "tensor.tcc"

#endif  // MPM_TENSOR_H_

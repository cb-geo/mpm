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

//! Compute Lode angle theta (cosine convention)
//! \param[in] stress Stress in Voigt notation where positive is tension
inline double compute_lode_angle(const Eigen::Matrix<double, 6, 1>& stress);

//! Compute derivative of p in terms of stress sigma
//! \param[in] stress Stress in Voigt notation where positive is tension
inline Eigen::Matrix<double, 6, 1> compute_dp_dsigma(
    const Eigen::Matrix<double, 6, 1>& stress);

//! Compute derivative of q in terms of stress sigma
//! \param[in] stress Stress in Voigt notation where positive is tension
inline Eigen::Matrix<double, 6, 1> compute_dq_dsigma(
    const Eigen::Matrix<double, 6, 1>& stress);

//! Compute derivative of J2 in terms of stress sigma
//! \param[in] stress Stress in Voigt notation where positive is tension
inline Eigen::Matrix<double, 6, 1> compute_dj2_dsigma(
    const Eigen::Matrix<double, 6, 1>& stress);

//! Compute derivative of J3 in terms of stress sigma
//! \param[in] stress Stress in Voigt notation where positive is tension
inline Eigen::Matrix<double, 6, 1> compute_dj3_dsigma(
    const Eigen::Matrix<double, 6, 1>& stress);

//! Compute derivative of Lode angle theta in terms of stress sigma
//! \param[in] stress Stress in Voigt notation where positive is tension
inline Eigen::Matrix<double, 6, 1> compute_dtheta_dsigma(
    const Eigen::Matrix<double, 6, 1>& stress);

}  // namespace tensor
}  // namespace mpm

#include "tensor.tcc"

#endif  // MPM_TENSOR_H_

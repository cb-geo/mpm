#ifndef MPM_MATERIAL_UTILITY_H_
#define MPM_MATERIAL_UTILITY_H_

#include <cmath>

#include "Eigen/Dense"

namespace mpm {
namespace material_utility {

//! Compute mean stress p
//! \param[in] stress Stress in Voigt notation where positive is tension
inline const double compute_mean_p(const Eigen::Matrix<double, 6, 1>& stress);

//! Compute deviatoric stress
//! \param[in] stress Stress in Voigt notation where positive is tension
inline const Eigen::Matrix<double, 6, 1> compute_deviatoric_stress(
    const Eigen::Matrix<double, 6, 1>& stress);

//! Compute J2 invariant
//! \param[in] stress Stress in Voigt notation where positive is tension
inline const double compute_j2(const Eigen::Matrix<double, 6, 1>& stress);

//! Compute J3 invariant
//! \param[in] stress Stress in Voigt notation where positive is tension
inline const double compute_j3(const Eigen::Matrix<double, 6, 1>& stress);

//! Compute deviatoric q
//! \param[in] stress Stress in Voigt notation where positive is tension
inline const double compute_deviatoric_q(
    const Eigen::Matrix<double, 6, 1>& stress);

//! Compute Lode angle theta (cosine convention)
//! \param[in] stress Stress in Voigt notation where positive is tension
inline const double compute_lode_angle(
    const Eigen::Matrix<double, 6, 1>& stress);

//! Compute derivative of p in terms of stress sigma
//! \param[in] stress Stress in Voigt notation where positive is tension
inline const Eigen::Matrix<double, 6, 1> compute_dp_dsigma(
    const Eigen::Matrix<double, 6, 1>& stress);

//! Compute derivative of q in terms of stress sigma
//! \param[in] stress Stress in Voigt notation where positive is tension
inline const Eigen::Matrix<double, 6, 1> compute_dq_dsigma(
    const Eigen::Matrix<double, 6, 1>& stress);

//! Compute derivative of J2 in terms of stress sigma
//! \param[in] stress Stress in Voigt notation where positive is tension
inline const Eigen::Matrix<double, 6, 1> compute_dj2_dsigma(
    const Eigen::Matrix<double, 6, 1>& stress);

//! Compute derivative of J3 in terms of stress sigma
//! \param[in] stress Stress in Voigt notation where positive is tension
inline const Eigen::Matrix<double, 6, 1> compute_dj3_dsigma(
    const Eigen::Matrix<double, 6, 1>& stress);

//! Compute derivative of Lode angle theta in terms of stress sigma
//! \param[in] stress Stress in Voigt notation where positive is tension
inline const Eigen::Matrix<double, 6, 1> compute_dtheta_dsigma(
    const Eigen::Matrix<double, 6, 1>& stress);

}  // namespace material_utility
}  // namespace mpm

#include "material_utility.tcc"

#endif  // MPM_MATERIAL_UTILITY_H_

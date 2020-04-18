#ifndef MPM_MATERIAL_UTILITY_H_
#define MPM_MATERIAL_UTILITY_H_

#include <cmath>

#include "data_types.h"

namespace mpm {
namespace materials {

//! Compute mean stress p (tension positive)
//! \param[in] stress Stress in Voigt notation where positive is tension
inline const double p(const Eigen::Matrix<double, 6, 1>& stress);

//! Compute deviatoric stress
//! \param[in] stress Stress in Voigt notation where positive is tension
inline const Eigen::Matrix<double, 6, 1> deviatoric_stress(
    const Eigen::Matrix<double, 6, 1>& stress);

//! Compute J2 invariant
//! \param[in] stress Stress in Voigt notation where positive is tension
inline const double j2(const Eigen::Matrix<double, 6, 1>& stress);

//! Compute J3 invariant
//! \param[in] stress Stress in Voigt notation where positive is tension
inline const double j3(const Eigen::Matrix<double, 6, 1>& stress);

//! Compute deviatoric q
//! \param[in] stress Stress in Voigt notation where positive is tension
inline const double q(const Eigen::Matrix<double, 6, 1>& stress);

//! Compute Lode angle theta (cosine convention)
//! \param[in] stress Stress in Voigt notation where positive is tension
inline const double lode_angle(const Eigen::Matrix<double, 6, 1>& stress);

//! Compute Lode angle theta (cosine convention)
//! \param[in] stress Stress in Voigt notation where positive is tension
//! \param[in] tolerance Default tolerance value specified by user
inline const double lode_angle(const Eigen::Matrix<double, 6, 1>& stress,
                               const double tolerance);

//! Compute derivative of p in terms of stress sigma
//! \param[in] stress Stress in Voigt notation where positive is tension
inline const Eigen::Matrix<double, 6, 1> dp_dsigma(
    const Eigen::Matrix<double, 6, 1>& stress);

//! Compute derivative of q in terms of stress sigma
//! \param[in] stress Stress in Voigt notation where positive is tension
inline const Eigen::Matrix<double, 6, 1> dq_dsigma(
    const Eigen::Matrix<double, 6, 1>& stress);

//! Compute derivative of J2 in terms of stress sigma
//! \param[in] stress Stress in Voigt notation where positive is tension
inline const Eigen::Matrix<double, 6, 1> dj2_dsigma(
    const Eigen::Matrix<double, 6, 1>& stress);

//! Compute derivative of J3 in terms of stress sigma
//! \param[in] stress Stress in Voigt notation where positive is tension
inline const Eigen::Matrix<double, 6, 1> dj3_dsigma(
    const Eigen::Matrix<double, 6, 1>& stress);

//! Compute derivative of Lode angle theta in terms of stress sigma
//! \param[in] stress Stress in Voigt notation where positive is tension
inline const Eigen::Matrix<double, 6, 1> dtheta_dsigma(
    const Eigen::Matrix<double, 6, 1>& stress);

//! Compute derivative of Lode angle theta in terms of stress sigma
//! \param[in] stress Stress in Voigt notation where positive is tension
//! \param[in] tolerance Default tolerance value specified by user
inline const Eigen::Matrix<double, 6, 1> dtheta_dsigma(
    const Eigen::Matrix<double, 6, 1>& stress, const double tolerance);

}  // namespace materials
}  // namespace mpm

#include "material_utility.tcc"

#endif  // MPM_MATERIAL_UTILITY_H_

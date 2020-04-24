#ifndef MPM_MATERIAL_UTILITY_H_
#define MPM_MATERIAL_UTILITY_H_

#include <cmath>

#include "data_types.h"

namespace mpm {
namespace materials {

//! Compute mean stress p (tension positive)
//! \param[in] stress Stress in Voigt notation where positive is tension
//! \retval mean stress p (tension positive)
inline double p(const Eigen::Matrix<double, 6, 1>& stress);

//! Compute deviatoric stress
//! \param[in] stress Stress in Voigt notation where positive is tension
//! \retval deviatoric stress
inline const Eigen::Matrix<double, 6, 1> deviatoric_stress(
    const Eigen::Matrix<double, 6, 1>& stress);

//! Compute J2 invariant
//! \param[in] stress Stress in Voigt notation where positive is tension
//! \retval J2 invariant
inline double j2(const Eigen::Matrix<double, 6, 1>& stress);

//! Compute J3 invariant
//! \param[in] stress Stress in Voigt notation where positive is tension
//! \retval J3 invariant
inline double j3(const Eigen::Matrix<double, 6, 1>& stress);

//! Compute deviatoric q
//! \param[in] stress Stress in Voigt notation where positive is tension
//! \retval deviatoric q
inline double q(const Eigen::Matrix<double, 6, 1>& stress);

//! Compute Lode angle theta (cosine convention)
//! \param[in] stress Stress in Voigt notation where positive is tension
//! \param[in] tolerance Default tolerance value specified by user
//! \retval Lode angle theta (cosine convention)
inline double lode_angle(
    const Eigen::Matrix<double, 6, 1>& stress,
    double tolerance = std::numeric_limits<double>::epsilon());

//! Compute derivative of p in terms of stress sigma
//! \param[in] stress Stress in Voigt notation where positive is tension
//! \retval derivative of p in terms of stress sigma
inline const Eigen::Matrix<double, 6, 1> dp_dsigma(
    const Eigen::Matrix<double, 6, 1>& stress);

//! Compute derivative of q in terms of stress sigma
//! \param[in] stress Stress in Voigt notation where positive is tension
//! \retval derivative of q in terms of stress sigma
inline const Eigen::Matrix<double, 6, 1> dq_dsigma(
    const Eigen::Matrix<double, 6, 1>& stress);

//! Compute derivative of J2 in terms of stress sigma
//! \param[in] stress Stress in Voigt notation where positive is tension
//! \retval derivative of J2 in terms of stress sigma
inline const Eigen::Matrix<double, 6, 1> dj2_dsigma(
    const Eigen::Matrix<double, 6, 1>& stress);

//! Compute derivative of J3 in terms of stress sigma
//! \param[in] stress Stress in Voigt notation where positive is tension
//! \retval derivative of J3 in terms of stress sigma
inline const Eigen::Matrix<double, 6, 1> dj3_dsigma(
    const Eigen::Matrix<double, 6, 1>& stress);

//! Compute derivative of Lode angle theta in terms of stress sigma
//! \param[in] stress Stress in Voigt notation where positive is tension
//! \param[in] tolerance Default tolerance value specified by user
//! \retval derivative of Lode angle theta in terms of stress sigma
inline const Eigen::Matrix<double, 6, 1> dtheta_dsigma(
    const Eigen::Matrix<double, 6, 1>& stress,
    double tolerance = std::numeric_limits<double>::epsilon());

//! Compute incremental dplastic_strain
//! \param[in] dstress Incremental stress in Voigt notation where positive is
//! tension \param[in] dstrain Incremental total strain in Voigt notation where
//! positive is tension \param[in] de Elastic stiffness matrix \param[in]
//! tolerance Default tolerance value specified by user \retval dplastic_strain
inline const Eigen::Matrix<double, 6, 1> dplastic_strain(
    const Eigen::Matrix<double, 6, 1>& dstress,
    const Eigen::Matrix<double, 6, 1>& dstrain,
    const Eigen::Matrix<double, 6, 6>& de);

//! Compute plastic deviatoric strain
//! \param[in] plastic_strain Plastic strain in Voigt notation where positive is
//! tension
//! \retval plastic deviatoric strain
inline double pdstrain(const Eigen::Matrix<double, 6, 1>& plastic_strain);

}  // namespace materials
}  // namespace mpm

#include "material_utility.tcc"

#endif  // MPM_MATERIAL_UTILITY_H_

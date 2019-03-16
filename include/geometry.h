#ifndef MPM_GEOMETRY_H_
#define MPM_GEOMETRY_H_

#include "Eigen/Dense"

namespace mpm {
namespace geometry {
//! Compute the Euler rotation matrix in orthogonal axis coordinate system
//! \param[in] angles Rotation angles wrt to the coordinate system
//! \tparam Tdim Dimension
template <int Tdim>
Eigen::Matrix<double, Tdim, Tdim> rotation_matrix(
    const Eigen::Matrix<double, Tdim, 1>& angles);

//! Compute the angle between two vectors in radians
//! \param[in] vector_a First vector
//! \param[in] vector_b Second vector
//! \tparam Tdim Dimension
template <int Tdim>
inline double angle_between_vectors(
    const Eigen::Matrix<double, Tdim, 1>& vector_a,
    const Eigen::Matrix<double, Tdim, 1>& vector_b);

//! Compute Euler angles with respect to the Cartesian coordinates
//! \param[in] new_axes New orthogonal coordinate systems
//! \retval euler_angles Euler Angles
//! \tparam Tdim Dimension
template <int Tdim>
inline Eigen::Matrix<double, Tdim, 1> euler_angles_cartesian(
    const Eigen::Matrix<double, Tdim, Tdim>& new_axes);
}  // namespace geometry
}  // namespace mpm

#include "geometry.tcc"

#endif  // MPM_GEOMETRY_H_

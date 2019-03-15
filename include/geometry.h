#ifndef MPM_GEOMETRY_H_
#define MPM_GEOMETRY_H_

#include "Eigen/Dense"

namespace mpm {
//! geometric manipulations / operations
namespace geometry {
//! Compute the Euler rotation matrix for an orthogonal axis
//! coordinate system
//! \param[in] angles Rotation angles depending on the dimension
//! \tparam Tdim Dimension
template <unsigned Tdim>
Eigen::Matrix<double, Tdim, Tdim> rotation_matrix(
    const Eigen::Matrix<double, Tdim, 1>& angles);

//! Compute the angle between two vectors in radians
//! \param[in] vector_a First vector
//! \param[in] vector_b Second vector
//! \tparam Tdim Dimension
template <unsigned Tdim>
double angle_between_vectors(const Eigen::Matrix<double, Tdim, 1>& vector_a,
                             const Eigen::Matrix<double, Tdim, 1>& vector_b);

//! Compute euler angles with respect to the Cartesian coordinates
//! \param[in] new_axes New orthogonal coordinate systems (2 vectors for 2D, 3
//! vectors for 3D)
//! \retval euler_angles Euler Angles (2 angles for 2D, 3
//! angles for 3D)
//! \tparam Tdim Dimension
template <unsigned Tdim>
Eigen::Matrix<double, Tdim, 1> euler_angles_cartesian(
    const Eigen::Matrix<double, Tdim, Tdim>& new_axes);
}  // namespace geometry
}  // namespace mpm

#endif  // MPM_GEOMETRY_H_

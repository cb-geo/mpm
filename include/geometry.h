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
template <int Tdim>
Eigen::Matrix<double, Tdim, Tdim> rotation_matrix(
    const Eigen::Matrix<double, Tdim, 1>& angles);

//! Compute the angle between two vectors in radians
//! \param[in] vector_a First vector
//! \param[in] vector_b Second vector
//! \tparam Tdim Dimension
template <int Tdim>
double angle_between_vectors(const Eigen::Matrix<double, Tdim, 1>& vector_a,
                             const Eigen::Matrix<double, Tdim, 1>& vector_b) {
  // angle between two vectors a and b = arccos( a dot b / ||a|| ||b||)
  return acos((vector_a.normalized()).dot((vector_b.normalized())));
}

//! Compute euler angles with respect to the Cartesian coordinates
//! \param[in] new_axes New orthogonal coordinate systems (2 vectors for 2D, 3
//! vectors for 3D)
//! \retval euler_angles Euler Angles (2 angles for 2D, 3
//! angles for 3D)
//! \tparam Tdim Dimension
template <int Tdim>
Eigen::Matrix<double, Tdim, 1> euler_angles_cartesian(
    const Eigen::Matrix<double, Tdim, Tdim>& new_axes) {

  // Make cartesian coordinate system
  Eigen::Matrix<double, Tdim, Tdim> original_axes;
  original_axes.setIdentity();

  // Compute line of nodes vector that bisects original x and y unit vector
  Eigen::Matrix<double, Tdim, 1> line_of_nodes =
      original_axes.col(0) + original_axes.col(1);
  line_of_nodes = line_of_nodes.normalized();

  // Make a vector of euler angles
  Eigen::Matrix<double, Tdim, 1> euler_angles;

  // Compute alpha
  euler_angles(0) = acos(
      (original_axes.col(0).normalized()).dot((line_of_nodes.normalized())));

  // Compute beta
  euler_angles(1) =
      acos((line_of_nodes.normalized()).dot((new_axes.col(0).normalized())));

  // Compute gamma
  if (Tdim == 3)
    euler_angles(2) = acos((new_axes.col(2).normalized())
                               .dot((original_axes.col(2).normalized())));

  return euler_angles;
}
}  // namespace geometry
}  // namespace mpm

#endif  // MPM_GEOMETRY_H_

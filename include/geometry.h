#ifndef MPM_GEOMETRY_H_
#define MPM_GEOMETRY_H_

#include "Eigen/Dense"

namespace mpm {

//! Geometry class
//! \brief Base class that performs geometric manipulations / operations
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Geometry {
 public:
  //! Constructor
  Geometry() = default;

  //! Compute the inverse Euler rotation matrix for an orthogonal axis
  //! coordinate system \param[in] angles Rotation angles depending on the
  //! dimension
  Eigen::Matrix<double, Tdim, Tdim> inverse_rotation_matrix(
      const Eigen::Matrix<double, Tdim, 1>& angles) const;

  //! Compute the angle between two vectors in radians
  //! \param[in] vector_a First vector
  //! \param[in] vector_b Second vector
  const double angle_between_vectors(
      const Eigen::Matrix<double, Tdim, 1>& vector_a,
      const Eigen::Matrix<double, Tdim, 1>& vector_b);
};
}  // namespace mpm

#include "geometry.tcc"

#endif  // MPM_GEOMETRY_H_

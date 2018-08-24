#ifndef MPM_GEOMETRY_H_
#define MPM_GEOMETRY_H_

#include <limits>
#include <memory>

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

  //! Return inverse rotation matrix for orthogonal axis coordinate system
  //! \param[in] angles Rotation angles depending on dimension
  //! \retval inverse of Euler rotation matrix R
  Eigen::Matrix<double, Tdim, Tdim> inverse_rotation_matrix(
      const Eigen::Matrix<double, Tdim, 1>& angles) const;

  //! Return angle between two vectors in radians
  //! \param[in] vector_a First vector
  //! \param[in] vecotr_b Second vector
  //! \retval anti-clockwise angle from vector_a to vector_b
  const double angle_between_two_vectors(
      const Eigen::Matrix<double, Tdim, 1>& vector_a,
      const Eigen::Matrix<double, Tdim, 1>& vector_b);

};  // Geometry class
}  // namespace mpm

#include "geometry.tcc"

#endif  // MPM_GEOMETRY_H_

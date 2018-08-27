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

  //! Return the inverse Euler rotation matrix for an orthogonal axis coordinate
  //! system \param[in] angles Rotation angles depending on the dimension
  Eigen::Matrix<double, Tdim, Tdim> inverse_rotation_matrix(
      const Eigen::Matrix<double, Tdim, 1>& angles) const;

  //! Return the angle between two vectors in radians
  //! \param[in] vector_a First vector
  //! \param[in] vecotr_b Second vector
  const double angle_between_vectors(
      const Eigen::Matrix<double, Tdim, 1>& vector_a,
      const Eigen::Matrix<double, Tdim, 1>& vector_b);

};  // Geometry class
}  // namespace mpm

#include "geometry.tcc"

#endif  // MPM_GEOMETRY_H_

#ifndef MPM_GEOMETRY_H_
#define MPM_GEOMETRY_H_

#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "Eigen/LU"

#include "logger.h"

namespace mpm {

//! Geometry class
//! \brief Base class that computes geometry mathematics
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Geometry {
 public:
  //! Constructor
  Geometry(){};

  //! Compute inverse of rotation matrix for orthogonal axis coordinate system
  //! \param[in] agnles Rotation angles depending on dimension
  //! \retval inverse of Euler rotation matrix R
  Eigen::MatrixXd compute_inverse_rotation_matrix(
      const Eigen::VectorXd& angles);

 protected:
  //! Logger
  std::unique_ptr<spdlog::logger> console_;

};  // Geometry class
}  // namespace mpm

#include "geometry.tcc"

#endif  // MPM_GEOMETRY_H_

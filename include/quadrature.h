#ifndef MPM_QUADRATURE_BASE_H_
#define MPM_QUADRATURE_BASE_H_

#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>

#include <Eigen/Dense>

//! MPM namespace
namespace mpm {

// Quadrature base class
//! \brief Base class for quadrature
//! \tparam Tdim Dimension
//! \tparam Tnquadratures number of quadratures
template <unsigned Tdim, unsigned Tnquadratures>
class QuadratureBase {
 public:
  //! Constructor
  //! Assign variables to zero
  QuadratureBase() {}

  //! Destructor
  virtual ~QuadratureBase() {}

  //! Return quadrature points
  //! \param[out] qpoints Quadrature points in local coordinates
  virtual Eigen::MatrixXd quadratures() = 0;

  //! Return weights
  //! \param[out] weights Weights for quadrature points
  virtual Eigen::VectorXd weights() = 0;

};

}  // namespace mpm

#endif  // MPM_QUADRATURE_BASE_H_
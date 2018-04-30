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
  QuadratureBase() {
    qpoints_.fill(std::numeric_limits<double>::quiet_NaN());
    weights_.resize(Tnquadratures);
    std::fill(weights_.begin(), weights_.end(),
              std::numeric_limits<double>::quiet_NaN());
  }

  //! Destructor
  virtual ~QuadratureBase() {}

  //! Return quadrature points
  //! \param[out] qpoints Quadrature points in local coordinates
  virtual Eigen::Matrix<double, Tnquadratures, Tdim> quadratures() { return qpoints_; }

  //! Return weights
  //! \param[out] weights Weights for quadrature points
  virtual std::vector<double> weights() { return weights_; }

 protected:
  Eigen::Matrix<double, Tnquadratures, Tdim> qpoints_;
  std::vector<double> weights_;
};

}  // namespace mpm

#endif  // MPM_QUADRATURE_BASE_H_
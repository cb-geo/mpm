#ifndef MPM_QUADRATURE_BASE_H_
#define MPM_QUADRATURE_BASE_H_

#include <algorithm>
#include <limits>
#include <vector>

#include <Eigen/Dense>

//! MPM namespace
namespace mpm {

// Quadrature base class
//! \brief Base class for quadrature
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Quadrature {
 public:
  //! Default constructor
  Quadrature() = default;

  //! Destructor
  virtual ~Quadrature() {}

  //! Return quadrature points
  //! \param[out] qpoints Quadrature points in local coordinates
  virtual Eigen::MatrixXd quadratures() const = 0;

  //! Return weights
  //! \param[out] weights Weights for quadrature points
  virtual Eigen::VectorXd weights() const = 0;
};

}  // namespace mpm

#endif  // MPM_QUADRATURE_BASE_H_

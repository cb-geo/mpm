#ifndef MPM_QUADRILATERAL_QUADRATURE_H_
#define MPM_QUADRILATERAL_QUADRATURE_H_

#include <exception>
#include <iostream>

#include <Eigen/Dense>

#include "quadrature.h"

//! MPM namespace
namespace mpm {

// Quadrilateral quadrature class derived from quadrature base class
//! \brief Quadrature (gauss points) for a quadrilateral  element
//! \tparam Tdim Dimension
//! \tparam Tnquadratures number of quadratures
template <unsigned Tdim, unsigned Tnquadratures>
class QuadrilateralQuadrature : public QuadratureBase<Tdim, Tnquadratures> {

 public:
  QuadrilateralQuadrature() : QuadratureBase<Tdim, Tnquadratures>() {
    static_assert(Tdim == 2, "Invalid dimension for a quadrilateral element");
    static_assert(
        ((Tnquadratures == 1) || (Tnquadratures == 4) || (Tnquadratures == 9)),
        "Invalid number of quadratures");
  }

  //! Return quadrature points
  //! \param[out] qpoints Quadrature points in local coordinates
  Eigen::MatrixXd quadratures();

  //! Return weights
  //! \param[out] weights Weights for quadrature points
  Eigen::VectorXd weights();
};

}  // namespace mpm

#include "quad_quadrature.tcc"

#endif  // MPM_QUADRILATERAL_QUADRATURE_H_

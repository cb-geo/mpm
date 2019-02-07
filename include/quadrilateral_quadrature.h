#ifndef MPM_QUADRILATERAL_QUADRATURE_H_
#define MPM_QUADRILATERAL_QUADRATURE_H_

#include <exception>

#include <Eigen/Dense>

#include "quadrature.h"

//! MPM namespace
namespace mpm {

// Quadrilateral quadrature class derived from Quadrature base class
//! \brief Quadrature (gauss points) for a quadrilateral  element
//! \tparam Tdim Dimension
//! \tparam Tnquadratures number of quadratures
template <unsigned Tdim, unsigned Tnquadratures>
class QuadrilateralQuadrature : public Quadrature<Tdim> {

 public:
  QuadrilateralQuadrature() : Quadrature<Tdim>() {
    static_assert(Tdim == 2, "Invalid dimension for a quadrilateral element");
    static_assert(
        ((Tnquadratures == 1) || (Tnquadratures == 4) || (Tnquadratures == 9)),
        "Invalid number of quadratures");
  }

  //! Return quadrature points
  //! \param[out] qpoints Quadrature points in local coordinates
  Eigen::MatrixXd quadratures() const override;

  //! Return weights
  //! \param[out] weights Weights for quadrature points
  Eigen::VectorXd weights() const override;
};

}  // namespace mpm

#include "quadrilateral_quadrature.tcc"

#endif  // MPM_QUADRILATERAL_QUADRATURE_H_

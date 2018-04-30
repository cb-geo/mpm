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
  QuadrilateralQuadrature();

 private:
  using QuadratureBase<Tdim, Tnquadratures>::qpoints_;
  using QuadratureBase<Tdim, Tnquadratures>::weights_;
};

}  // namespace mpm

#include "quad_quadrature.tcc"

#endif  // MPM_QUADRILATERAL_QUADRATURE_H_

#ifndef MPM_HEXAHEDRON_QUADRATURE_H_
#define MPM_HEXAHEDRON_QUADRATURE_H_

#include <exception>
#include <iostream>

#include <Eigen/Dense>

#include "quad_base.h"

//! MPM namespace
namespace mpm {

//! Hexahedron quadrature class derived from QuadratureBase class
//! \brief Quadrature (gauss points) for a hexahedron element
//! \tparam Tdim Dimension
//! \tparam Tnquadratures number of quadratures
template <unsigned Tdim, unsigned Tnquadratures>
class HexahedronQuadrature : public QuadratureBase<Tdim, Tnquadratures> {

 public:
  HexahedronQuadrature();

 private:
  using QuadratureBase<Tdim, Tnquadratures>::qpoints_;
  using QuadratureBase<Tdim, Tnquadratures>::weights_;
};

}  // namespace mpm

#include "hex_quad.tcc"

#endif  // MPM_HEXAHEDRON_QUADRATURE_H_

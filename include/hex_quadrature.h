#ifndef MPM_HEXAHEDRON_QUADRATURE_H_
#define MPM_HEXAHEDRON_QUADRATURE_H_

#include <exception>

#include <Eigen/Dense>

#include "quadrature.h"

//! MPM namespace
namespace mpm {

//! Hexahedron quadrature class derived from QuadratureBase class
//! \brief Quadrature (gauss points) for a hexahedron element
//! \tparam Tdim Dimension
//! \tparam Tnquadratures number of quadratures
template <unsigned Tdim, unsigned Tnquadratures>
class HexahedronQuadrature : public QuadratureBase<Tdim, Tnquadratures> {

 public:
  HexahedronQuadrature() : mpm::QuadratureBase<Tdim, Tnquadratures>() {
    static_assert(Tdim == 3, "Invalid dimension for a 3D hexahedron element");
    static_assert(
        (Tnquadratures == 1) || (Tnquadratures == 8) || (Tnquadratures == 27),
        "Invalid number of quadratures");
  }

  //! Return quadrature points
  //! \param[out] qpoints Quadrature points in local coordinates
  Eigen::MatrixXd quadratures() override;

  //! Return weights
  //! \param[out] weights Weights for quadrature points
  Eigen::VectorXd weights() override;
};

}  // namespace mpm

#include "hex_quadrature.tcc"

#endif  // MPM_HEXAHEDRON_QUADRATURE_H_

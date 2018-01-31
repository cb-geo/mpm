#ifndef MPM_HEXAHEDRONSHAPEFN_H_
#define MPM_HEXAHEDRONSHAPEFN_H_

#include <exception>
#include <iostream>

#include <Eigen/Dense>

#include "shapefn.h"

//! MPM namespace
namespace mpm {

//! Hexahedron shape function class derived from ShapeFn class
//! \brief Shape functions of a hexahedron element
//! \tparam Tdim Dimension
//! \tparam Tnfunctions Number of functions
template <unsigned Tdim, unsigned Tnfunctions>
class HexahedronShapeFn : public ShapeFn<Tdim> {

 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! constructor with number of shape functions
  HexahedronShapeFn() : mpm::ShapeFn<Tdim>() {
    static_assert(Tdim == 3, "Invalid dimension for a hexahedron element");
    static_assert((Tnfunctions == 8 || Tnfunctions == 20), 
            "Specified number of shape functions is not defined");    
  }

  //! Return number of functions
  unsigned nfunctions() const { return Tnfunctions; }

  //! Evaluate shape functions at given local coordinates
  //! \param[in] xi given local coordinates
  Eigen::VectorXd shapefn(const VectorDim& xi);

  //! Evaluate gradient of shape functions
  //! \param[in] xi given local coordinates
  Eigen::MatrixXd grad_shapefn(const VectorDim& xi);
};

}  // namespace mpm
#include "hex_shapefn.tcc"

#endif  // MPM_HEXAHEDRONSHAPEFN_H_

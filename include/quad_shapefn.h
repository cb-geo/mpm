#ifndef MPM_QUADRILATERALSHAPEFN_H_
#define MPM_QUADRILATERALSHAPEFN_H_

#include <exception>
#include <iostream>

#include <Eigen/Dense>

#include "shapefn.h"

namespace mpm {

//! Quadrilateral shape function class derived from ShapeFn class
//! \brief Shape functions of a quadrilateral element
//! \tparam Tdim Dimension
//! \tparam Tnfunctions Number of functions
template <unsigned Tdim, unsigned Tnfunctions>
class QuadrilateralShapeFn : public ShapeFn<Tdim> {

 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! constructor with number of shape functions
  QuadrilateralShapeFn() : mpm::ShapeFn<Tdim>() {
    static_assert(Tdim == 2, "Invalid dimension for a quadrilateral element");
    static_assert((Tnfunctions == 4 || Tnfunctions == 8 || Tnfunctions == 9),
            "Specified number of shape functions is not defined");
    shapefn_.resize(Tnfunctions, 1);
    grad_shapefn_.resize(Tnfunctions, Tdim);
  }

  //! Return number of functions
  unsigned nfunctions() const { return Tnfunctions; }

  //! Evaluate shape functions at given local coordinates
  //! \param[in] xi given local coordinates
  Eigen::VectorXd shapefn(const VectorDim& xi);

  //! Evaluate gradient of shape functions
  //! \param[in] xi given local coordinates
  Eigen::MatrixXd grad_shapefn(const VectorDim& xi);

 protected:
  // Shape functions
  using ShapeFn<Tdim>::shapefn_;
  // Gradient shape functions
  using ShapeFn<Tdim>::grad_shapefn_;
};

}  // namespace mpm
#include "quad_shapefn.tcc"

#endif  // MPM_QUADRILATERALSHAPEFN_H_

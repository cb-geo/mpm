#ifndef MPM_QUADRILATERALSHAPEFN_H_
#define MPM_QUADRILATERALSHAPEFN_H_

#include <exception>
#include <iostream>

#include <Eigen/Dense>

#include "shapefn.h"

namespace mpm {

//! Quadrilateral shape function class derived from ShapeFn class
//! \brief Shape functions of a quadrilateral element
//! \details 4-noded, 8-noded, and 9-noded quadrilateral element \n
//! Shape function, gradient shape function, B-matrix, indices \n
//! 4-node Quadrilateral Element \n
//! <pre>
//!
//! 3 0----------0 2
//!   |          |
//!   |          |
//!   |          |
//!   |          |
//! 0 0----------0 1
//!
//! </pre>
//! 8-node Quadrilateral Element
//! <pre>
//!
//!  3      6       2
//!   0-----0-----0
//!   |           |
//!   |           |
//! 7 0           0 5
//!   |           |
//!   |           |
//!   0-----0-----0
//! 0       4       1
//!
//! </pre>
//! 9-node Quadrilateral Element
//! <pre>
//!
//! 3       6       2
//!   0-----0-----0
//!   |           |
//!   |           |
//! 7 0   8 0     0 5
//!   |           |
//!   |           |
//!   0-----0-----0
//!  0      4       1
//!
//! </pre>
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
  }

  //! Return number of shape functions
  unsigned nfunctions() const override { return Tnfunctions; }

  //! Evaluate shape functions at given local coordinates
  //! \param[in] xi given local coordinates
  //! \retval shapefn Shape function of a given cell
  Eigen::VectorXd shapefn(const VectorDim& xi) const override;

  //! Evaluate gradient of shape functions
  //! \param[in] xi given local coordinates
  //! \retval grad_shapefn Gradient of shape function of a given cell
  Eigen::MatrixXd grad_shapefn(const VectorDim& xi) const override;

  //! Evaluate the B matrix at given local coordinates
  //! \param[in] xi given local coordinates
  //! \retval bmatrix B matrix
  std::vector<Eigen::MatrixXd> bmatrix(const VectorDim& xi) const override;

  //! Return the degree of shape function
  mpm::ShapeFnDegree degree() const override;

  //! Return nodal coordinates of a unit cell
  Eigen::MatrixXd unit_cell_coordinates() const override;

  //! Return the side indices of a cell to calculate the cell length
  //! \retval indices Outer-indices that form the sides of the cell
  Eigen::MatrixXi sides_indices() const override;

  //! Return the corner indices of a cell to calculate the cell volume
  //! \retval indices Outer-indices that form the cell
  Eigen::VectorXi corner_indices() const override;

  //! Return indices of a sub-tetrahedrons in a volume
  //! to check if a point is inside /outside of a hedron
  //! \retval indices Indices that form sub-tetrahedrons
  Eigen::MatrixXi inhedron_indices() const override;
};

}  // namespace mpm
#include "quad_shapefn.tcc"

#endif  // MPM_QUADRILATERALSHAPEFN_H_

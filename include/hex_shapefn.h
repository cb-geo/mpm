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
//! \details 8-noded and 20-noded hexahedron element \n
//! Shape function, gradient shape function, B-matrix, indices \n
//! 8-node (Trilinear) Hexahedron Element \n
//! <pre>
//!        3               2
//!          0_ _ _ _ _ _0
//!         /|           /|
//!        / |          / |
//!     7 0_ |_ _ _ _ _0 6|
//!       |  |         |  |
//!       |  |         |  |
//!       |  0_ _ _ _ _|_ 0
//!       | / 0        | / 1
//!       |/           |/
//!       0_ _ _ _ _ _ 0
//!     4               5
//!
//!
//! </pre>
//! 20-node (Serendipity) Hexahedron Element \n
//! <pre>
//!        3       13          2
//!          0_ _ _ 0 _ _ _  0
//!          /|             / |
//!      15 0 |         14 0  |
//!        /  0 9         /   |
//!     7 0_ _| _ 0 _ _ _ 0 6 0 11
//!       |   |   19     |    |
//!       |   |      8   |    |
//!       | 0 0_ _ _ 0 _ |_ _ 0  1
//!    17 0  /           0 18 /
//!       | 0 10         |  0 12
//!       |/             | /
//!       0_ _ _ 0 _ _ _ 0
//!     4        16         5
//!
//! </pre>
//!
//! 27-node (Triquadratic) Hexahedron Element \n
//! Check with GMSH \n
//! <pre>
//!          7           18             6
//!            0_ _ _ _ _ 0 _ _ _ _ _ 0
//!            /|                     /|
//!           / |                    / |
//!          /  |    25             /  |
//!      19 0   |     0         17 0   |
//!        /    |      23 0       /    |
//!       /  15 0                /     0 14
//!      /      |               /      |
//!  4  0_ _ _ _|_ 0 _ _ _ _ _ 0 5     |
//!     |       | 16           |       |
//!     |  0 24 |              |   0 22|
//!     |       |          8   |       |
//!     |     0 0_ _ _ _ _ 0_ _|_ _ _  0  1
//!     |      /               |      /
//!  17 0     /    0 25        0 18  /
//!     |    /                 |    /
//!     |10 0         0        |   0 12
//!     |  /         21        |  /
//!     | /                    | /
//!     |/                     |/
//!     0_ _ _ _ _ 0 _ _ _ _ _ 0
//!   4           16            5
//!
//!
//! </pre>
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

  //! Evaluate B matrix at given local coordinates
  //! \param[in] xi given local coordinates
  //! \retval bmatrix B matrix
  std::vector<Eigen::MatrixXd> bmatrix(const VectorDim& xi);

  //! Return nodal coordinates of a unit cell
  Eigen::MatrixXd unit_cell_coordinates() const;

  //! Return the side indices of a cell to calculate the cell length
  //! \retval indices Outer-indices that form the sides of the cell
  Eigen::MatrixXi sides_indices();

  //! Return the corner indices of a cell to calculate the cell volume
  //! \retval indices Outer-indices that form the cell
  Eigen::VectorXi corner_indices();

  //! Return indices of a sub-tetrahedrons in a volume
  //! to check if a point is inside /outside of a hedron
  //! \retval indices Indices that form sub-tetrahedrons
  Eigen::MatrixXi inhedron_indices();
};

}  // namespace mpm
#include "hex_shapefn.tcc"

#endif  // MPM_HEXAHEDRONSHAPEFN_H_

#ifndef MPM_QUADRILATERAL_ELEMENT_H_
#define MPM_QUADRILATERAL_ELEMENT_H_

#include <exception>
#include <map>

#include <Eigen/Dense>

#include "element.h"
#include "logger.h"

namespace mpm {

//! Quadrilateral element class derived from Element class
//! \brief Quadrilateral element
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
//! Face numbering for 4-node, 8-node and 9-node Quadrilateral Element \n
//! <pre>
//!
//!          F2
//!   3 0----------0 2
//!     |          |
//!  F3 |          | F1
//!     |          |
//!     |          |
//!   0 0----------0 1
//!          F0
//! </pre>
//!
//!
//! Namespace containing constants of face id

const std::map<unsigned, Eigen::VectorXi> face_indices_quadrilateral{

    {0, Eigen::Matrix<int, 2, 1>(0, 1)},
    {1, Eigen::Matrix<int, 2, 1>(1, 2)},
    {2, Eigen::Matrix<int, 2, 1>(2, 3)},
    {3, Eigen::Matrix<int, 2, 1>(3, 0)}

};
//!
//! \tparam Tdim Dimension
//! \tparam Tnfunctions Number of functions
template <unsigned Tdim, unsigned Tnfunctions>
class QuadrilateralElement : public Element<Tdim> {

 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! constructor with number of shape functions
  QuadrilateralElement() : mpm::Element<Tdim>() {
    static_assert(Tdim == 2, "Invalid dimension for a quadrilateral element");
    static_assert((Tnfunctions == 4 || Tnfunctions == 8 || Tnfunctions == 9),
                  "Specified number of shape functions is not defined");

    //! Logger
    std::string logger = "quadrilateral::<" + std::to_string(Tdim) + ", " +
                         std::to_string(Tnfunctions) + ">";
    console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
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

  //! Compute Jacobian
  //! \param[in] xi given local coordinates
  //! \param[in] nodal_coordinates Coordinates of nodes forming the cell
  //! \retval jacobian Jacobian matrix
  Eigen::Matrix<double, Tdim, Tdim> jacobian(
      const Eigen::Matrix<double, 2, 1>& xi,
      const Eigen::MatrixXd& nodal_coordinates) const override;

  //! Evaluate the B matrix at given local coordinates
  //! \param[in] xi given local coordinates
  //! \retval bmatrix B matrix
  std::vector<Eigen::MatrixXd> bmatrix(const VectorDim& xi) const override;

  //! Evaluate the B matrix at given local coordinates for a real cell
  //! \param[in] xi given local coordinates
  //! \param[in] nodal_coordinates Coordinates of nodes forming the cell
  //! \retval bmatrix B matrix
  std::vector<Eigen::MatrixXd> bmatrix(
      const VectorDim& xi,
      const Eigen::MatrixXd& nodal_coordinates) const override;

  //! Evaluate the mass matrix
  //! \param[in] xi_s Vector of local coordinates
  //! \retval mass_matrix mass matrix
  Eigen::MatrixXd mass_matrix(
      const std::vector<VectorDim>& xi_s) const override;

  //! Evaluate the Laplace matrix at given local coordinates for a real cell
  //! \param[in] xi_s Vector of local coordinates
  //! \param[in] nodal_coordinates Coordinates of nodes forming the cell
  //! \retval laplace_matrix Laplace matrix
  Eigen::MatrixXd laplace_matrix(
      const std::vector<VectorDim>& xi_s,
      const Eigen::MatrixXd& nodal_coordinates) const override;

  //! Return the degree of shape function
  mpm::ElementDegree degree() const override;

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

  //! Return indices of a face of an element
  //! \param[in] face_id given id of the face
  //! \retval indices Indices that make the face
  Eigen::VectorXi face_indices(unsigned face_id) const override;

 private:
  //! Logger
  std::unique_ptr<spdlog::logger> console_;
};

}  // namespace mpm
#include "quadrilateral_element.tcc"

#endif  // MPM_QUADRILATERAL_ELEMENT_H_

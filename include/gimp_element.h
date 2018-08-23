#ifndef MPM_GIMP_ELEMENT_H_
#define MPM_GIMP_ELEMENT_H_

#include <exception>

#include <Eigen/Dense>

#include "element.h"
#include "logger.h"

namespace mpm {

//!   13----------12----------11----------10
//!   |           |           |           |
//!   |           |           |           |
//!   |           |           |           |
//!   |           |           |           |
//!   |           |           |           |
//! 		        (-1, 1)	    (1,1)
//!   14----------3-----------2-----------9
//!   |           |           |           |
//!   |           |           |           |
//!   |           |   Point   |           |
//!   |           | location  |           |
//!   |           |           |           |
//!   15----------0-----------1-----------8
//!		         (-1,-1)	    (1,-1)
//!   |           |           |           |
//!   |           |           |           |
//!   |           |           |           |
//!   |           |           |           |
//!   |           |           |           |
//!   4-----------5-----------6-----------7
//! </pre>
//! \tparam Tdim Dimension
//! \tparam Tnfunctions Number of functions
template <unsigned Tdim, unsigned Tnfunctions>
class GimpElement : public Element<Tdim> {

 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! constructor with number of shape functions
  GimpElement() : mpm::Element<Tdim>() {
    static_assert(Tdim == 2, "Invalid dimension for a quadrilateral element");
    static_assert((Tnfunctions == 16),
                  "Specified number of shape functions is not defined");

    //! Logger
    std::string logger = "Gimp::<" + std::to_string(Tdim) + ", " +
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

  //! Return number of particles in center cell
  void number_of_particles(const unsigned nparticles) {
    nparticles_ = nparticles;
  }

 public:
  unsigned nparticles_;

 private:
  //! Logger
  std::unique_ptr<spdlog::logger> console_;
};

}  // namespace mpm
#include "gimp_element.tcc"

#endif  // MPM_QUADRILATERAL_ELEMENT_H_

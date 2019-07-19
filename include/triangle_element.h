#ifndef MPM_TRIANGLE_ELEMENT_H_
#define MPM_TRIANGLE_ELEMENT_H_

#include <algorithm>
#include <array>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "triangle_quadrature.h"
#include "triangle_shapefn.h"
#include "node_base.h"
#include "quadrature_base.h"
#include "shapefn_base.h"

#include "element.h"
#include "logger.h"

namespace mpm {

//! Triangle element class derived from Element class
//! \brief Triangle element
//! \details 3-noded and 6-noded triangle element \n
//! Shapte function, gradient shape function, B-matrix, indices \n
//! 3-node Triangle ELement \n
//! <pre>
//!   2 0
//!     |`\
//!     |  `\
//!     |    `\
//!     |      `\
//!     |        `\
//!   0 0----------0 1
//! </pre>
//! 6-node Triangle Element
//! <pre>
//!   2 0
//!     |`\
//!     |  `\
//!   5 0    `0 4
//!     |      `\
//!     |        `\
//!   0 0-----0----0 1
//!           3
//! </pre>
//!
//!
//! \tparam Tdim Dimension
//! \tparam Tnfunctions Number of functions
template <unsigned Tdim, unsigned Tnfunctions>
class TriangleElement : public Element<Tdim> {

 public:
  //! Define vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! constructor with number of shape functions
  TriangleElement() : mpm::Element<Tdim>() {
    static_assert(Tdim == 2, "Invalid dimension for a triangular element");
    static_assert((Tnfunctions == 3 || Tnfunctions == 6),
                  "Specified number of shape funcions is not defined");

    //! Logger
    std::string logger = "triangular::<" + std::to_string(Tdim) + ", " +
                         std::to_string(Tnfunctions) + ">";
    console_ = std::make_unique<spdlog::logger>(logger,mpm::stdout_sink);
  }

  //! Return number of shape functions
  unsigned nfunctions() const override { return Tnfunctions; }

  //! Evaluate shape functions at given local coordinates
  //! \param[in] xi given local coordinates
  //! \param[in] particle_size Particle size
  //! \param[in] deformation_gradient Deformation gradient
  //! \retval shapefn Shape function of a given cell
  Eigen::VectorXd shapefn(const VectorDim& xi, const VectorDim& particle_size,
                          const VectorDim& deformation_gradient) const_override;

  //! Evaluate local shape functions at given coordinates
  //! \param[in] xi given local coordinates
  //! \param[in] particle_size Particle size
  //! \param[in] deformation_gradient Deformation gradient
  //! \retval shapefn_local Shape function of a iven cell
  Eigen:VectorXd shapefn_local(
    const VectorDim& xi, const VectorDim& particle_size,
    const VectorDIm& deformation_gradient) const override;

  //! Evaluate gradient of shape functions
  //! \param[in] xi given local coordinates
  //! \param[in] particle_size Particle size
  //! \param[in] deformation_gradient Deformation gradient
  //! \retval grad_shapefn Gradient of shape function of a given cell
  Eigen::MatrixXd grad_shapefn(
    const VectorDim& xi, const VectorDim& particle_size,
    const VectorDim& deformation_gradient) const override;

  //! Compute Jacobian local
  //! \param[in] xi given local coordinates
  //! \param[in] nodal_coordinates Coordinates of nodes forming the cell
  //! \param[in] particle_size Particle size
  //! \param[in] deformation_gradient Deformation gradient
  //! \retval jacobian Jacobian matrix
  Eigen::Matrix<double, Tdim, Tdim> jacobian_local(
      const VectorDim& xi, const Eigen::MatrixXd& nodal_coordinates,
      const VectorDim& particle_size,
      const VectorDim& deformation_gradient) const override;

  //! Evaluate the B matrix at given local coordinates for a real cell
  //! \param[in] xi given local coordinates
  //! \param[in] nodal_coordinates Coordinates of nodes forming the cell
  //! \param[in] particle_size Particle size
  //! \param[in] deformation_gradient Deformation gradient
  //! \retval bmatrix B matrix
  std::vector<Eigen::MatrixXd> bmatrix(
      const VectorDim& xi, const Eigen::MatrixXd& nodal_coordinates,
      const VectorDim& particle_size,
      const VectorDim& deformation_gradient) const override;

  //! Evaluate the Ni Nj matrix
  //! \param[in] xi_s Vector of local coordinates
  //! \retval ni_nj_matrix Ni Nj matrix
  Eigen::MatrixXd ni_nj_matrix(
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

  //! Return the type of shape function
  mpm::ShapefnType shapefn_type() const override {
    return mpm::ShapefnType::NORMAL_MPM;
  }

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

  //! Return the number of faces in a triangle
  unsigned nfaces() const override { return 3; }

  //! Return unit element length
  double unit_element_length() const override { return 1.; }

  //! Return quadrature of the element
  std::shared_ptr<mpm::Quadrature<Tdim>> quadrature(
      unsigned nquadratures = 1) const override;
  
  //! Return the length of the triangle element
  //! \param[in] recompute Recompute the volume of the element
  //! \retval volume_ area of the element
  double volume(const bool& recompute) {
    // Recompute volume if volume is NaN or recompute is requested
    if (recompute || std::isnan(volume_)) {
      auto node1 = array_nodes_ptr_.at(0)->coordinates();
      auto node2 = array_nodes_ptr_.at(1)->coordinates();
      auto node3 = array_nodes_ptr_.at(2)->coordinates();
      // 2 * Area = (x2 * y3 - x3 * y2) - (x1 * y3 - x3 * y1) + (x1 * y2 - x2 * y1)
      volume_ = std::fabs(((node2.at(0) * node3.at(1)) - (node3.at(0) - node2.at(1))) -
                          ((node1.at(0) * node3.at(1)) - (node3.at(0) - node1.at(1))) +
                          ((node1.at(0) * node2.at(1)) - (node2.at(0) - node1.at(1))))/2.0;
    }
    return volume_;
  }
  
 private:
  //! Logger
  std::unique_ptr<spdlog::logger> console_;
};

} // namespace mpm
#include "triangle_element.tcc"

#endif  // MPM_TRIANGLE_ELEMENT_H_

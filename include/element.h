#ifndef MPM_ELEMENT_H_
#define MPM_ELEMENT_H_

#include <vector>

#include <Eigen/Dense>

namespace mpm {

// Degree of Element
enum ElementDegree { Linear = 1, Quadratic = 2 };

//! Base class of shape functions
//! \brief Base class that stores the information about shape functions
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Element {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Constructor
  //! Assign variables to zero
  Element(){};

  //! Destructor
  virtual ~Element() {}

  //! Return number of shape functions
  virtual unsigned nfunctions() const = 0;

  //! Evaluate shape functions at given local coordinates
  //! \param[in] xi given local coordinates
  virtual Eigen::VectorXd shapefn(const VectorDim& xi) const = 0;

  //! Evaluate gradient of shape functions
  //! \param[in] xi given local coordinates
  virtual Eigen::MatrixXd grad_shapefn(const VectorDim& xi) const = 0;

  //! Compute Jacobian
  //! \param[in] xi given local coordinates
  //! \param[in] nodal_coordinates Coordinates of nodes forming the cell
  //! \retval jacobian Jacobian matrix
  virtual Eigen::Matrix<double, Tdim, Tdim> jacobian(
      const VectorDim& xi, const Eigen::MatrixXd& nodal_coordinates) const = 0;

  //! Evaluate and return the B-matrix
  //! \param[in] xi given local coordinates
  //! \retval bmatrix B matrix
  virtual std::vector<Eigen::MatrixXd> bmatrix(const VectorDim& xi) const = 0;

  //! Evaluate the B matrix at given local coordinates for a real cell
  //! \param[in] xi given local coordinates
  //! \param[in] nodal_coordinates Coordinates of nodes forming the cell
  //! \retval bmatrix B matrix
  virtual std::vector<Eigen::MatrixXd> bmatrix(
      const VectorDim& xi, const Eigen::MatrixXd& nodal_coordinates) const = 0;

  //! Evaluate the mass matrix
  //! \param[in] xi_s Vector of local coordinates
  //! \retval mass_matrix mass matrix
  virtual Eigen::MatrixXd mass_matrix(
      const std::vector<VectorDim>& xi_s) const = 0;

  //! Evaluate the Laplace matrix at given local coordinates for a real cell
  //! \param[in] xi_s Vector of local coordinates
  //! \param[in] nodal_coordinates Coordinates of nodes forming the cell
  //! \retval laplace_matrix Laplace matrix
  virtual Eigen::MatrixXd laplace_matrix(
      const std::vector<VectorDim>& xi_s,
      const Eigen::MatrixXd& nodal_coordinates) const = 0;

  //! Return the degree of element
  virtual mpm::ElementDegree degree() const = 0;

  //! Return nodal coordinates of a unit cell
  virtual Eigen::MatrixXd unit_cell_coordinates() const = 0;

  //! Return the side indices of a cell to calculate the cell length
  //! \retval indices Outer-indices that form the sides of the cell
  virtual Eigen::MatrixXi sides_indices() const = 0;

  //! Return the corner indices of a cell to calculate the cell volume
  //! \retval indices Outer-indices that form the cell
  virtual Eigen::VectorXi corner_indices() const = 0;

  //! Return indices of a sub-tetrahedrons in a volume
  //! to check if a point is inside /outside of a hedron
  //! \retval indices Indices that form sub-tetrahedrons
  virtual Eigen::MatrixXi inhedron_indices() const = 0;

  //! Return indices of a face of an element
  //! \param[in] face_id given id of the face
  //! \retval indices Indices that make the face
  virtual Eigen::VectorXi face_indices(const unsigned face_id) const = 0;
};

}  // namespace mpm
#endif  // MPM_ELEMENTx_H_

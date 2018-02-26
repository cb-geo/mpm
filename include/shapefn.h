#ifndef MPM_SHAPEFN_H_
#define MPM_SHAPEFN_H_

#include <vector>

#include <Eigen/Dense>

namespace mpm {

//! Base class of shape functions
//! \brief Base class that stores the information about shape functions
//! \tparam Tdim Dimension
template <unsigned Tdim>
class ShapeFn {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Constructor
  //! Assign variables to zero
  ShapeFn(){};

  //! Destructor
  virtual ~ShapeFn() {}

  //! Return number of functions
  virtual unsigned nfunctions() const = 0;

  //! Evaluate shape functions at given local coordinates
  //! \param[in] xi given local coordinates
  virtual Eigen::VectorXd shapefn(const VectorDim& xi) = 0;

  //! Evaluate gradient of shape functions
  //! \param[in] xi given local coordinates
  virtual Eigen::MatrixXd grad_shapefn(const VectorDim& xi) = 0;

  //! Return the corner indices of a cell to calculate the cell volume
  //! \retval indices Outer-indices that form the cell
  virtual Eigen::VectorXi corner_indices() = 0;

  //! Return indices of a sub-tetrahedrons in a volume
  //! to check if a point is inside /outside of a hedron
  //! \retval indices Indices that form sub-tetrahedrons
  virtual Eigen::MatrixXi inhedron_indices() = 0;
};

}  // namespace mpm
#endif  // MPM_SHAPEFN_H_

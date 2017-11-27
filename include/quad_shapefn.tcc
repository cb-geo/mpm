// 4-node Quadrilateral Element
//! 3 0----------0 2
//!   |          |
//!   |          |
//!   |          |
//!   |          |
//! 0 0----------0 1

// 8-node Quadrilateral Element
//!  3      6       2
//!   0-----0-----0
//!   |           |
//!   |           |
//! 7 0           0 5
//!   |           |
//!   |           |
//!   0-----0-----0
//! 0       4       1

// 9-node Quadrilateral Element
//! 3       6       2
//!   0-----0-----0
//!   |           |
//!   |           |
//! 7 0   8 0     0 5
//!   |           |
//!   |           |
//!   0-----0-----0
//!  0      4       1

//! Return shape functions of a cell
//! \param[in] xi Coordinates of point of interest
//! \retval shapefn Shape function of a given cell
template <unsigned Tdim>
inline Eigen::MatrixXd mpm::QuadrilateralShapeFn<Tdim>::shapefn(
    const mpm::QuadrilateralShapeFn<Tdim>::VectorDim& xi) {
  switch (this->nfunctions_) {
    case 4:
      // 4-noded
      shapefn_.resize(4, 1);
      shapefn_(0) = 0.25 * (1 - xi(0)) * (1 - xi(1));
      shapefn_(1) = 0.25 * (1 + xi(0)) * (1 - xi(1));
      shapefn_(2) = 0.25 * (1 + xi(0)) * (1 + xi(1));
      shapefn_(3) = 0.25 * (1 - xi(0)) * (1 + xi(1));
      break;
    case 8:
      // 8-noded
      shapefn_.resize(8, 1);
      shapefn_(0) = -0.25 * (1. - xi(0)) * (1. - xi(1)) * (xi(0) + xi(1) + 1.);
      shapefn_(1) = 0.25 * (1. + xi(0)) * (1. - xi(1)) * (xi(0) - xi(1) - 1.);
      shapefn_(2) = 0.25 * (1. + xi(0)) * (1. + xi(1)) * (xi(0) + xi(1) - 1.);
      shapefn_(3) = -0.25 * (1. - xi(0)) * (1. + xi(1)) * (xi(0) - xi(1) + 1.);
      shapefn_(4) = 0.5 * (1. - (xi(0) * xi(0))) * (1. - xi(1));
      shapefn_(5) = 0.5 * (1. - (xi(1) * xi(1))) * (1. + xi(0));
      shapefn_(6) = 0.5 * (1. - (xi(0) * xi(0))) * (1. + xi(1));
      shapefn_(7) = 0.5 * (1. - (xi(1) * xi(1))) * (1. - xi(0));
      break;
    case 9:
      // 9-noded 
      shapefn_.resize(9, 1);
      shapefn_(0) = 0.25 * xi(0) * xi(1) * (xi(0) - 1.) * (xi(1) - 1.);
      shapefn_(1) = 0.25 * xi(0) * xi(1) * (xi(0) + 1.) * (xi(1) - 1.);
      shapefn_(2) = 0.25 * xi(0) * xi(1) * (xi(0) + 1.) * (xi(1) + 1.);
      shapefn_(3) = 0.25 * xi(0) * xi(1) * (xi(0) - 1.) * (xi(1) + 1.);
      shapefn_(4) = -0.5 * xi(1) * (xi(1) - 1.) * ((xi(0) * xi(0)) - 1.);
      shapefn_(5) = -0.5 * xi(0) * (xi(0) + 1.) * ((xi(1) * xi(1)) - 1.);
      shapefn_(6) = -0.5 * xi(1) * (xi(1) + 1.) * ((xi(0) * xi(0)) - 1.);
      shapefn_(7) = -0.5 * xi(0) * (xi(0) - 1.) * ((xi(1) * xi(1)) - 1.);
      shapefn_(8) = ((xi(0) * xi(0)) - 1.) * ((xi(1) * xi(1)) - 1.);
      break;
    default:
      // Throw error
      break;
  }
  return shapefn_;
}

//! Return gradient of shape functions of a cell
//! \param[in] xi Coordinates of point of interest
//! \retval grad_shapefn Gradient of shape function of a given cell
template <unsigned Tdim>
inline Eigen::MatrixXd mpm::QuadrilateralShapeFn<Tdim>::grad_shapefn(
    const mpm::QuadrilateralShapeFn<Tdim>::VectorDim& xi) {
  Eigen::MatrixXd grad_shapefn;
  switch (this->nfunctions_) {
    case 4:
      // 4-noded 
      grad_shapefn_.resize(4, 2);
      grad_shapefn_(0, 0) = -0.25 * (1 - xi(1));
      grad_shapefn_(1, 0) = 0.25 * (1 - xi(1));
      grad_shapefn_(2, 0) = 0.25 * (1 + xi(1));
      grad_shapefn_(3, 0) = -0.25 * (1 + xi(1));

      grad_shapefn_(0, 1) = -0.25 * (1 - xi(0));
      grad_shapefn_(1, 1) = -0.25 * (1 + xi(0));
      grad_shapefn_(2, 1) = 0.25 * (1 + xi(0));
      grad_shapefn_(3, 1) = 0.25 * (1 - xi(0));
      break;
    case 8:
      // 8-noded 
      grad_shapefn_.resize(8, 2);
      grad_shapefn_(0, 0) = 0.25 * (2. * xi(0) + xi(1)) * (1. - xi(1));
      grad_shapefn_(1, 0) = 0.25 * (2. * xi(0) - xi(1)) * (1. - xi(1));
      grad_shapefn_(2, 0) = 0.25 * (2. * xi(0) + xi(1)) * (1. + xi(1));
      grad_shapefn_(3, 0) = 0.25 * (2. * xi(0) - xi(1)) * (1. + xi(1));
      grad_shapefn_(4, 0) = -xi(0) * (1. - xi(1));
      grad_shapefn_(5, 0) = 0.5 * (1. - (xi(1) * xi(1)));
      grad_shapefn_(6, 0) = -xi(0) * (1. + xi(1));
      grad_shapefn_(7, 0) = -0.5 * (1. - (xi(1) * xi(1)));

      grad_shapefn_(0, 1) = 0.25 * (2. * xi(1) + xi(0)) * (1. - xi(0));
      grad_shapefn_(1, 1) = 0.25 * (2. * xi(1) - xi(0)) * (1. + xi(0));
      grad_shapefn_(2, 1) = 0.25 * (2. * xi(1) + xi(0)) * (1. + xi(0));
      grad_shapefn_(3, 1) = 0.25 * (2. * xi(1) - xi(0)) * (1. - xi(0));
      grad_shapefn_(4, 1) = -0.5 * (1. - (xi(0) * xi(0)));
      grad_shapefn_(5, 1) = -xi(1) * (1. + xi(0));
      grad_shapefn_(6, 1) = 0.5 * (1 - (xi(0) * xi(0)));
      grad_shapefn_(7, 1) = -xi(1) * (1. - xi(0));
      break;
    case 9:
      // 9-noded
      grad_shapefn_.resize(9, 2);
      grad_shapefn_(0, 0) = 0.25 * xi(1) * (xi(1) - 1.) * (2 * xi(0) - 1.);
      grad_shapefn_(1, 0) = 0.25 * xi(1) * (xi(1) - 1.) * (2 * xi(0) + 1.);
      grad_shapefn_(2, 0) = 0.25 * xi(1) * (xi(1) + 1.) * (2 * xi(0) + 1.);
      grad_shapefn_(3, 0) = 0.25 * xi(1) * (xi(1) + 1.) * (2 * xi(0) - 1.);
      grad_shapefn_(4, 0) = -xi(0) * xi(1) * (xi(1) - 1.);
      grad_shapefn_(5, 0) = -0.5 * (2. * xi(0) + 1.) * ((xi(1) * xi(1)) - 1.);
      grad_shapefn_(6, 0) = -xi(0) * xi(1) * (xi(1) + 1.);
      grad_shapefn_(7, 0) = -0.5 * (2. * xi(0) - 1.) * ((xi(1) * xi(1)) - 1.);
      grad_shapefn_(8, 0) = 2. * xi(0) * ((xi(1) * xi(1)) - 1.);
      grad_shapefn_(0, 1) = 0.25 * xi(0) * (xi(0) - 1.) * (2. * xi(1) - 1.);
      grad_shapefn_(1, 1) = 0.25 * xi(0) * (xi(0) + 1.) * (2. * xi(1) - 1.);
      grad_shapefn_(2, 1) = 0.25 * xi(0) * (xi(0) + 1.) * (2. * xi(1) + 1.);
      grad_shapefn_(3, 1) = 0.25 * xi(0) * (xi(0) - 1.) * (2. * xi(1) + 1.);
      grad_shapefn_(4, 1) = -0.5 * (2. * xi(1) - 1.) * ((xi(0) * xi(0)) - 1.);
      grad_shapefn_(5, 1) = -xi(0) * xi(1) * (xi(0) + 1.);
      grad_shapefn_(6, 1) = -0.5 * (2. * xi(1) + 1.) * ((xi(0) * xi(0)) - 1.);
      grad_shapefn_(7, 1) = -xi(0) * xi(1) * (xi(0) - 1.);
      grad_shapefn_(8, 1) = 2. * xi(1) * ((xi(0) * xi(0)) - 1.);
      break;
    default:
      // Throw error
      break;
  }
  return grad_shapefn_;
}

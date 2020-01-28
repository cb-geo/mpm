// 4-node Quadrilateral Element
//! 3 0----------0 2
//!   |          |
//!   |          |
//!   |          |
//!   |          |
//! 0 0----------0 1

//! Return shape functions of a 4-node Quadrilateral Element at a given local
//! coordinate, with particle size and deformation gradient
template <>
inline Eigen::VectorXd mpm::QuadrilateralElement<2, 4>::shapefn(
    const Eigen::Matrix<double, 2, 1>& xi,
    const Eigen::Matrix<double, 2, 1>& particle_size,
    const Eigen::Matrix<double, 2, 1>& deformation_gradient) const {
  Eigen::Matrix<double, 4, 1> shapefn;
  shapefn(0) = 0.25 * (1 - xi(0)) * (1 - xi(1));
  shapefn(1) = 0.25 * (1 + xi(0)) * (1 - xi(1));
  shapefn(2) = 0.25 * (1 + xi(0)) * (1 + xi(1));
  shapefn(3) = 0.25 * (1 - xi(0)) * (1 + xi(1));
  return shapefn;
}

//! Return gradient of shape functions of a 4-node Quadrilateral Element at a
//! given local coordinate, with particle size and deformation gradient
template <>
inline Eigen::MatrixXd mpm::QuadrilateralElement<2, 4>::grad_shapefn(
    const Eigen::Matrix<double, 2, 1>& xi,
    const Eigen::Matrix<double, 2, 1>& particle_size,
    const Eigen::Matrix<double, 2, 1>& deformation_gradient) const {
  Eigen::Matrix<double, 4, 2> grad_shapefn;
  grad_shapefn(0, 0) = -0.25 * (1 - xi(1));
  grad_shapefn(1, 0) = 0.25 * (1 - xi(1));
  grad_shapefn(2, 0) = 0.25 * (1 + xi(1));
  grad_shapefn(3, 0) = -0.25 * (1 + xi(1));

  grad_shapefn(0, 1) = -0.25 * (1 - xi(0));
  grad_shapefn(1, 1) = -0.25 * (1 + xi(0));
  grad_shapefn(2, 1) = 0.25 * (1 + xi(0));
  grad_shapefn(3, 1) = 0.25 * (1 - xi(0));

  return grad_shapefn;
}

//! Return nodal coordinates of a unit cell
template <>
inline Eigen::MatrixXd mpm::QuadrilateralElement<2, 4>::unit_cell_coordinates()
    const {
  // Coordinates of a unit cell
  Eigen::Matrix<double, 4, 2> unit_cell;
  // clang-format off
  unit_cell << -1., -1.,
                1., -1.,
                1.,  1.,
    // cppcheck-suppress *
               -1.,  1.;
  // clang-format on
  return unit_cell;
}

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

//! Return shape functions of a 8-node Quadrilateral Element at a given local
//! coordinate, with particle size and deformation gradient
template <>
inline Eigen::VectorXd mpm::QuadrilateralElement<2, 8>::shapefn(
    const Eigen::Matrix<double, 2, 1>& xi,
    const Eigen::Matrix<double, 2, 1>& particle_size,
    const Eigen::Matrix<double, 2, 1>& deformation_gradient) const {
  Eigen::Matrix<double, 8, 1> shapefn;
  shapefn(0) = -0.25 * (1. - xi(0)) * (1. - xi(1)) * (xi(0) + xi(1) + 1.);
  shapefn(1) = 0.25 * (1. + xi(0)) * (1. - xi(1)) * (xi(0) - xi(1) - 1.);
  shapefn(2) = 0.25 * (1. + xi(0)) * (1. + xi(1)) * (xi(0) + xi(1) - 1.);
  shapefn(3) = -0.25 * (1. - xi(0)) * (1. + xi(1)) * (xi(0) - xi(1) + 1.);
  shapefn(4) = 0.5 * (1. - (xi(0) * xi(0))) * (1. - xi(1));
  shapefn(5) = 0.5 * (1. - (xi(1) * xi(1))) * (1. + xi(0));
  shapefn(6) = 0.5 * (1. - (xi(0) * xi(0))) * (1. + xi(1));
  shapefn(7) = 0.5 * (1. - (xi(1) * xi(1))) * (1. - xi(0));
  return shapefn;
}

//! Return gradient of shape functions of a 8-node Quadrilateral Element at a
//! given local coordinate, with particle size and deformation gradient
template <>
inline Eigen::MatrixXd mpm::QuadrilateralElement<2, 8>::grad_shapefn(
    const Eigen::Matrix<double, 2, 1>& xi,
    const Eigen::Matrix<double, 2, 1>& particle_size,
    const Eigen::Matrix<double, 2, 1>& deformation_gradient) const {
  Eigen::Matrix<double, 8, 2> grad_shapefn;
  grad_shapefn(0, 0) = 0.25 * (2. * xi(0) + xi(1)) * (1. - xi(1));
  grad_shapefn(1, 0) = 0.25 * (2. * xi(0) - xi(1)) * (1. - xi(1));
  grad_shapefn(2, 0) = 0.25 * (2. * xi(0) + xi(1)) * (1. + xi(1));
  grad_shapefn(3, 0) = 0.25 * (2. * xi(0) - xi(1)) * (1. + xi(1));
  grad_shapefn(4, 0) = -xi(0) * (1. - xi(1));
  grad_shapefn(5, 0) = 0.5 * (1. - (xi(1) * xi(1)));
  grad_shapefn(6, 0) = -xi(0) * (1. + xi(1));
  grad_shapefn(7, 0) = -0.5 * (1. - (xi(1) * xi(1)));

  grad_shapefn(0, 1) = 0.25 * (2. * xi(1) + xi(0)) * (1. - xi(0));
  grad_shapefn(1, 1) = 0.25 * (2. * xi(1) - xi(0)) * (1. + xi(0));
  grad_shapefn(2, 1) = 0.25 * (2. * xi(1) + xi(0)) * (1. + xi(0));
  grad_shapefn(3, 1) = 0.25 * (2. * xi(1) - xi(0)) * (1. - xi(0));
  grad_shapefn(4, 1) = -0.5 * (1. - (xi(0) * xi(0)));
  grad_shapefn(5, 1) = -xi(1) * (1. + xi(0));
  grad_shapefn(6, 1) = 0.5 * (1 - (xi(0) * xi(0)));
  grad_shapefn(7, 1) = -xi(1) * (1. - xi(0));
  return grad_shapefn;
}

//! Return nodal coordinates of a unit cell
template <>
inline Eigen::MatrixXd mpm::QuadrilateralElement<2, 8>::unit_cell_coordinates()
    const {
  // Coordinates of a unit cell
  Eigen::Matrix<double, 8, 2> unit_cell;
  // clang-format off
  unit_cell << -1., -1.,
                1., -1.,
                1.,  1.,
               -1.,  1.,
                0., -1.,
                1.,  0.,
                0.,  1.,
  // cppcheck-suppress *
               -1.,  0.;
  // clang-format on
  return unit_cell;
}

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

//! Return shape functions of a 9-node Quadrilateral Element at a given local
//! coordinate, with particle size and deformation gradient
template <>
inline Eigen::VectorXd mpm::QuadrilateralElement<2, 9>::shapefn(
    const Eigen::Matrix<double, 2, 1>& xi,
    const Eigen::Matrix<double, 2, 1>& particle_size,
    const Eigen::Matrix<double, 2, 1>& deformation_gradient) const {
  Eigen::Matrix<double, 9, 1> shapefn;

  shapefn(0) = 0.25 * xi(0) * xi(1) * (xi(0) - 1.) * (xi(1) - 1.);
  shapefn(1) = 0.25 * xi(0) * xi(1) * (xi(0) + 1.) * (xi(1) - 1.);
  shapefn(2) = 0.25 * xi(0) * xi(1) * (xi(0) + 1.) * (xi(1) + 1.);
  shapefn(3) = 0.25 * xi(0) * xi(1) * (xi(0) - 1.) * (xi(1) + 1.);
  shapefn(4) = -0.5 * xi(1) * (xi(1) - 1.) * ((xi(0) * xi(0)) - 1.);
  shapefn(5) = -0.5 * xi(0) * (xi(0) + 1.) * ((xi(1) * xi(1)) - 1.);
  shapefn(6) = -0.5 * xi(1) * (xi(1) + 1.) * ((xi(0) * xi(0)) - 1.);
  shapefn(7) = -0.5 * xi(0) * (xi(0) - 1.) * ((xi(1) * xi(1)) - 1.);
  shapefn(8) = ((xi(0) * xi(0)) - 1.) * ((xi(1) * xi(1)) - 1.);

  return shapefn;
}

//! Return gradient of shape functions of a 9-node Quadrilateral Element at a
//! given local coordinate, with particle size and deformation gradient
template <>
inline Eigen::MatrixXd mpm::QuadrilateralElement<2, 9>::grad_shapefn(
    const Eigen::Matrix<double, 2, 1>& xi,
    const Eigen::Matrix<double, 2, 1>& particle_size,
    const Eigen::Matrix<double, 2, 1>& deformation_gradient) const {
  Eigen::Matrix<double, 9, 2> grad_shapefn;
  // 9-noded
  grad_shapefn(0, 0) = 0.25 * xi(1) * (xi(1) - 1.) * (2 * xi(0) - 1.);
  grad_shapefn(1, 0) = 0.25 * xi(1) * (xi(1) - 1.) * (2 * xi(0) + 1.);
  grad_shapefn(2, 0) = 0.25 * xi(1) * (xi(1) + 1.) * (2 * xi(0) + 1.);
  grad_shapefn(3, 0) = 0.25 * xi(1) * (xi(1) + 1.) * (2 * xi(0) - 1.);
  grad_shapefn(4, 0) = -xi(0) * xi(1) * (xi(1) - 1.);
  grad_shapefn(5, 0) = -0.5 * (2. * xi(0) + 1.) * ((xi(1) * xi(1)) - 1.);
  grad_shapefn(6, 0) = -xi(0) * xi(1) * (xi(1) + 1.);
  grad_shapefn(7, 0) = -0.5 * (2. * xi(0) - 1.) * ((xi(1) * xi(1)) - 1.);
  grad_shapefn(8, 0) = 2. * xi(0) * ((xi(1) * xi(1)) - 1.);
  grad_shapefn(0, 1) = 0.25 * xi(0) * (xi(0) - 1.) * (2. * xi(1) - 1.);
  grad_shapefn(1, 1) = 0.25 * xi(0) * (xi(0) + 1.) * (2. * xi(1) - 1.);
  grad_shapefn(2, 1) = 0.25 * xi(0) * (xi(0) + 1.) * (2. * xi(1) + 1.);
  grad_shapefn(3, 1) = 0.25 * xi(0) * (xi(0) - 1.) * (2. * xi(1) + 1.);
  grad_shapefn(4, 1) = -0.5 * (2. * xi(1) - 1.) * ((xi(0) * xi(0)) - 1.);
  grad_shapefn(5, 1) = -xi(0) * xi(1) * (xi(0) + 1.);
  grad_shapefn(6, 1) = -0.5 * (2. * xi(1) + 1.) * ((xi(0) * xi(0)) - 1.);
  grad_shapefn(7, 1) = -xi(0) * xi(1) * (xi(0) - 1.);
  grad_shapefn(8, 1) = 2. * xi(1) * ((xi(0) * xi(0)) - 1.);
  return grad_shapefn;
}

//! Return nodal coordinates of a unit cell
template <>
inline Eigen::MatrixXd mpm::QuadrilateralElement<2, 9>::unit_cell_coordinates()
    const {
  // Coordinates of a unit cell
  Eigen::Matrix<double, 9, 2> unit_cell;
  // clang-format off
  unit_cell << -1., -1.,
                1., -1.,
                1.,  1.,
               -1.,  1.,
                0., -1.,
                1.,  0.,
                0.,  1.,
               -1.,  0.,
  // cppcheck-suppress *
                0.,  0.;
  // clang-format on
  return unit_cell;
}

//! Return the degree of element
//! 4-noded quadrilateral
template <>
inline mpm::ElementDegree mpm::QuadrilateralElement<2, 4>::degree() const {
  return mpm::ElementDegree::Linear;
}

//! Return the degree of element
//! 8-noded quadrilateral
template <>
inline mpm::ElementDegree mpm::QuadrilateralElement<2, 8>::degree() const {
  return mpm::ElementDegree::Quadratic;
}

//! Return the degree of element
//! 9-noded quadrilateral
template <>
inline mpm::ElementDegree mpm::QuadrilateralElement<2, 9>::degree() const {
  return mpm::ElementDegree::Quadratic;
}

//! Return local shape functions of a Quadrilateral Element at a given local
//! coordinate, with particle size and deformation gradient
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::VectorXd
    mpm::QuadrilateralElement<Tdim, Tnfunctions>::shapefn_local(
        const VectorDim& xi, const VectorDim& particle_size,
        const VectorDim& deformation_gradient) const {
  return this->shapefn(xi, particle_size, deformation_gradient);
}

//! Compute Jacobian with particle size and deformation gradient
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::Matrix<double, Tdim, Tdim>
    mpm::QuadrilateralElement<Tdim, Tnfunctions>::jacobian(
        const VectorDim& xi, const Eigen::MatrixXd& nodal_coordinates,
        const VectorDim& particle_size,
        const VectorDim& deformation_gradient) const {

  // Get gradient shape functions
  const Eigen::MatrixXd grad_shapefn =
      mpm::QuadrilateralElement<Tdim, Tnfunctions>::grad_shapefn(
          xi, particle_size, deformation_gradient);

  try {
    // Check if matrices dimensions are correct
    if ((grad_shapefn.rows() != nodal_coordinates.rows()) ||
        (xi.size() != nodal_coordinates.cols()))
      throw std::runtime_error(
          "Jacobian calculation: Incorrect dimension of xi and "
          "nodal_coordinates");
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    return Eigen::Matrix<double, Tdim, Tdim>::Zero();
  }

  // Jacobian dx_i/dxi_j
  return (grad_shapefn.transpose() * nodal_coordinates);
}

//! Compute Jacobian local with particle size and deformation gradient
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::Matrix<double, Tdim, Tdim>
    mpm::QuadrilateralElement<Tdim, Tnfunctions>::jacobian_local(
        const VectorDim& xi, const Eigen::MatrixXd& nodal_coordinates,
        const VectorDim& particle_size,
        const VectorDim& deformation_gradient) const {
  // Jacobian dx_i/dxi_j
  return this->jacobian(xi, nodal_coordinates, particle_size,
                        deformation_gradient);
}

//! Compute Jacobian
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::MatrixXd mpm::QuadrilateralElement<Tdim, Tnfunctions>::dn_dx(
    const VectorDim& xi, const Eigen::MatrixXd& nodal_coordinates,
    const VectorDim& particle_size,
    const VectorDim& deformation_gradient) const {
  // Get gradient shape functions
  Eigen::MatrixXd grad_sf =
      this->grad_shapefn(xi, particle_size, deformation_gradient);

  // Jacobian dx_i/dxi_j
  Eigen::Matrix<double, Tdim, Tdim> jacobian =
      (grad_sf.transpose() * nodal_coordinates);

  // Gradient shapefn of the cell
  // dN/dx = [J]^-1 * dN/dxi
  return grad_sf * (jacobian.inverse()).transpose();
}

//! Return the B-matrix of a Quadrilateral Element at a given local
//! coordinate for a real cell
template <unsigned Tdim, unsigned Tnfunctions>
inline std::vector<Eigen::MatrixXd>
    mpm::QuadrilateralElement<Tdim, Tnfunctions>::bmatrix(
        const VectorDim& xi, const Eigen::MatrixXd& nodal_coordinates,
        const VectorDim& particle_size,
        const VectorDim& deformation_gradient) const {
  // Get gradient shape functions
  Eigen::MatrixXd grad_sf =
      this->grad_shapefn(xi, particle_size, deformation_gradient);

  // B-Matrix
  std::vector<Eigen::MatrixXd> bmatrix;
  bmatrix.reserve(Tnfunctions);

  try {
    // Check if matrices dimensions are correct
    if ((grad_sf.rows() != nodal_coordinates.rows()) ||
        (xi.rows() != nodal_coordinates.cols()))
      throw std::runtime_error(
          "BMatrix - Jacobian calculation: Incorrect dimension of xi and "
          "nodal_coordinates");
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    return bmatrix;
  }

  // Jacobian dx_i/dxi_j
  Eigen::Matrix<double, Tdim, Tdim> jacobian =
      (grad_sf.transpose() * nodal_coordinates);

  // Gradient shapefn of the cell
  // dN/dx = [J]^-1 * dN/dxi
  Eigen::MatrixXd grad_shapefn = grad_sf * (jacobian.inverse()).transpose();

  for (unsigned i = 0; i < Tnfunctions; ++i) {
    Eigen::Matrix<double, 3, Tdim> bi;
    // clang-format off
    bi(0, 0) = grad_shapefn(i, 0); bi(0, 1) = 0.;
    bi(1, 0) = 0.;                 bi(1, 1) = grad_shapefn(i, 1);
    bi(2, 0) = grad_shapefn(i, 1); bi(2, 1) = grad_shapefn(i, 0);
    bmatrix.push_back(bi);
    // clang-format on
  }
  return bmatrix;
}

//! Return ni_nj_matrix of a Hexahedron Element
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::MatrixXd
    mpm::QuadrilateralElement<Tdim, Tnfunctions>::ni_nj_matrix(
        const std::vector<VectorDim>& xi_s) const {
  // Ni Nj matrix
  Eigen::Matrix<double, Tnfunctions, Tnfunctions> ni_nj_matrix;
  ni_nj_matrix.setZero();
  for (const auto& xi : xi_s) {
    const Eigen::Matrix<double, Tnfunctions, 1> shape_fn =
        this->shapefn(xi, Eigen::Matrix<double, Tdim, 1>::Zero(),
                      Eigen::Matrix<double, Tdim, 1>::Zero());
    ni_nj_matrix += (shape_fn * shape_fn.transpose());
  }
  return ni_nj_matrix;
}

//! Return the laplace_matrix of a quadrilateral Element
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::MatrixXd
    mpm::QuadrilateralElement<Tdim, Tnfunctions>::laplace_matrix(
        const std::vector<VectorDim>& xi_s,
        const Eigen::MatrixXd& nodal_coordinates) const {

  try {
    // Check if matrices dimensions are correct
    if ((this->nfunctions() != nodal_coordinates.rows()) ||
        (xi_s.at(0).size() != nodal_coordinates.cols()))
      throw std::runtime_error(
          "Jacobian calculation: Incorrect dimension of xi & nodes");
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }

  // Laplace matrix
  Eigen::Matrix<double, Tnfunctions, Tnfunctions> laplace_matrix;
  laplace_matrix.setZero();
  for (const auto& xi : xi_s) {
    // Get gradient shape functions
    const Eigen::MatrixXd grad_sf =
        this->grad_shapefn(xi, Eigen::Matrix<double, Tdim, 1>::Zero(),
                           Eigen::Matrix<double, Tdim, 1>::Zero());

    // Jacobian dx_i/dxi_j
    const Eigen::Matrix<double, Tdim, Tdim> jacobian =
        (grad_sf.transpose() * nodal_coordinates);

    // Gradient shapefn of the cell
    // dN/dx = [J]^-1 * dN/dxi
    const Eigen::MatrixXd grad_shapefn = grad_sf * jacobian.inverse();

    laplace_matrix += (grad_shapefn * grad_shapefn.transpose());
  }
  return laplace_matrix;
}

//! Return the indices of a cell sides
//! \retval indices Sides that form the cell
//! \tparam Tdim Dimension
//! \tparam Tnfunctions Number of shape functions
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::MatrixXi
    mpm::QuadrilateralElement<Tdim, Tnfunctions>::sides_indices() const {
  Eigen::Matrix<int, 4, 2> indices;
  // clang-format off
  indices << 0, 1,
             1, 2,
             2, 3,
  // cppcheck-suppress *
             3, 0;
  // clang-format on
  return indices;
}

//! Return the corner indices of a cell to calculate the cell volume
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::VectorXi
    mpm::QuadrilateralElement<Tdim, Tnfunctions>::corner_indices() const {
  Eigen::Matrix<int, 4, 1> indices;
  // cppcheck-suppress *
  indices << 0, 1, 2, 3;
  return indices;
}

//! Return indices of a sub-tetrahedrons in a volume
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::MatrixXi
    mpm::QuadrilateralElement<Tdim, Tnfunctions>::inhedron_indices() const {
  Eigen::Matrix<int, 4, Tdim, Eigen::RowMajor> indices;

  // clang-format off
  indices << 0, 1,
             1, 2,
             2, 3,
  // cppcheck-suppress *
             3, 0;
  //clang-format on
  return indices;
}
  
//! Return indices of a face of the element
//! 4-noded quadrilateral
template <>
inline Eigen::VectorXi
    mpm::QuadrilateralElement<2, 4>::face_indices(unsigned face_id) const {
  
  //! Face ids and its associated nodal indices
  const std::map<unsigned, Eigen::Matrix<int, 2, 1>>
      face_indices_quadrilateral{{0, Eigen::Matrix<int, 2, 1>(0, 1)},
                                 {1, Eigen::Matrix<int, 2, 1>(1, 2)},
                                 {2, Eigen::Matrix<int, 2, 1>(2, 3)},
                                 {3, Eigen::Matrix<int, 2, 1>(3, 0)}}; 

  return face_indices_quadrilateral.at(face_id);
}

//! Return indices of a face of the element
//! 8-noded quadrilateral
template <>
inline Eigen::VectorXi
    mpm::QuadrilateralElement<2, 8>::face_indices(unsigned face_id) const {
  
  //! Face ids and its associated nodal indices
  const std::map<unsigned, Eigen::Matrix<int, 3, 1>>
      face_indices_quadrilateral{{0, Eigen::Matrix<int, 3, 1>(0, 1, 4)},
                                 {1, Eigen::Matrix<int, 3, 1>(1, 2, 5)},
                                 {2, Eigen::Matrix<int, 3, 1>(2, 3, 6)},
                                 {3, Eigen::Matrix<int, 3, 1>(3, 0, 7)}};

  return face_indices_quadrilateral.at(face_id);
}

//! Return indices of a face of the element
//! 9-noded quadrilateral
template <>
inline Eigen::VectorXi
    mpm::QuadrilateralElement<2, 9>::face_indices(unsigned face_id) const {
  
  //! Face ids and its associated nodal indices
  const std::map<unsigned, Eigen::Matrix<int, 3, 1>>
      face_indices_quadrilateral{{0, Eigen::Matrix<int, 3, 1>(0, 1, 4)},
                                 {1, Eigen::Matrix<int, 3, 1>(1, 2, 5)},
                                 {2, Eigen::Matrix<int, 3, 1>(2, 3, 6)},
                                 {3, Eigen::Matrix<int, 3, 1>(3, 0, 7)}};

  return face_indices_quadrilateral.at(face_id);
}


//! Return quadrature
template <unsigned Tdim, unsigned Tnfunctions>
inline std::shared_ptr<mpm::Quadrature<Tdim>>
    mpm::QuadrilateralElement<Tdim, Tnfunctions>::quadrature(
        unsigned nquadratures) const {
  switch (nquadratures) {
    case 1:
      return Factory<mpm::Quadrature<Tdim>>::instance()->create("QQ1");
      break;
    case 2:
      return Factory<mpm::Quadrature<Tdim>>::instance()->create("QQ2");
      break;
    case 3:
      return Factory<mpm::Quadrature<Tdim>>::instance()->create("QQ3");
      break;
    case 4:
      return Factory<mpm::Quadrature<Tdim>>::instance()->create("QQ4");
      break;
    default:
      return Factory<mpm::Quadrature<Tdim>>::instance()->create("QQ1");
      break;
  }
}

//! Compute volume
//! \param[in] nodal_coordinates Coordinates of nodes forming the cell
//! \retval volume Return the volume of cell
template <unsigned Tdim, unsigned Tnfunctions>
inline double mpm::QuadrilateralElement<Tdim, Tnfunctions>::compute_volume(
  const Eigen::MatrixXd& nodal_coordinates) const {
  //        b
  // 3 0---------0 2
  //   | \   q / |
  // a |   \  /  | c
  //   |   p \   |
  //   |  /    \ |
  // 0 0---------0 1
  //         d
  const double a = (nodal_coordinates.row(0) -
                    nodal_coordinates.row(3))
    .norm();
  const double b = (nodal_coordinates.row(2) -
                    nodal_coordinates.row(3))
    .norm();
  const double c = (nodal_coordinates.row(1) -
                    nodal_coordinates.row(2))
    .norm();
  const double d = (nodal_coordinates.row(0) -
                    nodal_coordinates.row(1))
    .norm();
  const double p = (nodal_coordinates.row(0) -
                    nodal_coordinates.row(2))
    .norm();
  const double q = (nodal_coordinates.row(1) -
                    nodal_coordinates.row(3))
    .norm();

  // K = 1/4 * sqrt ( 4p^2q^2 - (a^2 + c^2 - b^2 -d^2)^2)
  double volume =
    0.25 * std::sqrt(4 * p * p * q * q -
                     std::pow((a * a + c * c - b * b - d * d), 2.0));

  return volume;
}



//! Compute natural coordinates of a point (analytical)
template <>
inline bool mpm::QuadrilateralElement<2, 4>::isvalid_natural_coordinates_analytical() const { return true; }

//! Compute natural coordinates of a point (analytical)
template <>
inline bool mpm::QuadrilateralElement<2, 8>::isvalid_natural_coordinates_analytical() const { return false; }

//! Compute natural coordinates of a point (analytical)
template <>
inline bool mpm::QuadrilateralElement<2, 9>::isvalid_natural_coordinates_analytical() const { return false; }

//! Compute Natural coordinates of a point (analytical)
//! Analytical solution based on A consistent point-searching algorithm for
//! solution interpolation in unstructured meshes consisting of 4-node bilinear
//! quadrilateral elements - Zhao et al., 1999
template <>
inline Eigen::Matrix<double, 2, 1> mpm::QuadrilateralElement<2, 4>::natural_coordinates_analytical(
      const VectorDim& point,
      const Eigen::MatrixXd& nodal_coordinates) const {
  // Local point coordinates
  Eigen::Matrix<double, 2, 1> xi;
  xi.fill(std::numeric_limits<double>::max());


  const double xa = point(0);
  const double ya = point(1);

  if (nodal_coordinates.rows() == 0 || nodal_coordinates.cols() == 0)
    throw std::runtime_error("Nodal coordinates matrix is empty, cell is probably not initialized");

  const double x1 = nodal_coordinates(0, 0);
  const double y1 = nodal_coordinates(0, 1);
  const double x2 = nodal_coordinates(1, 0);
  const double y2 = nodal_coordinates(1, 1);
  const double x3 = nodal_coordinates(2, 0);
  const double y3 = nodal_coordinates(2, 1);
  const double x4 = nodal_coordinates(3, 0);
  const double y4 = nodal_coordinates(3, 1);

  const double a1 = x1 + x2 + x3 + x4;
  const double a2 = -x1 + x2 + x3 - x4;
  const double a3 = -x1 - x2 + x3 + x4;
  const double a4 = x1 - x2 + x3 - x4;

  const double b1 = y1 + y2 + y3 + y4;
  const double b2 = -y1 + y2 + y3 - y4;
  const double b3 = -y1 - y2 + y3 + y4;
  const double b4 = y1 - y2 + y3 - y4;

  const double c1 = 4. * xa - a1;
  const double c2 = 4. * ya - b1;

  // General solution of xi and eta based on solving with Sympy
  // a2 * xi(0) + a3 * xi(1) + a4 * xi(0) * xi(1) = 4 x_a - a1
  // b2 * xi(0) + b3 * xi(1) + b4 * xi(0) * xi(1) = 4 y_a - b1
  const double u1 =
    (-a1 * b4 + a4 * b1 - 4 * a4 * ya + 4 * b4 * xa +
     (-a3 * b4 + a4 * b3) *
     ((-a1 * b4 + a2 * b3 - a3 * b2 + a4 * b1 - 4 * a4 * ya +
       4 * b4 * xa) /
      (2 * (a3 * b4 - a4 * b3)) -
      (std::sqrt(a1 * a1 * b4 * b4 - 2 * a1 * a2 * b3 * b4 -
                 2 * a1 * a3 * b2 * b4 - 2 * a1 * a4 * b1 * b4 +
                 4 * a1 * a4 * b2 * b3 + 8 * a1 * a4 * b4 * ya -
                 8 * a1 * b4 * b4 * xa + a2 * a2 * b3 * b3 +
                 4 * a2 * a3 * b1 * b4 - 2 * a2 * a3 * b2 * b3 -
                 16 * a2 * a3 * b4 * ya - 2 * a2 * a4 * b1 * b3 +
                 8 * a2 * a4 * b3 * ya + 8 * a2 * b3 * b4 * xa +
                 a3 * a3 * b2 * b2 - 2 * a3 * a4 * b1 * b2 +
                 8 * a3 * a4 * b2 * ya + 8 * a3 * b2 * b4 * xa +
                 a4 * a4 * b1 * b1 - 8 * a4 * a4 * b1 * ya +
                 16 * a4 * a4 * ya * ya + 8 * a4 * b1 * b4 * xa -
                 16 * a4 * b2 * b3 * xa - 32 * a4 * b4 * xa * ya +
                 16 * b4 * b4 * xa * xa)) /
      (2 * (a3 * b4 - a4 * b3)))) /
    (a2 * b4 - a4 * b2);

  const double u2 =
    (-a1 * b4 + a2 * b3 - a3 * b2 + a4 * b1 - 4 * a4 * ya + 4 * b4 * xa) /
    (2 * (a3 * b4 - a4 * b3)) -
    (std::sqrt(
      a1 * a1 * b4 * b4 - 2 * a1 * a2 * b3 * b4 - 2 * a1 * a3 * b2 * b4 -
      2 * a1 * a4 * b1 * b4 + 4 * a1 * a4 * b2 * b3 +
      8 * a1 * a4 * b4 * ya - 8 * a1 * b4 * b4 * xa + a2 * a2 * b3 * b3 +
      4 * a2 * a3 * b1 * b4 - 2 * a2 * a3 * b2 * b3 -
      16 * a2 * a3 * b4 * ya - 2 * a2 * a4 * b1 * b3 +
      8 * a2 * a4 * b3 * ya + 8 * a2 * b3 * b4 * xa + a3 * a3 * b2 * b2 -
      2 * a3 * a4 * b1 * b2 + 8 * a3 * a4 * b2 * ya +
      8 * a3 * b2 * b4 * xa + a4 * a4 * b1 * b1 - 8 * a4 * a4 * b1 * ya +
      16 * a4 * a4 * ya * ya + 8 * a4 * b1 * b4 * xa -
      16 * a4 * b2 * b3 * xa - 32 * a4 * b4 * xa * ya +
      16 * b4 * b4 * xa * xa)) /
    (2 * (a3 * b4 - a4 * b3));

  // Second solution of a quadratic equation
  const double v1 =
    (-a1 * b4 + a4 * b1 - 4 * a4 * ya + 4 * b4 * xa +
     (-a3 * b4 + a4 * b3) *
     ((-a1 * b4 + a2 * b3 - a3 * b2 + a4 * b1 - 4 * a4 * ya +
       4 * b4 * xa) /
      (2 * (a3 * b4 - a4 * b3)) +
      (std::sqrt(a1 * a1 * b4 * b4 - 2 * a1 * a2 * b3 * b4 -
                 2 * a1 * a3 * b2 * b4 - 2 * a1 * a4 * b1 * b4 +
                 4 * a1 * a4 * b2 * b3 + 8 * a1 * a4 * b4 * ya -
                 8 * a1 * b4 * b4 * xa + a2 * a2 * b3 * b3 +
                 4 * a2 * a3 * b1 * b4 - 2 * a2 * a3 * b2 * b3 -
                 16 * a2 * a3 * b4 * ya - 2 * a2 * a4 * b1 * b3 +
                 8 * a2 * a4 * b3 * ya + 8 * a2 * b3 * b4 * xa +
                 a3 * a3 * b2 * b2 - 2 * a3 * a4 * b1 * b2 +
                 8 * a3 * a4 * b2 * ya + 8 * a3 * b2 * b4 * xa +
                 a4 * a4 * b1 * b1 - 8 * a4 * a4 * b1 * ya +
                 16 * a4 * a4 * ya * ya + 8 * a4 * b1 * b4 * xa -
                 16 * a4 * b2 * b3 * xa - 32 * a4 * b4 * xa * ya +
                 16 * b4 * b4 * xa * xa)) /
      (2 * (a3 * b4 - a4 * b3)))) /
    (a2 * b4 - a4 * b2);

  const double v2 =
    (-a1 * b4 + a2 * b3 - a3 * b2 + a4 * b1 - 4 * a4 * ya + 4 * b4 * xa) /
    (2 * (a3 * b4 - a4 * b3)) +
    (std::sqrt(
      a1 * a1 * b4 * b4 - 2 * a1 * a2 * b3 * b4 - 2 * a1 * a3 * b2 * b4 -
      2 * a1 * a4 * b1 * b4 + 4 * a1 * a4 * b2 * b3 +
      8 * a1 * a4 * b4 * ya - 8 * a1 * b4 * b4 * xa + a2 * a2 * b3 * b3 +
      4 * a2 * a3 * b1 * b4 - 2 * a2 * a3 * b2 * b3 -
      16 * a2 * a3 * b4 * ya - 2 * a2 * a4 * b1 * b3 +
      8 * a2 * a4 * b3 * ya + 8 * a2 * b3 * b4 * xa + a3 * a3 * b2 * b2 -
      2 * a3 * a4 * b1 * b2 + 8 * a3 * a4 * b2 * ya +
      8 * a3 * b2 * b4 * xa + a4 * a4 * b1 * b1 - 8 * a4 * a4 * b1 * ya +
      16 * a4 * a4 * ya * ya + 8 * a4 * b1 * b4 * xa -
      16 * a4 * b2 * b3 * xa - 32 * a4 * b4 * xa * ya +
      16 * b4 * b4 * xa * xa)) /
    (2 * (a3 * b4 - a4 * b3));

  // Choosing a quadratic solution
  if (u1 >= -1. && u1 <= 1. && u2 >= -1. && u2 <= 1.) {
    xi(0) = u1;
    xi(1) = u2;
    return xi;
  } else if (v1 >= -1. && v1 <= 1. && v2 >= -1. && v2 <= 1.) {
    xi(0) = v1;
    xi(1) = v2;
    return xi;
  }

  // Case1: a4 == 0 and b4 != 0: Eq 10
  if (a4 == 0 && b4 == 0) {
    xi(0) = (b3 * c1 - a3 * c2) / (a2 * b3 - a3 * b2);
    xi(1) = (-b2 * c1 + a2 * c2) / (a2 * b3 - a3 * b2);
  } else if (a4 == 0 && b4 != 0) {  // Case 2: Eq 11
    if (a2 == 0 && a3 != 0) {       // Case 2.1 Eq 12
      xi(1) = c1 / a3;
      xi(0) = (c2 - b3 * xi(1)) / (b2 + b4 * xi(1));
    } else if (a2 != 0 && a3 == 0) {  // Case 2.2 Eq 13
      xi(0) = c1 / a2;
      xi(1) = (c2 - b2 * xi(0)) / (b3 + b4 * xi(0));
    } else {  // a2 != 0 && a3 != 0 // Case 2.3 Eq 14
      const double aa = b4 * a3 / a2;
      const double bb = ((b2 * a3 - b4 * c1) / a2) - b3;
      const double cc = -(b2 * c1 / a2) + c2;
      // There are two possible solutions
      const double u2 = (-bb + std::sqrt(bb * bb - 4 * aa * cc)) / (2 * aa);
      const double u1 = (c1 - a3 * u2) / a2;
      // Second solution of a quadratic equation
      const double v2 = (-bb - std::sqrt(bb * bb - 4 * aa * cc)) / (2 * aa);
      const double v1 = (c1 - a3 * v2) / a2;
      if (u1 >= -1. && u1 <= 1. && u2 >= -1. && u2 <= 1.) {
        xi(0) = u1;
        xi(1) = u2;
      } else if (v1 >= -1. && v1 <= 1. && v2 >= -1. && v2 <= 1.) {
        xi(0) = v1;
        xi(1) = v2;
      }
    }
  } else if (a4 != 0 && b4 == 0) {  // Case 3 Eq 16
    if (b2 == 0 && b3 != 0) {       // Case 3.1 Eq 17
      xi(1) = c2 / b3;
      xi(0) = (c1 - a3 * xi(1)) / (a2 + a4 * xi(1));
    } else if (b2 != 0 && b3 == 0) {  // Case 3.2 Eq 18
      xi(0) = c2 / b2;
      xi(1) = (c1 - a2 * xi(0)) / (a3 + a4 * xi(0));
    } else {  // b2 != 0 && b3 != 0  // Case 3.3 Eq 19
      const double aa = a4 * b3 / b2;
      const double bb = ((a2 * b3 - a4 * c2) / b2) - b3;
      const double cc = -(a2 * c2 / b2) + c1;
      // There are two possible solutions
      const double u2 = (-bb + std::sqrt((bb * bb - 4 * aa * cc))) / (2 * aa);
      const double u1 = (c2 - b3 * u2) / b2;
      // Second solution of a quadratic equation
      const double v2 = (-bb - std::sqrt((bb * bb - 4 * aa * cc))) / (2 * aa);
      const double v1 = (c2 - b3 * v2) / b2;
      if (u1 >= -1. && u1 <= 1. && u2 >= -1. && u2 <= 1.) {
        xi(0) = u1;
        xi(1) = u2;
      } else if (v1 >= -1. && v1 <= 1. && v2 >= -1. && v2 <= 1.) {
        xi(0) = v1;
        xi(1) = v2;
      }
    }
  } else {  // a4 != 0 && b4 != 0   // Case 4 Eq 21
    const double a2s = a2 / a4;
    const double a3s = a3 / a4;

    const double b2s = b2 / b4;
    const double b3s = b3 / b4;

    const double c1s = c1 / a4;
    const double c2s = c2 / b4;

    if ((a2s - b2s) == 0) {  // Case 4.1 Eq 25
      xi(1) = (c1s - c2s) / (a3s - b3s);
      xi(0) = (c1s - a3s * xi(1)) / (a2s + xi(1));
    } else {
      const double alpha = (c1s - c2s) / (a2s - b2s);
      const double beta = (a3s - b3s) / (a2s - b2s);

      if (beta == 0) {  // Case 4.2a Eq 28
        xi(0) = alpha;
        xi(1) = (c1s - a2s * xi(0)) / (a3s + xi(0));
      } else {  // Case 4.2b Eq 29
        // There are two possible solutions
        const double u2 = (-(a2s * beta + a3s - alpha) +
                           std::sqrt((a2s * beta + a3s - alpha) *
                                     (a2s * beta + a3s - alpha) -
                                     (4 * beta * (c1s - a2s * alpha)))) /
          (2. * beta);
        const double u1 = alpha - beta * u2;
        // Second solution of a quadratic equation
        const double v2 = (-(a2s * beta + a3s - alpha) -
                           std::sqrt((a2s * beta + a3s - alpha) *
                                     (a2s * beta + a3s - alpha) -
                                     (4 * beta * (c1s - a2s * alpha)))) /
          (2. * beta);
        const double v1 = alpha - beta * v2;
        if (u1 >= -1. && u1 <= 1. && u2 >= -1. && u2 <= 1.) {
          xi(0) = u1;
          xi(1) = u2;
        } else if (v1 >= -1. && v1 <= 1. && v2 >= -1. && v2 <= 1.) {
          xi(0) = v1;
          xi(1) = v2;
        }
      }
    }
  }
  return xi;
}

//! Compute natural coordinates of a point (analytical)
template <>
inline Eigen::Matrix<double, 2, 1> mpm::QuadrilateralElement<2, 8>::natural_coordinates_analytical(
      const VectorDim& point,
      const Eigen::MatrixXd& nodal_coordinates) const {
  // Local point coordinates
  Eigen::Matrix<double, 2, 1> xi;
  xi.fill(std::numeric_limits<double>::max());
  throw std::runtime_error("Analytical solution for Quad<2, 8> has not been implemented");
  return xi;
}

//! Compute natural coordinates of a point (analytical)
template <>
inline Eigen::Matrix<double, 2, 1> mpm::QuadrilateralElement<2, 9>::natural_coordinates_analytical(
      const VectorDim& point,
      const Eigen::MatrixXd& nodal_coordinates) const {
  // Local point coordinates
  Eigen::Matrix<double, 2, 1> xi;
  xi.fill(std::numeric_limits<double>::max());
  throw std::runtime_error("Analytical solution for Quad<2, 9> has not been implemented");
  return xi;
}

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
  Eigen::MatrixXd grad_shapefn = grad_sf * jacobian.inverse();

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

//! Return mass_matrix of a Hexahedron Element
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::MatrixXd
    mpm::QuadrilateralElement<Tdim, Tnfunctions>::mass_matrix(
        const std::vector<VectorDim>& xi_s) const {
  // Mass matrix
  Eigen::Matrix<double, Tnfunctions, Tnfunctions> mass_matrix;
  mass_matrix.setZero();
  for (const auto& xi : xi_s) {
    const Eigen::Matrix<double, Tnfunctions, 1> shape_fn =
        this->shapefn(xi, Eigen::Matrix<double, Tdim, 1>::Zero(),
                      Eigen::Matrix<double, Tdim, 1>::Zero());
    mass_matrix += (shape_fn * shape_fn.transpose());
  }
  return mass_matrix;
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
             3, 0;
  // clang-format on
  return indices;
}

//! Return the corner indices of a cell to calculate the cell volume
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::VectorXi
    mpm::QuadrilateralElement<Tdim, Tnfunctions>::corner_indices() const {
  Eigen::Matrix<int, 4, 1> indices;
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
    default:
      return Factory<mpm::Quadrature<Tdim>>::instance()->create("QQ1");
      break;
  }
}

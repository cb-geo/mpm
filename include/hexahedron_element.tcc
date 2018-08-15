// 8-node (Trilinear) Hexahedron Element
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

//! Return shape function of a 8-noded hexahedron
//! \param[in] xi Coordinates of point of interest
//! \retval shapefn Shape function of a given cell
template <>
inline Eigen::VectorXd mpm::HexahedronElement<3, 8>::shapefn(
    const Eigen::Matrix<double, 3, 1>& xi) const {
  // 8-noded
  Eigen::Matrix<double, 8, 1> shapefn;
  shapefn(0) = 0.125 * (1 - xi(0)) * (1 - xi(1)) * (1 - xi(2));
  shapefn(1) = 0.125 * (1 + xi(0)) * (1 - xi(1)) * (1 - xi(2));
  shapefn(2) = 0.125 * (1 + xi(0)) * (1 + xi(1)) * (1 - xi(2));
  shapefn(3) = 0.125 * (1 - xi(0)) * (1 + xi(1)) * (1 - xi(2));
  shapefn(4) = 0.125 * (1 - xi(0)) * (1 - xi(1)) * (1 + xi(2));
  shapefn(5) = 0.125 * (1 + xi(0)) * (1 - xi(1)) * (1 + xi(2));
  shapefn(6) = 0.125 * (1 + xi(0)) * (1 + xi(1)) * (1 + xi(2));
  shapefn(7) = 0.125 * (1 - xi(0)) * (1 + xi(1)) * (1 + xi(2));
  return shapefn;
}

//! Return gradient of shape functions of a 8-noded hexahedron
//! \param[in] xi Coordinates of point of interest
//! \retval grad_shapefn Gradient of shape function of a given cell
template <>
inline Eigen::MatrixXd mpm::HexahedronElement<3, 8>::grad_shapefn(
    const Eigen::Matrix<double, 3, 1>& xi) const {
  Eigen::Matrix<double, 8, 3> grad_shapefn;
  grad_shapefn(0, 0) = -0.125 * (1 - xi(1)) * (1 - xi(2));
  grad_shapefn(1, 0) = 0.125 * (1 - xi(1)) * (1 - xi(2));
  grad_shapefn(2, 0) = 0.125 * (1 + xi(1)) * (1 - xi(2));
  grad_shapefn(3, 0) = -0.125 * (1 + xi(1)) * (1 - xi(2));
  grad_shapefn(4, 0) = -0.125 * (1 - xi(1)) * (1 + xi(2));
  grad_shapefn(5, 0) = 0.125 * (1 - xi(1)) * (1 + xi(2));
  grad_shapefn(6, 0) = 0.125 * (1 + xi(1)) * (1 + xi(2));
  grad_shapefn(7, 0) = -0.125 * (1 + xi(1)) * (1 + xi(2));

  grad_shapefn(0, 1) = -0.125 * (1 - xi(0)) * (1 - xi(2));
  grad_shapefn(1, 1) = -0.125 * (1 + xi(0)) * (1 - xi(2));
  grad_shapefn(2, 1) = 0.125 * (1 + xi(0)) * (1 - xi(2));
  grad_shapefn(3, 1) = 0.125 * (1 - xi(0)) * (1 - xi(2));
  grad_shapefn(4, 1) = -0.125 * (1 - xi(0)) * (1 + xi(2));
  grad_shapefn(5, 1) = -0.125 * (1 + xi(0)) * (1 + xi(2));
  grad_shapefn(6, 1) = 0.125 * (1 + xi(0)) * (1 + xi(2));
  grad_shapefn(7, 1) = 0.125 * (1 - xi(0)) * (1 + xi(2));

  grad_shapefn(0, 2) = -0.125 * (1 - xi(0)) * (1 - xi(1));
  grad_shapefn(1, 2) = -0.125 * (1 + xi(0)) * (1 - xi(1));
  grad_shapefn(2, 2) = -0.125 * (1 + xi(0)) * (1 + xi(1));
  grad_shapefn(3, 2) = -0.125 * (1 - xi(0)) * (1 + xi(1));
  grad_shapefn(4, 2) = 0.125 * (1 - xi(0)) * (1 - xi(1));
  grad_shapefn(5, 2) = 0.125 * (1 + xi(0)) * (1 - xi(1));
  grad_shapefn(6, 2) = 0.125 * (1 + xi(0)) * (1 + xi(1));
  grad_shapefn(7, 2) = 0.125 * (1 - xi(0)) * (1 + xi(1));
  return grad_shapefn;
}

//! Return nodal coordinates of a unit cell
template <>
inline Eigen::MatrixXd mpm::HexahedronElement<3, 8>::unit_cell_coordinates()
    const {
  // Coordinates of a unit cell
  Eigen::Matrix<double, 8, 3> unit_cell;
  // clang-format off
  unit_cell << -1., -1., -1.,
                1., -1., -1.,
                1.,  1., -1.,
               -1.,  1., -1.,
               -1., -1.,  1.,
                1., -1.,  1.,
                1.,  1.,  1.,
               -1.,  1.,  1.;
  // clang-format on
  return unit_cell;
}

// 20-node (Serendipity) Hexahedron Element
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

//! Return shape function of a 20-noded hexahedron
//! \param[in] xi Coordinates of point of interest
//! \retval shapefn Shape function of a given cell
template <>
inline Eigen::VectorXd mpm::HexahedronElement<3, 20>::shapefn(
    const Eigen::Matrix<double, 3, 1>& xi) const {
  Eigen::Matrix<double, 20, 1> shapefn;
  shapefn(0) = -0.125 * (1 - xi(0)) * (1 - xi(1)) * (1 - xi(2)) *
               (2 + xi(0) + xi(1) + xi(2));
  shapefn(1) = -0.125 * (1 + xi(0)) * (1 - xi(1)) * (1 - xi(2)) *
               (2 - xi(0) + xi(1) + xi(2));
  shapefn(2) = -0.125 * (1 + xi(0)) * (1 + xi(1)) * (1 - xi(2)) *
               (2 - xi(0) - xi(1) + xi(2));
  shapefn(3) = -0.125 * (1 - xi(0)) * (1 + xi(1)) * (1 - xi(2)) *
               (2 + xi(0) - xi(1) + xi(2));
  shapefn(4) = -0.125 * (1 - xi(0)) * (1 - xi(1)) * (1 + xi(2)) *
               (2 + xi(0) + xi(1) - xi(2));
  shapefn(5) = -0.125 * (1 + xi(0)) * (1 - xi(1)) * (1 + xi(2)) *
               (2 - xi(0) + xi(1) - xi(2));
  shapefn(6) = -0.125 * (1 + xi(0)) * (1 + xi(1)) * (1 + xi(2)) *
               (2 - xi(0) - xi(1) - xi(2));
  shapefn(7) = -0.125 * (1 - xi(0)) * (1 + xi(1)) * (1 + xi(2)) *
               (2 + xi(0) - xi(1) - xi(2));

  shapefn(8) = 0.25 * (1 - xi(0) * xi(0)) * (1 - xi(1)) * (1 - xi(2));
  shapefn(11) = 0.25 * (1 - xi(1) * xi(1)) * (1 + xi(0)) * (1 - xi(2));
  shapefn(13) = 0.25 * (1 - xi(0) * xi(0)) * (1 + xi(1)) * (1 - xi(2));
  shapefn(9) = 0.25 * (1 - xi(1) * xi(1)) * (1 - xi(0)) * (1 - xi(2));
  shapefn(10) = 0.25 * (1 - xi(2) * xi(2)) * (1 - xi(0)) * (1 - xi(1));
  shapefn(12) = 0.25 * (1 - xi(2) * xi(2)) * (1 + xi(0)) * (1 - xi(1));
  shapefn(14) = 0.25 * (1 - xi(2) * xi(2)) * (1 + xi(0)) * (1 + xi(1));
  shapefn(15) = 0.25 * (1 - xi(2) * xi(2)) * (1 - xi(0)) * (1 + xi(1));
  shapefn(16) = 0.25 * (1 - xi(0) * xi(0)) * (1 - xi(1)) * (1 + xi(2));
  shapefn(18) = 0.25 * (1 - xi(1) * xi(1)) * (1 + xi(0)) * (1 + xi(2));
  shapefn(19) = 0.25 * (1 - xi(0) * xi(0)) * (1 + xi(1)) * (1 + xi(2));
  shapefn(17) = 0.25 * (1 - xi(1) * xi(1)) * (1 - xi(0)) * (1 + xi(2));
  return shapefn;
}

//! Return gradient of shape functions of a 20-noded hexahedron
//! \param[in] xi Coordinates of point of interest
//! \retval grad_shapefn Gradient of shape function of a given cell
template <>
inline Eigen::MatrixXd mpm::HexahedronElement<3, 20>::grad_shapefn(
    const Eigen::Matrix<double, 3, 1>& xi) const {
  Eigen::Matrix<double, 20, 3> grad_shapefn;

  grad_shapefn(0, 0) =
      0.125 * (2 * xi(0) + xi(1) + xi(2) + 1) * (1 - xi(1)) * (1 - xi(2));
  grad_shapefn(1, 0) =
      -0.125 * (-2 * xi(0) + xi(1) + xi(2) + 1) * (1 - xi(1)) * (1 - xi(2));
  grad_shapefn(2, 0) =
      -0.125 * (-2 * xi(0) - xi(1) + xi(2) + 1) * (1 + xi(1)) * (1 - xi(2));
  grad_shapefn(3, 0) =
      0.125 * (2 * xi(0) - xi(1) + xi(2) + 1) * (1 + xi(1)) * (1 - xi(2));
  grad_shapefn(4, 0) =
      0.125 * (2 * xi(0) + xi(1) - xi(2) + 1) * (1 - xi(1)) * (1 + xi(2));
  grad_shapefn(5, 0) =
      -0.125 * (-2 * xi(0) + xi(1) - xi(2) + 1) * (1 - xi(1)) * (1 + xi(2));
  grad_shapefn(6, 0) =
      -0.125 * (-2 * xi(0) - xi(1) - xi(2) + 1) * (1 + xi(1)) * (1 + xi(2));
  grad_shapefn(7, 0) =
      0.125 * (2 * xi(0) - xi(1) - xi(2) + 1) * (1 + xi(1)) * (1 + xi(2));
  grad_shapefn(8, 0) = -0.5 * xi(0) * (1 - xi(1)) * (1 - xi(2));
  grad_shapefn(9, 0) = -0.25 * (1 - xi(1) * xi(1)) * (1 - xi(2));
  grad_shapefn(10, 0) = -0.25 * (1 - xi(2) * xi(2)) * (1 - xi(1));
  grad_shapefn(11, 0) = 0.25 * (1 - xi(1) * xi(1)) * (1 - xi(2));
  grad_shapefn(12, 0) = 0.25 * (1 - xi(2) * xi(2)) * (1 - xi(1));
  grad_shapefn(13, 0) = -0.5 * xi(0) * (1 + xi(1)) * (1 - xi(2));
  grad_shapefn(14, 0) = 0.25 * (1 - xi(2) * xi(2)) * (1 + xi(1));
  grad_shapefn(15, 0) = -0.25 * (1 - xi(2) * xi(2)) * (1 + xi(1));
  grad_shapefn(16, 0) = -0.5 * xi(0) * (1 - xi(1)) * (1 + xi(2));
  grad_shapefn(18, 0) = 0.25 * (1 - xi(1) * xi(1)) * (1 + xi(2));
  grad_shapefn(19, 0) = -0.5 * xi(0) * (1 + xi(1)) * (1 + xi(2));
  grad_shapefn(17, 0) = -0.25 * (1 - xi(1) * xi(1)) * (1 + xi(2));

  grad_shapefn(0, 1) =
      0.125 * (xi(0) + 2 * xi(1) + xi(2) + 1) * (1 - xi(0)) * (1 - xi(2));
  grad_shapefn(1, 1) =
      0.125 * (-xi(0) + 2 * xi(1) + xi(2) + 1) * (1 + xi(0)) * (1 - xi(2));
  grad_shapefn(2, 1) =
      -0.125 * (-xi(0) - 2 * xi(1) + xi(2) + 1) * (1 + xi(0)) * (1 - xi(2));
  grad_shapefn(3, 1) =
      -0.125 * (xi(0) - 2 * xi(1) + xi(2) + 1) * (1 - xi(1)) * (1 - xi(2));
  grad_shapefn(4, 1) =
      0.125 * (xi(0) + 2 * xi(1) - xi(2) + 1) * (1 - xi(0)) * (1 + xi(2));
  grad_shapefn(5, 1) =
      0.125 * (-xi(0) + 2 * xi(1) - xi(2) + 1) * (1 + xi(0)) * (1 + xi(2));
  grad_shapefn(6, 1) =
      -0.125 * (-xi(0) - 2 * xi(1) - xi(2) + 1) * (1 + xi(0)) * (1 + xi(2));
  grad_shapefn(7, 1) =
      -0.125 * (xi(0) - 2 * xi(1) - xi(2) + 1) * (1 - xi(1)) * (1 + xi(2));
  grad_shapefn(8, 1) = -0.25 * (1 - xi(0) * xi(0)) * (1 - xi(2));
  grad_shapefn(9, 1) = -0.5 * xi(1) * (1 - xi(0)) * (1 - xi(2));
  grad_shapefn(10, 1) = -0.25 * (1 - xi(2) * xi(2)) * (1 - xi(0));
  grad_shapefn(11, 1) = -0.5 * xi(1) * (1 + xi(0)) * (1 - xi(2));
  grad_shapefn(12, 1) = -0.25 * (1 - xi(2) * xi(2)) * (1 + xi(0));
  grad_shapefn(13, 1) = 0.25 * (1 - xi(0) * xi(0)) * (1 - xi(2));
  grad_shapefn(14, 1) = 0.25 * (1 - xi(2) * xi(2)) * (1 + xi(0));
  grad_shapefn(15, 1) = 0.25 * (1 - xi(2) * xi(2)) * (1 - xi(0));
  grad_shapefn(16, 1) = -0.25 * (1 - xi(0) * xi(0)) * (1 + xi(2));
  grad_shapefn(18, 1) = -0.5 * xi(1) * (1 + xi(0)) * (1 + xi(2));
  grad_shapefn(19, 1) = 0.25 * (1 - xi(0) * xi(0)) * (1 + xi(2));
  grad_shapefn(17, 1) = -0.5 * xi(1) * (1 - xi(0)) * (1 + xi(2));

  grad_shapefn(0, 2) =
      0.125 * (xi(0) + xi(1) + 2 * xi(2) + 1) * (1 - xi(0)) * (1 - xi(1));
  grad_shapefn(1, 2) =
      0.125 * (-xi(0) + xi(1) + 2 * xi(2) + 1) * (1 + xi(0)) * (1 - xi(1));
  grad_shapefn(2, 2) =
      0.125 * (-xi(0) - xi(1) + 2 * xi(2) + 1) * (1 + xi(0)) * (1 + xi(1));
  grad_shapefn(3, 2) =
      0.125 * (xi(0) - xi(1) + 2 * xi(2) + 1) * (1 - xi(1)) * (1 + xi(1));
  grad_shapefn(4, 2) =
      -0.125 * (xi(0) + xi(1) - 2 * xi(2) + 1) * (1 - xi(0)) * (1 - xi(1));
  grad_shapefn(5, 2) =
      -0.125 * (-xi(0) + xi(1) - 2 * xi(2) + 1) * (1 + xi(0)) * (1 - xi(1));
  grad_shapefn(6, 2) =
      -0.125 * (-xi(0) - xi(1) - 2 * xi(2) + 1) * (1 + xi(0)) * (1 + xi(1));
  grad_shapefn(7, 2) =
      -0.125 * (xi(0) - xi(1) - 2 * xi(2) + 1) * (1 - xi(1)) * (1 + xi(1));
  grad_shapefn(8, 2) = -0.25 * (1 - xi(0) * xi(0)) * (1 - xi(1));
  grad_shapefn(9, 2) = -0.25 * (1 - xi(1) * xi(1)) * (1 - xi(0));
  grad_shapefn(10, 2) = -0.5 * xi(2) * (1 - xi(0)) * (1 - xi(1));
  grad_shapefn(11, 2) = -0.25 * (1 - xi(1) * xi(1)) * (1 + xi(0));
  grad_shapefn(12, 2) = -0.5 * xi(2) * (1 + xi(0)) * (1 - xi(1));
  grad_shapefn(13, 2) = -0.25 * (1 - xi(0) * xi(0)) * (1 + xi(1));
  grad_shapefn(14, 2) = -0.5 * xi(2) * (1 + xi(0)) * (1 + xi(1));
  grad_shapefn(15, 2) = -0.5 * xi(2) * (1 - xi(0)) * (1 + xi(1));
  grad_shapefn(16, 2) = 0.25 * (1 - xi(0) * xi(0)) * (1 - xi(1));
  grad_shapefn(18, 2) = 0.25 * (1 - xi(1) * xi(1)) * (1 + xi(0));
  grad_shapefn(19, 2) = 0.25 * (1 - xi(0) * xi(0)) * (1 + xi(1));
  grad_shapefn(17, 2) = 0.25 * (1 - xi(1) * xi(1)) * (1 - xi(0));
  return grad_shapefn;
}

//! Compute Jacobian
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::Matrix<double, Tdim, Tdim>
    mpm::HexahedronElement<Tdim, Tnfunctions>::jacobian(
        const Eigen::Matrix<double, 3, 1>& xi,
        const Eigen::MatrixXd& nodal_coordinates) const {
  // Get gradient shape functions
  const Eigen::MatrixXd grad_shapefn = this->grad_shapefn(xi);
  try {
    // Check if dimensions are correct
    if ((grad_shapefn.rows() != nodal_coordinates.rows()) ||
        (xi.size() != nodal_coordinates.cols()))
      throw std::runtime_error(
          "Jacobian calculation: Incorrect dimension of xi and "
          "nodal_coordinates");
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }

  // Jacobian
  return (grad_shapefn.transpose() * nodal_coordinates);
}

//! Return B-matrix of a Hexahedron Element
template <unsigned Tdim, unsigned Tnfunctions>
inline std::vector<Eigen::MatrixXd>
    mpm::HexahedronElement<Tdim, Tnfunctions>::bmatrix(
        const VectorDim& xi) const {
  // Get gradient shape functions
  Eigen::MatrixXd grad_shapefn = this->grad_shapefn(xi);

  // B-Matrix
  std::vector<Eigen::MatrixXd> bmatrix;
  bmatrix.reserve(Tnfunctions);

  for (unsigned i = 0; i < Tnfunctions; ++i) {
    // clang-format off
    Eigen::Matrix<double, 6, Tdim> bi;
    bi(0, 0) = grad_shapefn(i, 0); bi(0, 1) = 0.;                 bi(0, 2) = 0.;
    bi(1, 0) = 0.;                 bi(1, 1) = grad_shapefn(i, 1); bi(1, 2) = 0.;
    bi(2, 0) = 0.;                 bi(2, 1) = 0.;                 bi(2, 2) = grad_shapefn(i, 2);
    bi(3, 0) = grad_shapefn(i, 1); bi(3, 1) = grad_shapefn(i, 0); bi(3, 2) = 0.;
    bi(4, 0) = 0.;                 bi(4, 1) = grad_shapefn(i, 2); bi(4, 2) = grad_shapefn(i, 1);
    bi(5, 0) = grad_shapefn(i, 2); bi(5, 1) = 0.;                 bi(5, 2) = grad_shapefn(i, 0);
    bmatrix.push_back(bi);
  }
  return bmatrix;
}

//! Return B-matrix of a Hexahedron Element at a given local
//! coordinate for a real cell
template <unsigned Tdim, unsigned Tnfunctions>
inline std::vector<Eigen::MatrixXd>
    mpm::HexahedronElement<Tdim, Tnfunctions>::bmatrix(
        const VectorDim& xi, const Eigen::MatrixXd& nodal_coordinates) const {
  // Get gradient shape functions
  Eigen::MatrixXd grad_sf = this->grad_shapefn(xi);

  try {
    // Check if matrices dimensions are correct
    if ((grad_sf.rows() != nodal_coordinates.rows()) ||
      (xi.size() != nodal_coordinates.cols()))
    throw std::runtime_error(
        "BMatrix - Jacobian calculation: Incorrect dimension of xi and "
        "nodal_coordinates");
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }

  // Jacobian dx_i/dxi_j
  Eigen::Matrix<double, Tdim, Tdim> jacobian =
      (grad_sf.transpose() * nodal_coordinates);

  // Gradient shapefn of the cell
  // dN/dx = [J]^-1 * dN/dxi
  Eigen::MatrixXd grad_shapefn = grad_sf * jacobian.inverse();
  
  // B-Matrix
  std::vector<Eigen::MatrixXd> bmatrix;
  bmatrix.reserve(Tnfunctions);

  for (unsigned i = 0; i < Tnfunctions; ++i) {
    // clang-format off
    Eigen::Matrix<double, 6, Tdim> bi;
    bi(0, 0) = grad_shapefn(i, 0); bi(0, 1) = 0.;                 bi(0, 2) = 0.;
    bi(1, 0) = 0.;                 bi(1, 1) = grad_shapefn(i, 1); bi(1, 2) = 0.;
    bi(2, 0) = 0.;                 bi(2, 1) = 0.;                 bi(2, 2) = grad_shapefn(i, 2);
    bi(3, 0) = grad_shapefn(i, 1); bi(3, 1) = grad_shapefn(i, 0); bi(3, 2) = 0.;
    bi(4, 0) = 0.;                 bi(4, 1) = grad_shapefn(i, 2); bi(4, 2) = grad_shapefn(i, 1);
    bi(5, 0) = grad_shapefn(i, 2); bi(5, 1) = 0.;                 bi(5, 2) = grad_shapefn(i, 0);
    bmatrix.push_back(bi);
  }
  return bmatrix;
}

//! Return mass_matrix of a Hexahedron Element
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::MatrixXd
    mpm::HexahedronElement<Tdim, Tnfunctions>::mass_matrix(
      const std::vector<VectorDim>& xi_s) const {
  // Mass matrix
  Eigen::Matrix<double, Tnfunctions, Tnfunctions>  mass_matrix;
  mass_matrix.setZero();
  for (const auto& xi : xi_s) {
    const Eigen::Matrix<double, Tnfunctions, 1> shape_fn = this->shapefn(xi);
    mass_matrix += (shape_fn * shape_fn.transpose());
  }
  return mass_matrix;
}

//! Return the degree of element
//! 8-noded hexahedron
template <>
inline  mpm::ElementDegree
mpm::HexahedronElement<3, 8>::degree() const {
  return mpm::ElementDegree::Linear;
}

//! Return the degree of shape function
//! 8-noded hexahedron
template <>
inline  mpm::ElementDegree
mpm::HexahedronElement<3, 20>::degree() const {
  return mpm::ElementDegree::Quadratic;
}

//! Return nodal coordinates of a unit cell
template <>
inline Eigen::MatrixXd mpm::HexahedronElement<3, 20>::unit_cell_coordinates() const {
  // Coordinates of a unit cell
  Eigen::Matrix<double, 20, 3> unit_cell;
  // clang-format off
  unit_cell << -1., -1., -1.,
                1., -1., -1.,
                1.,  1., -1.,
               -1.,  1., -1.,
               -1., -1.,  1.,
                1., -1.,  1.,
                1.,  1.,  1.,
               -1.,  1.,  1.,
                0., -1., -1.,
               -1.,  0., -1.,
               -1., -1.,  0.,
                1.,  0., -1.,
                1., -1.,  0.,
                0.,  1., -1.,
                1.,  1.,  0.,
               -1.,  1.,  0.,
                0., -1.,  1.,
               -1.,  0.,  1.,
                1.,  0.,  1.,
                0.,  1.,  1.;
  // clang-format on
  return unit_cell;
}

// 27-node (Triquadratic) Hexahedron Element
//! Check with GMSH
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

//! Return the indices of a cell sides
//! \retval indices Sides that form the cell
//! \tparam Tdim Dimension
//! \tparam Tnfunctions Number of shape functions
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::MatrixXi
    mpm::HexahedronElement<Tdim, Tnfunctions>::sides_indices() const {
  Eigen::Matrix<int, 12, 2> indices;
  // clang-format off
  indices << 0, 1,
             1, 2,
             2, 3,
             3, 0,
             4, 5,
             5, 6,
             6, 7,
             7, 4,
             0, 4,
             1, 5,
             2, 6,
             3, 7;
  // clang-format on
  return indices;
}

//! Return the corner indices of a cell to calculate the cell volume
//! \retval indices Outer-indices that form the cell
//! \tparam Tdim Dimension
//! \tparam Tnfunctions Number of shape functions
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::VectorXi
    mpm::HexahedronElement<Tdim, Tnfunctions>::corner_indices() const {
  Eigen::Matrix<int, 8, 1> indices;
  indices << 0, 1, 2, 3, 4, 5, 6, 7;
  return indices;
}

//! Return indices of a sub-tetrahedrons in a volume
//! to check if a point is inside /outside of a hedron
//! \retval indices Indices that form sub-tetrahedrons
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::MatrixXi
    mpm::HexahedronElement<Tdim, Tnfunctions>::inhedron_indices() const {
  Eigen::Matrix<int, 12, Tdim, Eigen::RowMajor> indices;

  // clang-format off
  indices << 0, 5, 4,
             0, 1, 5,
             3, 6, 7,
             3, 2, 6,
             2, 1, 6,
             6, 1, 5,
             7, 6, 5,
             5, 4, 7,
             7, 4, 0,
             7, 0, 3,
             3, 0, 1,
             3, 1, 2;
  //clang-format on
  return indices;
}
  

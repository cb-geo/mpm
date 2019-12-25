// 3-node Triangle Element
//!   2 0
//!     |`\
//!     |  `\
//!     |    `\
//!     |      `\
//!     |        `\
//!   0 0----------0 1

//! Return shape functions of a 3-node Triangle Element at a given local
//! coordinate, with particle size and deformation gradient
template <>
inline Eigen::VectorXd mpm::TriangleElement<2, 3>::shapefn(
    const Eigen::Matrix<double, 2, 1>& xi,
    const Eigen::Matrix<double, 2, 1>& particle_size,
    const Eigen::Matrix<double, 2, 1>& deformation_gradient) const {
  Eigen::Matrix<double, 3, 1> shapefn;
  shapefn(0) = 1 - (xi(0) + xi(1));
  shapefn(1) = xi(0);
  shapefn(2) = xi(1);
  return shapefn;
}

//! Return gradient of shape functions of a 3-node Triangle Element at a
//! given local coordinate, with particle size and deformation gradient
template <>
inline Eigen::MatrixXd mpm::TriangleElement<2, 3>::grad_shapefn(
    const Eigen::Matrix<double, 2, 1>& xi,
    const Eigen::Matrix<double, 2, 1>& particle_size,
    const Eigen::Matrix<double, 2, 1>& deformation_gradient) const {
  Eigen::Matrix<double, 3, 2> grad_shapefn;

  grad_shapefn(0, 0) = -1.;
  grad_shapefn(1, 0) = 1.;
  grad_shapefn(2, 0) = 0.;
  grad_shapefn(0, 1) = -1.;
  grad_shapefn(1, 1) = 0.;
  grad_shapefn(2, 1) = 1.;
  return grad_shapefn;
}

//! Return nodal coordinates of a unit cell
template <>
inline Eigen::MatrixXd mpm::TriangleElement<2, 3>::unit_cell_coordinates()
    const {
  // Coordinates of a unit cell
  Eigen::Matrix<double, 3, 2> unit_cell;
  // clang-format off
  unit_cell << 0., 0.,
               1., 0.,
    // cppcheck-suppress *
               0., 1.;
  // clang-format on
  return unit_cell;
}

// 6-node Triangle Element
//!   2 0
//!     |`\
//!     |  `\
//!   5 0    `0 4
//!     |      `\
//!     |        `\
//!   0 0-----0----0 1
//!           3

//! Return shape functions of a 6-node Triangle Element at a given local
//! coordinate, with particle size and deformation gradient
template <>
inline Eigen::VectorXd mpm::TriangleElement<2, 6>::shapefn(
    const Eigen::Matrix<double, 2, 1>& xi,
    const Eigen::Matrix<double, 2, 1>& particle_size,
    const Eigen::Matrix<double, 2, 1>& deformation_gradient) const {
  Eigen::Matrix<double, 6, 1> shapefn;
  shapefn(0) = (1. - xi(0) - xi(1)) * (1. - 2. * xi(0) - 2. * xi(1));
  shapefn(1) = xi(0) * (2. * xi(0) - 1.);
  shapefn(2) = xi(1) * (2. * xi(1) - 1.);
  shapefn(3) = 4. * xi(0) * (1. - xi(0) - xi(1));
  shapefn(4) = 4. * xi(0) * xi(1);
  shapefn(5) = 4. * xi(1) * (1. - xi(0) - xi(1));
  return shapefn;
}

//! Return gradient of shape functions of a 6-node Triangle Element at a
//! given local coordinate, with particle size and deformation gradient
template <>
inline Eigen::MatrixXd mpm::TriangleElement<2, 6>::grad_shapefn(
    const Eigen::Matrix<double, 2, 1>& xi,
    const Eigen::Matrix<double, 2, 1>& particle_size,
    const Eigen::Matrix<double, 2, 1>& deformation_gradient) const {
  Eigen::Matrix<double, 6, 2> grad_shapefn;
  grad_shapefn(0, 0) = 4. * xi(0) + 4. * xi(1) - 3.;
  grad_shapefn(1, 0) = 4. * xi(0) - 1.;
  grad_shapefn(2, 0) = 0;
  grad_shapefn(3, 0) = 4. - 8. * xi(0) - 4. * xi(1);
  grad_shapefn(4, 0) = 4. * xi(1);
  grad_shapefn(5, 0) = -4. * xi(1);

  grad_shapefn(0, 1) = 4. * xi(0) + 4. * xi(1) - 3.;
  grad_shapefn(1, 1) = 0;
  grad_shapefn(2, 1) = 4. * xi(1) - 1.;
  grad_shapefn(3, 1) = -4. * xi(0);
  grad_shapefn(4, 1) = 4. * xi(0);
  grad_shapefn(5, 1) = 4. - 4. * xi(0) - 8. * xi(1);
  return grad_shapefn;
}

//! Return nodal coordinates of a unit cell
template <>
inline Eigen::MatrixXd mpm::TriangleElement<2, 6>::unit_cell_coordinates()
    const {
  // Coordinates of a unit cell
  Eigen::Matrix<double, 6, 2> unit_cell;
  // clang-format off
  unit_cell << 0. , 0. ,
               1. , 0. ,
               0. , 1. ,
               0.5, 0. ,
               0.5, 0.5,
    // cppcheck-suppress *
               0. , 0.5;
  // clang-format on
  return unit_cell;
}

//! Return the degree of element
//! 3-noded triangle
template <>
inline mpm::ElementDegree mpm::TriangleElement<2, 3>::degree() const {
  return mpm::ElementDegree::Linear;
}

//! Return the degree of element
//! 6-noded triangle
template <>
inline mpm::ElementDegree mpm::TriangleElement<2, 6>::degree() const {
  return mpm::ElementDegree::Quadratic;
}

//! Return local shape functions of a Triangle Element at a given local
//! coordinate, with particle size and deformation gradient
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::VectorXd mpm::TriangleElement<Tdim, Tnfunctions>::shapefn_local(
    const VectorDim& xi, const VectorDim& particle_size,
    const VectorDim& deformation_gradient) const {
  return this->shapefn(xi, particle_size, deformation_gradient);
}

//! Compute Jacobian with particle size and deformation gradient
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::Matrix<double, Tdim, Tdim>
    mpm::TriangleElement<Tdim, Tnfunctions>::jacobian(
        const VectorDim& xi, const Eigen::MatrixXd& nodal_coordinates,
        const VectorDim& particle_size,
        const VectorDim& deformation_gradient) const {

  // Get gradient shape functions
  const Eigen::MatrixXd grad_shapefn =
      mpm::TriangleElement<Tdim, Tnfunctions>::grad_shapefn(
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
    mpm::TriangleElement<Tdim, Tnfunctions>::jacobian_local(
        const VectorDim& xi, const Eigen::MatrixXd& nodal_coordinates,
        const VectorDim& particle_size,
        const VectorDim& deformation_gradient) const {
  // Jacobian dx_i/dxi_j
  return this->jacobian(xi, nodal_coordinates, particle_size,
                        deformation_gradient);
}

//! Compute Jacobian
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::MatrixXd mpm::TriangleElement<Tdim, Tnfunctions>::dn_dx(
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

//! Return the B-matrix of a Triangle Element at a given local
//! coordinate for a real cell
template <unsigned Tdim, unsigned Tnfunctions>
inline std::vector<Eigen::MatrixXd>
    mpm::TriangleElement<Tdim, Tnfunctions>::bmatrix(
        const VectorDim& xi, const Eigen::MatrixXd& nodal_coordinates,
        const VectorDim& particle_size,
        const VectorDim& deformation_gradient) const {
  // Get gradient shape functions
  Eigen::MatrixXd grad_sf =
      this->grad_shapefn(xi, particle_size, deformation_gradient);

  // B-matrix
  std::vector<Eigen::MatrixXd> bmatrix;
  bmatrix.reserve(Tnfunctions);

  try {
    // Check if matrces dimensions are correct
    if ((grad_sf.rows() != nodal_coordinates.rows()) ||
        (xi.rows() != nodal_coordinates.cols()))
      throw std::runtime_error(
          "Bmatrix - Jacobian calculation: Incorrect dimension of xi and "
          "nodal_coordinates");
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    return bmatrix;
  }

  // Jacobian dx_i/dxi_j
  Eigen::Matrix<double, Tdim, Tdim> jacobian =
      (grad_sf.transpose() * nodal_coordinates);

  // Gradient shapefn of the cell
  // dN/dx = [J]^-1 * dN/dxi_j
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

//! Return ni_nj_matrix of a Triangle Element
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::MatrixXd mpm::TriangleElement<Tdim, Tnfunctions>::ni_nj_matrix(
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

//! Return the laplace_matrix of a Triangle Element
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::MatrixXd mpm::TriangleElement<Tdim, Tnfunctions>::laplace_matrix(
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
inline Eigen::MatrixXi mpm::TriangleElement<Tdim, Tnfunctions>::sides_indices()
    const {
  Eigen::Matrix<int, 3, 2> indices;
  // clang-format off
  indices << 0, 1,
             1, 2,
    // cppcheck-suppress *
             2, 0;
  // clang-format on
  return indices;
}

//! Return the corner indices of a cell to calculate the cell volume
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::VectorXi mpm::TriangleElement<Tdim, Tnfunctions>::corner_indices()
    const {
  Eigen::Matrix<int, 3, 1> indices;
  // cppcheck-suppress *
  indices << 0, 1, 2;
  return indices;
}

//! Return indices of a sub-tetrahedrons in a volume
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::MatrixXi
    mpm::TriangleElement<Tdim, Tnfunctions>::inhedron_indices() const {
  Eigen::Matrix<int, 3, Tdim, Eigen::RowMajor> indices;

  // clang-format off
  indices << 0, 1,
             1, 2,
    // cppcheck-suppress *
             2, 0;
  //clang-format on
  return indices;
}

//! Return indices of a face of the element
//! 3-noded triangle
template <>
inline Eigen::VectorXi
    mpm::TriangleElement<2, 3>::face_indices(unsigned face_id) const {
  
  //! Face ids and its associated nodal indices
  const std::map<unsigned, Eigen::Matrix<int, 2, 1>>
      face_indices_triangle{{0, Eigen::Matrix<int, 2, 1>(0, 1)},
                            {1, Eigen::Matrix<int, 2, 1>(1, 2)},
                            {2, Eigen::Matrix<int, 2, 1>(2, 0)}}; 

  return face_indices_triangle.at(face_id);
}

//! Return indices of a face of the element
//! 6-noded triangle
template <>
inline Eigen::VectorXi
    mpm::TriangleElement<2, 6>::face_indices(unsigned face_id) const {
  
  //! Face ids and its associated nodal indices
  const std::map<unsigned, Eigen::Matrix<int, 3, 1>>
      face_indices_triangle{{0, Eigen::Matrix<int, 3, 1>(0, 1, 3)},
                            {1, Eigen::Matrix<int, 3, 1>(1, 2, 4)},
                            {2, Eigen::Matrix<int, 3, 1>(2, 0, 5)}};

  return face_indices_triangle.at(face_id);
}

//! Return quadrature
template <unsigned Tdim, unsigned Tnfunctions>
inline std::shared_ptr<mpm::Quadrature<Tdim>>
    mpm::TriangleElement<Tdim, Tnfunctions>::quadrature(
        unsigned nquadratures) const {
  switch (nquadratures) {
    case 1:
      return Factory<mpm::Quadrature<Tdim>>::instance()->create("QT1");
      break;
    case 2:
      return Factory<mpm::Quadrature<Tdim>>::instance()->create("QT2");
      break;
    default:
      return Factory<mpm::Quadrature<Tdim>>::instance()->create("QT1");
      break;
  }
}

//! Compute volume
//! \param[in] nodal_coordinates Coordinates of nodes forming the cell
//! \retval volume Return the volume of cell
template <unsigned Tdim, unsigned Tnfunctions>
inline double mpm::TriangleElement<Tdim, Tnfunctions>::compute_volume(
  const Eigen::MatrixXd& nodal_coordinates) const {
  //   2 0
  //     |`\
  //     |  `\
  //     |    `\
  //     |      `\
  //     |        `\
  //   0 0----------0 1
  //
  auto node0 = nodal_coordinates.row(0);
  auto node1 = nodal_coordinates.row(1);
  auto node2 = nodal_coordinates.row(2);

  // Area = 0.5 * [ (x1 * y2 - x2 * y1)
  //              - (x0 * y2 - x2 * y0)
  //              + (x0 * y1 - x1 * y0) ]
  double volume_ = std::fabs(((node1(0) * node2(1)) - (node2(0) * node1(1))) -
                             ((node0(0) * node2(1)) - (node2(0) * node0(1))) +
                             ((node0(0) * node1(1)) - (node1(0) * node0(1)))) *
    0.5;
      
  return volume_;
}

//! Compute natural coordinates of a point (analytical)
template <>
inline bool mpm::TriangleElement<2, 3>::isvalid_natural_coordinates_analytical() const { return true; }

//! Compute natural coordinates of a point (analytical)
template <>
inline bool mpm::TriangleElement<2, 6>::isvalid_natural_coordinates_analytical() const { return false; }

//! Compute Natural coordinates of a point (analytical)
template <>
inline Eigen::Matrix<double, 2, 1> mpm::TriangleElement<2, 3>::natural_coordinates_analytical(
      const VectorDim& point,
      const Eigen::MatrixXd& nodal_coordinates) const {
  // Local point coordinates
  Eigen::Matrix<double, 2, 1> xi;
  xi.fill(std::numeric_limits<double>::max());

  // initialize cartesian coordinates of point of interest and vertices of the cell
  const double xa = point(0);
  const double ya = point(1);
  const double x1 = nodal_coordinates(0, 0);
  const double y1 = nodal_coordinates(0, 1);
  const double x2 = nodal_coordinates(1, 0);
  const double y2 = nodal_coordinates(1, 1);
  const double x3 = nodal_coordinates(2, 0);
  const double y3 = nodal_coordinates(2, 1);
  
  // calculate matrix determinant of T
  const double T = (x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3);
  
  // determine the barycentric coordinates lambda1, lambda2 and lambda3
  const double lambda1 = ((y2 - y3) * (xa - x3) + (x3 - x2) * (ya - y3)) / T;
  const double lambda2 = ((y3 - y1) * (xa - x3) + (x1 - x3) * (ya - y3)) / T;
  const double lambda3 = 1 - lambda1 - lambda2;
  
  xi(0) = lambda2;
  xi(1) = lambda3;

  return xi;
}


//! Compute natural coordinates of a point (analytical)
template <>
inline Eigen::Matrix<double, 2, 1> mpm::TriangleElement<2, 6>::natural_coordinates_analytical(
      const VectorDim& point,
      const Eigen::MatrixXd& nodal_coordinates) const {
  // Local point coordinates
  Eigen::Matrix<double, 2, 1> xi;
  xi.fill(std::numeric_limits<double>::max());
  throw std::runtime_error("Analytical solution for Triangle<2, 6> has not been implemented");
  return xi;
}

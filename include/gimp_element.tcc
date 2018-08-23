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
//!   |           |   particle   |           |
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

//! Return nodal coordinates of a unit cell
template <>
inline Eigen::MatrixXd mpm::GimpElement<2, 16>::unit_cell_coordinates() const {
  // Coordinates of a unit cell
  Eigen::Matrix<double, 16, 2> unit_cell;
  // clang-format off
  unit_cell << -1., -1.,
                1., -1.,
               -1.,  1.,
                1.,  1.,
               -3., -3.,
               -1., -3.,
                1., -3.,
                3., -3.,
                3., -1.,
                3.,  1.,
                3.,  3.,
                1.,  3.,
               -1.,  3.,
               -3.,  3.,
               -3.,  1.,
               -3., -1.;

  // clang-format on
  return unit_cell;
}

//! Return shape functions of a 16-node Quadrilateral GIMP Element at a given
//! local coordinate
template <>
inline Eigen::VectorXd mpm::GimpElement<2, 16>::shapefn(
    const Eigen::Matrix<double, 2, 1>& xi) const {

  Eigen::Matrix<double, 16, 1> shapefn;

  //! Matrix to store local node coordinates
  Eigen::Matrix<double, 16, 2> local_node = this->unit_cell_coordinates();

  //! length and volume of element in local coordinate (TODO: double check)
  const double element_length = 2;
  const double element_volume = 2 * 2;

  //! particle length = element length / number of particles in cell.
  const double particle_length = sqrt(element_volume / nparticles_);

  //! local node x & y valyes
  double xni, yni;
  const double numnodes = 16;

  //! local particle - local node
  double xpxi, ypyi;

  //! local shapefunctions for x and y coordinates
  double sxni, syni;

  for (unsigned n = 0; n < numnodes; n++) {
    xni = local_node(n, 0);
    yni = local_node(n, 1);
    xpxi = xi(0) - xni;  // local particle x - local node x
    ypyi = xi(1) - yni;  // local particle x - local node x

    //! x direction
    if (xpxi >= element_length + particle_length) {
      sxni = 0;
    } else if (xpxi > -element_length - particle_length &&
               xpxi <= -element_length + particle_length) {
      sxni = ((element_length + particle_length + xpxi) *
              (element_length + particle_length + xpxi)) /
             4 * element_length * particle_length;
    } else if (xpxi > -element_length + particle_length &&
               xpxi <= -particle_length) {
      sxni = 1 + (xpxi / element_length);
    } else if (xpxi > -particle_length && xpxi <= particle_length) {
      sxni = 1 - (((xpxi * xpxi) + (particle_length * particle_length)) / 2 *
                  element_length * particle_length);
    } else if (xpxi > particle_length &&
               xpxi <= element_length - particle_length) {
      sxni = 1 - (xpxi / element_length);
    } else if (xpxi > element_length - particle_length &&
               xpxi <= element_length + particle_length) {
      sxni = ((element_length + particle_length + xpxi) *
              (element_length + particle_length + xpxi)) /
             4 * element_length * particle_length;
    }
    if (ypyi >= element_length + particle_length) {
      syni = 0;
    } else if (ypyi > -element_length - particle_length &&
               ypyi <= -element_length + particle_length) {
      syni = ((element_length + particle_length + ypyi) *
              (element_length + particle_length + ypyi)) /
             4 * element_length * particle_length;
    } else if (ypyi > -element_length + particle_length &&
               ypyi <= -particle_length) {
      syni = 1 + (ypyi / element_length);
    } else if (ypyi > -particle_length && ypyi <= particle_length) {
      syni = 1 - (((ypyi * ypyi) + (particle_length * particle_length)) / 2 *
                  element_length * particle_length);
    } else if (ypyi > particle_length &&
               ypyi <= element_length - particle_length) {
      syni = 1 - (ypyi / element_length);
    } else if (ypyi > element_length - particle_length &&
               ypyi <= element_length + particle_length) {
      syni = ((element_length + particle_length + ypyi) *
              (element_length + particle_length + ypyi)) /
             4 * element_length * particle_length;
    }
    shapefn(n) = sxni * syni;
  }

  return shapefn;
}

//! Return gradient of shape functions of a 9-node Quadrilateral Element at a
//! given local coordinate
template <>
inline Eigen::MatrixXd mpm::GimpElement<2, 16>::grad_shapefn(
    const Eigen::Matrix<double, 2, 1>& xi) const {
  Eigen::Matrix<double, 16, 2> grad_shapefn;

  //! Matrix to store local node coordinates
  Eigen::Matrix<double, 16, 2> local_node = this->unit_cell_coordinates();

  //! length and volume of element in local coordinate (TODO: double check)
  const double element_length = 2;
  const double element_volume = 2 * 2;

  //! particle length = element length / number of particles in cell.
  const double particle_length = sqrt(element_volume / nparticles_);

  //! local node x & y valyes
  double xni, yni;
  const double numnodes = 16;

  //! local particle - local node
  double xpxi, ypyi;

  //! local shapefunctions for x and y coordinates
  double sxni, syni;

  //! shape function grad
  double dxni, dyni;

  for (unsigned n = 0; n < numnodes; n++) {
    xni = local_node(n, 0);
    yni = local_node(n, 1);
    xpxi = xi(0) - xni;  // local particle x - local node x
    ypyi = xi(1) - yni;  // local particle x - local node x

    //! x direction
    if (xpxi >= element_length + particle_length) {
      sxni = 0;
      dxni = 0;
    } else if (xpxi > -element_length - particle_length &&
               xpxi <= -element_length + particle_length) {
      sxni = ((element_length + particle_length + xpxi) *
              (element_length + particle_length + xpxi)) /
             4 * element_length * particle_length;
      dxni = (element_length + particle_length + xpxi) / 2 * element_length *
             particle_length;
    } else if (xpxi > -element_length + particle_length &&
               xpxi <= -particle_length) {
      sxni = 1 + (xpxi / element_length);
      dxni = 1 / element_length;
    } else if (xpxi > -particle_length && xpxi <= particle_length) {
      sxni = 1 - (((xpxi * xpxi) + (particle_length * particle_length)) / 2 *
                  element_length * particle_length);
      dxni = -(xpxi / element_length * particle_length);
    } else if (xpxi > particle_length &&
               xpxi <= element_length - particle_length) {
      sxni = 1 - (xpxi / element_length);
      dxni = -(1 / element_length);
    } else if (xpxi > element_length - particle_length &&
               xpxi <= element_length + particle_length) {
      sxni = ((element_length + particle_length + xpxi) *
              (element_length + particle_length + xpxi)) /
             4 * element_length * particle_length;
      dxni = -((element_length + particle_length + xpxi) / 2 * element_length *
               particle_length);
    }
    //! y direction
    if (ypyi >= element_length + particle_length) {
      syni = 0;
      dyni = 0;
    } else if (ypyi > -element_length - particle_length &&
               ypyi <= -element_length + particle_length) {
      syni = ((element_length + particle_length + ypyi) *
              (element_length + particle_length + ypyi)) /
             4 * element_length * particle_length;
      dyni = (element_length + particle_length + ypyi) / 2 * element_length *
             particle_length;
    } else if (ypyi > -element_length + particle_length &&
               ypyi <= -particle_length) {
      syni = 1 + (ypyi / element_length);
      dyni = 1 / element_length;
    } else if (ypyi > -particle_length && ypyi <= particle_length) {
      syni = 1 - (((ypyi * ypyi) + (particle_length * particle_length)) / 2 *
                  element_length * particle_length);
      dyni = -(ypyi / element_length * particle_length);
    } else if (ypyi > particle_length &&
               ypyi <= element_length - particle_length) {
      syni = 1 - (ypyi / element_length);
      dyni = -(1 / element_length);
    } else if (ypyi > element_length - particle_length &&
               ypyi <= element_length + particle_length) {
      syni = ((element_length + particle_length + ypyi) *
              (element_length + particle_length + ypyi)) /
             4 * element_length * particle_length;
      dyni = -((element_length + particle_length + ypyi) / 2 * element_length *
               particle_length);
    }
    grad_shapefn(1, n) = dxni * syni;
    grad_shapefn(0, n) = dyni * sxni;
  }
  return grad_shapefn;
}
//! Return the degree of element
//! 9-noded quadrilateral
template <>
inline mpm::ElementDegree mpm::GimpElement<2, 16>::degree() const {
  return mpm::ElementDegree::Quadratic;
}

///! Compute Jacobian
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::Matrix<double, Tdim, Tdim>
    mpm::GimpElement<Tdim, Tnfunctions>::jacobian(
        const Eigen::Matrix<double, 2, 1>& xi,
        const Eigen::MatrixXd& nodal_coordinates) const {
  // Get gradient shape functions
  const Eigen::MatrixXd grad_shapefn = this->grad_shapefn(xi);

  try {
    // Check if matrices dimensions are correct
    if ((grad_shapefn.rows() != nodal_coordinates.rows()) ||
        (xi.size() != nodal_coordinates.cols()))
      throw std::runtime_error(
          "Jacobian calculation: Incorrect dimension of xi and "
          "nodal_coordinates");
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }

  // Jacobian dx_i/dxi_j
  return (grad_shapefn.transpose() * nodal_coordinates);
}

//! Return the B-matrix of a Quadrilateral Element at a given local
//! coordinate
template <unsigned Tdim, unsigned Tnfunctions>
inline std::vector<Eigen::MatrixXd>
    mpm::GimpElement<Tdim, Tnfunctions>::bmatrix(const VectorDim& xi) const {
  // Get gradient shape functions
  Eigen::MatrixXd grad_shapefn = this->grad_shapefn(xi);

  // B-Matrix
  std::vector<Eigen::MatrixXd> bmatrix;
  bmatrix.reserve(Tnfunctions);

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

//! Return the B-matrix of a Quadrilateral Element at a given local
//! coordinate for a real cell
template <unsigned Tdim, unsigned Tnfunctions>
inline std::vector<Eigen::MatrixXd>
    mpm::GimpElement<Tdim, Tnfunctions>::bmatrix(
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
inline Eigen::MatrixXd mpm::GimpElement<Tdim, Tnfunctions>::mass_matrix(
    const std::vector<VectorDim>& xi_s) const {
  // Mass matrix
  Eigen::Matrix<double, Tnfunctions, Tnfunctions> mass_matrix;
  mass_matrix.setZero();
  for (const auto& xi : xi_s) {
    const Eigen::Matrix<double, Tnfunctions, 1> shape_fn = this->shapefn(xi);
    mass_matrix += (shape_fn * shape_fn.transpose());
  }
  return mass_matrix;
}

//! Return the laplace_matrix of a quadrilateral Element
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::MatrixXd mpm::GimpElement<Tdim, Tnfunctions>::laplace_matrix(
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
    const Eigen::MatrixXd grad_sf = this->grad_shapefn(xi);

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
inline Eigen::MatrixXi mpm::GimpElement<Tdim, Tnfunctions>::sides_indices()
    const {
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
inline Eigen::VectorXi mpm::GimpElement<Tdim, Tnfunctions>::corner_indices()
    const {
  Eigen::Matrix<int, 4, 1> indices;
  indices << 0, 1, 2, 3;
  return indices;
}

//! Return indices of a sub-tetrahedrons in a volume
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::MatrixXi mpm::GimpElement<Tdim, Tnfunctions>::inhedron_indices()
    const {
  Eigen::Matrix<int, 4, Tdim, Eigen::RowMajor> indices;

  // clang-format off
  indices << 0, 1,
             1, 2,
             2, 3,
             3, 0;
  //clang-format on
  return indices;
}

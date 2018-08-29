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
//!   |           |particle   |           |
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

// RENAME TO QUAD GIMP
//! Return shape functions of a 16-node Quadrilateral GIMP Element at a given
//! local coordinate
template <>
inline Eigen::VectorXd
    mpm::QuadrilateralGIMPElement::QuadrilateralElement<2, 16>::shapefn(
        const Eigen::Matrix<double, 2, 1>& xi,
        const unsigned& number_of_particles,
        const Eigen::Matrix<double, 2, 1>& deformation_gradient) const {

  Eigen::Matrix<double, 16, 1> shapefn;

  //! Matrix to store local node coordinates
  Eigen::Matrix<double, 16, 2> local_node = this->unit_cell_coordinates();

  //! length and volume of element in local coordinate (TODO: double check)
  const double element_length = 2;
  const double element_volume = 4;  // 2*2

  //! particle length = element length / number of particles in cell.
  const double particle_length =
      std::sqrt(element_volume / number_of_particles);

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

//! Return gradient of shape functions of a 16-node Quadrilateral Element at a
//! given local coordinate
template <>
inline Eigen::MatrixXd
    mpm::QuadrilateralGIMPElement::QuadrilateralElement<2, 16>::grad_shapefn(
        const Eigen::Matrix<double, 2, 1>& xi,
        const unsigned& number_of_particles,
        const Eigen::Matrix<double, 2, 1>& deformation_gradient) const {
  Eigen::Matrix<double, 16, 2> grad_shapefn;

  //! Matrix to store local node coordinates
  Eigen::Matrix<double, 16, 2> local_node = this->unit_cell_coordinates();

  //! length and volume of element in local coordinate (TODO: double check)
  const double element_length = 2;
  const double element_volume = 2 * 2;

  //! particle length = element length / number of particles in cell.
  const double particle_length = sqrt(element_volume / number_of_particles);

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
    ypyi = xi(1) - yni;  // local particle y - local node y

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

//! Compute Jacobian with particle size and deformation gradient
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::Matrix<double, Tdim, Tdim> mpm::QuadrilateralGIMPElement::
    QuadrilateralElement<Tdim, Tnfunctions>::jacobian(
        const VectorDim& xi, const Eigen::MatrixXd& nodal_coordinates,
        const unsigned& number_of_particles,
        const VectorDim& deformation_gradient) const {
  // Get gradient shape functions
  const Eigen::MatrixXd grad_shapefn =
      this->grad_shapefn(xi, number_of_particles, deformation_gradient);

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

//! Return the B-matrix of a Quadrilateral Element at a given local
//! coordinate for a real cell
template <unsigned Tdim, unsigned Tnfunctions>
inline std::vector<Eigen::MatrixXd> mpm::QuadrilateralGIMPElement::
    QuadrilateralElement<Tdim, Tnfunctions>::bmatrix(
        const VectorDim& xi, const Eigen::MatrixXd& nodal_coordinates,
        const unsigned& number_of_particles,
        const VectorDim& deformation_gradient) const {
  // Get gradient shape functions
  const Eigen::MatrixXd grad_sf =
      this->grad_shapefn(xi, number_of_particles, deformation_gradient);

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

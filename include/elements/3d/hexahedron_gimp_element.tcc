// Return natural nodal coordinates
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::MatrixXd mpm::HexahedronGIMPElement<
    Tdim, Tnfunctions>::natural_nodal_coordinates() const {
  //! Natural coordinates of nodes
  const Eigen::Matrix<double, Tnfunctions, Tdim> local_nodes =
      (Eigen::Matrix<double, Tnfunctions, Tdim>() << -1., 1., -1, 1., 1., -1,
       1., 1., 1, -1., 1., 1, -1., -1., -1, 1., -1., -1, 1., -1., 1, -1., -1.,
       1, -3., 3, -3, -1., 3, -3, 1., 3, -3, 3., 3, -3, -3., 3, -1, -1., 3, -1,
       1., 3, -1, 3., 3, -1, -3., 3, 1, -1., 3, 1, 1., 3, 1, 3., 3, 1, -3., 3,
       3, -1., 3, 3, 1., 3, 3, 3., 3, 3, -3., 1., -3, -1., 1., -3, 1., 1., -3,
       3., 1., -3, -3., 1., -1, 3., 1., -1, -3., 1., 1, 3., 1., 1, -3., 1., 3,
       -1., 1., 3, 1., 1., 3, 3., 1., 3, -3., -1., -3, -1., -1., -3, 1., -1.,
       -3, 3., -1., -3, -3., -1., -1, 3., -1., -1, -3., -1., 1, 3., -1., 1, -3.,
       -1., 3, -1., -1., 3, 1., -1., 3, 3., -1., 3, -3., -3., -3, -1., -3., -3,
       1., -3., -3, 3., -3., -3, -3., -3., -1, -1., -3., -1, 1., -3., -1, 3.,
       -3., -1, -3., -3., 1, -1., -3., 1, 1., -3., 1, 3., -3., 1, -3., -3., 3,
       -1., -3., 3, 1., -3., 3, 3., -3., 3)
          .finished();
  return local_nodes;
}

//! Return shape functions of a 64-node Hexahedron GIMP Element at a given
//! local coordinate
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::VectorXd mpm::HexahedronGIMPElement<Tdim, Tnfunctions>::shapefn(
    const Eigen::Matrix<double, Tdim, 1>& xi,
    const Eigen::Matrix<double, Tdim, 1>& particle_size,
    const Eigen::Matrix<double, Tdim, 1>& deformation_gradient) const {

  //! length of element in local coordinate
  const double element_length = 2.;
  //! Natural nodal coordinates
  const Eigen::Matrix<double, Tnfunctions, Tdim> local_nodes =
      this->natural_nodal_coordinates();
  //! To store shape functions
  Eigen::Matrix<double, Tnfunctions, 1> shapefn;

  try {
    //! loop to iterate over nodes
    for (unsigned n = 0; n < Tnfunctions; ++n) {
      //! local shape function in current plane (x, y or z)
      Eigen::Matrix<double, Tdim, 1> sni;
      //! loop to iterate over dimensions
      for (unsigned i = 0; i < Tdim; ++i) {
        double lp = particle_size(i) * 0.5;
        double ni = local_nodes(n, i);
        double npni = xi(i) - ni;  // local particle  - local node
        //! Conditional shape function statement see: Bardenhagen 2004
        if (npni <= (-element_length - lp)) {
          sni(i) = 0.;
        } else if ((-element_length - lp) < npni &&
                   npni <= (-element_length + lp)) {
          sni(i) = std::pow(element_length + lp + npni, 2.) /
                   (4. * (element_length * lp));
        } else if ((-element_length + lp) < npni && npni <= -lp) {
          sni(i) = 1. + (npni / element_length);
        } else if (-lp < npni && npni <= lp) {
          sni(i) =
              1. - (((npni * npni) + (lp * lp)) / (2. * element_length * lp));
        } else if (lp < npni && npni <= (element_length - lp)) {
          sni(i) = 1. - (npni / element_length);
        } else if ((element_length - lp) < npni &&
                   npni <= element_length + lp) {
          sni(i) = std::pow(element_length + lp - npni, 2.) /
                   (4. * element_length * lp);
        } else if ((element_length + lp) < npni) {
          sni(i) = 0.;
        } else {
          throw std::runtime_error(
              "GIMP shapefn: Point location outside area of influence");
        }
      }
      shapefn(n) =
          sni(0) * sni(1) * sni(2);  // See: Pruijn, N.S., 2016. Eq(4.30)
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    return shapefn;
  }
  return shapefn;
}

//! Return gradient of shape functions of a 64-node Hexahedron Element at a
//! given local coordinate
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::MatrixXd
    mpm::HexahedronGIMPElement<Tdim, Tnfunctions>::grad_shapefn(
        const Eigen::Matrix<double, Tdim, 1>& xi,
        const Eigen::Matrix<double, Tdim, 1>& particle_size,
        const Eigen::Matrix<double, Tdim, 1>& deformation_gradient) const {

  //! length of element in local coordinate
  const double element_length = 2.;
  //! Natural nodal coordinates
  const Eigen::Matrix<double, Tnfunctions, Tdim> local_nodes =
      this->natural_nodal_coordinates();
  //! To store grad shape functions
  Eigen::Matrix<double, Tnfunctions, Tdim> grad_shapefn;
  try {
    //! loop to iterate over nodes
    for (unsigned n = 0; n < Tnfunctions; ++n) {
      //! local shape function in current plane (x, y or z)
      Eigen::Matrix<double, Tdim, 1> sni;
      //! local grad shape function in current plane (x, y or z)
      Eigen::Matrix<double, Tdim, 1> dni;
      //! loop to iterate over dimensions
      for (unsigned i = 0; i < Tdim; ++i) {
        double lp = particle_size(i) * 0.5;
        double ni = local_nodes(n, i);
        double npni = xi(i) - ni;  // local particle  - local node
        //! Conditional shape function statement
        // see: Pruijn, N.S., 2016. Eq(4.30)
        if (npni <= (-element_length - lp)) {
          sni(i) = 0.;
          dni(i) = 0.;
        } else if ((-element_length - lp) < npni &&
                   npni <= (-element_length + lp)) {

          sni(i) = std::pow(element_length + lp + npni, 2.) /
                   (4. * (element_length * lp));
          dni(i) = (element_length + lp + npni) / (2. * element_length * lp);
        } else if ((-element_length + lp) < npni && npni <= -lp) {
          sni(i) = 1. + (npni / element_length);
          dni(i) = 1. / element_length;
        } else if (-lp < npni && npni <= lp) {
          sni(i) =
              1. - (((npni * npni) + (lp * lp)) / (2. * element_length * lp));
          dni(i) = -(npni / (element_length * lp));
        } else if (lp < npni && npni <= (element_length - lp)) {
          sni(i) = 1. - (npni / element_length);
          dni(i) = -(1. / element_length);
        } else if ((element_length - lp) < npni &&
                   npni <= (element_length + lp)) {
          sni(i) = std::pow(element_length + lp - npni, 2.) /
                   (4. * element_length * lp);
          dni(i) = -((element_length + lp - npni) / (2. * element_length * lp));
        } else if ((element_length + lp) < npni) {
          sni(i) = 0.;
          dni(i) = 0.;
        } else {
          throw std::runtime_error(
              "GIMP grad shapefn: Point location outside area of influence");
        }
      }
      // see: Pruijn, N.S., 2016. Eq(4.32)
      grad_shapefn(n, 0) = dni(0) * sni(1) * sni(2);
      grad_shapefn(n, 1) = sni(0) * dni(1) * sni(2);
      grad_shapefn(n, 2) = sni(0) * sni(1) * dni(2);
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    return grad_shapefn;
  }
  return grad_shapefn;
}

//! Return local shape functions of a GIMP Hexahedron Element at a given
//! Return local shape functions of a Hexahedron Element at a given local
//! coordinate, with particle size and deformation gradient
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::VectorXd
    mpm::HexahedronGIMPElement<Tdim, Tnfunctions>::shapefn_local(
        const Eigen::Matrix<double, Tdim, 1>& xi,
        const Eigen::Matrix<double, Tdim, 1>& particle_size,
        const Eigen::Matrix<double, Tdim, 1>& deformation_gradient) const {
  return mpm::HexahedronElement<Tdim, 8>::shapefn(xi, particle_size,
                                                  deformation_gradient);
}

//! Compute Jacobian
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::Matrix<double, Tdim, Tdim>
    mpm::HexahedronGIMPElement<Tdim, Tnfunctions>::jacobian(
        const Eigen::Matrix<double, 3, 1>& xi,
        const Eigen::MatrixXd& nodal_coordinates,
        const Eigen::Matrix<double, 3, 1>& particle_size,
        const Eigen::Matrix<double, 3, 1>& deformation_gradient) const {
  // Get gradient shape functions
  const Eigen::MatrixXd grad_shapefn =
      this->grad_shapefn(xi, particle_size, deformation_gradient);
  try {
    // Check if dimensions are correct
    if ((grad_shapefn.rows() != nodal_coordinates.rows()) ||
        (xi.size() != nodal_coordinates.cols()))
      throw std::runtime_error(
          "Jacobian calculation: Incorrect dimension of xi and "
          "nodal_coordinates");
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    return Eigen::Matrix<double, Tdim, Tdim>::Zero();
  }

  // Jacobian
  return (grad_shapefn.transpose() * nodal_coordinates);
}
//! Compute Jacobian local with particle size and deformation gradient
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::Matrix<double, Tdim, Tdim>
    mpm::HexahedronGIMPElement<Tdim, Tnfunctions>::jacobian_local(
        const Eigen::Matrix<double, 3, 1>& xi,
        const Eigen::MatrixXd& nodal_coordinates,
        const Eigen::Matrix<double, 3, 1>& particle_size,
        const Eigen::Matrix<double, 3, 1>& deformation_gradient) const {
  // Jacobian dx_i/dxi_j
  return this->jacobian(xi, nodal_coordinates, particle_size,
                        deformation_gradient);
}
//! Compute Bmatrix
template <unsigned Tdim, unsigned Tnfunctions>
inline std::vector<Eigen::MatrixXd>
    mpm::HexahedronGIMPElement<Tdim, Tnfunctions>::bmatrix(
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
    // clang-format off
    Eigen::Matrix<double, 6, Tdim> bi;
    bi(0, 0) = grad_shapefn(i, 0); bi(0, 1) = 0.;                 bi(0, 2) = 0.;
    bi(1, 0) = 0.;                 bi(1, 1) = grad_shapefn(i, 1); bi(1, 2) = 0.;
    bi(2, 0) = 0.;                 bi(2, 1) = 0.;                 bi(2, 2) = grad_shapefn(i, 2);
    bi(3, 0) = grad_shapefn(i, 1); bi(3, 1) = grad_shapefn(i, 0); bi(3, 2) = 0.;
    bi(4, 0) = 0.;                 bi(4, 1) = grad_shapefn(i, 2); bi(4, 2) = grad_shapefn(i, 1);
    bi(5, 0) = grad_shapefn(i, 2); bi(5, 1) = 0.;                 bi(5, 2) = grad_shapefn(i, 0);
    // clang-format on
    bmatrix.push_back(bi);
  }
  return bmatrix;
}

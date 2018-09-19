// Return natural nodal coordinates
template <>
inline Eigen::MatrixXd
    mpm::QuadrilateralGIMPElement<2, 16>::natural_nodal_coordinates() const {
  //! Natural coordinates of nodes
  Eigen::Matrix<double, 16, 2> local_nodes;
  // clang-format off
  (local_nodes << -1., -1.,
                  1., -1.,
                  1.,  1.,
                 -1.,  1.,
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
                 -3., -1.).finished();
  // clang-format on
  return local_nodes;
}

//! Return shape functions of a 16-node Quadrilateral GIMP Element at a given
//! local coordinate
template <>
inline Eigen::VectorXd mpm::QuadrilateralGIMPElement<2, 16>::shapefn(
    const Eigen::Matrix<double, 2, 1>& xi, const VectorDim& particle_size,
    const Eigen::Matrix<double, 2, 1>& deformation_gradient) const {

  //! Dimension
  const unsigned Dim = 2;
  //! Nodes in GIMP function
  const unsigned Nfunctions = 16;

  //! length of element in local coordinate
  const double element_length = 2.;

  //! Natural nodal coordinates
  const Eigen::Matrix<double, 16, 2> local_nodes =
      this->natural_nodal_coordinates();
  //! To store shape functions
  Eigen::Matrix<double, 16, 1> shapefn;
  //! loop to iterate over nodes
  for (unsigned n = 0; n < Nfunctions; ++n) {
    //! local shape function in current plane (x, y or z)
    Eigen::Matrix<double, 2, 1> sni;
    //! loop to iterate over dimensions
    for (unsigned i = 0; i < Dim; ++i) {
      double ni = local_nodes(n, i);
      double npni = xi(i) - ni;  // local particle  - local node

      //! Conditional shape function statement see: Bardenhagen 2004
      if (npni <= (-element_length - particle_size(i))) {
        sni(i) = 0.;
      } else if ((-element_length - particle_size(i)) < npni &&
                 npni <= (-element_length + particle_size(i))) {
        sni(i) = std::pow(element_length + particle_size(i) + npni, 2.) /
                 (4. * (element_length * particle_size(i)));
      } else if ((-element_length + particle_size(i)) < npni &&
                 npni <= -particle_size(i)) {
        sni(i) = 1. + (npni / element_length);
      } else if (-particle_size(i) < npni && npni <= particle_size(i)) {
        sni(i) = 1. - (((npni * npni) + (particle_size(i) * particle_size(i))) /
                       (2. * element_length * particle_size(i)));
      } else if (particle_size(i) < npni &&
                 npni <= element_length - particle_size(i)) {
        sni(i) = 1. - (npni / element_length);
      } else if ((element_length - particle_size(i)) < npni &&
                 npni <= element_length + particle_size(i)) {
        sni(i) = std::pow(element_length + particle_size(i) - npni, 2.) /
                 (4. * element_length * particle_size(i));
      } else if ((element_length + particle_size(i)) < npni) {
        sni(i) = 0.;
      }
    }
    shapefn(n) = sni(0) * sni(1);
  }

  return shapefn;
}

//! Return gradient of shape functions of a 16-node Quadrilateral Element at a
//! given local coordinate
template <>
inline Eigen::MatrixXd mpm::QuadrilateralGIMPElement<2, 16>::grad_shapefn(
    const Eigen::Matrix<double, 2, 1>& xi, const VectorDim& particle_size,
    const Eigen::Matrix<double, 2, 1>& deformation_gradient) const {

  //! Dimension
  const unsigned Dim = 2;
  //! Nodes in GIMP function
  const unsigned Nfunctions = 16;

  //! Natural nodal coordinates
  const Eigen::Matrix<double, 16, 2> local_nodes =
      this->natural_nodal_coordinates();
  //! To store grad shape functions
  Eigen::Matrix<double, 16, 2> grad_shapefn;
  //! length of element in local coordinate
  const double element_length = 2.;
  //! loop to iterate over nodes
  for (unsigned n = 0; n < Nfunctions; ++n) {
    //! local shape function in current plane (x, y or z)
    Eigen::Matrix<double, 2, 1> sni;
    //! local grad shape function in current plane (x, y or z)
    Eigen::Matrix<double, 2, 1> dni;
    //! loop to iterate over dimensions
    for (unsigned i = 0; i < Dim; ++i) {
      double ni = local_nodes(n, i);
      double npni = xi(i) - ni;  // local particle  - local node
      //! Conditional shape function statement see: Bardenhagen 2004
      if (npni <= (-element_length - particle_size(i))) {
        sni(i) = 0.;
        dni(i) = 0.;

      } else if ((-element_length - particle_size(i)) < npni &&
                 npni <= (-element_length + particle_size(i))) {

        sni(i) = std::pow(element_length + particle_size(i) + npni, 2.) /
                 (4. * (element_length * particle_size(i)));
        dni(i) = (element_length + particle_size(i) + npni) /
                 (2. * element_length * particle_size(i));
      } else if ((-element_length + particle_size(i)) < npni &&
                 npni <= -particle_size(i)) {
        sni(i) = 1. + (npni / element_length);
        dni(i) = 1. / element_length;
      } else if (-particle_size(i) < npni && npni <= particle_size(i)) {
        sni(i) = 1. - (((npni * npni) + (particle_size(i) * particle_size(i))) /
                       (2. * element_length * particle_size(i)));
        dni(i) = -(npni / (element_length * particle_size(i)));
      } else if (particle_size(i) < npni &&
                 npni <= element_length - particle_size(i)) {
        sni(i) = 1. - (npni / element_length);
        dni(i) = -(1. / element_length);
      } else if ((element_length - particle_size(i)) < npni &&
                 npni <= element_length + particle_size(i)) {
        sni(i) = std::pow(element_length + particle_size(i) - npni, 2.) /
                 (4. * element_length * particle_size(i));
        dni(i) = -((element_length + particle_size(i) - npni) /
                   (2. * element_length * particle_size(i)));
      } else if ((element_length + particle_size(i)) < npni) {
        sni(i) = 0.;
        dni(i) = 0.;
      }
    }
    grad_shapefn(n, 0) = dni(0) * sni(1);
    grad_shapefn(n, 1) = dni(1) * sni(0);
  }
  return grad_shapefn;
}

//! Return the B-matrix of a Quadrilateral Element at a given local
//! coordinate for a real cell
template <unsigned Tdim, unsigned Tnfunctions>
inline std::vector<Eigen::MatrixXd>
    mpm::QuadrilateralGIMPElement<Tdim, Tnfunctions>::bmatrix(
        const VectorDim& xi, const Eigen::MatrixXd& nodal_coordinates,
        const VectorDim& particle_size,
        const VectorDim& deformation_gradient) const {
  return this->mpm::QuadrilateralElement<Tdim, Tnfunctions>::bmatrix(
      xi, nodal_coordinates);
}

//! Compute Jacobian with particle size and deformation gradient
template <unsigned Tdim, unsigned Tnfunctions>
inline Eigen::Matrix<double, Tdim, Tdim>
    mpm::QuadrilateralGIMPElement<Tdim, Tnfunctions>::jacobian(
        const VectorDim& xi, const Eigen::MatrixXd& nodal_coordinates,
        const VectorDim& particle_size,
        const VectorDim& deformation_gradient) const {
  return this->mpm::QuadrilateralElement<Tdim, Tnfunctions>::jacobian(
      xi, nodal_coordinates);
}

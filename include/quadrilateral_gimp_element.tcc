//!   13----------12----------11----------10
//!   |           |           |           |
//!   |           |           |           |
//!   |           |           |           |
//!   |        (-1, 1)      (1,1)         |
//!   14----------3-----------2-----------9
//!   |           |           |           |
//!   |           |particle   |           |
//!   |           | location  |           |
//!   |           |           |           |
//!   15----------0-----------1-----------8
//!   |        (-1,-1)      (1,-1)        |
//!   |           |           |           |
//!   |           |           |           |
//!   |           |           |           |
//!   4-----------5-----------6-----------7

template <>
inline Eigen::MatrixXd
    mpm::QuadrilateralGIMPElement<2, 16>::local_node_coordinates() const {
  // Coordinates of a unit cell
  //! Matrix to store local node coordinates
  Eigen::Matrix<double, 16, 2> local_node;
  // clang-format off
  local_node << -1., -1.,
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
  return local_node;
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

  //! Matrix to store local node coordinates
  const Eigen::Matrix<double, 16, 2> local_node =
      this->local_node_coordinates();
  //! To store shape functions
  Eigen::Matrix<double, 16, 1> shapefn;
  //! loop to iterate over nodes
  for (unsigned n = 0; n < Nfunctions; ++n) {
    //! local shape function in current plane (x, y or z)
    Eigen::Matrix<double, 2, 1> sni;
    //! loop to iterate over dimensions
    for (unsigned i = 0; i < Dim; ++i) {
      double ni = local_node(n, i);
      double npni = xi(i) - ni;  // local particle  - local node

      //! Conditional shape function statement see: Bardenhagen 2004
      if (npni >= element_length + particle_size(i)) {
        sni(i) = 0;
      } else if (npni > -element_length - particle_size(i) &&
                 npni <= -element_length + particle_size(i)) {
        sni(i) = ((element_length + particle_size(i) + npni) *
                  (element_length + particle_size(i) + npni)) /
                 4 * element_length * particle_size(i);
      } else if (npni > -element_length + particle_size(i) &&
                 npni <= -particle_size(i)) {
        sni(i) = 1 + (npni / element_length);
      } else if (npni > -particle_size(i) && npni <= particle_size(i)) {
        sni(i) = 1 - (((npni * npni) + (particle_size(i) * particle_size(i))) /
                      2 * element_length * particle_size(i));
      } else if (npni > particle_size(i) &&
                 npni <= element_length - particle_size(i)) {
        sni(i) = 1 - (npni / element_length);
      } else if (npni > element_length - particle_size(i) &&
                 npni <= element_length + particle_size(i)) {
        sni(i) = ((element_length + particle_size(i) + npni) *
                  (element_length + particle_size(i) + npni)) /
                 4 * element_length * particle_size(i);
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

  //! Matrix to store local node coordinates
  const Eigen::Matrix<double, 16, 2> local_node =
      this->local_node_coordinates();
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
      double ni = local_node(n, i);
      double npni = xi(i) - ni;  // local particle  - local node

      //! Conditional grad shape function statement see: Zhang 2016
      if (npni >= element_length + particle_size(i)) {
        sni(i) = 0;
        dni(i) = 0;
      } else if (npni > -element_length - particle_size(i) &&
                 npni <= -element_length + particle_size(i)) {
        sni(i) = ((element_length + particle_size(i) + npni) *
                  (element_length + particle_size(i) + npni)) /
                 4 * element_length * particle_size(i);
        dni(i) = (element_length + particle_size(i) + npni) / 2 *
                 element_length * particle_size(i);
      } else if (npni > -element_length + particle_size(i) &&
                 npni <= -particle_size(i)) {
        sni(i) = 1 + (npni / element_length);
        dni(i) = 1 / element_length;
      } else if (npni > -particle_size(i) && npni <= particle_size(i)) {
        sni(i) = 1 - (((npni * npni) + (particle_size(i) * particle_size(i))) /
                      2 * element_length * particle_size(i));
        dni(i) = -(npni / element_length * particle_size(i));
      } else if (npni > particle_size(i) &&
                 npni <= element_length - particle_size(i)) {
        sni(i) = 1 - (npni / element_length);
        dni(i) = -(1 / element_length);
      } else if (npni > element_length - particle_size(i) &&
                 npni <= element_length + particle_size(i)) {
        sni(i) = ((element_length + particle_size(i) + npni) *
                  (element_length + particle_size(i) + npni)) /
                 4 * element_length * particle_size(i);
        dni(i) = -((element_length + particle_size(i) + npni) / 2 *
                   element_length * particle_size(i));
      }
    }
    grad_shapefn(1, n) = dni(0) * sni(1);
    grad_shapefn(0, n) = dni(1) * sni(0);
  }
  return grad_shapefn;
}

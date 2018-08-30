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
    mpm::QuadrilateralGIMPElement<2, 16>::QuadrilateralElement<2, 4>::shapefn(
        const Eigen::Matrix<double, 2, 1>& xi, const unsigned& particle_size,
        const Eigen::Matrix<double, 2, 1>& deformation_gradient) const {

  Eigen::Matrix<double, 16, 1> shapefn;

  //! Matrix to store local node coordinates
  Eigen::Matrix<double, 16, 2> local_node;  //= this->unit_cell_coordinates();

  //! length and volume of element in local coordinate (TODO: double check)
  const double element_length = 2;
  const double element_volume = 4;  // 2*2

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
    if (xpxi >= element_length + particle_size) {
      sxni = 0;
    } else if (xpxi > -element_length - particle_size &&
               xpxi <= -element_length + particle_size) {
      sxni = ((element_length + particle_size + xpxi) *
              (element_length + particle_size + xpxi)) /
             4 * element_length * particle_size;
    } else if (xpxi > -element_length + particle_size &&
               xpxi <= -particle_size) {
      sxni = 1 + (xpxi / element_length);
    } else if (xpxi > -particle_size && xpxi <= particle_size) {
      sxni = 1 - (((xpxi * xpxi) + (particle_size * particle_size)) / 2 *
                  element_length * particle_size);
    } else if (xpxi > particle_size && xpxi <= element_length - particle_size) {
      sxni = 1 - (xpxi / element_length);
    } else if (xpxi > element_length - particle_size &&
               xpxi <= element_length + particle_size) {
      sxni = ((element_length + particle_size + xpxi) *
              (element_length + particle_size + xpxi)) /
             4 * element_length * particle_size;
    }
    if (ypyi >= element_length + particle_size) {
      syni = 0;
    } else if (ypyi > -element_length - particle_size &&
               ypyi <= -element_length + particle_size) {
      syni = ((element_length + particle_size + ypyi) *
              (element_length + particle_size + ypyi)) /
             4 * element_length * particle_size;
    } else if (ypyi > -element_length + particle_size &&
               ypyi <= -particle_size) {
      syni = 1 + (ypyi / element_length);
    } else if (ypyi > -particle_size && ypyi <= particle_size) {
      syni = 1 - (((ypyi * ypyi) + (particle_size * particle_size)) / 2 *
                  element_length * particle_size);
    } else if (ypyi > particle_size && ypyi <= element_length - particle_size) {
      syni = 1 - (ypyi / element_length);
    } else if (ypyi > element_length - particle_size &&
               ypyi <= element_length + particle_size) {
      syni = ((element_length + particle_size + ypyi) *
              (element_length + particle_size + ypyi)) /
             4 * element_length * particle_size;
    }
    shapefn(n) = sxni * syni;
  }

  return shapefn;
}

//! Return gradient of shape functions of a 16-node Quadrilateral Element at a
//! given local coordinate
template <>
inline Eigen::MatrixXd mpm::QuadrilateralGIMPElement<2, 16>::
    QuadrilateralElement<2, 4>::grad_shapefn(
        const Eigen::Matrix<double, 2, 1>& xi, const unsigned& particle_size,
        const Eigen::Matrix<double, 2, 1>& deformation_gradient) const {
  Eigen::Matrix<double, 16, 2> grad_shapefn;

  //! Matrix to store local node coordinates
  Eigen::Matrix<double, 16, 2> local_node;  // = TODO:
                                            // this->unit_cell_coordinates();

  //! length and volume of element in local coordinate (TODO: double check)
  const double element_length = 2;
  const double element_volume = 2 * 2;

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
    if (xpxi >= element_length + particle_size) {
      sxni = 0;
      dxni = 0;
    } else if (xpxi > -element_length - particle_size &&
               xpxi <= -element_length + particle_size) {
      sxni = ((element_length + particle_size + xpxi) *
              (element_length + particle_size + xpxi)) /
             4 * element_length * particle_size;
      dxni = (element_length + particle_size + xpxi) / 2 * element_length *
             particle_size;
    } else if (xpxi > -element_length + particle_size &&
               xpxi <= -particle_size) {
      sxni = 1 + (xpxi / element_length);
      dxni = 1 / element_length;
    } else if (xpxi > -particle_size && xpxi <= particle_size) {
      sxni = 1 - (((xpxi * xpxi) + (particle_size * particle_size)) / 2 *
                  element_length * particle_size);
      dxni = -(xpxi / element_length * particle_size);
    } else if (xpxi > particle_size && xpxi <= element_length - particle_size) {
      sxni = 1 - (xpxi / element_length);
      dxni = -(1 / element_length);
    } else if (xpxi > element_length - particle_size &&
               xpxi <= element_length + particle_size) {
      sxni = ((element_length + particle_size + xpxi) *
              (element_length + particle_size + xpxi)) /
             4 * element_length * particle_size;
      dxni = -((element_length + particle_size + xpxi) / 2 * element_length *
               particle_size);
    }
    //! y direction
    if (ypyi >= element_length + particle_size) {
      syni = 0;
      dyni = 0;
    } else if (ypyi > -element_length - particle_size &&
               ypyi <= -element_length + particle_size) {
      syni = ((element_length + particle_size + ypyi) *
              (element_length + particle_size + ypyi)) /
             4 * element_length * particle_size;
      dyni = (element_length + particle_size + ypyi) / 2 * element_length *
             particle_size;
    } else if (ypyi > -element_length + particle_size &&
               ypyi <= -particle_size) {
      syni = 1 + (ypyi / element_length);
      dyni = 1 / element_length;
    } else if (ypyi > -particle_size && ypyi <= particle_size) {
      syni = 1 - (((ypyi * ypyi) + (particle_size * particle_size)) / 2 *
                  element_length * particle_size);
      dyni = -(ypyi / element_length * particle_size);
    } else if (ypyi > particle_size && ypyi <= element_length - particle_size) {
      syni = 1 - (ypyi / element_length);
      dyni = -(1 / element_length);
    } else if (ypyi > element_length - particle_size &&
               ypyi <= element_length + particle_size) {
      syni = ((element_length + particle_size + ypyi) *
              (element_length + particle_size + ypyi)) /
             4 * element_length * particle_size;
      dyni = -((element_length + particle_size + ypyi) / 2 * element_length *
               particle_size);
    }
    grad_shapefn(1, n) = dxni * syni;
    grad_shapefn(0, n) = dyni * sxni;
  }
  return grad_shapefn;
}

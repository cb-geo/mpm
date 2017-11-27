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

// 27-node (Triquadratic) Hexahedron Element
//! Check with GMSH
//!          7           18             6
//!            0_ _ _ _ _ 0 _ _ _ _ _ 0
//!            /|                     /|
//!           / |                    / |
//!          /  |    25             /  |
//!      19 0   |     0        17  0 17|
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

//! Return shape function of nodes
//! \param[in] xi Coordinates of point of interest
//! \retval shapefn Shape function of a given cell
//! \tparam Tdim Dimension
template <unsigned Tdim>
inline Eigen::MatrixXd mpm::HexahedronShapeFn<Tdim>::shapefn(
    const mpm::HexahedronShapeFn<Tdim>::VectorDim& xi) {
  switch (this->nfunctions_) {
    case 8:
      // 8-noded
      shapefn_.resize(8, 1);
      shapefn_(0) = 0.125 * (1 - xi(0)) * (1 - xi(1)) * (1 - xi(2));
      shapefn_(1) = 0.125 * (1 + xi(0)) * (1 - xi(1)) * (1 - xi(2));
      shapefn_(2) = 0.125 * (1 + xi(0)) * (1 + xi(1)) * (1 - xi(2));
      shapefn_(3) = 0.125 * (1 - xi(0)) * (1 + xi(1)) * (1 - xi(2));
      shapefn_(4) = 0.125 * (1 - xi(0)) * (1 - xi(1)) * (1 + xi(2));
      shapefn_(5) = 0.125 * (1 + xi(0)) * (1 - xi(1)) * (1 + xi(2));
      shapefn_(6) = 0.125 * (1 + xi(0)) * (1 + xi(1)) * (1 + xi(2));
      shapefn_(7) = 0.125 * (1 - xi(0)) * (1 + xi(1)) * (1 + xi(2));
      break;
    case 20:
      // 20-noded
      shapefn_.resize(20, 1);
      shapefn_(0) = -0.125 * (1 - xi(0)) * (1 - xi(1)) * (1 - xi(2)) *
                    (2 + xi(0) + xi(1) + xi(2));
      shapefn_(1) = -0.125 * (1 + xi(0)) * (1 - xi(1)) * (1 - xi(2)) *
                    (2 - xi(0) + xi(1) + xi(2));
      shapefn_(2) = -0.125 * (1 + xi(0)) * (1 + xi(1)) * (1 - xi(2)) *
                    (2 - xi(0) - xi(1) + xi(2));
      shapefn_(3) = -0.125 * (1 - xi(0)) * (1 + xi(1)) * (1 - xi(2)) *
                    (2 + xi(0) - xi(1) + xi(2));
      shapefn_(4) = -0.125 * (1 - xi(0)) * (1 - xi(1)) * (1 + xi(2)) *
                    (2 + xi(0) + xi(1) - xi(2));
      shapefn_(5) = -0.125 * (1 + xi(0)) * (1 - xi(1)) * (1 + xi(2)) *
                    (2 - xi(0) + xi(1) - xi(2));
      shapefn_(6) = -0.125 * (1 + xi(0)) * (1 + xi(1)) * (1 + xi(2)) *
                    (2 - xi(0) - xi(1) - xi(2));
      shapefn_(7) = -0.125 * (1 - xi(0)) * (1 + xi(1)) * (1 + xi(2)) *
                    (2 + xi(0) - xi(1) - xi(2));

      shapefn_(8) = 0.25 * (1 - xi(0) * xi(0)) * (1 - xi(1)) * (1 - xi(2));
      shapefn_(11) = 0.25 * (1 - xi(1) * xi(1)) * (1 + xi(0)) * (1 - xi(2));
      shapefn_(13) = 0.25 * (1 - xi(0) * xi(0)) * (1 + xi(1)) * (1 - xi(2));
      shapefn_(9) = 0.25 * (1 - xi(1) * xi(1)) * (1 - xi(0)) * (1 - xi(2));
      shapefn_(10) = 0.25 * (1 - xi(2) * xi(2)) * (1 - xi(0)) * (1 - xi(1));
      shapefn_(12) = 0.25 * (1 - xi(2) * xi(2)) * (1 + xi(0)) * (1 - xi(1));
      shapefn_(14) = 0.25 * (1 - xi(2) * xi(2)) * (1 + xi(0)) * (1 + xi(1));
      shapefn_(15) = 0.25 * (1 - xi(2) * xi(2)) * (1 - xi(0)) * (1 + xi(1));
      shapefn_(16) = 0.25 * (1 - xi(0) * xi(0)) * (1 - xi(1)) * (1 + xi(2));
      shapefn_(18) = 0.25 * (1 - xi(1) * xi(1)) * (1 + xi(0)) * (1 + xi(2));
      shapefn_(19) = 0.25 * (1 - xi(0) * xi(0)) * (1 + xi(1)) * (1 + xi(2));
      shapefn_(17) = 0.25 * (1 - xi(1) * xi(1)) * (1 - xi(0)) * (1 + xi(2));
      break;
    default:
      // Throw error
      break;
  }
  return shapefn_;
}

//! Return gradient of shape functions of a cell
//! \param[in] xi Coordinates of point of interest
//! \retval grad_shapefn Gradient of shape function of a given cell
//! \tparam Tdim Dimension
template <unsigned Tdim>
inline Eigen::MatrixXd mpm::HexahedronShapeFn<Tdim>::grad_shapefn(
    const mpm::HexahedronShapeFn<Tdim>::VectorDim& xi) {
  Eigen::MatrixXd grad_shapefn;
  switch (this->nfunctions_) {
    case 8:
      // 8-noded
      grad_shapefn_.resize(8, 3);
      grad_shapefn_(0, 0) = -0.125 * (1 - xi(1)) * (1 - xi(2));
      grad_shapefn_(1, 0) = 0.125 * (1 - xi(1)) * (1 - xi(2));
      grad_shapefn_(2, 0) = 0.125 * (1 + xi(1)) * (1 - xi(2));
      grad_shapefn_(3, 0) = -0.125 * (1 + xi(1)) * (1 - xi(2));
      grad_shapefn_(4, 0) = -0.125 * (1 - xi(1)) * (1 + xi(2));
      grad_shapefn_(5, 0) = 0.125 * (1 - xi(1)) * (1 + xi(2));
      grad_shapefn_(6, 0) = 0.125 * (1 + xi(1)) * (1 + xi(2));
      grad_shapefn_(7, 0) = -0.125 * (1 + xi(1)) * (1 + xi(2));

      grad_shapefn_(0, 1) = -0.125 * (1 - xi(0)) * (1 - xi(2));
      grad_shapefn_(1, 1) = -0.125 * (1 + xi(0)) * (1 - xi(2));
      grad_shapefn_(2, 1) = 0.125 * (1 + xi(0)) * (1 - xi(2));
      grad_shapefn_(3, 1) = 0.125 * (1 - xi(0)) * (1 - xi(2));
      grad_shapefn_(4, 1) = -0.125 * (1 - xi(0)) * (1 + xi(2));
      grad_shapefn_(5, 1) = -0.125 * (1 + xi(0)) * (1 + xi(2));
      grad_shapefn_(6, 1) = 0.125 * (1 + xi(0)) * (1 + xi(2));
      grad_shapefn_(7, 1) = 0.125 * (1 - xi(0)) * (1 + xi(2));

      grad_shapefn_(0, 2) = -0.125 * (1 - xi(0)) * (1 - xi(1));
      grad_shapefn_(1, 2) = -0.125 * (1 + xi(0)) * (1 - xi(1));
      grad_shapefn_(2, 2) = -0.125 * (1 + xi(0)) * (1 + xi(1));
      grad_shapefn_(3, 2) = -0.125 * (1 - xi(0)) * (1 + xi(1));
      grad_shapefn_(4, 2) = 0.125 * (1 - xi(0)) * (1 - xi(1));
      grad_shapefn_(5, 2) = 0.125 * (1 + xi(0)) * (1 - xi(1));
      grad_shapefn_(6, 2) = 0.125 * (1 + xi(0)) * (1 + xi(1));
      grad_shapefn_(7, 2) = 0.125 * (1 - xi(0)) * (1 + xi(1));
      return grad_shapefn_;
      break;
    case 20:
      // 20-noded
      grad_shapefn_.resize(20, 3);

      grad_shapefn_(0, 0) =
          0.125 * (2 * xi(0) + xi(1) + xi(2) + 1) * (1 - xi(1)) * (1 - xi(2));
      grad_shapefn_(1, 0) =
          -0.125 * (-2 * xi(0) + xi(1) + xi(2) + 1) * (1 - xi(1)) * (1 - xi(2));
      grad_shapefn_(2, 0) =
          -0.125 * (-2 * xi(0) - xi(1) + xi(2) + 1) * (1 + xi(1)) * (1 - xi(2));
      grad_shapefn_(3, 0) =
          0.125 * (2 * xi(0) - xi(1) + xi(2) + 1) * (1 + xi(1)) * (1 - xi(2));
      grad_shapefn_(4, 0) =
          0.125 * (2 * xi(0) + xi(1) - xi(2) + 1) * (1 - xi(1)) * (1 + xi(2));
      grad_shapefn_(5, 0) =
          -0.125 * (-2 * xi(0) + xi(1) - xi(2) + 1) * (1 - xi(1)) * (1 + xi(2));
      grad_shapefn_(6, 0) =
          -0.125 * (-2 * xi(0) - xi(1) - xi(2) + 1) * (1 + xi(1)) * (1 + xi(2));
      grad_shapefn_(7, 0) =
          0.125 * (2 * xi(0) - xi(1) - xi(2) + 1) * (1 + xi(1)) * (1 + xi(2));
      grad_shapefn_(8, 0) = -0.5 * xi(0) * (1 - xi(1)) * (1 - xi(2));
      grad_shapefn_(11, 0) = 0.25 * (1 - xi(1) * xi(1)) * (1 - xi(2));
      grad_shapefn_(13, 0) = -0.5 * xi(0) * (1 + xi(1)) * (1 - xi(2));
      grad_shapefn_(9, 0) = -0.25 * (1 - xi(1) * xi(1)) * (1 - xi(2));
      grad_shapefn_(10, 0) = -0.25 * (1 - xi(2) * xi(2)) * (1 - xi(1));
      grad_shapefn_(12, 0) = 0.25 * (1 - xi(2) * xi(2)) * (1 - xi(1));
      grad_shapefn_(14, 0) = 0.25 * (1 - xi(2) * xi(2)) * (1 + xi(1));
      grad_shapefn_(15, 0) = -0.25 * (1 - xi(2) * xi(2)) * (1 + xi(1));
      grad_shapefn_(16, 0) = -0.5 * xi(0) * (1 - xi(1)) * (1 + xi(2));
      grad_shapefn_(18, 0) = 0.25 * (1 - xi(1) * xi(1)) * (1 + xi(2));
      grad_shapefn_(19, 0) = -0.5 * xi(0) * (1 + xi(1)) * (1 + xi(2));
      grad_shapefn_(17, 0) = -0.25 * (1 - xi(1) * xi(1)) * (1 + xi(2));

      grad_shapefn_(0, 1) =
          0.125 * (xi(0) + 2 * xi(1) + xi(2) + 1) * (1 - xi(0)) * (1 - xi(2));
      grad_shapefn_(1, 1) =
          0.125 * (-xi(0) + 2 * xi(1) + xi(2) + 1) * (1 + xi(0)) * (1 - xi(2));
      grad_shapefn_(2, 1) =
          -0.125 * (-xi(0) - 2 * xi(1) + xi(2) + 1) * (1 + xi(0)) * (1 - xi(2));
      grad_shapefn_(3, 1) =
          -0.125 * (xi(0) - 2 * xi(1) + xi(2) + 1) * (1 - xi(1)) * (1 - xi(2));
      grad_shapefn_(4, 1) =
          0.125 * (xi(0) + 2 * xi(1) - xi(2) + 1) * (1 - xi(0)) * (1 + xi(2));
      grad_shapefn_(5, 1) =
          0.125 * (-xi(0) + 2 * xi(1) - xi(2) + 1) * (1 + xi(0)) * (1 + xi(2));
      grad_shapefn_(6, 1) =
          -0.125 * (-xi(0) - 2 * xi(1) - xi(2) + 1) * (1 + xi(0)) * (1 + xi(2));
      grad_shapefn_(7, 1) =
          -0.125 * (xi(0) - 2 * xi(1) - xi(2) + 1) * (1 - xi(1)) * (1 + xi(2));
      grad_shapefn_(8, 1) = -0.25 * (1 - xi(0) * xi(0)) * (1 - xi(2));
      grad_shapefn_(11, 1) = -0.5 * xi(1) * (1 + xi(0)) * (1 - xi(2));
      grad_shapefn_(13, 1) = 0.25 * (1 - xi(0) * xi(0)) * (1 - xi(2));
      grad_shapefn_(9, 1) = -0.5 * xi(1) * (1 - xi(0)) * (1 - xi(2));
      grad_shapefn_(10, 1) = -0.25 * (1 - xi(2) * xi(2)) * (1 - xi(0));
      grad_shapefn_(12, 1) = -0.25 * (1 - xi(2) * xi(2)) * (1 + xi(0));
      grad_shapefn_(14, 1) = 0.25 * (1 - xi(2) * xi(2)) * (1 + xi(0));
      grad_shapefn_(15, 1) = 0.25 * (1 - xi(2) * xi(2)) * (1 - xi(0));
      grad_shapefn_(16, 1) = -0.25 * (1 - xi(0) * xi(0)) * (1 + xi(2));
      grad_shapefn_(18, 1) = -0.5 * xi(1) * (1 + xi(0)) * (1 + xi(2));
      grad_shapefn_(19, 1) = 0.25 * (1 - xi(0) * xi(0)) * (1 + xi(2));
      grad_shapefn_(17, 1) = -0.5 * xi(1) * (1 - xi(0)) * (1 + xi(2));

      grad_shapefn_(0, 2) =
          0.125 * (xi(0) + xi(1) + 2 * xi(2) + 1) * (1 - xi(0)) * (1 - xi(1));
      grad_shapefn_(1, 2) =
          0.125 * (-xi(0) + xi(1) + 2 * xi(2) + 1) * (1 + xi(0)) * (1 - xi(1));
      grad_shapefn_(2, 2) =
          0.125 * (-xi(0) - xi(1) + 2 * xi(2) + 1) * (1 + xi(0)) * (1 + xi(1));
      grad_shapefn_(3, 2) =
          0.125 * (xi(0) - xi(1) + 2 * xi(2) + 1) * (1 - xi(1)) * (1 + xi(1));
      grad_shapefn_(4, 2) =
          -0.125 * (xi(0) + xi(1) - 2 * xi(2) + 1) * (1 - xi(0)) * (1 - xi(1));
      grad_shapefn_(5, 2) =
          -0.125 * (-xi(0) + xi(1) - 2 * xi(2) + 1) * (1 + xi(0)) * (1 - xi(1));
      grad_shapefn_(6, 2) =
          -0.125 * (-xi(0) - xi(1) - 2 * xi(2) + 1) * (1 + xi(0)) * (1 + xi(1));
      grad_shapefn_(7, 2) =
          -0.125 * (xi(0) - xi(1) - 2 * xi(2) + 1) * (1 - xi(1)) * (1 + xi(1));
      grad_shapefn_(8, 2) = -0.25 * (1 - xi(0) * xi(0)) * (1 - xi(1));
      grad_shapefn_(11, 2) = -0.25 * (1 - xi(1) * xi(1)) * (1 + xi(0));
      grad_shapefn_(13, 2) = -0.25 * (1 - xi(0) * xi(0)) * (1 + xi(1));
      grad_shapefn_(9, 2) = -0.25 * (1 - xi(1) * xi(1)) * (1 - xi(0));
      grad_shapefn_(10, 2) = -0.5 * xi(2) * (1 - xi(0)) * (1 - xi(1));
      grad_shapefn_(12, 2) = -0.5 * xi(2) * (1 + xi(0)) * (1 - xi(1));
      grad_shapefn_(14, 2) = -0.5 * xi(2) * (1 + xi(0)) * (1 + xi(1));
      grad_shapefn_(15, 2) = -0.5 * xi(2) * (1 - xi(0)) * (1 + xi(1));
      grad_shapefn_(16, 2) = 0.25 * (1 - xi(0) * xi(0)) * (1 - xi(1));
      grad_shapefn_(18, 2) = 0.25 * (1 - xi(1) * xi(1)) * (1 + xi(0));
      grad_shapefn_(19, 2) = 0.25 * (1 - xi(0) * xi(0)) * (1 + xi(1));
      grad_shapefn_(17, 2) = 0.25 * (1 - xi(1) * xi(1)) * (1 - xi(0));
      break;
    default:
      // Throw error
      break;
  }
  return grad_shapefn_;
}

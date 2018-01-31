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
inline Eigen::VectorXd mpm::HexahedronShapeFn<3, 8>::shapefn(
    const Eigen::Matrix<double, 3, 1>& xi) {
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
inline Eigen::MatrixXd mpm::HexahedronShapeFn<3, 8>::grad_shapefn(
    const Eigen::Matrix<double, 3, 1>& xi) {
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
inline Eigen::VectorXd mpm::HexahedronShapeFn<3, 20>::shapefn(
    const Eigen::Matrix<double, 3, 1>& xi) {
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
inline Eigen::MatrixXd mpm::HexahedronShapeFn<3, 20>::grad_shapefn(
    const Eigen::Matrix<double, 3, 1>& xi) {
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
  grad_shapefn(11, 0) = 0.25 * (1 - xi(1) * xi(1)) * (1 - xi(2));
  grad_shapefn(13, 0) = -0.5 * xi(0) * (1 + xi(1)) * (1 - xi(2));
  grad_shapefn(9, 0) = -0.25 * (1 - xi(1) * xi(1)) * (1 - xi(2));
  grad_shapefn(10, 0) = -0.25 * (1 - xi(2) * xi(2)) * (1 - xi(1));
  grad_shapefn(12, 0) = 0.25 * (1 - xi(2) * xi(2)) * (1 - xi(1));
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
  grad_shapefn(11, 1) = -0.5 * xi(1) * (1 + xi(0)) * (1 - xi(2));
  grad_shapefn(13, 1) = 0.25 * (1 - xi(0) * xi(0)) * (1 - xi(2));
  grad_shapefn(9, 1) = -0.5 * xi(1) * (1 - xi(0)) * (1 - xi(2));
  grad_shapefn(10, 1) = -0.25 * (1 - xi(2) * xi(2)) * (1 - xi(0));
  grad_shapefn(12, 1) = -0.25 * (1 - xi(2) * xi(2)) * (1 + xi(0));
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
  grad_shapefn(11, 2) = -0.25 * (1 - xi(1) * xi(1)) * (1 + xi(0));
  grad_shapefn(13, 2) = -0.25 * (1 - xi(0) * xi(0)) * (1 + xi(1));
  grad_shapefn(9, 2) = -0.25 * (1 - xi(1) * xi(1)) * (1 - xi(0));
  grad_shapefn(10, 2) = -0.5 * xi(2) * (1 - xi(0)) * (1 - xi(1));
  grad_shapefn(12, 2) = -0.5 * xi(2) * (1 + xi(0)) * (1 - xi(1));
  grad_shapefn(14, 2) = -0.5 * xi(2) * (1 + xi(0)) * (1 + xi(1));
  grad_shapefn(15, 2) = -0.5 * xi(2) * (1 - xi(0)) * (1 + xi(1));
  grad_shapefn(16, 2) = 0.25 * (1 - xi(0) * xi(0)) * (1 - xi(1));
  grad_shapefn(18, 2) = 0.25 * (1 - xi(1) * xi(1)) * (1 + xi(0));
  grad_shapefn(19, 2) = 0.25 * (1 - xi(0) * xi(0)) * (1 + xi(1));
  grad_shapefn(17, 2) = 0.25 * (1 - xi(1) * xi(1)) * (1 - xi(0));
  return grad_shapefn;
}

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


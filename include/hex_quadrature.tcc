// Getting the quadratures for Tnquadratures = 1
template <>
inline Eigen::MatrixXd mpm::HexahedronQuadrature<3, 1>::quadratures() {
  Eigen::Matrix<double, 1, 3> quadratures;
  quadratures(0, 0) = 0.;
  quadratures(0, 1) = 0.;
  quadratures(0, 2) = 0.;

  return quadratures;
}

// Getting the weights for Tnquadratures = 1
template <>
inline Eigen::VectorXd mpm::HexahedronQuadrature<3, 1>::weights() {
  Eigen::VectorXd weights(1);
  weights(0) = 8.;

  return weights;
}

// Getting the quadratures for Tnquadratures = 8
template <>
inline Eigen::MatrixXd mpm::HexahedronQuadrature<3, 8>::quadratures() {
  Eigen::Matrix<double, 8, 3> quadratures;
  quadratures(0, 0) = -1. / std::sqrt(3.);
  quadratures(0, 1) = -1. / std::sqrt(3.);
  quadratures(0, 2) = -1. / std::sqrt(3.);

  quadratures(1, 0) = 1. / std::sqrt(3.);
  quadratures(1, 1) = -1. / std::sqrt(3.);
  quadratures(1, 2) = -1. / std::sqrt(3.);

  quadratures(2, 0) = 1. / std::sqrt(3.);
  quadratures(2, 1) = 1. / std::sqrt(3.);
  quadratures(2, 2) = -1. / std::sqrt(3.);

  quadratures(3, 0) = -1. / std::sqrt(3.);
  quadratures(3, 1) = 1. / std::sqrt(3.);
  quadratures(3, 2) = -1. / std::sqrt(3.);

  quadratures(4, 0) = -1. / std::sqrt(3.);
  quadratures(4, 1) = -1. / std::sqrt(3.);
  quadratures(4, 2) = 1. / std::sqrt(3.);

  quadratures(5, 0) = 1. / std::sqrt(3.);
  quadratures(5, 1) = -1. / std::sqrt(3.);
  quadratures(5, 2) = 1. / std::sqrt(3.);

  quadratures(6, 0) = 1. / std::sqrt(3.);
  quadratures(6, 1) = 1. / std::sqrt(3.);
  quadratures(6, 2) = 1. / std::sqrt(3.);

  quadratures(7, 0) = -1. / std::sqrt(3.);
  quadratures(7, 1) = 1. / std::sqrt(3.);
  quadratures(7, 2) = 1. / std::sqrt(3.);

  return quadratures;
}

// Getting the weights for Tnquadratures = 8
template <>
inline Eigen::VectorXd mpm::HexahedronQuadrature<3, 8>::weights() {
  Eigen::VectorXd weights(8);
  weights(0) = 1.;
  weights(1) = 1.;
  weights(2) = 1.;
  weights(3) = 1.;
  weights(4) = 1.;
  weights(5) = 1.;
  weights(6) = 1.;
  weights(7) = 1.;

  return weights;
}

// Getting the quadratures for Tnquadratures = 27
template <>
inline Eigen::MatrixXd mpm::HexahedronQuadrature<3, 27>::quadratures() {
  Eigen::Matrix<double, 27, 3> quadratures;
  const double qpoint_a = std::sqrt(3. / 5.);

  // Along xi(0) = -1 plane
  // Line xi(1) = -1
  quadratures(0, 0) = -qpoint_a;
  quadratures(0, 1) = -qpoint_a;
  quadratures(0, 2) = -qpoint_a;

  quadratures(1, 0) = -qpoint_a;
  quadratures(1, 1) = -qpoint_a;
  quadratures(1, 2) = 0.;

  quadratures(2, 0) = -qpoint_a;
  quadratures(2, 1) = -qpoint_a;
  quadratures(2, 2) = qpoint_a;

  // Line xi(1) = 0
  quadratures(3, 0) = -qpoint_a;
  quadratures(3, 1) = 0.;
  quadratures(3, 2) = -qpoint_a;

  quadratures(4, 0) = -qpoint_a;
  quadratures(4, 1) = 0.;
  quadratures(4, 2) = 0.;

  quadratures(5, 0) = -qpoint_a;
  quadratures(5, 1) = 0.;
  quadratures(5, 2) = qpoint_a;

  // Line xi(1) = +1
  quadratures(6, 0) = -qpoint_a;
  quadratures(6, 1) = qpoint_a;
  quadratures(6, 2) = -qpoint_a;

  quadratures(7, 0) = -qpoint_a;
  quadratures(7, 1) = qpoint_a;
  quadratures(7, 2) = 0.;

  quadratures(8, 0) = -qpoint_a;
  quadratures(8, 1) = qpoint_a;
  quadratures(8, 2) = qpoint_a;

  // Along xi(0) = 0 plane
  // Line xi(1) = -1
  quadratures(9, 0) = 0.;
  quadratures(9, 1) = -qpoint_a;
  quadratures(9, 2) = -qpoint_a;

  quadratures(10, 0) = 0.;
  quadratures(10, 1) = -qpoint_a;
  quadratures(10, 2) = 0.;

  quadratures(11, 0) = 0.;
  quadratures(11, 1) = -qpoint_a;
  quadratures(11, 2) = qpoint_a;

  // Line xi(1) = 0
  quadratures(12, 0) = 0.;
  quadratures(12, 1) = 0.;
  quadratures(12, 2) = -qpoint_a;

  quadratures(13, 0) = 0.;
  quadratures(13, 1) = 0.;
  quadratures(13, 2) = 0.;

  quadratures(14, 0) = 0.;
  quadratures(14, 1) = 0.;
  quadratures(14, 2) = qpoint_a;

  // Line xi(1) = +1
  quadratures(15, 0) = 0.;
  quadratures(15, 1) = qpoint_a;
  quadratures(15, 2) = -qpoint_a;

  quadratures(16, 0) = 0.;
  quadratures(16, 1) = qpoint_a;
  quadratures(16, 2) = 0.;

  quadratures(17, 0) = 0.;
  quadratures(17, 1) = qpoint_a;
  quadratures(17, 2) = qpoint_a;

  // Along xi(0) = +1 plane
  // Line xi(1) = -1
  quadratures(18, 0) = qpoint_a;
  quadratures(18, 1) = -qpoint_a;
  quadratures(18, 2) = -qpoint_a;

  quadratures(19, 0) = qpoint_a;
  quadratures(19, 1) = -qpoint_a;
  quadratures(19, 2) = 0.;

  quadratures(20, 0) = qpoint_a;
  quadratures(20, 1) = -qpoint_a;
  quadratures(20, 2) = qpoint_a;

  // Line xi(1) = 0
  quadratures(21, 0) = qpoint_a;
  quadratures(21, 1) = 0.;
  quadratures(21, 2) = -qpoint_a;

  quadratures(22, 0) = qpoint_a;
  quadratures(22, 1) = 0.;
  quadratures(22, 2) = 0.;

  quadratures(23, 0) = qpoint_a;
  quadratures(23, 1) = 0.;
  quadratures(23, 2) = qpoint_a;

  // Line xi(1) = +1
  quadratures(24, 0) = qpoint_a;
  quadratures(24, 1) = qpoint_a;
  quadratures(24, 2) = -qpoint_a;

  quadratures(25, 0) = qpoint_a;
  quadratures(25, 1) = qpoint_a;
  quadratures(25, 2) = 0.;

  quadratures(26, 0) = qpoint_a;
  quadratures(26, 1) = qpoint_a;
  quadratures(26, 2) = qpoint_a;

  return quadratures;
}

// Getting the weights for Tnquadratures = 27
template <>
inline Eigen::VectorXd mpm::HexahedronQuadrature<3, 27>::weights() {
  Eigen::VectorXd weights(27);
  const double weight_b = 5. / 9.;
  const double weight_c = 8. / 9.;

  weights(0) = weight_b * weight_b * weight_b;
  weights(1) = weight_b * weight_b * weight_c;
  weights(2) = weight_b * weight_b * weight_b;
  weights(3) = weight_b * weight_c * weight_b;
  weights(4) = weight_b * weight_c * weight_c;
  weights(5) = weight_b * weight_c * weight_b;
  weights(6) = weight_b * weight_b * weight_b;
  weights(7) = weight_b * weight_b * weight_c;
  weights(8) = weight_b * weight_b * weight_b;

  weights(9) = weight_c * weight_b * weight_b;
  weights(10) = weight_c * weight_b * weight_c;
  weights(11) = weight_c * weight_b * weight_b;
  weights(12) = weight_c * weight_c * weight_b;
  weights(13) = weight_c * weight_c * weight_c;
  weights(14) = weight_c * weight_c * weight_b;
  weights(15) = weight_c * weight_b * weight_b;
  weights(16) = weight_c * weight_b * weight_c;
  weights(17) = weight_c * weight_b * weight_b;

  weights(18) = weight_b * weight_b * weight_b;
  weights(19) = weight_b * weight_b * weight_c;
  weights(20) = weight_b * weight_b * weight_b;
  weights(21) = weight_b * weight_c * weight_b;
  weights(22) = weight_b * weight_c * weight_c;
  weights(23) = weight_b * weight_c * weight_b;
  weights(24) = weight_b * weight_b * weight_b;
  weights(25) = weight_b * weight_b * weight_c;
  weights(26) = weight_b * weight_b * weight_b;

  return weights;
}
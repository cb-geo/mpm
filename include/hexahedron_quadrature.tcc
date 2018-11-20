// Getting the quadratures for Tnquadratures = 1
template <>
inline Eigen::MatrixXd mpm::HexahedronQuadrature<3, 1>::quadratures() const {
  Eigen::Matrix<double, 3, 1> quadratures;
  quadratures(0, 0) = 0.;
  quadratures(1, 0) = 0.;
  quadratures(2, 0) = 0.;

  return quadratures;
}

// Getting the weights for Tnquadratures = 1
template <>
inline Eigen::VectorXd mpm::HexahedronQuadrature<3, 1>::weights() const {
  Eigen::VectorXd weights(1);
  weights(0) = 8.;

  return weights;
}

// Getting the quadratures for Tnquadratures = 8
template <>
inline Eigen::MatrixXd mpm::HexahedronQuadrature<3, 8>::quadratures() const {
  Eigen::Matrix<double, 3, 8> quadratures;
  quadratures(0, 0) = -1. / std::sqrt(3.);
  quadratures(1, 0) = -1. / std::sqrt(3.);
  quadratures(2, 0) = -1. / std::sqrt(3.);

  quadratures(0, 1) = 1. / std::sqrt(3.);
  quadratures(1, 1) = -1. / std::sqrt(3.);
  quadratures(2, 1) = -1. / std::sqrt(3.);

  quadratures(0, 2) = 1. / std::sqrt(3.);
  quadratures(1, 2) = 1. / std::sqrt(3.);
  quadratures(2, 2) = -1. / std::sqrt(3.);

  quadratures(0, 3) = -1. / std::sqrt(3.);
  quadratures(1, 3) = 1. / std::sqrt(3.);
  quadratures(2, 3) = -1. / std::sqrt(3.);

  quadratures(0, 4) = -1. / std::sqrt(3.);
  quadratures(1, 4) = -1. / std::sqrt(3.);
  quadratures(2, 4) = 1. / std::sqrt(3.);

  quadratures(0, 5) = 1. / std::sqrt(3.);
  quadratures(1, 5) = -1. / std::sqrt(3.);
  quadratures(2, 5) = 1. / std::sqrt(3.);

  quadratures(0, 6) = 1. / std::sqrt(3.);
  quadratures(1, 6) = 1. / std::sqrt(3.);
  quadratures(2, 6) = 1. / std::sqrt(3.);

  quadratures(0, 7) = -1. / std::sqrt(3.);
  quadratures(1, 7) = 1. / std::sqrt(3.);
  quadratures(2, 7) = 1. / std::sqrt(3.);

  return quadratures;
}

// Getting the weights for Tnquadratures = 8
template <>
inline Eigen::VectorXd mpm::HexahedronQuadrature<3, 8>::weights() const {
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
inline Eigen::MatrixXd mpm::HexahedronQuadrature<3, 27>::quadratures() const {
  Eigen::Matrix<double, 3, 27> quadratures;
  const double qpoint_a = std::sqrt(3. / 5.);

  // Along xi(0) = -1 plane
  // Line xi(1) = -1
  quadratures(0, 0) = -qpoint_a;
  quadratures(1, 0) = -qpoint_a;
  quadratures(2, 0) = -qpoint_a;

  quadratures(0, 1) = -qpoint_a;
  quadratures(1, 1) = -qpoint_a;
  quadratures(2, 1) = 0.;

  quadratures(0, 2) = -qpoint_a;
  quadratures(1, 2) = -qpoint_a;
  quadratures(2, 2) = qpoint_a;

  // Line xi(1) = 0
  quadratures(0, 3) = -qpoint_a;
  quadratures(1, 3) = 0.;
  quadratures(2, 3) = -qpoint_a;

  quadratures(0, 4) = -qpoint_a;
  quadratures(1, 4) = 0.;
  quadratures(2, 4) = 0.;

  quadratures(0, 5) = -qpoint_a;
  quadratures(1, 5) = 0.;
  quadratures(2, 5) = qpoint_a;

  // Line xi(1) = +1
  quadratures(0, 6) = -qpoint_a;
  quadratures(1, 6) = qpoint_a;
  quadratures(2, 6) = -qpoint_a;

  quadratures(0, 7) = -qpoint_a;
  quadratures(1, 7) = qpoint_a;
  quadratures(2, 7) = 0.;

  quadratures(0, 8) = -qpoint_a;
  quadratures(1, 8) = qpoint_a;
  quadratures(2, 8) = qpoint_a;

  // Along xi(0) = 0 plane
  // Line xi(1) = -1
  quadratures(0, 9) = 0.;
  quadratures(1, 9) = -qpoint_a;
  quadratures(2, 9) = -qpoint_a;

  quadratures(0, 10) = 0.;
  quadratures(1, 10) = -qpoint_a;
  quadratures(2, 10) = 0.;

  quadratures(0, 11) = 0.;
  quadratures(1, 11) = -qpoint_a;
  quadratures(2, 11) = qpoint_a;

  // Line xi(1) = 0
  quadratures(0, 12) = 0.;
  quadratures(1, 12) = 0.;
  quadratures(2, 12) = -qpoint_a;

  quadratures(0, 13) = 0.;
  quadratures(1, 13) = 0.;
  quadratures(2, 13) = 0.;

  quadratures(0, 14) = 0.;
  quadratures(1, 14) = 0.;
  quadratures(2, 14) = qpoint_a;

  // Line xi(1) = +1
  quadratures(0, 15) = 0.;
  quadratures(1, 15) = qpoint_a;
  quadratures(2, 15) = -qpoint_a;

  quadratures(0, 16) = 0.;
  quadratures(1, 16) = qpoint_a;
  quadratures(2, 16) = 0.;

  quadratures(0, 17) = 0.;
  quadratures(1, 17) = qpoint_a;
  quadratures(2, 17) = qpoint_a;

  // Along xi(0) = +1 plane
  // Line xi(1) = -1
  quadratures(0, 18) = qpoint_a;
  quadratures(1, 18) = -qpoint_a;
  quadratures(2, 18) = -qpoint_a;

  quadratures(0, 19) = qpoint_a;
  quadratures(1, 19) = -qpoint_a;
  quadratures(2, 19) = 0.;

  quadratures(0, 20) = qpoint_a;
  quadratures(1, 20) = -qpoint_a;
  quadratures(2, 20) = qpoint_a;

  // Line xi(1) = 0
  quadratures(0, 21) = qpoint_a;
  quadratures(1, 21) = 0.;
  quadratures(2, 21) = -qpoint_a;

  quadratures(0, 22) = qpoint_a;
  quadratures(1, 22) = 0.;
  quadratures(2, 22) = 0.;

  quadratures(0, 23) = qpoint_a;
  quadratures(1, 23) = 0.;
  quadratures(2, 23) = qpoint_a;

  // Line xi(1) = +1
  quadratures(0, 24) = qpoint_a;
  quadratures(1, 24) = qpoint_a;
  quadratures(2, 24) = -qpoint_a;

  quadratures(0, 25) = qpoint_a;
  quadratures(1, 25) = qpoint_a;
  quadratures(2, 25) = 0.;

  quadratures(0, 26) = qpoint_a;
  quadratures(1, 26) = qpoint_a;
  quadratures(2, 26) = qpoint_a;

  return quadratures;
}

// Getting the weights for Tnquadratures = 27
template <>
inline Eigen::VectorXd mpm::HexahedronQuadrature<3, 27>::weights() const {
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

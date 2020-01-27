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

// Getting the quadratures for Tnquadratures = 64
template <>
inline Eigen::MatrixXd mpm::HexahedronQuadrature<3, 64>::quadratures() const {
  Eigen::Matrix<double, 3, 64> quadratures;

  const double val_sqrt_6by5 = std::sqrt(6. / 5.);
  const double val_sqrt_0_340 = std::sqrt(3. / 7. - 2. / 7. * val_sqrt_6by5);
  const double val_sqrt_0_861 = std::sqrt(3. / 7. + 2. / 7. * val_sqrt_6by5);

  // Along xi(2) = -val_sqrt_0_861 plane
  quadratures(0, 0) = -val_sqrt_0_861;
  quadratures(1, 0) = -val_sqrt_0_861;
  quadratures(2, 0) = -val_sqrt_0_861;

  quadratures(0, 1) = -val_sqrt_0_861;
  quadratures(1, 1) = -val_sqrt_0_340;
  quadratures(2, 1) = -val_sqrt_0_861;

  quadratures(0, 2) = -val_sqrt_0_861;
  quadratures(1, 2) = val_sqrt_0_340;
  quadratures(2, 2) = -val_sqrt_0_861;

  quadratures(0, 3) = -val_sqrt_0_861;
  quadratures(1, 3) = val_sqrt_0_861;
  quadratures(2, 3) = -val_sqrt_0_861;

  quadratures(0, 4) = -val_sqrt_0_340;
  quadratures(1, 4) = -val_sqrt_0_861;
  quadratures(2, 4) = -val_sqrt_0_861;

  quadratures(0, 5) = -val_sqrt_0_340;
  quadratures(1, 5) = -val_sqrt_0_340;
  quadratures(2, 5) = -val_sqrt_0_861;

  quadratures(0, 6) = -val_sqrt_0_340;
  quadratures(1, 6) = val_sqrt_0_340;
  quadratures(2, 6) = -val_sqrt_0_861;

  quadratures(0, 7) = -val_sqrt_0_340;
  quadratures(1, 7) = val_sqrt_0_861;
  quadratures(2, 7) = -val_sqrt_0_861;

  quadratures(0, 8) = val_sqrt_0_340;
  quadratures(1, 8) = -val_sqrt_0_861;
  quadratures(2, 8) = -val_sqrt_0_861;

  quadratures(0, 9) = val_sqrt_0_340;
  quadratures(1, 9) = -val_sqrt_0_340;
  quadratures(2, 9) = -val_sqrt_0_861;

  quadratures(0, 10) = val_sqrt_0_340;
  quadratures(1, 10) = val_sqrt_0_340;
  quadratures(2, 10) = -val_sqrt_0_861;

  quadratures(0, 11) = val_sqrt_0_340;
  quadratures(1, 11) = val_sqrt_0_861;
  quadratures(2, 11) = -val_sqrt_0_861;

  quadratures(0, 12) = val_sqrt_0_861;
  quadratures(1, 12) = -val_sqrt_0_861;
  quadratures(2, 12) = -val_sqrt_0_861;

  quadratures(0, 13) = val_sqrt_0_861;
  quadratures(1, 13) = -val_sqrt_0_340;
  quadratures(2, 13) = -val_sqrt_0_861;

  quadratures(0, 14) = val_sqrt_0_861;
  quadratures(1, 14) = val_sqrt_0_340;
  quadratures(2, 14) = -val_sqrt_0_861;

  quadratures(0, 15) = val_sqrt_0_861;
  quadratures(1, 15) = val_sqrt_0_861;
  quadratures(2, 15) = -val_sqrt_0_861;

  // Along xi(2) = -val_sqrt_0_340 plane
  quadratures(0, 16) = -val_sqrt_0_861;
  quadratures(1, 16) = -val_sqrt_0_861;
  quadratures(2, 16) = -val_sqrt_0_340;

  quadratures(0, 17) = -val_sqrt_0_861;
  quadratures(1, 17) = -val_sqrt_0_340;
  quadratures(2, 17) = -val_sqrt_0_340;

  quadratures(0, 18) = -val_sqrt_0_861;
  quadratures(1, 18) = val_sqrt_0_340;
  quadratures(2, 18) = -val_sqrt_0_340;

  quadratures(0, 19) = -val_sqrt_0_861;
  quadratures(1, 19) = val_sqrt_0_861;
  quadratures(2, 19) = -val_sqrt_0_340;

  quadratures(0, 20) = -val_sqrt_0_340;
  quadratures(1, 20) = -val_sqrt_0_861;
  quadratures(2, 20) = -val_sqrt_0_340;

  quadratures(0, 21) = -val_sqrt_0_340;
  quadratures(1, 21) = -val_sqrt_0_340;
  quadratures(2, 21) = -val_sqrt_0_340;

  quadratures(0, 22) = -val_sqrt_0_340;
  quadratures(1, 22) = val_sqrt_0_340;
  quadratures(2, 22) = -val_sqrt_0_340;

  quadratures(0, 23) = -val_sqrt_0_340;
  quadratures(1, 23) = val_sqrt_0_861;
  quadratures(2, 23) = -val_sqrt_0_340;

  quadratures(0, 24) = val_sqrt_0_340;
  quadratures(1, 24) = -val_sqrt_0_861;
  quadratures(2, 24) = -val_sqrt_0_340;

  quadratures(0, 25) = val_sqrt_0_340;
  quadratures(1, 25) = -val_sqrt_0_340;
  quadratures(2, 25) = -val_sqrt_0_340;

  quadratures(0, 26) = val_sqrt_0_340;
  quadratures(1, 26) = val_sqrt_0_340;
  quadratures(2, 26) = -val_sqrt_0_340;

  quadratures(0, 27) = val_sqrt_0_340;
  quadratures(1, 27) = val_sqrt_0_861;
  quadratures(2, 27) = -val_sqrt_0_340;

  quadratures(0, 28) = val_sqrt_0_861;
  quadratures(1, 28) = -val_sqrt_0_861;
  quadratures(2, 28) = -val_sqrt_0_340;

  quadratures(0, 29) = val_sqrt_0_861;
  quadratures(1, 29) = -val_sqrt_0_340;
  quadratures(2, 29) = -val_sqrt_0_340;

  quadratures(0, 30) = val_sqrt_0_861;
  quadratures(1, 30) = val_sqrt_0_340;
  quadratures(2, 30) = -val_sqrt_0_340;

  quadratures(0, 31) = val_sqrt_0_861;
  quadratures(1, 31) = val_sqrt_0_861;
  quadratures(2, 31) = -val_sqrt_0_340;

  // Along xi(2) = val_sqrt_0_340 plane
  quadratures(0, 32) = -val_sqrt_0_861;
  quadratures(1, 32) = -val_sqrt_0_861;
  quadratures(2, 32) = val_sqrt_0_340;

  quadratures(0, 33) = -val_sqrt_0_861;
  quadratures(1, 33) = -val_sqrt_0_340;
  quadratures(2, 33) = val_sqrt_0_340;

  quadratures(0, 34) = -val_sqrt_0_861;
  quadratures(1, 34) = val_sqrt_0_340;
  quadratures(2, 34) = val_sqrt_0_340;

  quadratures(0, 35) = -val_sqrt_0_861;
  quadratures(1, 35) = val_sqrt_0_861;
  quadratures(2, 35) = val_sqrt_0_340;

  quadratures(0, 36) = -val_sqrt_0_340;
  quadratures(1, 36) = -val_sqrt_0_861;
  quadratures(2, 36) = val_sqrt_0_340;

  quadratures(0, 37) = -val_sqrt_0_340;
  quadratures(1, 37) = -val_sqrt_0_340;
  quadratures(2, 37) = val_sqrt_0_340;

  quadratures(0, 38) = -val_sqrt_0_340;
  quadratures(1, 38) = val_sqrt_0_340;
  quadratures(2, 38) = val_sqrt_0_340;

  quadratures(0, 39) = -val_sqrt_0_340;
  quadratures(1, 39) = val_sqrt_0_861;
  quadratures(2, 39) = val_sqrt_0_340;

  quadratures(0, 40) = val_sqrt_0_340;
  quadratures(1, 40) = -val_sqrt_0_861;
  quadratures(2, 40) = val_sqrt_0_340;

  quadratures(0, 41) = val_sqrt_0_340;
  quadratures(1, 41) = -val_sqrt_0_340;
  quadratures(2, 41) = val_sqrt_0_340;

  quadratures(0, 42) = val_sqrt_0_340;
  quadratures(1, 42) = val_sqrt_0_340;
  quadratures(2, 42) = val_sqrt_0_340;

  quadratures(0, 43) = val_sqrt_0_340;
  quadratures(1, 43) = val_sqrt_0_861;
  quadratures(2, 43) = val_sqrt_0_340;

  quadratures(0, 44) = val_sqrt_0_861;
  quadratures(1, 44) = -val_sqrt_0_861;
  quadratures(2, 44) = val_sqrt_0_340;

  quadratures(0, 45) = val_sqrt_0_861;
  quadratures(1, 45) = -val_sqrt_0_340;
  quadratures(2, 45) = val_sqrt_0_340;

  quadratures(0, 46) = val_sqrt_0_861;
  quadratures(1, 46) = val_sqrt_0_340;
  quadratures(2, 46) = val_sqrt_0_340;

  quadratures(0, 47) = val_sqrt_0_861;
  quadratures(1, 47) = val_sqrt_0_861;
  quadratures(2, 47) = val_sqrt_0_340;

  // Along xi(2) = val_sqrt_0_861 plane
  quadratures(0, 48) = -val_sqrt_0_861;
  quadratures(1, 48) = -val_sqrt_0_861;
  quadratures(2, 48) = val_sqrt_0_861;

  quadratures(0, 49) = -val_sqrt_0_861;
  quadratures(1, 49) = -val_sqrt_0_340;
  quadratures(2, 49) = val_sqrt_0_861;

  quadratures(0, 50) = -val_sqrt_0_861;
  quadratures(1, 50) = val_sqrt_0_340;
  quadratures(2, 50) = val_sqrt_0_861;

  quadratures(0, 51) = -val_sqrt_0_861;
  quadratures(1, 51) = val_sqrt_0_861;
  quadratures(2, 51) = val_sqrt_0_861;

  quadratures(0, 52) = -val_sqrt_0_340;
  quadratures(1, 52) = -val_sqrt_0_861;
  quadratures(2, 52) = val_sqrt_0_861;

  quadratures(0, 53) = -val_sqrt_0_340;
  quadratures(1, 53) = -val_sqrt_0_340;
  quadratures(2, 53) = val_sqrt_0_861;

  quadratures(0, 54) = -val_sqrt_0_340;
  quadratures(1, 54) = val_sqrt_0_340;
  quadratures(2, 54) = val_sqrt_0_861;

  quadratures(0, 55) = -val_sqrt_0_340;
  quadratures(1, 55) = val_sqrt_0_861;
  quadratures(2, 55) = val_sqrt_0_861;

  quadratures(0, 56) = val_sqrt_0_340;
  quadratures(1, 56) = -val_sqrt_0_861;
  quadratures(2, 56) = val_sqrt_0_861;

  quadratures(0, 57) = val_sqrt_0_340;
  quadratures(1, 57) = -val_sqrt_0_340;
  quadratures(2, 57) = val_sqrt_0_861;

  quadratures(0, 58) = val_sqrt_0_340;
  quadratures(1, 58) = val_sqrt_0_340;
  quadratures(2, 58) = val_sqrt_0_861;

  quadratures(0, 59) = val_sqrt_0_340;
  quadratures(1, 59) = val_sqrt_0_861;
  quadratures(2, 59) = val_sqrt_0_861;

  quadratures(0, 60) = val_sqrt_0_861;
  quadratures(1, 60) = -val_sqrt_0_861;
  quadratures(2, 60) = val_sqrt_0_861;

  quadratures(0, 61) = val_sqrt_0_861;
  quadratures(1, 61) = -val_sqrt_0_340;
  quadratures(2, 61) = val_sqrt_0_861;

  quadratures(0, 62) = val_sqrt_0_861;
  quadratures(1, 62) = val_sqrt_0_340;
  quadratures(2, 62) = val_sqrt_0_861;

  quadratures(0, 63) = val_sqrt_0_861;
  quadratures(1, 63) = val_sqrt_0_861;
  quadratures(2, 63) = val_sqrt_0_861;

  return quadratures;
}

// Getting the weights for Tnquadratures = 64
template <>
inline Eigen::VectorXd mpm::HexahedronQuadrature<3, 64>::weights() const {
  Eigen::VectorXd weights(64);
  const double val_sqrt_30 = std::sqrt(30.);
  const double val_0_652 = (18 + val_sqrt_30) / 36.;  // Corresponds to 0.340
  const double val_0_348 = (18 - val_sqrt_30) / 36.;  // Corresponds to 0.861

  weights(0) = val_0_348 * val_0_348 * val_0_348;
  weights(1) = val_0_348 * val_0_348 * val_0_652;
  weights(2) = val_0_348 * val_0_348 * val_0_652;
  weights(3) = val_0_348 * val_0_348 * val_0_348;
  weights(4) = val_0_348 * val_0_652 * val_0_348;
  weights(5) = val_0_348 * val_0_652 * val_0_652;
  weights(6) = val_0_348 * val_0_652 * val_0_652;
  weights(7) = val_0_348 * val_0_652 * val_0_348;
  weights(8) = val_0_348 * val_0_652 * val_0_348;
  weights(9) = val_0_348 * val_0_652 * val_0_652;
  weights(10) = val_0_348 * val_0_652 * val_0_652;
  weights(11) = val_0_348 * val_0_652 * val_0_348;
  weights(12) = val_0_348 * val_0_348 * val_0_348;
  weights(13) = val_0_348 * val_0_652 * val_0_348;
  weights(14) = val_0_348 * val_0_652 * val_0_348;
  weights(15) = val_0_348 * val_0_348 * val_0_348;

  weights(16) = val_0_652 * val_0_348 * val_0_348;
  weights(17) = val_0_652 * val_0_348 * val_0_652;
  weights(18) = val_0_652 * val_0_348 * val_0_652;
  weights(19) = val_0_652 * val_0_348 * val_0_348;
  weights(20) = val_0_652 * val_0_652 * val_0_348;
  weights(21) = val_0_652 * val_0_652 * val_0_652;
  weights(22) = val_0_652 * val_0_652 * val_0_652;
  weights(23) = val_0_652 * val_0_652 * val_0_348;
  weights(24) = val_0_652 * val_0_652 * val_0_348;
  weights(25) = val_0_652 * val_0_652 * val_0_652;
  weights(26) = val_0_652 * val_0_652 * val_0_652;
  weights(27) = val_0_652 * val_0_652 * val_0_348;
  weights(28) = val_0_652 * val_0_348 * val_0_348;
  weights(29) = val_0_652 * val_0_652 * val_0_348;
  weights(30) = val_0_652 * val_0_652 * val_0_348;
  weights(31) = val_0_652 * val_0_348 * val_0_348;

  weights(32) = val_0_652 * val_0_348 * val_0_348;
  weights(33) = val_0_652 * val_0_348 * val_0_652;
  weights(34) = val_0_652 * val_0_348 * val_0_652;
  weights(35) = val_0_652 * val_0_348 * val_0_348;
  weights(36) = val_0_652 * val_0_652 * val_0_348;
  weights(37) = val_0_652 * val_0_652 * val_0_652;
  weights(38) = val_0_652 * val_0_652 * val_0_652;
  weights(39) = val_0_652 * val_0_652 * val_0_348;
  weights(40) = val_0_652 * val_0_652 * val_0_348;
  weights(41) = val_0_652 * val_0_652 * val_0_652;
  weights(42) = val_0_652 * val_0_652 * val_0_652;
  weights(43) = val_0_652 * val_0_652 * val_0_348;
  weights(44) = val_0_652 * val_0_348 * val_0_348;
  weights(45) = val_0_652 * val_0_652 * val_0_348;
  weights(46) = val_0_652 * val_0_652 * val_0_348;
  weights(47) = val_0_652 * val_0_348 * val_0_348;

  weights(48) = val_0_348 * val_0_348 * val_0_348;
  weights(49) = val_0_348 * val_0_348 * val_0_652;
  weights(50) = val_0_348 * val_0_348 * val_0_652;
  weights(51) = val_0_348 * val_0_348 * val_0_348;
  weights(52) = val_0_348 * val_0_652 * val_0_348;
  weights(53) = val_0_348 * val_0_652 * val_0_652;
  weights(54) = val_0_348 * val_0_652 * val_0_652;
  weights(55) = val_0_348 * val_0_652 * val_0_348;
  weights(56) = val_0_348 * val_0_652 * val_0_348;
  weights(57) = val_0_348 * val_0_652 * val_0_652;
  weights(58) = val_0_348 * val_0_652 * val_0_652;
  weights(59) = val_0_348 * val_0_652 * val_0_348;
  weights(60) = val_0_348 * val_0_348 * val_0_348;
  weights(61) = val_0_348 * val_0_652 * val_0_348;
  weights(62) = val_0_348 * val_0_652 * val_0_348;
  weights(63) = val_0_348 * val_0_348 * val_0_348;

  return weights;
}

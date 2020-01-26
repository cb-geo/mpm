// Getting the quadratures for Tnquadratures = 1
template <>
inline Eigen::MatrixXd mpm::QuadrilateralQuadrature<2, 1>::quadratures() const {
  Eigen::Matrix<double, 2, 1> quadratures;
  quadratures(0, 0) = 0.;
  quadratures(1, 0) = 0.;

  return quadratures;
}

// Getting the weights for Tnquadratures = 1
template <>
inline Eigen::VectorXd mpm::QuadrilateralQuadrature<2, 1>::weights() const {
  Eigen::VectorXd weights(1);
  weights(0) = 4.;

  return weights;
}

// Getting the quadratures for Tnquadratures = 4
template <>
inline Eigen::MatrixXd mpm::QuadrilateralQuadrature<2, 4>::quadratures() const {
  Eigen::Matrix<double, 2, 4> quadratures;
  const double val_1_by_sqrt3 = 1. / std::sqrt(3.);

  quadratures(0, 0) = -val_1_by_sqrt3;
  quadratures(1, 0) = -val_1_by_sqrt3;
  quadratures(0, 1) = val_1_by_sqrt3;
  quadratures(1, 1) = -val_1_by_sqrt3;
  quadratures(0, 2) = val_1_by_sqrt3;
  quadratures(1, 2) = val_1_by_sqrt3;
  quadratures(0, 3) = -val_1_by_sqrt3;
  quadratures(1, 3) = val_1_by_sqrt3;

  return quadratures;
}

// Getting the weights for Tnquadratures = 4
template <>
inline Eigen::VectorXd mpm::QuadrilateralQuadrature<2, 4>::weights() const {
  Eigen::VectorXd weights(4);
  weights(0) = 1.;
  weights(1) = 1.;
  weights(2) = 1.;
  weights(3) = 1.;

  return weights;
}

// Getting the quadratures for Tnquadratures = 9
template <>
inline Eigen::MatrixXd mpm::QuadrilateralQuadrature<2, 9>::quadratures() const {
  Eigen::Matrix<double, 2, 9> quadratures;
  const double val_sqrt_3by5 = std::sqrt(3. / 5.);

  quadratures(0, 0) = -val_sqrt_3by5;
  quadratures(1, 0) = -val_sqrt_3by5;

  quadratures(0, 1) = val_sqrt_3by5;
  quadratures(1, 1) = -val_sqrt_3by5;

  quadratures(0, 2) = val_sqrt_3by5;
  quadratures(1, 2) = val_sqrt_3by5;

  quadratures(0, 3) = -val_sqrt_3by5;
  quadratures(1, 3) = val_sqrt_3by5;

  quadratures(0, 4) = 0;
  quadratures(1, 4) = -val_sqrt_3by5;

  quadratures(0, 5) = val_sqrt_3by5;
  quadratures(1, 5) = 0;

  quadratures(0, 6) = 0;
  quadratures(1, 6) = val_sqrt_3by5;

  quadratures(0, 7) = -val_sqrt_3by5;
  quadratures(1, 7) = 0;

  quadratures(0, 8) = 0;
  quadratures(1, 8) = 0;

  return quadratures;
}

// Getting the weights for Tnquadratures = 9
template <>
inline Eigen::VectorXd mpm::QuadrilateralQuadrature<2, 9>::weights() const {
  Eigen::VectorXd weights(9);
  const double val_25_by_81 = 25. / 81.;
  const double val_40_by_81 = 40. / 81.;
  const double val_64_by_81 = 64. / 81.;

  weights(0) = val_25_by_81;
  weights(1) = val_25_by_81;
  weights(2) = val_25_by_81;
  weights(3) = val_25_by_81;
  weights(4) = val_40_by_81;
  weights(5) = val_40_by_81;
  weights(6) = val_40_by_81;
  weights(7) = val_40_by_81;
  weights(8) = val_64_by_81;

  return weights;
}

// Getting the quadratures for Tnquadratures = 16
template <>
inline Eigen::MatrixXd mpm::QuadrilateralQuadrature<2, 16>::quadratures() const {
  Eigen::Matrix<double, 2, 16> quadratures;
  const double val_sqrt_6by5 = std::sqrt(6. / 5.);
  const double val_sqrt_0_340 = std::sqrt(3. / 7. - 2. / 7. * val_sqrt_6by5);
  const double val_sqrt_0_861 = std::sqrt(3. / 7. + 2. / 7. * val_sqrt_6by5);

  quadratures(0, 0) = -val_sqrt_0_861;
  quadratures(1, 0) = -val_sqrt_0_861;

  quadratures(0, 1) = -val_sqrt_0_861;
  quadratures(1, 1) = -val_sqrt_0_340;

  quadratures(0, 2) = -val_sqrt_0_861;
  quadratures(1, 2) = val_sqrt_0_340;

  quadratures(0, 3) = -val_sqrt_0_861;
  quadratures(1, 3) = val_sqrt_0_861;

  quadratures(0, 4) = -val_sqrt_0_340;
  quadratures(1, 4) = -val_sqrt_0_861;

  quadratures(0, 5) = -val_sqrt_0_340;
  quadratures(1, 5) = -val_sqrt_0_340;

  quadratures(0, 6) = -val_sqrt_0_340;
  quadratures(1, 6) = val_sqrt_0_340;

  quadratures(0, 7) = -val_sqrt_0_340;
  quadratures(1, 7) = val_sqrt_0_861;

  quadratures(0, 8) = val_sqrt_0_340;
  quadratures(1, 8) = -val_sqrt_0_861;

  quadratures(0, 9) = val_sqrt_0_340;
  quadratures(1, 9) = -val_sqrt_0_340;

  quadratures(0, 10) = val_sqrt_0_340;
  quadratures(1, 10) = val_sqrt_0_340;

  quadratures(0, 11) = val_sqrt_0_340;
  quadratures(1, 11) = val_sqrt_0_861;

  quadratures(0, 12) = val_sqrt_0_861;
  quadratures(1, 12) = -val_sqrt_0_861;

  quadratures(0, 13) = val_sqrt_0_861;
  quadratures(1, 13) = -val_sqrt_0_340;

  quadratures(0, 14) = val_sqrt_0_861;
  quadratures(1, 14) = val_sqrt_0_340;

  quadratures(0, 15) = val_sqrt_0_861;
  quadratures(1, 15) = val_sqrt_0_861;

  return quadratures;
}

// Getting the weights for Tnquadratures = 16
template <>
inline Eigen::VectorXd mpm::QuadrilateralQuadrature<2, 16>::weights() const {
  Eigen::VectorXd weights(16);
  const double val_sqrt_30 = std::sqrt(30.);
  const double val_0_652 = (18 + val_sqrt_30) / 36.; // Corresponds to 0.340
  const double val_0_348 = (18 - val_sqrt_30) / 36.; // Corresponds to 0.861
  const double val_0_121 = val_0_348 * val_0_348;
  const double val_0_227 = val_0_348 * val_0_652;
  const double val_0_425 = val_0_652 * val_0_652;

  weights(0) = val_0_121;
  weights(1) = val_0_227;
  weights(2) = val_0_227;
  weights(3) = val_0_121;
  weights(4) = val_0_227;
  weights(5) = val_0_425;
  weights(6) = val_0_425;
  weights(7) = val_0_227;
  weights(8) = val_0_227;
  weights(9) = val_0_425;
  weights(10) = val_0_425;
  weights(11) = val_0_227;
  weights(12) = val_0_121;
  weights(13) = val_0_227;
  weights(14) = val_0_227;
  weights(15) = val_0_121;

  return weights;
}

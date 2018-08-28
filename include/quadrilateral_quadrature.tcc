// Getting the quadratures for Tnquadratures = 1
template <>
inline Eigen::MatrixXd mpm::QuadrilateralQuadrature<2, 1>::quadratures() {
  Eigen::Matrix<double, 1, 2> quadratures;
  quadratures(0, 0) = 0.;
  quadratures(0, 1) = 0.;

  return quadratures;
}

// Getting the weights for Tnquadratures = 1
template <>
inline Eigen::VectorXd mpm::QuadrilateralQuadrature<2, 1>::weights() {
  Eigen::VectorXd weights(1);
  weights(0) = 4.;

  return weights;
}

// Getting the quadratures for Tnquadratures = 4
template <>
inline Eigen::MatrixXd mpm::QuadrilateralQuadrature<2, 4>::quadratures() {
  Eigen::Matrix<double, 4, 2> quadratures;
  const double val_1_by_sqrt3 = 1. / std::sqrt(3.);

  quadratures(0, 0) = -val_1_by_sqrt3;
  quadratures(0, 1) = -val_1_by_sqrt3;
  quadratures(1, 0) = val_1_by_sqrt3;
  quadratures(1, 1) = -val_1_by_sqrt3;
  quadratures(2, 0) = val_1_by_sqrt3;
  quadratures(2, 1) = val_1_by_sqrt3;
  quadratures(3, 0) = -val_1_by_sqrt3;
  quadratures(3, 1) = val_1_by_sqrt3;

  return quadratures;
}

// Getting the weights for Tnquadratures = 4
template <>
inline Eigen::VectorXd mpm::QuadrilateralQuadrature<2, 4>::weights() {
  Eigen::VectorXd weights(4);
  weights(0) = 1.;
  weights(1) = 1.;
  weights(2) = 1.;
  weights(3) = 1.;

  return weights;
}

// Getting the quadratures for Tnquadratures = 9
template <>
inline Eigen::MatrixXd mpm::QuadrilateralQuadrature<2, 9>::quadratures() {
  Eigen::Matrix<double, 9, 2> quadratures;
  const double val_sqrt_3by5 = std::sqrt(3. / 5.);

  quadratures(0, 0) = -val_sqrt_3by5;
  quadratures(0, 1) = -val_sqrt_3by5;

  quadratures(1, 0) = val_sqrt_3by5;
  quadratures(1, 1) = -val_sqrt_3by5;

  quadratures(2, 0) = val_sqrt_3by5;
  quadratures(2, 1) = val_sqrt_3by5;

  quadratures(3, 0) = -val_sqrt_3by5;
  quadratures(3, 1) = val_sqrt_3by5;

  quadratures(4, 0) = 0;
  quadratures(4, 1) = -val_sqrt_3by5;

  quadratures(5, 0) = val_sqrt_3by5;
  quadratures(5, 1) = 0;

  quadratures(6, 0) = 0;
  quadratures(6, 1) = val_sqrt_3by5;

  quadratures(7, 0) = -val_sqrt_3by5;
  quadratures(7, 1) = 0;

  quadratures(8, 0) = 0;
  quadratures(8, 1) = 0;

  return quadratures;
}

// Getting the weights for Tnquadratures = 9
template <>
inline Eigen::VectorXd mpm::QuadrilateralQuadrature<2, 9>::weights() {
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
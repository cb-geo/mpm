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

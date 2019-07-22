// Getting the quadratures for Tnquadratures = 1
template <>
inline Eigen::MatrixXd mpm::TriangleQuadrature<2, 1>::quadratures() const {
  Eigen::Matrix<double, 2, 1> quadratures;
  quadratures(0, 0) = 1. / 3.;
  quadratures(1, 0) = 1. / 3.;

  return quadratures;
}

// Getting the weights for Tnquadratures = 1
template <>
inline Eigen::VectorXd mpm::QuadrilateralQuadrature<2, 1>::weights() const {
  Eigen::VectorXd weights(1);
  weights(0) = 0.5;

  return weights;
}

// Getting the quadratures for Tnquadratures = 3
template <>
inline Eigen::MatrixXd mpm::QuadrilateralQuadrature<2, 3>::quadratures() const {
  Eigen::Matrix<double, 2, 3> quadratures;
  const double val_1_by_sqrt3 = 1. / std::sqrt(3.);

  quadratures(0, 0) = 1. / 6.;
  quadratures(1, 0) = 1. / 6.;
  quadratures(0, 1) = 2. / 3.;
  quadratures(1, 1) = 1. / 6.;
  quadratures(0, 2) = 1. / 6.;
  quadratures(1, 2) = 2. / 3.;

  return quadratures;
}

// Getting the weights for Tnquadratures = 3
template <>
inline Eigen::VectorXd mpm::QuadrilateralQuadrature<2, 3>::weights() const {
  Eigen::VectorXd weights(3);
  weights(0) = 1. / 3.;
  weights(1) = 1. / 3.;
  weights(2) = 1. / 3.;

  return weights;
}

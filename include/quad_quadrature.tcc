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
    const double Val_1_by_sqrt3 = 1. / std::sqrt(3.);

    quadratures(0, 0) = -Val_1_by_sqrt3;
    quadratures(0, 1) = -Val_1_by_sqrt3;
    quadratures(1, 0) = Val_1_by_sqrt3;
    quadratures(1, 1) = -Val_1_by_sqrt3;
    quadratures(2, 0) = Val_1_by_sqrt3;
    quadratures(2, 1) = Val_1_by_sqrt3;
    quadratures(3, 0) = -Val_1_by_sqrt3;
    quadratures(3, 1) = Val_1_by_sqrt3;
    
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
    const double Val_sqrt_3by5 = std::sqrt(3. / 5.);
    
    quadratures(0, 0) = -Val_sqrt_3by5;
    quadratures(0, 1) = -Val_sqrt_3by5;

    quadratures(1, 0) = Val_sqrt_3by5;
    quadratures(1, 1) = -Val_sqrt_3by5;

    quadratures(2, 0) = Val_sqrt_3by5;
    quadratures(2, 1) = Val_sqrt_3by5;

    quadratures(3, 0) = -Val_sqrt_3by5;
    quadratures(3, 1) = Val_sqrt_3by5;

    quadratures(4, 0) = 0;
    quadratures(4, 1) = -Val_sqrt_3by5;

    quadratures(5, 0) = Val_sqrt_3by5;
    quadratures(5, 1) = 0;

    quadratures(6, 0) = 0;
    quadratures(6, 1) = Val_sqrt_3by5;

    quadratures(7, 0) = -Val_sqrt_3by5;
    quadratures(7, 1) = 0;

    quadratures(8, 0) = 0;
    quadratures(8, 1) = 0;
    
  return quadratures;
}

// Getting the weights for Tnquadratures = 9
template <>
inline Eigen::VectorXd mpm::QuadrilateralQuadrature<2, 9>::weights() {
  Eigen::VectorXd weights(9);
    const double Val_25_by_81 = 25. / 81.;
    const double Val_40_by_81 = 40. / 81.;
    const double Val_64_by_81 = 64. / 81.;

    weights(0) = Val_25_by_81;
    weights(1) = Val_25_by_81;
    weights(2) = Val_25_by_81;
    weights(3) = Val_25_by_81;
    weights(4) = Val_40_by_81;
    weights(5) = Val_40_by_81;
    weights(6) = Val_40_by_81;
    weights(7) = Val_40_by_81;
    weights(8) = Val_64_by_81;

  return weights;
}
// Getting the quadratures for Tnquadratures = 1
template <>
inline Eigen::MatrixXd mpm::QuadrilateralRandomQuadrature<2, 1>::quadratures() const {

  Eigen::Matrix<double, 2, 1> quadratures;

  for (unsigned i = 0; i < 1; ++i) {
    // Generate random quadratures in the interval (-0.99, +0.99)
    quadratures(0, i) = 1.98 * (double(std::rand()) / (double(RAND_MAX) + 1.0)) - 0.99;
    quadratures(1, i) = 1.98 * (double(std::rand()) / (double(RAND_MAX) + 1.0)) - 0.99;
  }

  return quadratures;
}

// Getting the quadratures for Tnquadratures = 4
template <>
inline Eigen::MatrixXd mpm::QuadrilateralRandomQuadrature<2, 4>::quadratures() const {

  Eigen::Matrix<double, 2, 4> quadratures;

  for (unsigned i = 0; i < 4; ++i) {
    // Generate random quadratures in the interval (-0.99, +0.99)
    quadratures(0, i) = 1.98 * (double(std::rand()) / (double(RAND_MAX) + 1.0)) - 0.99;
    quadratures(1, i) = 1.98 * (double(std::rand()) / (double(RAND_MAX) + 1.0)) - 0.99;
  }

  return quadratures;
}

// Getting the quadratures for Tnquadratures = 9
template <>
inline Eigen::MatrixXd mpm::QuadrilateralRandomQuadrature<2, 9>::quadratures() const {

  Eigen::Matrix<double, 2, 9> quadratures;

  for (unsigned i = 0; i < 9; ++i) {
    // Generate random quadratures in the interval (-0.99, +0.99)
    quadratures(0, i) = 1.98 * (double(std::rand()) / (double(RAND_MAX) + 1.0)) - 0.99;
    quadratures(1, i) = 1.98 * (double(std::rand()) / (double(RAND_MAX) + 1.0)) - 0.99;
  }

  return quadratures;
}

// Getting the quadratures for Tnquadratures = 16
template <>
inline Eigen::MatrixXd mpm::QuadrilateralRandomQuadrature<2, 16>::quadratures() const {

  Eigen::Matrix<double, 2, 16> quadratures;

  for (unsigned i = 0; i < 16; ++i) {
    // Generate random quadratures in the interval (-0.99, +0.99)
    quadratures(0, i) = 1.98 * (double(std::rand()) / (double(RAND_MAX) + 1.0)) - 0.99;
    quadratures(1, i) = 1.98 * (double(std::rand()) / (double(RAND_MAX) + 1.0)) - 0.99;
  }

  return quadratures;
}

// Getting the quadratures for Tnquadratures = 25
template <>
inline Eigen::MatrixXd mpm::QuadrilateralRandomQuadrature<2, 25>::quadratures() const {

  Eigen::Matrix<double, 2, 25> quadratures;

  for (unsigned i = 0; i < 25; ++i) {
    // Generate random quadratures in the interval (-0.99, +0.99)
    quadratures(0, i) = 1.98 * (double(std::rand()) / (double(RAND_MAX) + 1.0)) - 0.99;
    quadratures(1, i) = 1.98 * (double(std::rand()) / (double(RAND_MAX) + 1.0)) - 0.99;
  }

  return quadratures;
}

// Getting the weights for Tnquadratures = 1
template <>
inline Eigen::VectorXd mpm::QuadrilateralRandomQuadrature<2, 1>::weights() const {
 
  Eigen::Matrix<double, 1, 1> weights;

  for (unsigned i = 0; i < 1; ++i) {
    // Compute weights (TO DO)
    weights(i) = 1.;
  }

  return weights;
}

// Getting the weights for Tnquadratures = 4
template <>
inline Eigen::VectorXd mpm::QuadrilateralRandomQuadrature<2, 4>::weights() const {
 
  Eigen::Matrix<double, 4, 1> weights;

  for (unsigned i = 0; i < 4; ++i) {
    // Compute weights (TO DO)
    weights(i) = 1.;
  }

  return weights;
}

// Getting the weights for Tnquadratures = 9
template <>
inline Eigen::VectorXd mpm::QuadrilateralRandomQuadrature<2, 9>::weights() const {
 
  Eigen::Matrix<double, 9, 1> weights;

  for (unsigned i = 0; i < 9; ++i) {
    // Compute weights (TO DO)
    weights(i) = 1.;
  }

  return weights;
}

// Getting the weights for Tnquadratures = 16
template <>
inline Eigen::VectorXd mpm::QuadrilateralRandomQuadrature<2, 16>::weights() const {
 
  Eigen::Matrix<double, 16, 1> weights;

  for (unsigned i = 0; i < 16; ++i) {
    // Compute weights (TO DO)
    weights(i) = 1.;
  }

  return weights;
}

// Getting the weights for Tnquadratures = 25
template <>
inline Eigen::VectorXd mpm::QuadrilateralRandomQuadrature<2, 25>::weights() const {
 
  Eigen::Matrix<double, 25, 1> weights;

  for (unsigned i = 0; i < 25; ++i) {
    // Compute weights (TO DO)
    weights(i) = 1.;
  }

  return weights;
}
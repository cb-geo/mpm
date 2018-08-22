//! Compute inverse of 2d rotation matrix for orthogonal axis coordinate system
template <>
inline Eigen::MatrixXd mpm::Geometry<2>::compute_inverse_rotation_matrix(
    const Eigen::VectorXd& angles) {

  // Get cos and sin of angles
  const double cos_alpha = cos(angles(0));
  const double sin_alpha = sin(angles(0));
  const double cos_beta = cos(angles(1));
  const double sin_beta = sin(angles(1));

  Eigen::Matrix<double, 2, 2> rotation_matrix;

  // clang-format off
  rotation_matrix << cos_alpha*cos_beta - sin_alpha*sin_beta,  -cos_alpha*sin_beta - sin_alpha*cos_beta,
                     sin_alpha*cos_beta + cos_alpha*sin_beta,  -sin_alpha*sin_beta + cos_alpha*cos_beta;                                 
  // clang-format on            

  // Invert rotation matrix
  Eigen::Matrix<double, 2, 2> inverse_rotation_matrix;
  inverse_rotation_matrix = rotation_matrix.inverse();

  return inverse_rotation_matrix;
}

//! Compute inverse of 3d rotation matrix for orthogonal axis coordinate system
template <>
inline Eigen::MatrixXd mpm::Geometry<3>::compute_inverse_rotation_matrix(const 
    Eigen::VectorXd& angles) {

  // Get cos and sin of angles
  const double cos_alpha = cos(angles(0));
  const double sin_alpha = sin(angles(0));
  const double cos_beta = cos(angles(1));
  const double sin_beta = sin(angles(1));
  const double cos_gamma = cos(angles(2));
  const double sin_gamma = sin(angles(2));

  Eigen::Matrix<double, 3, 3> rotation_matrix;

  // clang-format off
  rotation_matrix << cos_alpha*cos_beta - sin_alpha*cos_gamma*sin_beta,  -cos_alpha*sin_beta - sin_alpha*cos_gamma*cos_beta,   sin_gamma*sin_alpha,
                     sin_alpha*cos_beta + cos_alpha*cos_gamma*sin_beta,  -sin_alpha*sin_beta + cos_alpha*cos_gamma*cos_beta,  -sin_gamma*cos_alpha,
                     sin_gamma*sin_beta,                                  sin_gamma*cos_beta,                                  cos_gamma;
  // clang-format on

  // Invert rotation matrix
  Eigen::Matrix<double, 3, 3> inverse_rotation_matrix;
  inverse_rotation_matrix = rotation_matrix.inverse();

  return inverse_rotation_matrix;
}
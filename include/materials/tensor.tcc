//! Compute mean stress p
inline double mpm::tensor::compute_mean_p(
    const Eigen::Matrix<double, 6, 1>& stress) {
  return (-1. / 3. * (stress(0) + stress(1) + stress(2)));
}

//! Compute deviatoric stress
inline Eigen::Matrix<double, 6, 1> mpm::tensor::compute_deviatoric_stress(
    const Eigen::Matrix<double, 6, 1>& stress) {

  // Compute mean_p in tension positive
  double mean_p = -1.0 * compute_mean_p(stress);

  // Compute deviatoric by subtracting volumetric part
  Eigen::Matrix<double, 6, 1> dev_stress = stress;
  for (unsigned i = 0; i < 3; ++i) dev_stress(i) -= mean_p;

  return dev_stress;
}

//! Compute J2 invariant
inline double mpm::tensor::compute_j2(
    const Eigen::Matrix<double, 6, 1>& stress) {

  double j2 = (std::pow((stress(0) - stress(1)), 2) +
               std::pow((stress(1) - stress(2)), 2) +
               std::pow((stress(0) - stress(2)), 2)) /
                  6.0 +
              std::pow(stress(3), 2) + std::pow(stress(4), 2) +
              std::pow(stress(5), 2);

  return j2;
}

//! Compute J3 invariant
inline double mpm::tensor::compute_j3(
    const Eigen::Matrix<double, 6, 1>& stress) {

  // Compute deviatoric stress
  Eigen::Matrix<double, 6, 1> dev_stress = compute_deviatoric_stress(stress);

  // Compute J3
  double j3 = (dev_stress(0) * dev_stress(1) * dev_stress(2)) -
              (dev_stress(2) * std::pow(dev_stress(3), 2)) +
              ((2 * dev_stress(3) * dev_stress(4) * dev_stress(5)) -
               (dev_stress(0) * std::pow(dev_stress(4), 2)) -
               (dev_stress(1) * std::pow(dev_stress(5), 2)));

  return j3;
}

//! Compute deviatoric q
inline double mpm::tensor::compute_deviatoric_q(
    const Eigen::Matrix<double, 6, 1>& stress) {

  // Compute J2 from
  double j2 = compute_j2(stress);

  double deviatoric_q = std::sqrt(3 * j2);

  return deviatoric_q;
}

//! Compute Lode angle
inline double mpm::tensor::compute_lode_angle(
    const Eigen::Matrix<double, 6, 1>& stress) {

  // Compute j2 and j3
  double j2 = compute_j2(stress);
  double j3 = compute_j3(stress);

  // Compute Lode angle value
  double lode_angle_val = 0.0;
  if (abs(j2) > 1.0E-6) {
    lode_angle_val = (3. * std::sqrt(3.) / 2.) * (j3 / std::pow(j2, 1.5));
  }
  if (lode_angle_val > 1.0) lode_angle_val = 1.0;
  if (lode_angle_val < -1.0) lode_angle_val = -1.0;

  // Compute Lode angle (sin convention)
  double lode_angle = (1. / 3.) * asin(lode_angle_val);
  if (lode_angle > M_PI / 6.) lode_angle = M_PI / 6.;
  if (lode_angle < -M_PI / 6.) lode_angle = -M_PI / 6.;
}
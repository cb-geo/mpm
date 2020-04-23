//! Compute mean stress p (tension positive)
inline double mpm::materials::p(const Eigen::Matrix<double, 6, 1>& stress) {
  // Compute and return mean p
  return (1. / 3. * (stress(0) + stress(1) + stress(2)));
}

//! Compute deviatoric stress
inline const Eigen::Matrix<double, 6, 1> mpm::materials::deviatoric_stress(
    const Eigen::Matrix<double, 6, 1>& stress) {

  // Compute mean p
  const double p = mpm::materials::p(stress);

  // Compute deviatoric by subtracting volumetric part
  Eigen::Matrix<double, 6, 1> deviatoric_stress = stress;
  for (unsigned i = 0; i < 3; ++i) deviatoric_stress(i) -= p;

  return deviatoric_stress;
}

//! Compute J2 invariant
inline double mpm::materials::j2(const Eigen::Matrix<double, 6, 1>& stress) {

  const double j2 = (std::pow((stress(0) - stress(1)), 2) +
                     std::pow((stress(1) - stress(2)), 2) +
                     std::pow((stress(0) - stress(2)), 2)) /
                        6.0 +
                    std::pow(stress(3), 2) + std::pow(stress(4), 2) +
                    std::pow(stress(5), 2);

  return j2;
}

//! Compute J3 invariant
inline double mpm::materials::j3(const Eigen::Matrix<double, 6, 1>& stress) {

  // Compute deviatoric stress
  Eigen::Matrix<double, 6, 1> deviatoric_stress =
      mpm::materials::deviatoric_stress(stress);

  // Compute J3
  const double j3 =
      (deviatoric_stress(0) * deviatoric_stress(1) * deviatoric_stress(2)) -
      (deviatoric_stress(2) * std::pow(deviatoric_stress(3), 2)) +
      ((2 * deviatoric_stress(3) * deviatoric_stress(4) *
        deviatoric_stress(5)) -
       (deviatoric_stress(0) * std::pow(deviatoric_stress(4), 2)) -
       (deviatoric_stress(1) * std::pow(deviatoric_stress(5), 2)));

  return j3;
}

//! Compute deviatoric q
inline double mpm::materials::q(const Eigen::Matrix<double, 6, 1>& stress) {

  // Compute J2 from
  const double j2 = mpm::materials::j2(stress);

  // Compute and return q
  return (std::sqrt(3 * j2));
}

//! Compute Lode angle
inline double mpm::materials::lode_angle(
    const Eigen::Matrix<double, 6, 1>& stress, const double tolerance) {

  // Compute j2 and j3
  const double j2 = mpm::materials::j2(stress);
  const double j3 = mpm::materials::j3(stress);

  // Compute Lode angle value
  double lode_angle_val = 0.0;
  if (std::abs(j2) > tolerance) {
    lode_angle_val = (3. * std::sqrt(3.) / 2.) * (j3 / std::pow(j2, 1.5));
  }
  if (lode_angle_val > 1.0) lode_angle_val = 1.0;
  if (lode_angle_val < -1.0) lode_angle_val = -1.0;

  // Compute and return Lode angle (cos convention, between 0 and pi/3)
  return (1. / 3.) * acos(lode_angle_val);
}

//! Compute derivative of p in terms of stress sigma
inline const Eigen::Matrix<double, 6, 1> mpm::materials::dp_dsigma(
    const Eigen::Matrix<double, 6, 1>& stress) {

  Eigen::Matrix<double, 6, 1> dp_dsigma = Eigen::Matrix<double, 6, 1>::Zero();
  dp_dsigma(0) = 1. / 3.;
  dp_dsigma(1) = 1. / 3.;
  dp_dsigma(2) = 1. / 3.;

  return dp_dsigma;
}

//! Compute derivative of q in terms of stress sigma
inline const Eigen::Matrix<double, 6, 1> mpm::materials::dq_dsigma(
    const Eigen::Matrix<double, 6, 1>& stress) {

  // Compute q
  const double q = mpm::materials::q(stress);

  // Compute deviatoric stress
  const auto deviatoric_stress = mpm::materials::deviatoric_stress(stress);

  // Compute dq / dsigma
  Eigen::Matrix<double, 6, 1> dq_dsigma = Eigen::Matrix<double, 6, 1>::Zero();
  if (std::abs(q) > std::numeric_limits<double>::epsilon()) {
    dq_dsigma(0) = 3. / (2. * q) * deviatoric_stress(0);
    dq_dsigma(1) = 3. / (2. * q) * deviatoric_stress(1);
    dq_dsigma(2) = 3. / (2. * q) * deviatoric_stress(2);
    dq_dsigma(3) = 3. / q * deviatoric_stress(3);
    dq_dsigma(4) = 3. / q * deviatoric_stress(4);
    dq_dsigma(5) = 3. / q * deviatoric_stress(5);
  }

  return dq_dsigma;
}

//! Compute derivative of J2 in terms of stress sigma
inline const Eigen::Matrix<double, 6, 1> mpm::materials::dj2_dsigma(
    const Eigen::Matrix<double, 6, 1>& stress) {

  // Compute deviatoric stress and put them as jd2 / dsigma
  auto dj2_dsigma = mpm::materials::deviatoric_stress(stress);
  dj2_dsigma(3) *= 2.0;
  dj2_dsigma(4) *= 2.0;
  dj2_dsigma(5) *= 2.0;

  return dj2_dsigma;
}

//! Compute derivative of J3 in terms of stress sigma
inline const Eigen::Matrix<double, 6, 1> mpm::materials::dj3_dsigma(
    const Eigen::Matrix<double, 6, 1>& stress) {

  // Compute J2
  const double j2 = mpm::materials::j2(stress);

  // Compute deviatoric stress
  const auto deviatoric_stress = mpm::materials::deviatoric_stress(stress);

  // Compute dj3 / dsigma
  Eigen::Matrix<double, 3, 1> dev1;
  dev1(0) = deviatoric_stress(0);
  dev1(1) = deviatoric_stress(3);
  dev1(2) = deviatoric_stress(5);
  Eigen::Matrix<double, 3, 1> dev2;
  dev2(0) = deviatoric_stress(3);
  dev2(1) = deviatoric_stress(1);
  dev2(2) = deviatoric_stress(4);
  Eigen::Matrix<double, 3, 1> dev3;
  dev3(0) = deviatoric_stress(5);
  dev3(1) = deviatoric_stress(4);
  dev3(2) = deviatoric_stress(2);

  Eigen::Matrix<double, 6, 1> dj3_dsigma = Eigen::Matrix<double, 6, 1>::Zero();
  dj3_dsigma(0) = dev1.dot(dev1) - (2. / 3.) * j2;
  dj3_dsigma(1) = dev2.dot(dev2) - (2. / 3.) * j2;
  dj3_dsigma(2) = dev3.dot(dev3) - (2. / 3.) * j2;
  dj3_dsigma(3) = 2.0 * dev1.dot(dev2);
  dj3_dsigma(4) = 2.0 * dev2.dot(dev3);
  dj3_dsigma(5) = 2.0 * dev1.dot(dev3);

  return dj3_dsigma;
}

//! Compute derivative of Lode angle theta in terms of stress sigma
inline const Eigen::Matrix<double, 6, 1> mpm::materials::dtheta_dsigma(
    const Eigen::Matrix<double, 6, 1>& stress, const double tolerance) {

  // Compute J2
  const double j2 = mpm::materials::j2(stress);

  // Compute J3
  const double j3 = mpm::materials::j3(stress);

  // Compute dj2_dsigma
  const auto dj2_dsigma = mpm::materials::dj2_dsigma(stress);

  // Compute dj3_dsigma
  const auto dj3_dsigma = mpm::materials::dj3_dsigma(stress);

  // Define derivatives of R in terms of J2 and J3
  double dr_dj2 = -9.0 / 4.0 * sqrt(3.0) * j3;
  double dr_dj3 = 3.0 / 2.0 * sqrt(3.0);

  // Compute derivative of theta in terms of R
  double dtheta_dr = -1.0 / 3.0;

  // Update when J2 is non zero
  if (std::abs(j2) > tolerance) {
    // Declare R defined as R = cos(3 theta)
    double r = j3 / 2.0 * std::pow(j2 / 3.0, -1.5);
    // Update derivatives of R
    dr_dj2 *= std::pow(j2, -2.5);
    dr_dj3 *= std::pow(j2, -1.5);
    // Update derivative of theta in terms of R, check for sqrt of zero
    const double factor =
        (std::abs(1 - r * r) < tolerance) ? tolerance : (1 - r * r);
    dtheta_dr = -1.0 / (3.0 * sqrt(factor));
  }

  // Compute and return dtheta / dsigma
  return (dtheta_dr * ((dr_dj2 * dj2_dsigma) + (dr_dj3 * dj3_dsigma)));
}

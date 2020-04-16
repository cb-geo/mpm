//! Compute mean stress p
inline const double mpm::material_utility::compute_mean_p(
    const Eigen::Matrix<double, 6, 1>& stress) {
  // Compute mean p in compresion positive
  return (-1. / 3. * (stress(0) + stress(1) + stress(2)));
}

//! Compute deviatoric stress
inline const Eigen::Matrix<double, 6, 1>
    mpm::material_utility::compute_deviatoric_stress(
        const Eigen::Matrix<double, 6, 1>& stress) {

  // Compute mean p in tension positive
  const double mean_p = -1.0 * compute_mean_p(stress);

  // Compute deviatoric by subtracting volumetric part
  Eigen::Matrix<double, 6, 1> dev_stress = stress;
  for (unsigned i = 0; i < 3; ++i) dev_stress(i) -= mean_p;

  return dev_stress;
}

//! Compute J2 invariant
inline const double mpm::material_utility::compute_j2(
    const Eigen::Matrix<double, 6, 1>& stress) {

  const double j2 = (std::pow((stress(0) - stress(1)), 2) +
                     std::pow((stress(1) - stress(2)), 2) +
                     std::pow((stress(0) - stress(2)), 2)) /
                        6.0 +
                    std::pow(stress(3), 2) + std::pow(stress(4), 2) +
                    std::pow(stress(5), 2);

  return j2;
}

//! Compute J3 invariant
inline const double mpm::material_utility::compute_j3(
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
inline const double mpm::material_utility::compute_deviatoric_q(
    const Eigen::Matrix<double, 6, 1>& stress) {

  // Compute J2 from
  const double j2 = compute_j2(stress);

  const double deviatoric_q = std::sqrt(3 * j2);

  return deviatoric_q;
}

//! Compute Lode angle
inline const double mpm::material_utility::compute_lode_angle(
    const Eigen::Matrix<double, 6, 1>& stress) {

  // Compute j2 and j3
  const double j2 = compute_j2(stress);
  const double j3 = compute_j3(stress);

  // Compute Lode angle value
  double lode_angle_val = 0.0;
  if (abs(j2) > 1.0E-6) {
    lode_angle_val = (3. * std::sqrt(3.) / 2.) * (j3 / std::pow(j2, 1.5));
  }
  if (lode_angle_val > 1.0) lode_angle_val = 1.0;
  if (lode_angle_val < -1.0) lode_angle_val = -1.0;

  // Compute Lode angle (sin convention)
  double lode_angle = (1. / 3.) * acos(lode_angle_val);
  if (lode_angle > M_PI / 3.) lode_angle = M_PI / 3.;
  if (lode_angle < -0.) lode_angle = 0.;

  return lode_angle;
}

//! Compute derivative of p in terms of stress sigma
inline const Eigen::Matrix<double, 6, 1>
    mpm::material_utility::compute_dp_dsigma(
        const Eigen::Matrix<double, 6, 1>& stress) {

  Eigen::Matrix<double, 6, 1> dp_dsigma = Eigen::Matrix<double, 6, 1>::Zero();
  dp_dsigma(0) = 1. / 3.;
  dp_dsigma(1) = 1. / 3.;
  dp_dsigma(2) = 1. / 3.;

  return dp_dsigma;
}

//! Compute derivative of q in terms of stress sigma
inline const Eigen::Matrix<double, 6, 1>
    mpm::material_utility::compute_dq_dsigma(
        const Eigen::Matrix<double, 6, 1>& stress) {

  // Compute q
  const double deviatoric_q = compute_deviatoric_q(stress);

  // Compute deviatoric stress
  Eigen::Matrix<double, 6, 1> dev_stress = compute_deviatoric_stress(stress);

  // Compute dq / dsigma
  Eigen::Matrix<double, 6, 1> dq_dsigma = Eigen::Matrix<double, 6, 1>::Zero();
  if (abs(deviatoric_q) > 1.E-6) {
    dq_dsigma(0) = 3. / 2. / deviatoric_q * dev_stress(0);
    dq_dsigma(1) = 3. / 2. / deviatoric_q * dev_stress(1);
    dq_dsigma(2) = 3. / 2. / deviatoric_q * dev_stress(2);
    dq_dsigma(3) = 3. / deviatoric_q * dev_stress(3);
    dq_dsigma(4) = 3. / deviatoric_q * dev_stress(4);
    dq_dsigma(5) = 3. / deviatoric_q * dev_stress(5);
  }

  return dq_dsigma;
}

//! Compute derivative of J2 in terms of stress sigma
inline const Eigen::Matrix<double, 6, 1>
    mpm::material_utility::compute_dj2_dsigma(
        const Eigen::Matrix<double, 6, 1>& stress) {

  // Compute deviatoric stress
  Eigen::Matrix<double, 6, 1> dev_stress = compute_deviatoric_stress(stress);

  // Compute dj2 / dsigma
  Eigen::Matrix<double, 6, 1> dj2_dsigma = dev_stress;
  dj2_dsigma(3) *= 2.0;
  dj2_dsigma(4) *= 2.0;
  dj2_dsigma(5) *= 2.0;

  return dj2_dsigma;
}

//! Compute derivative of J3 in terms of stress sigma
inline const Eigen::Matrix<double, 6, 1>
    mpm::material_utility::compute_dj3_dsigma(
        const Eigen::Matrix<double, 6, 1>& stress) {

  // Compute J2
  const double j2 = compute_j2(stress);

  // Compute deviatoric stress
  const Eigen::Matrix<double, 6, 1> dev_stress =
      compute_deviatoric_stress(stress);

  // Compute dj3 / dsigma
  Eigen::Matrix<double, 3, 1> dev1;
  dev1(0) = dev_stress(0);
  dev1(1) = dev_stress(3);
  dev1(2) = dev_stress(5);
  Eigen::Matrix<double, 3, 1> dev2;
  dev2(0) = dev_stress(3);
  dev2(1) = dev_stress(1);
  dev2(2) = dev_stress(4);
  Eigen::Matrix<double, 3, 1> dev3;
  dev3(0) = dev_stress(5);
  dev3(1) = dev_stress(4);
  dev3(2) = dev_stress(2);

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
inline const Eigen::Matrix<double, 6, 1>
    mpm::material_utility::compute_dtheta_dsigma(
        const Eigen::Matrix<double, 6, 1>& stress) {

  // Compute J2
  const double j2 = compute_j2(stress);

  // Compute J3
  const double j3 = compute_j3(stress);

  // Compute dj2_dsigma
  const Eigen::Matrix<double, 6, 1> dj2_dsigma = compute_dj2_dsigma(stress);

  // Compute dj3_dsigma
  const Eigen::Matrix<double, 6, 1> dj3_dsigma = compute_dj3_dsigma(stress);

  // Declare R as zero to avoid division by zero J2
  // R is defined as R = cos(3 theta)
  double R = 0.0;

  // Define derivatives of R in terms of J2 and J3
  double dR_dj2 = -9.0 / 4.0 * sqrt(3.0) * j3;
  double dR_dj3 = 3.0 / 2.0 * sqrt(3.0);

  // Compute derivative of theta in terms of R
  double dtheta_dR = -1.0 / 3.0;

  // Update when J2 is non zero
  if (abs(j2) > 1.0E-6) {
    // Update R
    R = j3 / 2.0 * std::pow(j2 / 3.0, -1.5);
    // Update derivatives of R
    dR_dj2 *= std::pow(j2, -2.5);
    dR_dj3 *= std::pow(j2, -1.5);
    // Update derivative of theta in terms of R, check for sqrt of zero
    if (abs(1 - R * R) < 1.0E-6) {
      dtheta_dR = -1.0 / 3.0 / sqrt(1.0E-6);
    } else {
      dtheta_dR = -1.0 / 3.0 / sqrt(1 - R * R);
    }
  }

  // Compute dtheta / dsigma
  const Eigen::Matrix<double, 6, 1> dtheta_dsigma =
      dtheta_dR * ((dR_dj2 * dj2_dsigma) + (dR_dj3 * dj3_dsigma));

  return dtheta_dsigma;
}

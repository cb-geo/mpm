//! Constructor with id and material properties
template <unsigned Tdim>
mpm::MohrCoulomb<Tdim>::MohrCoulomb(unsigned id,
                                    const Json& material_properties)
    : Material<Tdim>(id, material_properties) {
  try {
    density_ = material_properties["density"].template get<double>();
    youngs_modulus_ =
        material_properties["youngs_modulus"].template get<double>();
    poisson_ratio_ =
        material_properties["poisson_ratio"].template get<double>();
    friction_ = material_properties["friction"].template get<double>();
    dilation_ = material_properties["dilation"].template get<double>();
    cohesion_ = material_properties["cohesion"].template get<double>();
    residual_friction_ =
        material_properties["residual_friction"].template get<double>();
    residual_dilation_ =
        material_properties["residual_dilation"].template get<double>();
    residual_cohesion_ =
        material_properties["residual_cohesion"].template get<double>();
    peak_epds_ = material_properties["peak_epds"].template get<double>();
    crit_epds_ = material_properties["crit_epds"].template get<double>();
    tension_cutoff_ =
        material_properties["tension_cutoff"].template get<double>();
    porosity_ = material_properties["porosity"].template get<double>();
    properties_ = material_properties;
    // Calculate bulk modulus
    bulk_modulus_ = youngs_modulus_ / (3.0 * (1. - 2. * poisson_ratio_));
    shear_modulus_ = youngs_modulus_ / (2.0 * (1 + poisson_ratio_));

    // Set elastic tensor
    this->compute_elastic_tensor();

  } catch (std::exception& except) {
    console_->error("Material parameter not set: {}\n", except.what());
  }
}

//! Return elastic tensor
template <unsigned Tdim>
bool mpm::MohrCoulomb<Tdim>::compute_elastic_tensor() {
  // Shear modulus
  const double G = youngs_modulus_ / (2.0 * (1. + poisson_ratio_));

  const double a1 = bulk_modulus_ + (4.0 / 3.0) * G;
  const double a2 = bulk_modulus_ - (2.0 / 3.0) * G;

  // clang-format off
  // compute elasticityTensor
  de_(0,0)=a1;    de_(0,1)=a2;    de_(0,2)=a2;    de_(0,3)=0;    de_(0,4)=0;    de_(0,5)=0;
  de_(1,0)=a2;    de_(1,1)=a1;    de_(1,2)=a2;    de_(1,3)=0;    de_(1,4)=0;    de_(1,5)=0;
  de_(2,0)=a2;    de_(2,1)=a2;    de_(2,2)=a1;    de_(2,3)=0;    de_(2,4)=0;    de_(2,5)=0;
  de_(3,0)= 0;    de_(3,1)= 0;    de_(3,2)= 0;    de_(3,3)=G;    de_(3,4)=0;    de_(3,5)=0;
  de_(4,0)= 0;    de_(4,1)= 0;    de_(4,2)= 0;    de_(4,3)=0;    de_(4,4)=G;    de_(4,5)=0;
  de_(5,0)= 0;    de_(5,1)= 0;    de_(5,2)= 0;    de_(5,3)=0;    de_(5,4)=0;    de_(5,5)=G;
  // clang-format on
  return true;
}

//! Return j2, j3, rho and theta
template <unsigned Tdim>
void mpm::MohrCoulomb<Tdim>::compute_rho_theta(const Vector6d& stress) {

  const double ONETHIRDPI = 1.047197551;
  // Mean stress
  double mean_p = (stress(0) + stress(1) + stress(2)) / 3.;

  // Deviatoric stress components
  Vector6d dev_stress = Vector6d::Zero();
  dev_stress(0) = stress(0) - mean_p;
  dev_stress(1) = stress(1) - mean_p;
  dev_stress(2) = stress(2) - mean_p;
  dev_stress(3) = stress(3);
  if (Tdim == 3) {
    dev_stress(4) = stress(4);
    dev_stress(5) = stress(5);
  }

  // Second invariant J2
  j2_ = (pow((stress(0) - stress(1)), 2) + pow((stress(1) - stress(2)), 2) +
         pow((stress(0) - stress(2)), 2)) /
            6.0 +
        pow(stress(3), 2);
  if (Tdim == 3) j2_ += pow(stress(4), 2) + pow(stress(5), 2);

  // Third invariant J3
  j3_ = (dev_stress(0) * dev_stress(1) * dev_stress(2)) -
        (dev_stress(2) * pow(dev_stress(3), 2));
  if (Tdim == 3)
    j3_ += ((2 * dev_stress(3) * dev_stress(4) * dev_stress(5)) -
            (dev_stress(0) * pow(dev_stress(4), 2) -
             dev_stress(1) * pow(dev_stress(5), 2)));

  // Theta value
  double theta_val = 0.;
  if (fabs(j2_) > 0.0) theta_val = (3. * sqrt(3.) / 2.) * (j3_ / pow(j2_, 1.5));
  if (theta_val > 0.99) theta_val = 1.0;
  if (theta_val < -0.99) theta_val = -1.0;

  // Theta
  theta_ = (1. / 3.) * acos(theta_val);
  if (theta_ > ONETHIRDPI) theta_ = ONETHIRDPI;
  if (theta_ < 0.0) theta_ = 0.;

  rho_ = sqrt(2. * j2_);
}

//! Compute dF/dSigma and dP/dSigma
template <unsigned Tdim>
void mpm::MohrCoulomb<Tdim>::compute_df_dp(const Vector6d& stress,
                                           Vector6d& df_dsigma_,
                                           Vector6d& dp_dsigma_) {
  const double ONETHIRDPI = 1.047197551;

  // Mean stress
  double mean_p = (stress(0) + stress(1) + stress(2)) / 3.0;
  if (mean_p >= 0.0) mean_p = 1.0;

  // Deviatoric stress components
  Vector6d dev_stress = Vector6d::Zero();
  dev_stress(0) = stress(0) - mean_p;
  dev_stress(1) = stress(1) - mean_p;
  dev_stress(2) = stress(2) - mean_p;
  dev_stress(3) = stress(3);
  if (Tdim == 3) {
    dev_stress(4) = stress(4);
    dev_stress(5) = stress(5);
  }

  // dF / dEpsilon
  double df_depsilon = tan(phi_) / sqrt(3.);

  //  dF / dRho
  double df_drho =
      sqrt(3. / 2.) * ((sin(theta_ + ONETHIRDPI) / (sqrt(3.) * cos(phi_))) +
                       (cos(theta_ + ONETHIRDPI) * tan(phi_) / 3.));
  // dF / dTheta
  double df_dtheta = sqrt(3. / 2.) * rho_ *
                     ((cos(theta_ + ONETHIRDPI) / (sqrt(3.) * cos(phi_))) -
                      (sin(theta_ + ONETHIRDPI) * tan(phi_) / 3.));

  // compute dEpsilon / dSigma
  Vector6d depsilon_dsigma = Vector6d::Zero();
  depsilon_dsigma(0) = depsilon_dsigma(1) = depsilon_dsigma(2) = 1. / sqrt(3.);

  // compute dRho / dSigma
  Vector6d drho_dsigma = Vector6d::Zero();
  double multiplier = 1.;
  if (fabs(rho_) > 0.) multiplier = 1. / rho_;
  drho_dsigma = multiplier * dev_stress;
  if (Tdim == 2) drho_dsigma(4) = drho_dsigma(5) = 0.;

  // compute dTheta / dSigma
  Vector6d dtheta_dsigma = Vector6d::Zero();
  double r_val = 0.;
  if (fabs(j2_) > 1.E-22) r_val = (3. * sqrt(3.) / 2.) * (j3_ / pow(j2_, 1.5));

  // dTheta / dr
  double divider = 1 - (r_val * r_val);
  if (divider <= 0.) divider = 0.001;
  const double dtheta_dr = -1 / (3. * sqrt(divider));

  // dR/dJ2
  double dr_dj2 = (-9 * sqrt(3.) / 4.) * j3_;
  if (fabs(j2_) > 1.E-22) dr_dj2 = dr_dj2 / pow(j2_, 2.5);

  // dR/dJ3
  double dr_dj3 = 1.5 * sqrt(3.);
  if (fabs(j2_) > 1.E-22) dr_dj3 = dr_dj3 / pow(j2_, 1.5);

  Vector6d dj2_dsigma = dev_stress;
  Vector6d dj3_dsigma = Vector6d::Zero();

  Eigen::Matrix<double, 3, 1> dev1, dev2, dev3;
  dev1(0) = dev_stress(0);
  dev1(1) = dev_stress(3);
  dev1(2) = dev_stress(5);
  dev2(0) = dev_stress(3);
  dev2(1) = dev_stress(1);
  dev2(2) = dev_stress(4);
  dev3(0) = dev_stress(5);
  dev3(1) = dev_stress(4);
  dev3(2) = dev_stress(2);
  dj3_dsigma(0) = dev1.dot(dev1) - (2. / 3.) * j2_;
  dj3_dsigma(1) = dev2.dot(dev2) - (2. / 3.) * j2_;
  dj3_dsigma(2) = dev3.dot(dev3) - (2. / 3.) * j2_;
  dj3_dsigma(3) = dev1.dot(dev2);
  if (Tdim == 3) {
    dj3_dsigma(4) = dev2.dot(dev3);
    dj3_dsigma(5) = dev1.dot(dev3);
  }
  dtheta_dsigma = dtheta_dr * ((dr_dj2 * dj2_dsigma) + (dr_dj3 * dj3_dsigma));
  if (Tdim == 2) dtheta_dsigma(4) = dtheta_dsigma(5) = 0.;

  df_dsigma_ = (df_depsilon * depsilon_dsigma) + (df_drho * drho_dsigma) +
               (df_dtheta * dtheta_dsigma);
  if (Tdim == 2) df_dsigma_(4) = df_dsigma_(5) = 0.;

  // compute dP/dSigma
  const double r_mc = (3. - sin(phi_)) / (6 * cos(phi_));

  double e_val = (3. - sin(phi_)) / (3. + sin(phi_));
  if ((e_val - 0.5) < 0.) e_val = 0.501;
  if ((e_val - 1.) > 0.) e_val = 1.0;

  double sqpart = (4. * (1 - e_val * e_val) * pow(cos(theta_), 2)) +
                  (5 * e_val * e_val) - (4. * e_val);
  if (sqpart < 0.) sqpart = 0.00001;

  double r_mw_den = (2. * (1 - e_val * e_val) * cos(theta_)) +
                    ((2. * e_val - 1) * sqrt(sqpart));
  if (fabs(r_mw_den) < 1.E-22) r_mw_den = 0.001;

  const double r_mw_num = (4. * (1. - e_val * e_val) * pow(cos(theta_), 2)) +
                          pow((2. * e_val - 1.), 2);
  const double r_mw = (r_mw_num / r_mw_den) * r_mc;

  const double xi = 0.1;
  double omega =
      pow((xi * c_ * tan(psi_)), 2) + pow((r_mw * sqrt(3. / 2.) * rho_), 2);
  if (omega < 1.E-22) omega = 0.001;

  const double L = r_mw_num;
  const double M = r_mw_den;

  // dL/dTheta
  const double dl_dtheta =
      -8. * (1. - e_val * e_val) * cos(theta_) * sin(theta_);
  const double dm_dtheta = (-2. * (1. - e_val * e_val) * sin(theta_)) +
                           (0.5 * (2. * e_val - 1.) * dl_dtheta) / sqrt(sqpart);
  const double drmw_dtheta = ((M * dl_dtheta) - (L * dm_dtheta)) / (M * M);

  const double dp_depsilon = tan(psi_) / sqrt(3.);
  const double dp_drho = 3. * rho_ * r_mw * r_mw / (2. * sqrt(omega));
  const double dp_dtheta =
      (3. * rho_ * rho_ * r_mw * r_mc * drmw_dtheta) / (2. * sqrt(omega));

  dp_dsigma_ = (dp_depsilon * depsilon_dsigma) + (dp_drho * drho_dsigma) +
               (dp_dtheta * dtheta_dsigma);

  // compute softening part
  double dphi_dp_strain = 0.;
  double dc_dp_strain = 0.;
  // if(epds_ > epds_peak_ && epds_ < epds_crit_) {
  //   dphi_dp_strain = (phi_resd_ - phi_) / (epds_crit_ - epds_peak_);
  //   dc_dp_strain = (c_resd_ - c_) / (epds_crit_ - epds_peak_);
  // }
  // double  epsilon = (1. / sqrt(3.)) * (stress(0) + stress(1) + stress(2));
  // double df_dphi = sqrt(3./2.) * rho_ * ((sin(phi_) * sin(theta_ +
  // ONETHIRDPI) / (sqrt(3.)*cos(phi_) * cos(phi_))) + (cos(theta_
  // +ONETHIRDPI)/(3.*cos(phi_)*cos(phi_))) ) + (epsilon /
  // (sqrt(3.)*cos(phi_)*cos(phi_)));

  // double df_dc = -1.;
  // softening_ = (-1.) * ((dF_dPhi*dphi_dp_strain) + (df_dc*dc_dp_strain)) *
  // dp_drho;
}

//! Compute stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::MohrCoulomb<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr) {

  const double PI = std::atan(1.0) * 4.;
  const double ONETHIRDPI = PI / 3.;

  // Friction and dilation in radians
  const double phi_max = friction_ * PI / 180.;
  const double psi_max = dilation_ * PI / 180.;
  const double c_max = cohesion_;
  const double phi_min = residual_friction_ * PI / 180.;
  const double psi_min = residual_dilation_ * PI / 180.;
  const double c_min = residual_cohesion_;

  // Current MC parameters using a linear softening rule
  double epds = 0;
  if ((peak_epds_ - epds) >= 0.) {
    phi_ = phi_max;
    psi_ = psi_max;
    c_ = c_max;
  } else if ((epds - peak_epds_) > 0. && (crit_epds_ - epds) > 0.) {
    phi_ = phi_min + ((phi_max - phi_min) * (epds - crit_epds_) /
                      (peak_epds_ - crit_epds_));
    psi_ = psi_min + ((psi_max - psi_min) * (epds - crit_epds_) /
                      (peak_epds_ - crit_epds_));
    c_ = c_min +
         ((c_max - c_min) * (epds - crit_epds_) / (peak_epds_ - crit_epds_));
  } else if ((epds - crit_epds_) >= 0.) {
    phi_ = phi_max;
    psi_ = psi_max;
    c_ = c_max;
  }

  // Yield function for the current stress state
  this->compute_rho_theta(stress);
  double epsilon = (1. / sqrt(3.)) * (stress(0) + stress(1) + stress(2));
  double yield_func = sqrt(3. / 2.) * rho_ *
                          ((sin(theta_ + ONETHIRDPI) / (sqrt(3.) * cos(phi_))) +
                           (cos(theta_ + ONETHIRDPI) * tan(phi_) / 3.)) +
                      (epsilon / 3.) * tan(phi_) - c_;
  bool yield_state;
  if (yield_func > 1.E-22)
    yield_state = true;
  else
    yield_state = false;

  // compute plastic multiplier from the current stress state
  Vector6d dF_dSigma, dP_dSigma;
  this->softening_ = 0;
  this->compute_df_dp(stress, dF_dSigma, dP_dSigma);
  double lambda = dF_dSigma.dot(this->de_ * dstrain) /
                  ((dF_dSigma.dot(this->de_ * dP_dSigma)) + softening_);
  if (!yield_state) lambda = 0.;

  // compute the trial stress
  // sigma_trial = sigma + De * dstrain
  j2_ = j3_ = rho_ = theta_ = 0.;
  Vector6d trial_stress = stress + (this->de_ * dstrain);
  this->compute_rho_theta(trial_stress);
  epsilon =
      (1. / sqrt(3.)) * (trial_stress(0) + trial_stress(1) + trial_stress(2));
  double yield_func_trial =
      sqrt(3. / 2.) * rho_ *
          ((sin(theta_ + ONETHIRDPI) / (sqrt(3.) * cos(phi_))) +
           (cos(theta_ + ONETHIRDPI) * tan(phi_) / 3.)) +
      (epsilon / 3.) * tan(phi_) - c_;

  bool yield_state_trial;
  if (yield_func_trial > 1.E-22)
    yield_state_trial = true;
  else
    yield_state_trial = false;

  Vector6d dF_dSigma_trial, dP_dSigma_trial;
  softening_ = 0;
  this->compute_df_dp(trial_stress, dF_dSigma_trial, dP_dSigma_trial);
  double lambda_trial =
      yield_func_trial /
      ((dF_dSigma_trial.transpose() * de_).dot(dP_dSigma_trial.transpose()) +
       softening_);

  double p_multiplier;
  if (yield_state)
    p_multiplier = lambda;
  else if (!yield_state) {
    if (yield_state_trial)
      p_multiplier = lambda_trial;
    else if (!yield_state_trial)
      p_multiplier = 0.;
  }

  // update stress (plastic correction)
  Vector6d stress_update = trial_stress - (p_multiplier * de_ * dP_dSigma);
  // compute plastic deviatoric strain
  Vector6d dstress = stress - stress_update;
  Vector6d dpstrain = dstrain - (de_.inverse()) * dstress;
  if (Tdim == 2) dpstrain(4) = dpstrain(5) = 0.;
  // PDS_ += dPstrain;

  // compute equivalent plastic deviatoric strain
  // double epds_inc =
  // (2./3.)*sqrt(0.5*((dPstrain(0)-dPstrain(1))*(dPstrain(0)-dPstrain(1)) +
  // (dPstrain(1)-dPstrain(2))*(dPstrain(1)-dPstrain(2)) +
  // (dPstrain(0)-dPstrain(2))*(dPstrain(0)-dPstrain(2)))
  // + 3.*((dPstrain(3)*dPstrain(3)) + (dPstrain(4)*dPstrain(4)) +
  // (dPstrain(5)*dPstrain(5)))); epds_ += epds_inc;

  return stress_update;
}

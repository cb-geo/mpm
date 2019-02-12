//! Constructor with id and material properties
template <unsigned Tdim>
mpm::MohrCoulomb<Tdim>::MohrCoulomb(unsigned id,
                                    const Json& material_properties)
    : Material<Tdim>(id, material_properties) {
  try {
    // Elastic parameters
    youngs_modulus_ =
        material_properties["youngs_modulus"].template get<double>();
    poisson_ratio_ =
        material_properties["poisson_ratio"].template get<double>();

    // Peak friction, dilation and cohesion
    phi_max_ =
        material_properties["friction"].template get<double>() * PI / 180.;
    psi_max_ =
        material_properties["dilation"].template get<double>() * PI / 180.;
    cohesion_max_ = material_properties["cohesion"].template get<double>();

    // Residual friction, dilation and cohesion
    phi_residual_ =
        material_properties["residual_friction"].template get<double>() * PI /
        180.;
    psi_residual_ =
        material_properties["residual_dilation"].template get<double>() * PI /
        180.;
    cohesion_residual_ =
        material_properties["residual_cohesion"].template get<double>();

    // plastic deviatoric strain
    peak_epds_ = material_properties["peak_epds"].template get<double>();
    crit_epds_ = material_properties["critical_epds"].template get<double>();

    tension_cutoff_ =
        material_properties["tension_cutoff"].template get<double>();

    density_ = material_properties["density"].template get<double>();
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

//! Initialise state variables
//! equivalent plastic deviatoric strain (epds)
template <unsigned Tdim>
bool mpm::MohrCoulomb<Tdim>::initialise_state_variables(
    std::map<std::string, double>* state_vars) {
  bool status = false;

  // Equivalent plastic deviatoric strain
  status = state_vars->insert(std::make_pair("epds", 0.)).second;
  if (!status) {
    (*state_vars)["epds"] = 0.;
    status = true;
  }

  // Plastic deviatoric strain components
  status = state_vars->insert(std::make_pair("pds0", 0.)).second;
  if (!status) {
    (*state_vars)["pds0"] = 0.;
    status = true;
  }
  status = state_vars->insert(std::make_pair("pds1", 0.)).second;
  if (!status) {
    (*state_vars)["pds1"] = 0.;
    status = true;
  }
  status = state_vars->insert(std::make_pair("pds2", 0.)).second;
  if (!status) {
    (*state_vars)["pds2"] = 0.;
    status = true;
  }
  status = state_vars->insert(std::make_pair("pds3", 0.)).second;
  if (!status) {
    (*state_vars)["pds3"] = 0.;
    status = true;
  }
  status = state_vars->insert(std::make_pair("pds4", 0.)).second;
  if (!status) {
    (*state_vars)["pds4"] = 0.;
    status = true;
  }
  status = state_vars->insert(std::make_pair("pds5", 0.)).second;
  if (!status) {
    (*state_vars)["pds5"] = 0.;
    status = true;
  }
  return status;
}

//! Return elastic tensor
template <unsigned Tdim>
bool mpm::MohrCoulomb<Tdim>::compute_elastic_tensor() {
  // Shear modulus
  const double G = youngs_modulus_ / (2.0 * (1. + poisson_ratio_));
  const double a1 = bulk_modulus_ + (4.0 / 3.0) * G;
  const double a2 = bulk_modulus_ - (2.0 / 3.0) * G;

  // compute elastic stiffness matrix
  // clang-format off
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
bool mpm::MohrCoulomb<Tdim>::compute_rho_theta(const Vector6d& stress,
                                               double* j2, double* j3,
                                               double* rho, double* theta) {
  double mean_p = (stress(0) + stress(1) + stress(2)) / 3.;
  Vector6d dev_stress = Vector6d::Zero();
  dev_stress(0) = stress(0) - mean_p;
  dev_stress(1) = stress(1) - mean_p;
  dev_stress(2) = stress(2) - mean_p;
  dev_stress(3) = stress(3);
  if (Tdim == 3) {
    dev_stress(4) = stress(4);
    dev_stress(5) = stress(5);
  }

  // compute J2
  (*j2) = (pow((stress(0) - stress(1)), 2) + pow((stress(1) - stress(2)), 2) +
           pow((stress(0) - stress(2)), 2)) /
              6.0 +
          pow(stress(3), 2);
  if (Tdim == 3) (*j2) += pow(stress(4), 2) + pow(stress(5), 2);

  // compute J3
  (*j3) = (dev_stress(0) * dev_stress(1) * dev_stress(2)) -
          (dev_stress(2) * pow(dev_stress(3), 2));
  if (Tdim == 3)
    (*j3) += ((2 * dev_stress(3) * dev_stress(4) * dev_stress(5)) -
              (dev_stress(0) * pow(dev_stress(4), 2) -
               dev_stress(1) * pow(dev_stress(5), 2)));

  // compute theta value
  double theta_val = 0.;
  if (fabs(*j2) > 0.0)
    theta_val = (3. * sqrt(3.) / 2.) * ((*j3) / pow((*j2), 1.5));
  if (theta_val > 0.99) theta_val = 1.0;
  if (theta_val < -0.99) theta_val = -1.0;

  (*theta) = (1. / 3.) * acos(theta_val);
  if ((*theta) > PI / 3.) (*theta) = PI / 3.;
  if ((*theta) < 0.0) (*theta) = 0.;

  (*rho) = sqrt(2 * (*j2));
  return true;
}

//! Compute yield function in tension and shear
template <unsigned Tdim>
Eigen::Matrix<double, 2, 1> mpm::MohrCoulomb<Tdim>::compute_yield(
    double epsilon, double rho, double theta, double phi, double cohesion) {
  Eigen::Matrix<double, 2, 1> yield_function;

  // Tension
  yield_function(0) =
      sqrt(2. / 3.) * cos(theta) * rho + epsilon / sqrt(3.) - tension_cutoff_;

  // Shear
  yield_function(1) = sqrt(3. / 2.) * rho *
                          ((sin(theta + PI / 3.) / (sqrt(3.) * cos(phi))) +
                           (cos(theta + PI / 3.) * tan(phi) / 3.)) +
                      (epsilon / sqrt(3.)) * tan(phi) - cohesion;

  return yield_function;
}

//! Check the yield state and return the value of yield function
template <unsigned Tdim>
int mpm::MohrCoulomb<Tdim>::check_yield(
    const Eigen::Matrix<double, 2, 1>& yield_function, double epsilon,
    double rho, double theta, double phi, double cohesion) {

  const double yield_tension = yield_function(0);

  const double yield_shear = yield_function(1);

  // Yield type 0: elastic, 1: tension failure, 2: shear failure
  int yield_type = 0;
  // Check for tension or shear
  if (yield_tension > 1.E-22 && yield_shear > 1.E-22) {
    double n_phi = (1. + sin(phi)) / (1. - sin(phi));
    double sigma_p = tension_cutoff_ * n_phi - 2. * cohesion * sqrt(n_phi);
    double alpha_p = sqrt(1. + n_phi * n_phi) + n_phi;
    // Compute the shear-tension edge
    double h = sqrt(2. / 3.) * cos(theta) * rho + epsilon / sqrt(3.) -
               tension_cutoff_ +
               alpha_p * (sqrt(2. / 3.) * cos(theta - 4. * PI / 3.) * rho +
                          epsilon / sqrt(3.) - sigma_p);
    // Tension
    if (h > 1.E-22) yield_type = 1;
    // Shear
    else
      yield_type = 2;
  }

  // Shear
  if (yield_tension < 1.E-22 && yield_shear > 1.E-22) yield_type = 2;

  // Tension
  if (yield_tension > 1.E-22 && yield_shear < 1.E-22) yield_type = 1;

  return yield_type;
}

//! Compute dF/dSigma and dP/dSigma
template <unsigned Tdim>
void mpm::MohrCoulomb<Tdim>::compute_df_dp(
    int yield_type, double j2, double j3, double rho, double theta, double phi,
    double psi, double cohesion, const Vector6d& stress, double epds,
    Vector6d* df_dsigma, Vector6d* dp_dsigma, double* softening,
    const ParticleBase<Tdim>* ptr) {

  // mean stress
  double mean_p = (stress(0) + stress(1) + stress(2)) / 3.0;

  // deviatoric stress
  Vector6d dev_stress = Vector6d::Zero();
  dev_stress(0) = stress(0) - mean_p;
  dev_stress(1) = stress(1) - mean_p;
  dev_stress(2) = stress(2) - mean_p;
  dev_stress(3) = stress(3);
  if (Tdim == 3) {
    dev_stress(4) = stress(4);
    dev_stress(5) = stress(5);
  }

  // Compute dF / dEpsilon,  dF / dRho, dF / dTheta
  double df_depsilon, df_drho, df_dtheta;
  // Values in tension yield
  if (yield_type == 1) {
    df_depsilon = 1. / sqrt(3.);
    df_drho = sqrt(2. / 3.) * cos(theta);
    df_dtheta = -sqrt(2. / 3.) * rho * sin(theta);
  }
  // Values in shear yield
  else {
    df_depsilon = tan(phi) / sqrt(3.);
    df_drho = sqrt(3. / 2.) * ((sin(theta + PI / 3.) / (sqrt(3.) * cos(phi))) +
                               (cos(theta + PI / 3.) * tan(phi) / 3.));
    df_dtheta = sqrt(3. / 2.) * rho *
                ((cos(theta + PI / 3.) / (sqrt(3.) * cos(phi))) -
                 (sin(theta + PI / 3.) * tan(phi) / 3.));
  }

  // Compute dEpsilon / dSigma (the same in two yield types)
  Vector6d depsilon_dsigma = Vector6d::Zero();
  depsilon_dsigma(0) = 1. / sqrt(3.);
  depsilon_dsigma(1) = 1. / sqrt(3.);
  depsilon_dsigma(2) = 1. / sqrt(3.);

  // Compute dRho / dSigma (the same in two yield types)
  Vector6d drho_dsigma = Vector6d::Zero();
  double multiplier = 1.;
  if (fabs(rho) > 0.) multiplier = 1. / rho;
  drho_dsigma = multiplier * dev_stress;
  if (Tdim == 2) {
    drho_dsigma(4) = 0.;
    drho_dsigma(5) = 0.;
  }

  // Compute dTheta / dSigma (the same in two yield types)
  double r_val = 0.;
  if (fabs(j2) > 1.E-22) r_val = (3. * sqrt(3.) / 2.) * (j3 / pow(j2, 1.5));

  double divider = 1 - (r_val * r_val);
  if (divider <= 0.) divider = 1.0E-3;
  double dtheta_dr = -1 / (3. * sqrt(divider));

  double dr_dj2 = (-9 * sqrt(3.) / 4.) * j3;
  if (fabs(j2) > 1.0E-22) dr_dj2 = dr_dj2 / pow(j2, 2.5);

  double dr_dj3 = 1.5 * sqrt(3.);
  if (fabs(j2) > 1.0E-22) dr_dj3 = dr_dj3 / pow(j2, 1.5);

  Vector6d dj2_dsigma = dev_stress;
  dj2_dsigma(3) = dev_stress(3);

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

  Vector6d dj3_dsigma = Vector6d::Zero();
  dj3_dsigma(0) = dev1.dot(dev1) - (2. / 3.) * j2;
  dj3_dsigma(1) = dev2.dot(dev2) - (2. / 3.) * j2;
  dj3_dsigma(2) = dev3.dot(dev3) - (2. / 3.) * j2;
  dj3_dsigma(3) = dev1.dot(dev2);
  if (Tdim == 3) {
    dj3_dsigma(4) = dev2.dot(dev3);
    dj3_dsigma(5) = dev1.dot(dev3);
  }

  // compute dtheta / dsigma (the same in two yield types)
  Vector6d dtheta_dsigma = Vector6d::Zero();
  dtheta_dsigma = dtheta_dr * ((dr_dj2 * dj2_dsigma) + (dr_dj3 * dj3_dsigma));
  if (Tdim == 2) {
    dtheta_dsigma(4) = 0.;
    dtheta_dsigma(5) = 0.;
  }

  // Compute dF/dSigma
  (*df_dsigma) = (df_depsilon * depsilon_dsigma) + (df_drho * drho_dsigma) +
                 (df_dtheta * dtheta_dsigma);
  if (Tdim == 2) {
    (*df_dsigma)(4) = 0.;
    (*df_dsigma)(5) = 0.;
  }

  // compute dp/dsigma, dp/dj
  double dp_dj = 0.;
  if (yield_type == 1) {
    double et_value = 0.6;
    double xit = 0.1;
    double rt_den =
        2. * (1 - et_value * et_value) * cos(theta) +
        (2. * et_value - 1) *
            sqrt(4. * (1 - et_value * et_value) * cos(theta) * cos(theta) +
                 5. * et_value * et_value - 4. * et_value);
    double rt_num = 4. * (1 - et_value * et_value) * cos(theta) * cos(theta) +
                    (2. * et_value - 1) * (2. * et_value - 1);
    double rt = rt_num / (3. * rt_den);
    double dp_drt = 1.5 * rho * rho * rt /
                    sqrt(xit * xit * tension_cutoff_ * tension_cutoff_ +
                         1.5 * rt * rt * rho * rho);
    double dp_drho = 1.5 * rho * rt * rt /
                     sqrt(xit * xit * tension_cutoff_ * tension_cutoff_ +
                          1.5 * rt * rt * rho * rho);
    double dp_depsilon = 1. / sqrt(3.);
    double drtden_dtheta =
        -2. * (1 - et_value * et_value) * sin(theta) -
        (2. * et_value - 1) * 4. * (1 - et_value * et_value) * cos(theta) *
            sin(theta) /
            sqrt(4. * (1 - et_value * et_value) * cos(theta) * cos(theta) +
                 5. * et_value * et_value - 4. * et_value);
    double drtnum_dtheta =
        -8. * (1 - et_value * et_value) * cos(theta) * sin(theta);
    double drt_dtheta = (drtnum_dtheta * rt_den - drtden_dtheta * rt_num) /
                        (3. * rt_den * rt_den);
    // compute the value of dp/dsigma and dp/dj in tension yield
    (*dp_dsigma) = (dp_depsilon * depsilon_dsigma) + (dp_drho * drho_dsigma) +
                   (dp_drt * drt_dtheta * dtheta_dsigma);
    dp_dj = dp_drho * sqrt(2.);
  } else {
    double r_mc = (3. - sin(phi)) / (6 * cos(phi));
    double e_val = (3. - sin(phi)) / (3. + sin(phi));
    if ((e_val - 0.5) < 0.) e_val = 0.501;
    if ((e_val - 1.) > 0.) e_val = 1.0;
    double sqpart = (4. * (1 - e_val * e_val) * pow(cos(theta), 2)) +
                    (5 * e_val * e_val) - (4. * e_val);
    if (sqpart < 0.) sqpart = 0.00001;
    double r_mw_den = (2. * (1 - e_val * e_val) * cos(theta)) +
                      ((2. * e_val - 1) * sqrt(sqpart));
    if (fabs(r_mw_den) < 1.e-22) r_mw_den = 0.001;
    double r_mw_num = (4. * (1. - e_val * e_val) * pow(cos(theta), 2)) +
                      pow((2. * e_val - 1.), 2);
    double r_mw = (r_mw_num / r_mw_den) * r_mc;
    const double xi = 0.1;
    double omega = pow((xi * cohesion * tan(psi)), 2) +
                   pow((r_mw * sqrt(3. / 2.) * rho), 2);
    if (omega < 1.e-22) omega = 0.001;
    double l = r_mw_num;
    double m = r_mw_den;
    double dl_dtheta = -8. * (1. - e_val * e_val) * cos(theta) * sin(theta);
    double dm_dtheta = (-2. * (1. - e_val * e_val) * sin(theta)) +
                       (0.5 * (2. * e_val - 1.) * dl_dtheta) / sqrt(sqpart);
    double drmw_dtheta = ((m * dl_dtheta) - (l * dm_dtheta)) / (m * m);
    double dp_depsilon = tan(psi) / sqrt(3.);
    double dp_drho = 3. * rho * r_mw * r_mw / (2. * sqrt(omega));
    double dp_dtheta =
        (3. * rho * rho * r_mw * r_mc * drmw_dtheta) / (2. * sqrt(omega));
    // compute the value of dp/dsigma and dp/dj in shear yield
    (*dp_dsigma) = (dp_depsilon * depsilon_dsigma) + (dp_drho * drho_dsigma) +
                   (dp_dtheta * dtheta_dsigma);
    dp_dj = dp_drho * sqrt(2.);
  }

  // compute softening part
  double dphi_dpstrain = 0.;
  double dc_dpstrain = 0.;
  if (yield_type == 2 && epds > peak_epds_ && epds < crit_epds_) {
    dphi_dpstrain = (phi_residual_ - phi) / (crit_epds_ - peak_epds_);
    dc_dpstrain = (cohesion_residual_ - cohesion) / (crit_epds_ - peak_epds_);
  }
  double epsilon = (1. / sqrt(3.)) * (stress(0) + stress(1) + stress(2));
  double df_dphi = sqrt(3. / 2.) * rho *
                       ((sin(phi) * sin(theta + PI / 3.) /
                         (sqrt(3.) * cos(phi) * cos(phi))) +
                        (cos(theta + PI / 3.) / (3. * cos(phi) * cos(phi)))) +
                   (epsilon / (sqrt(3.) * cos(phi) * cos(phi)));
  double df_dc = -1.;
  (*softening) =
      (-1.) * ((df_dphi * dphi_dpstrain) + (df_dc * dc_dpstrain)) * dp_dj;
}

//! Compute stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::MohrCoulomb<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr, std::map<std::string, double>* state_vars) {

  // Current MC parameters using a linear softening rule
  // plastic deviatoric strain
  Eigen::Matrix<double, 6, 1> plastic_deviatoric_strain;
  plastic_deviatoric_strain(0) = (*state_vars).at("pds0");
  plastic_deviatoric_strain(1) = (*state_vars).at("pds1");
  plastic_deviatoric_strain(2) = (*state_vars).at("pds2");
  plastic_deviatoric_strain(3) = (*state_vars).at("pds3");
  plastic_deviatoric_strain(4) = (*state_vars).at("pds4");
  plastic_deviatoric_strain(5) = (*state_vars).at("pds5");

  const double epds =
      (2. / 3.) *
      sqrt(3. / 2. *
               (plastic_deviatoric_strain(0) * plastic_deviatoric_strain(0) +
                plastic_deviatoric_strain(1) * plastic_deviatoric_strain(1) +
                plastic_deviatoric_strain(2) * plastic_deviatoric_strain(2)) +
           3. * (plastic_deviatoric_strain(3) * plastic_deviatoric_strain(3) +
                 plastic_deviatoric_strain(4) * plastic_deviatoric_strain(4) +
                 plastic_deviatoric_strain(5) * plastic_deviatoric_strain(5)));
  Vector6d trial_stress = stress + (this->de_ * dstrain);

  double j2 = 0;
  double j3 = 0;
  double rho = 0;
  double theta = 0;
  this->compute_rho_theta(stress, &j2, &j3, &rho, &theta);

  double epsilon = (1. / sqrt(3.)) * (stress(0) + stress(1) + stress(2));

  // Friction, dilation and cohesion
  double phi = phi_max_;
  double psi = psi_max_;
  double cohesion = cohesion_max_;

  Eigen::Matrix<double, 2, 1> yield_function =
      this->compute_yield(epsilon, rho, theta, phi, cohesion);

  int yield_type =
      this->check_yield(yield_function, epsilon, rho, theta, phi, cohesion);
  
  if (yield_type != 1 && (epds - peak_epds_) > 0. && (crit_epds_ - epds) > 0.) {
    phi = phi_residual_ + ((phi_max_ - phi_residual_) * (epds - crit_epds_) /
                            (peak_epds_ - crit_epds_));
    psi = psi_residual_ + ((psi_max_ - psi_residual_) * (epds - crit_epds_) /
                            (peak_epds_ - crit_epds_));
    cohesion =
        cohesion_residual_ + ((cohesion_max_ - cohesion_residual_) *
                              (epds - crit_epds_) / (peak_epds_ - crit_epds_));
  } else if (yield_type != 1 && (epds - crit_epds_) >= 0.) {
    phi = phi_residual_;
    psi = psi_residual_;
    cohesion = cohesion_residual_;
  }

  // Yield function for the current stress state
  yield_function = this->compute_yield(epsilon, rho, theta, phi, cohesion);
  yield_type = this->check_yield(yield_function, epsilon, rho, theta, phi, cohesion);
  double softening = 0.;

  // Compute plastic multiplier from the current stress state
  Vector6d df_dsigma = Vector6d::Zero();
  Vector6d dp_dsigma = Vector6d::Zero();
  this->compute_df_dp(yield_type, j2, j3, rho, theta, phi, psi, cohesion,
                      stress, epds, &df_dsigma, &dp_dsigma, &softening, ptr);
  // Check the epds
  if ((*state_vars)["epds"] < peak_epds_ && epds > peak_epds_) softening = 0;
  double lambda = df_dsigma.dot(this->de_ * dstrain) /
                  ((df_dsigma.dot(this->de_ * dp_dsigma)) + softening);

  // Compute the trial stress
  j2 = j3 = rho = theta = 0.;
  this->compute_rho_theta(trial_stress, &j2, &j3, &rho, &theta);

  epsilon =
      (1. / sqrt(3.)) * (trial_stress(0) + trial_stress(1) + trial_stress(2));

  // Trial yield function
  Eigen::Matrix<double, 2, 1> yield_function_trial =
      this->compute_yield(epsilon, rho, theta, phi, cohesion);
  int yield_type_trial =
      this->check_yield(yield_function_trial, epsilon, rho, theta, phi, cohesion);

  double softening_trial = 0.;
  Vector6d df_dsigma_trial = Vector6d::Zero();
  Vector6d dp_dsigma_trial = Vector6d::Zero();
  this->compute_df_dp(yield_type_trial, j2, j3, rho, theta, phi, psi, cohesion,
                      trial_stress, epds, &df_dsigma_trial, &dp_dsigma_trial,
                      &softening_trial, ptr);

  // Check the epds
  if ((*state_vars)["epds"] < peak_epds_ && epds > peak_epds_)
    softening_trial = 0;
  double lambda_trial = 0.;

  if (yield_type_trial == 1)
    lambda_trial =
        yield_function_trial(0) /
        ((df_dsigma_trial.transpose() * de_).dot(dp_dsigma_trial.transpose()) +
         softening_trial);
  else
    lambda_trial =
        yield_function_trial(1) /
        ((df_dsigma_trial.transpose() * de_).dot(dp_dsigma_trial.transpose()) +
         softening_trial);
  // Compute the correction stress
  double p_multiplier = 0.;
  Vector6d dp_dsigma_final = Vector6d::Zero();
  bool update_pds = false;
  // check if it is a load process or an unload process
  double dyfun_t = yield_function_trial(0) - yield_function(0);
  double dyfun_s = yield_function_trial(1) - yield_function(1);
  if (yield_type != 0) {
    if (yield_type == 1 && dyfun_t > 0) {
      p_multiplier = lambda;
      dp_dsigma_final = dp_dsigma;
    } else if (yield_type == 2 && dyfun_s > 0) {
      p_multiplier = lambda;
      dp_dsigma_final = dp_dsigma;
      update_pds = true;
    }
  } else if (yield_type_trial != 0) {
    p_multiplier = lambda_trial;
    dp_dsigma_final = dp_dsigma_trial;
    if(yield_type_trial == 2) update_pds = true;
  }
  // update stress (plastic correction)
  Vector6d updated_stress =
      trial_stress - (p_multiplier * this->de_ * dp_dsigma_final);
  // compute plastic deviatoric strain
  Vector6d dstress = updated_stress - stress;
  Vector6d dpstrain = dstrain - (this->de_.inverse()) * dstress;
  if (Tdim == 2) dpstrain(4) = dpstrain(5) = 0.;

  // Record the epds
  (*state_vars)["epds"] = epds;

  // Update plastic deviatoric strain if it is a load process
  if (update_pds) {
    Vector6d dp_ds = dpstrain;
    double dpvstrain = p_multiplier * (dp_dsigma_final(0) + dp_dsigma_final(1) +
                                       dp_dsigma_final(2));
    dp_ds(0) -= (1. / 3. * dpvstrain);
    dp_ds(1) -= (1. / 3. * dpvstrain);
    dp_ds(2) -= (1. / 3. * dpvstrain);
    dp_ds(3) = 0.5 * dpstrain(3);
    dp_ds(4) = 0.5 * dpstrain(4);
    dp_ds(5) = 0.5 * dpstrain(5);
    plastic_deviatoric_strain += dp_ds;
    (*state_vars).at("pds0") = plastic_deviatoric_strain(0);
    (*state_vars).at("pds1") = plastic_deviatoric_strain(1);
    (*state_vars).at("pds2") = plastic_deviatoric_strain(2);
    (*state_vars).at("pds3") = plastic_deviatoric_strain(3);
    (*state_vars).at("pds4") = plastic_deviatoric_strain(4);
    (*state_vars).at("pds5") = plastic_deviatoric_strain(5);
  }
  return updated_stress;
}

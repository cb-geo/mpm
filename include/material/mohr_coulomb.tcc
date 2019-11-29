//! Constructor with id and material properties
template <unsigned Tdim>
mpm::MohrCoulomb<Tdim>::MohrCoulomb(unsigned id,
                                    const Json& material_properties)
    : Material<Tdim>(id, material_properties) {
  try {
    // General parameters
    // Density
    density_ = material_properties["density"].template get<double>();
    // Young's modulus
    youngs_modulus_ =
        material_properties["youngs_modulus"].template get<double>();
    // Poisson ratio
    poisson_ratio_ =
        material_properties["poisson_ratio"].template get<double>();
    // Softening status
    softening_ = material_properties["softening"].template get<bool>();
    // Peak friction, dilation and cohesion
    phi_peak_ =
        material_properties["friction"].template get<double>() * M_PI / 180.;
    psi_peak_ =
        material_properties["dilation"].template get<double>() * M_PI / 180.;
    cohesion_peak_ = material_properties["cohesion"].template get<double>();
    // Residual friction, dilation and cohesion
    phi_residual_ =
        material_properties["residual_friction"].template get<double>() * M_PI /
        180.;
    psi_residual_ =
        material_properties["residual_dilation"].template get<double>() * M_PI /
        180.;
    cohesion_residual_ =
        material_properties["residual_cohesion"].template get<double>();
    // Peak equivalent plastic deviatoric strain
    epds_peak_ = material_properties["peak_epds"].template get<double>();
    // Residual equivalent plastic deviatoric strain
    epds_residual_ =
        material_properties["critical_epds"].template get<double>();
    // Tensile strength
    tension_cutoff_ =
        material_properties["tension_cutoff"].template get<double>();
    // Properties
    properties_ = material_properties;
    // Bulk modulus
    bulk_modulus_ = youngs_modulus_ / (3.0 * (1. - 2. * poisson_ratio_));
    // Shear modulus
    shear_modulus_ = youngs_modulus_ / (2.0 * (1 + poisson_ratio_));
    // Set elastic tensor
    this->compute_elastic_tensor();
  } catch (std::exception& except) {
    console_->error("Material parameter not set: {}\n", except.what());
  }
}

//! Initialise state variables
template <unsigned Tdim>
mpm::dense_map mpm::MohrCoulomb<Tdim>::initialise_state_variables() {
  mpm::dense_map state_vars = {// MC parameters
                               // Friction (phi)
                               {"phi", this->phi_peak_},
                               // Dilation (psi)
                               {"psi", this->psi_peak_},
                               // Cohesion
                               {"cohesion", this->cohesion_peak_},
                               // Stress invariants
                               // j3
                               {"j3", 0.},
                               // Epsilon
                               {"epsilon", 0.},
                               // Rho
                               {"rho", 0.},
                               // Theta
                               {"theta", 0.},
                               // Plastic strain
                               // Equivalent plastic deviatoric strain
                               {"epds", 0.},
                               // Plastic strain components
                               {"pstrain0", 0.},
                               {"pstrain1", 0.},
                               {"pstrain2", 0.},
                               {"pstrain3", 0.},
                               {"pstrain4", 0.},
                               {"pstrain5", 0.}};

  return state_vars;
}

//! Compute elastic tensor
template <unsigned Tdim>
bool mpm::MohrCoulomb<Tdim>::compute_elastic_tensor() {
  // Shear modulus
  const double G = shear_modulus_;
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

//! Compute stress invariants
template <unsigned Tdim>
bool mpm::MohrCoulomb<Tdim>::compute_stress_invariants(
    const Vector6d& stress, mpm::dense_map* state_vars) {
  // Compute the mean pressure
  const double mean_p = (stress(0) + stress(1) + stress(2)) / 3.;
  // Compute the deviatoric stress
  Vector6d dev_stress = stress;
  for (unsigned i = 0; i < 3; ++i) dev_stress(i) -= mean_p;
  // Compute J2
  double j2 = (std::pow((stress(0) - stress(1)), 2) +
               std::pow((stress(1) - stress(2)), 2) +
               std::pow((stress(0) - stress(2)), 2)) /
                  6.0 +
              std::pow(stress(3), 2) + std::pow(stress(4), 2) +
              std::pow(stress(5), 2);
  // Compute J3
  (*state_vars).at("j3") = dev_stress(0) * dev_stress(1) * dev_stress(2) -
                           dev_stress(2) * std::pow(dev_stress(3), 2) +
                           2 * dev_stress(3) * dev_stress(4) * dev_stress(5) -
                           dev_stress(0) * std::pow(dev_stress(4), 2) -
                           dev_stress(1) * std::pow(dev_stress(5), 2);
  // Compute theta value (Lode angle)
  double theta_val = 0.;
  if (fabs(j2) > 0.0)
    theta_val =
        (1.5 * std::sqrt(3.)) * ((*state_vars).at("j3") / std::pow(j2, 1.5));
  // Check theta value
  if (theta_val > 1.0) theta_val = 1.0;
  if (theta_val < -1.0) theta_val = -1.0;
  // Compute theta
  (*state_vars)["theta"] = (1. / 3.) * acos(theta_val);
  // Check theta
  if ((*state_vars).at("theta") > M_PI / 3.) (*state_vars)["theta"] = M_PI / 3.;
  if ((*state_vars).at("theta") < 0.0) (*state_vars)["theta"] = 0.;
  // Compute rho
  (*state_vars)["rho"] = std::sqrt(2. * (j2));
  // Compute epsilon
  (*state_vars)["epsilon"] =
      (1. / std::sqrt(3.)) * (stress(0) + stress(1) + stress(2));

  return true;
}

//! Compute yield function and yield state
template <unsigned Tdim>
typename mpm::FailureState mpm::MohrCoulomb<Tdim>::compute_yield_state(
    Eigen::Matrix<double, 2, 1>* yield_function,
    const mpm::dense_map& state_vars) {
  // Tolerance for yield function
  const double Tolerance = -1E-1;
  // Get stress invariants
  const double& epsilon = state_vars.at("epsilon");
  const double& rho = state_vars.at("rho");
  const double& theta = state_vars.at("theta");
  // Get MC parameters
  const double& phi = state_vars.at("phi");
  const double& cohesion = state_vars.at("cohesion");
  // Compute yield functions (tension & shear)
  // Tension
  (*yield_function)(0) = std::sqrt(2. / 3.) * cos(theta) * rho +
                         epsilon / std::sqrt(3.) - tension_cutoff_;
  // Shear
  (*yield_function)(1) =
      std::sqrt(1.5) * rho *
          ((sin(theta + M_PI / 3.) / (std::sqrt(3.) * cos(phi))) +
           (cos(theta + M_PI / 3.) * tan(phi) / 3.)) +
      (epsilon / std::sqrt(3.)) * tan(phi) - cohesion;
  // Initialise yield status (0: elastic, 1: tension failure, 2: shear failure)
  auto yield_type = FailureState::Elastic;
  // Check for tension and shear
  if ((*yield_function)(0) > Tolerance && (*yield_function)(1) > Tolerance) {
    // Compute tension and shear edge parameters
    double n_phi = (1. + sin(phi)) / (1. - sin(phi));
    double sigma_p = tension_cutoff_ * n_phi - 2. * cohesion * std::sqrt(n_phi);
    double alpha_p = std::sqrt(1. + n_phi * n_phi) + n_phi;
    // Compute the shear-tension edge
    double h =
        (*yield_function)(0) +
        alpha_p * (std::sqrt(2. / 3.) * cos(theta - 4. * M_PI / 3.) * rho +
                   epsilon / std::sqrt(3.) - sigma_p);
    // Tension
    if (h > 1.E-22) yield_type = FailureState::Tensile;
    // Shear
    else
      yield_type = FailureState::Shear;
  }
  // Shear failure
  if ((*yield_function)(0) < Tolerance && (*yield_function)(1) > Tolerance)
    yield_type = FailureState::Shear;
  // Tension failure
  if ((*yield_function)(0) > Tolerance && (*yield_function)(1) < Tolerance)
    yield_type = FailureState::Tensile;

  return yield_type;
}

//! Compute dF/dSigma and dP/dSigma
template <unsigned Tdim>
void mpm::MohrCoulomb<Tdim>::compute_df_dp(mpm::FailureState yield_type,
                                           const mpm::dense_map* state_vars,
                                           const Vector6d& stress,
                                           Vector6d* df_dsigma,
                                           Vector6d* dp_dsigma,
                                           double* softening) {
  // Get stress invariants
  const double& j3 = (*state_vars).at("j3");
  const double& rho = (*state_vars).at("rho");
  const double& theta = (*state_vars).at("theta");
  // Get MC parameters
  const double& phi = (*state_vars).at("phi");
  const double& psi = (*state_vars).at("psi");
  const double& cohesion = (*state_vars).at("cohesion");
  // Get equivalent plastic deviatoric strain
  const double& epds = (*state_vars).at("epds");
  // Compute the mean stress
  double mean_p = (stress(0) + stress(1) + stress(2)) / 3.0;
  // Compute deviatoric stress
  Vector6d dev_stress = stress;
  for (unsigned i = 0; i < 3; ++i) dev_stress(i) -= mean_p;
  // Compute dF / dEpsilon,  dF / dRho, dF / dTheta
  double df_depsilon, df_drho, df_dtheta;
  // Values in tension yield
  if (yield_type == FailureState::Tensile) {
    df_depsilon = 1. / std::sqrt(3.);
    df_drho = std::sqrt(2. / 3.) * cos(theta);
    df_dtheta = -std::sqrt(2. / 3.) * rho * sin(theta);
  }
  // Values in shear yield / elastic
  else {
    df_depsilon = tan(phi) / std::sqrt(3.);
    df_drho = std::sqrt(1.5) *
              ((sin(theta + M_PI / 3.) / (std::sqrt(3.) * cos(phi))) +
               (cos(theta + M_PI / 3.) * tan(phi) / 3.));
    df_dtheta = std::sqrt(1.5) * rho *
                ((cos(theta + M_PI / 3.) / (std::sqrt(3.) * cos(phi))) -
                 (sin(theta + M_PI / 3.) * tan(phi) / 3.));
  }
  // Compute dEpsilon / dSigma
  Vector6d depsilon_dsigma = Vector6d::Zero();
  depsilon_dsigma(0) = depsilon_dsigma(1) = depsilon_dsigma(2) =
      1. / std::sqrt(3.);
  // Compute dRho / dSigma
  Vector6d drho_dsigma = Vector6d::Zero();
  double multiplier = 1.;
  if (fabs(rho) > 0.) multiplier = 1. / rho;
  drho_dsigma = multiplier * dev_stress;
  if (Tdim == 2) {
    drho_dsigma(4) = 0.;
    drho_dsigma(5) = 0.;
  }
  // Compute dTheta / dSigma
  // Compute r
  double r_val = 0.;
  if (fabs(rho) > 1.E-22)
    r_val = (3. * std::sqrt(6.)) * (j3 / std::pow(rho, 3));
  // Compute dTheta / dr
  double divider = 1 - (r_val * r_val);
  if (divider <= 0.) divider = 1.E-3;
  double dtheta_dr = -1 / (3. * std::sqrt(divider));
  // Compute dr / dJ2
  double dr_dj2 = (-9 * std::sqrt(6.)) * j3;
  if (fabs(rho) > 1.E-22) dr_dj2 /= std::pow(rho, 5);
  // Compute dr / dJ3
  double dr_dj3 = 3. * std::sqrt(6.);
  if (fabs(rho) > 1.E-22) dr_dj3 /= std::pow(rho, 3);
  // Compute dJ2 / dSigma
  Vector6d dj2_dsigma = dev_stress;
  // Compute dJ3 / dSigma
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
  dj3_dsigma(0) = dev1.dot(dev1) - (1. / 3.) * rho * rho;
  dj3_dsigma(1) = dev2.dot(dev2) - (1. / 3.) * rho * rho;
  dj3_dsigma(2) = dev3.dot(dev3) - (1. / 3.) * rho * rho;
  dj3_dsigma(3) = dev1.dot(dev2);
  dj3_dsigma(4) = dev2.dot(dev3);
  dj3_dsigma(5) = dev1.dot(dev3);
  // Compute dtheta / dsigma
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
  // Initialise dp/dq
  double dp_dq = 0.;
  // Compute dp/dsigma and dp/dj in tension yield
  if (yield_type == FailureState::Tensile) {
    // Define deviatoric eccentricity
    double et_value = 0.6;
    // Define meridional eccentricity
    double xit = 0.1;
    // Compute Rt
    double sqpart = 4. * (1 - et_value * et_value) * cos(theta) * cos(theta) +
                    5. * et_value * et_value - 4. * et_value;
    if (sqpart < 1.E-22) sqpart = 1.E-5;
    double rt_den = 2. * (1 - et_value * et_value) * cos(theta) +
                    (2. * et_value - 1) * std::sqrt(sqpart);
    double rt_num = 4. * (1 - et_value * et_value) * cos(theta) * cos(theta) +
                    (2. * et_value - 1) * (2. * et_value - 1);
    if (fabs(rt_den) < 1.E-22) rt_den = 1.E-5;
    double rt = rt_num / (3. * rt_den);
    // Compute dP/dRt
    double dp_drt = 1.5 * rho * rho * rt /
                    std::sqrt(xit * xit * tension_cutoff_ * tension_cutoff_ +
                              1.5 * rt * rt * rho * rho);
    // Compute dP/dRho
    double dp_drho = 1.5 * rho * rt * rt /
                     std::sqrt(xit * xit * tension_cutoff_ * tension_cutoff_ +
                               1.5 * rt * rt * rho * rho);
    // Compute dP/dEpsilon
    double dp_depsilon = 1. / std::sqrt(3.);
    // Compute dRt/dThera
    double drtden_dtheta =
        -2. * (1 - et_value * et_value) * sin(theta) -
        (2. * et_value - 1) * 4. * (1 - et_value * et_value) * cos(theta) *
            sin(theta) /
            std::sqrt(4. * (1 - et_value * et_value) * cos(theta) * cos(theta) +
                      5. * et_value * et_value - 4. * et_value);
    double drtnum_dtheta =
        -8. * (1 - et_value * et_value) * cos(theta) * sin(theta);
    double drt_dtheta = (drtnum_dtheta * rt_den - drtden_dtheta * rt_num) /
                        (3. * rt_den * rt_den);
    // Compute dP/dSigma
    (*dp_dsigma) = (dp_depsilon * depsilon_dsigma) + (dp_drho * drho_dsigma) +
                   (dp_drt * drt_dtheta * dtheta_dsigma);
    // Compute dP/dJ
    dp_dq = dp_drho * std::sqrt(2. / 3.);
  }
  // Compute dp/dsigma and dp/dj in shear yield
  else {
    // Compute Rmc
    double r_mc = (3. - sin(phi)) / (6 * cos(phi));
    // Compute deviatoric eccentricity
    double e_val = (3. - sin(phi)) / (3. + sin(phi));
    if (e_val <= 0.5) e_val = 0.5 + 1.E-10;
    if (e_val > 1.) e_val = 1.;
    // Compute Rmw
    double sqpart = (4. * (1 - e_val * e_val) * std::pow(cos(theta), 2)) +
                    (5 * e_val * e_val) - (4. * e_val);
    if (sqpart < 1.E-22) sqpart = 1.E-5;
    double m = (2. * (1 - e_val * e_val) * cos(theta)) +
               ((2. * e_val - 1) * std::sqrt(sqpart));
    if (fabs(m) < 1.E-22) m = 1.E-5;
    double l = (4. * (1. - e_val * e_val) * std::pow(cos(theta), 2)) +
               std::pow((2. * e_val - 1.), 2);
    double r_mw = (l / m) * r_mc;
    // Initialise meridional eccentricity
    double xi = 0.1;
    double omega = std::pow((xi * cohesion_peak_ * tan(psi)), 2) +
                   std::pow((r_mw * std::sqrt(1.5) * rho), 2);
    if (omega < 1.E-22) omega = 1.E-5;
    double dl_dtheta = -8. * (1. - e_val * e_val) * cos(theta) * sin(theta);
    double dm_dtheta =
        (-2. * (1. - e_val * e_val) * sin(theta)) +
        (0.5 * (2. * e_val - 1.) * dl_dtheta) / std::sqrt(sqpart);
    double drmw_dtheta = ((m * dl_dtheta) - (l * dm_dtheta)) / (m * m);
    double dp_depsilon = tan(psi) / std::sqrt(3.);
    double dp_drho = 3. * rho * r_mw * r_mw / (2. * std::sqrt(omega));
    double dp_dtheta =
        (3. * rho * rho * r_mw * r_mc * drmw_dtheta) / (2. * std::sqrt(omega));
    // compute the value of dp/dsigma and dp/dj in shear yield
    (*dp_dsigma) = (dp_depsilon * depsilon_dsigma) + (dp_drho * drho_dsigma) +
                   (dp_dtheta * dtheta_dsigma);
    dp_dq = dp_drho * std::sqrt(2. / 3.);
  }
  // Compute softening part
  double dphi_dpstrain = 0.;
  double dc_dpstrain = 0.;
  (*softening) = 0.;
  if (softening_ && epds > epds_peak_ && epds < epds_residual_) {
    // Compute dPhi/dPstrain
    dphi_dpstrain = (phi_residual_ - phi_peak_) / (epds_residual_ - epds_peak_);
    // Compute dc/dPstrain
    dc_dpstrain =
        (cohesion_residual_ - cohesion_peak_) / (epds_residual_ - epds_peak_);
    // Compute dF/dPstrain
    double df_dphi =
        std::sqrt(1.5) * rho *
            ((sin(phi) * sin(theta + M_PI / 3.) /
              (std::sqrt(3.) * cos(phi) * cos(phi))) +
             (cos(theta + M_PI / 3.) / (3. * cos(phi) * cos(phi)))) +
        (mean_p / (cos(phi) * cos(phi)));
    double df_dc = -1.;
    (*softening) =
        (-1.) * ((df_dphi * dphi_dpstrain) + (df_dc * dc_dpstrain)) * dp_dq;
  }
}

//! Compute stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::MohrCoulomb<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr, mpm::dense_map* state_vars) {
  // Get equivalent plastic deviatoric strain
  const double& epds = (*state_vars).at("epds");
  // Update MC parameters using a linear softening rule
  if (softening_ && (epds - epds_peak_) > 0. && (epds_residual_ - epds) > 0.) {
    (*state_vars)["phi"] =
        phi_residual_ + ((phi_peak_ - phi_residual_) * (epds - epds_residual_) /
                         (epds_peak_ - epds_residual_));
    (*state_vars)["psi"] =
        psi_residual_ + ((psi_peak_ - psi_residual_) * (epds - epds_residual_) /
                         (epds_peak_ - epds_residual_));
    (*state_vars)["cohesion"] =
        cohesion_residual_ +
        ((cohesion_peak_ - cohesion_residual_) * (epds - epds_residual_) /
         (epds_peak_ - epds_residual_));
  }
  if (softening_ && (epds - epds_residual_) >= 0.) {
    (*state_vars)["phi"] = phi_residual_;
    (*state_vars)["psi"] = psi_residual_;
    (*state_vars)["cohesion"] = cohesion_residual_;
  }
  //-------------------------------------------------------------------------
  // Elastic-predictor stage: compute the trial stress
  Vector6d trial_stress = stress + (this->de_ * dstrain);
  // Compute stress invariants based on trial stress
  this->compute_stress_invariants(trial_stress, state_vars);
  // Compute yield function based on the trial stress
  Eigen::Matrix<double, 2, 1> yield_function_trial;
  auto yield_type_trial =
      this->compute_yield_state(&yield_function_trial, (*state_vars));
  // Return the updated stress in elastic state
  if (yield_type_trial == FailureState::Elastic) return trial_stress;
  //-------------------------------------------------------------------------
  // Plastic-corrector stage: correct the stress back to the yield surface
  // Define tolerance of yield function
  double Tolerance = 1E-1;
  // Compute plastic multiplier based on trial stress (Lambda trial)
  double softening_trial = 0.;
  Vector6d df_dsigma_trial = Vector6d::Zero();
  Vector6d dp_dsigma_trial = Vector6d::Zero();
  this->compute_df_dp(yield_type_trial, state_vars, trial_stress,
                      &df_dsigma_trial, &dp_dsigma_trial, &softening_trial);
  double yield_trial = 0.;
  if (yield_type_trial == FailureState::Tensile)
    yield_trial = yield_function_trial(0);
  if (yield_type_trial == FailureState::Shear)
    yield_trial = yield_function_trial(1);
  double lambda_trial =
      yield_trial /
      ((df_dsigma_trial.transpose() * de_).dot(dp_dsigma_trial.transpose()) +
       softening_trial);
  // Compute stress invariants based on stress input
  this->compute_stress_invariants(stress, state_vars);
  // Compute yield function based on stress input
  Eigen::Matrix<double, 2, 1> yield_function;
  auto yield_type = this->compute_yield_state(&yield_function, (*state_vars));
  // Initialise value of yield function based on stress
  double yield{std::numeric_limits<double>::max()};
  if (yield_type == FailureState::Tensile) yield = yield_function(0);
  if (yield_type == FailureState::Shear) yield = yield_function(1);
  // Compute plastic multiplier based on stress input (Lambda)
  double softening = 0.;
  Vector6d df_dsigma = Vector6d::Zero();
  Vector6d dp_dsigma = Vector6d::Zero();
  this->compute_df_dp(yield_type, state_vars, stress, &df_dsigma, &dp_dsigma,
                      &softening);
  double lambda =
      ((df_dsigma.transpose() * this->de_).dot(dstrain)) /
      (((df_dsigma.transpose() * this->de_).dot(dp_dsigma)) + softening);
  // Compute the correction stress
  double p_multiplier = 0.;
  Vector6d dp_dsigma_final = Vector6d::Zero();
  // Correction stress based on stress
  if (fabs(yield) < Tolerance) {
    p_multiplier = lambda;
    dp_dsigma_final = dp_dsigma;
  }
  // Correction stress based on trial stress
  else {
    p_multiplier = lambda_trial;
    dp_dsigma_final = dp_dsigma_trial;
  }
  // Correct stress back to the yield surface
  Vector6d updated_stress =
      trial_stress - (p_multiplier * this->de_ * dp_dsigma_final);

  // Check the update stress
  // Compute stress invariants based on updated stress
  this->compute_stress_invariants(updated_stress, state_vars);
  // Compute yield function based on updated stress
  yield_type_trial =
      this->compute_yield_state(&yield_function_trial, (*state_vars));
  // Define the maximum iteration step
  int itr_max = 100;
  // Initialise counter of iteration step
  int itr = 0;
  // Correct the stress again
  while ((yield_function_trial(0) > Tolerance ||
          yield_function_trial(1) > Tolerance) &&
         itr < itr_max) {
    // Compute plastic multiplier based on updated stress
    this->compute_df_dp(yield_type_trial, state_vars, updated_stress,
                        &df_dsigma_trial, &dp_dsigma_trial, &softening_trial);
    if (yield_type_trial == FailureState::Tensile)
      yield_trial = yield_function_trial(0);
    if (yield_type_trial == FailureState::Shear)
      yield_trial = yield_function_trial(1);
    // Compute plastic multiplier based on updated stress
    lambda_trial =
        yield_trial /
        ((df_dsigma_trial.transpose() * de_).dot(dp_dsigma_trial.transpose()) +
         softening_trial);
    // Correct stress back to the yield surface
    updated_stress -= (lambda_trial * this->de_ * dp_dsigma_trial);
    // Compute stress invariants based on updated stress
    this->compute_stress_invariants(updated_stress, state_vars);
    // Compute yield function based on updated stress
    yield_type_trial =
        this->compute_yield_state(&yield_function_trial, (*state_vars));
    // Count the iteration step
    itr++;
  }
  // Compute incremental of plastic strain
  Vector6d dstress = updated_stress - stress;
  Vector6d dpstrain = dstrain - (this->de_.inverse()) * dstress;
  if (Tdim == 2) dpstrain(4) = dpstrain(5) = 0.;
  // Update plastic strain
  (*state_vars).at("pstrain0") += dpstrain(0);
  (*state_vars).at("pstrain1") += dpstrain(1);
  (*state_vars).at("pstrain2") += dpstrain(2);
  (*state_vars).at("pstrain3") += dpstrain(3);
  (*state_vars).at("pstrain4") += dpstrain(4);
  (*state_vars).at("pstrain5") += dpstrain(5);
  // Update equivalent plastic deviatoric strain
  (*state_vars).at("epds") = std::sqrt(
      2. / 9. *
          (std::pow(
               ((*state_vars).at("pstrain0") - (*state_vars).at("pstrain1")),
               2.) +
           std::pow(
               ((*state_vars).at("pstrain1") - (*state_vars).at("pstrain2")),
               2.) +
           std::pow(
               ((*state_vars).at("pstrain2") - (*state_vars).at("pstrain0")),
               2.)) +
      1. / 3. *
          (std::pow(((*state_vars).at("pstrain3")), 2.) +
           std::pow(((*state_vars).at("pstrain4")), 2.) +
           std::pow(((*state_vars).at("pstrain5")), 2.)));

  return updated_stress;
}
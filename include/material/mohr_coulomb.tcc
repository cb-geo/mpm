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
    friction_angle_ = material_properties["friction"].template get<double>();
    dilation_angle_ = material_properties["dilation"].template get<double>();
    cohesion_ = material_properties["cohesion"].template get<double>();
    residual_friction_angle_ =
        material_properties["residual_friction"].template get<double>();
    residual_dilation_angle_ =
        material_properties["residual_dilation"].template get<double>();
    residual_cohesion_ =
        material_properties["residual_cohesion"].template get<double>();
    peak_epds_ = material_properties["peak_epds"].template get<double>();
    crit_epds_ = material_properties["critical_epds"].template get<double>();
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

//! Initialise state variables
//! equivalent plastic deviatoric strain (epds)
template <unsigned Tdim>
bool mpm::MohrCoulomb<Tdim>::initialise_state_variables(
    std::map<std::string, double>* state_vars) {
  bool status = false;
  status = state_vars->insert(std::make_pair("epds0", 0.)).second;
  if (!status) {
    (*state_vars)["epds0"] = 0.;
    status = true;
  }
  status = state_vars->insert(std::make_pair("epds1", 0.)).second;
  if (!status) {
    (*state_vars)["epds1"] = 0.;
    status = true;
  }
  status = state_vars->insert(std::make_pair("epds2", 0.)).second;
  if (!status) {
    (*state_vars)["epds2"] = 0.;
    status = true;
  }
  status = state_vars->insert(std::make_pair("epds3", 0.)).second;
  if (!status) {
    (*state_vars)["epds3"] = 0.;
    status = true;
  }
  status = state_vars->insert(std::make_pair("epds4", 0.)).second;
  if (!status) {
    (*state_vars)["epds4"] = 0.;
    status = true;
  }
  status = state_vars->insert(std::make_pair("epds5", 0.)).second;
  if (!status) {
    (*state_vars)["epds5"] = 0.;
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
bool mpm::MohrCoulomb<Tdim>::compute_rho_theta(const Vector6d stress,
                                               double& _j2, double& _j3,
                                               double& _rho, double& _theta) {
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
  // compute j2
  _j2 = (pow((stress(0) - stress(1)), 2) + pow((stress(1) - stress(2)), 2) +
         pow((stress(0) - stress(2)), 2)) /
            6.0 +
        pow(stress(3), 2);
  if (Tdim == 3) _j2 += pow(stress(4), 2) + pow(stress(5), 2);
  // compute j3
  _j3 = (dev_stress(0) * dev_stress(1) * dev_stress(2)) -
        (dev_stress(2) * pow(dev_stress(3), 2));
  if (Tdim == 3)
    _j3 += ((2 * dev_stress(3) * dev_stress(4) * dev_stress(5)) -
            (dev_stress(0) * pow(dev_stress(4), 2) -
             dev_stress(1) * pow(dev_stress(5), 2)));
  // compute theta value
  double theta_val = 0.;
  if (fabs(_j2) > 0.0) theta_val = (3. * sqrt(3.) / 2.) * (_j3 / pow(_j2, 1.5));
  if (theta_val > 0.99) theta_val = 1.0;
  if (theta_val < -0.99) theta_val = -1.0;
  _theta = (1. / 3.) * acos(theta_val);
  if (_theta > PI / 3.) _theta = PI / 3.;
  if (_theta < 0.0) _theta = 0.;

  _rho = sqrt(2 * _j2);
  return true;
}

//! Check the yield state and return the value of yield function
template <unsigned Tdim>
int mpm::MohrCoulomb<Tdim>::check_yield(
    Eigen::Matrix<double, 2, 1>& _yield_function, const double _epsilon,
    const double _rho, const double _theta) {
  double yield_tension = sqrt(2. / 3.) * cos(_theta) * _rho +
                         _epsilon / sqrt(3.) - tension_cutoff_;
  double yield_shear = sqrt(3. / 2.) * _rho *
                           ((sin(_theta + PI / 3.) / (sqrt(3.) * cos(phi_))) +
                            (cos(_theta + PI / 3.) * tan(phi_) / 3.)) +
                       (_epsilon / sqrt(3.)) * tan(phi_) - c_;
  int _yield_type = 0;
  _yield_function(0) = yield_tension;
  _yield_function(1) = yield_shear;
  // 0 for elastic, 1 for tension failure, 2 for shear failure
  if (yield_tension > 1.E-22 && yield_shear > 1.E-22) {
    // Compute the shear-tension edge
    double Nphi = (1. + sin(phi_)) / (1. - sin(phi_));
    double SigmaP = tension_cutoff_ * Nphi - 2. * c_ * sqrt(Nphi);
    double alphaP = sqrt(1. + Nphi * Nphi) + Nphi;
    double h = sqrt(2. / 3.) * cos(_theta) * _rho + _epsilon / sqrt(3.) -
               tension_cutoff_ +
               alphaP * (sqrt(2. / 3.) * cos(_theta - 4. * PI / 3.) * _rho +
                         _epsilon / sqrt(3.) - SigmaP);
    if (h > 1.E-22)
      _yield_type = 1;
    else
      _yield_type = 2;
  } else if (yield_shear > 1.E-22 && yield_tension < 1.E-22)
    _yield_type = 2;
  else if (yield_tension > 1.E-22 && yield_shear < 1.E-22)
    _yield_type = 1;
  return _yield_type;
}

//! Return dF/dSigma and dP/dSigma
template <unsigned Tdim>
bool mpm::MohrCoulomb<Tdim>::compute_df_dp(
    const int _yield_type, const double _j2, const double _j3,
    const double _rho, const double _theta, const Vector6d stress,
    Vector6d& _dF_dSigma, Vector6d& _dP_dSigma, const double _epds,
    double& _softening, const ParticleBase<Tdim>* ptr) {
  // Compute deviatoric stress
  double mean_p = (stress(0) + stress(1) + stress(2)) / 3.0;
  // if (mean_p >= 0.0)
  // mean_p = 1.0;
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
  double dF_dEpsilon, dF_dRho, dF_dTheta;
  // Values in tension yield
  if (_yield_type == 1) {
    dF_dEpsilon = 1. / sqrt(3.);
    dF_dRho = sqrt(2. / 3.) * cos(_theta);
    dF_dTheta = -sqrt(2. / 3.) * _rho * sin(_theta);
  }
  // Values in shear yield
  else {
    dF_dEpsilon = tan(phi_) / sqrt(3.);
    dF_dRho =
        sqrt(3. / 2.) * ((sin(_theta + PI / 3.) / (sqrt(3.) * cos(phi_))) +
                         (cos(_theta + PI / 3.) * tan(phi_) / 3.));
    dF_dTheta = sqrt(3. / 2.) * _rho *
                ((cos(_theta + PI / 3.) / (sqrt(3.) * cos(phi_))) -
                 (sin(_theta + PI / 3.) * tan(phi_) / 3.));
  }

  Vector6d dEpsilon_dSigma, dRho_dSigma, dTheta_dSigma;
  dEpsilon_dSigma = dRho_dSigma = dTheta_dSigma = Vector6d::Zero();
  // Compute dEpsilon / dSigma (the same in two yield types)
  dEpsilon_dSigma(0) = dEpsilon_dSigma(1) = dEpsilon_dSigma(2) = 1. / sqrt(3.);
  // Compute dRho / dSigma (the same in two yield types)
  double multiplier = 1.;
  if (fabs(_rho) > 0.) multiplier = 1. / _rho;
  dRho_dSigma = multiplier * dev_stress;
  if (Tdim == 2) dRho_dSigma(4) = dRho_dSigma(5) = 0.;
  // Compute dTheta / dSigma (the same in two yield types)
  double r_val = 0.;
  if (fabs(_j2) > 1.E-22) r_val = (3. * sqrt(3.) / 2.) * (_j3 / pow(_j2, 1.5));
  double devider = 1 - (r_val * r_val);
  if (devider <= 0.) devider = 0.001;
  double dTheta_dR = -1 / (3. * sqrt(devider));
  double dR_dJ2 = (-9 * sqrt(3.) / 4.) * _j3;
  if (fabs(_j2) > 1.E-22) dR_dJ2 = dR_dJ2 / pow(_j2, 2.5);
  double dR_dJ3 = 1.5 * sqrt(3.);
  if (fabs(_j2) > 1.E-22) dR_dJ3 = dR_dJ3 / pow(_j2, 1.5);
  Vector6d dJ2_dSigma = dev_stress;
  dJ2_dSigma(3) = dev_stress(3);
  Vector6d dJ3_dSigma = Vector6d::Zero();
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
  dJ3_dSigma(0) = dev1.dot(dev1) - (2. / 3.) * _j2;
  dJ3_dSigma(1) = dev2.dot(dev2) - (2. / 3.) * _j2;
  dJ3_dSigma(2) = dev3.dot(dev3) - (2. / 3.) * _j2;
  dJ3_dSigma(3) = dev1.dot(dev2);
  if (Tdim == 3) {
    dJ3_dSigma(4) = dev2.dot(dev3);
    dJ3_dSigma(5) = dev1.dot(dev3);
  }
  dTheta_dSigma = dTheta_dR * ((dR_dJ2 * dJ2_dSigma) + (dR_dJ3 * dJ3_dSigma));
  if (Tdim == 2) dTheta_dSigma(4) = dTheta_dSigma(5) = 0.;

  // Compute dF/dSigma
  _dF_dSigma = (dF_dEpsilon * dEpsilon_dSigma) + (dF_dRho * dRho_dSigma) +
               (dF_dTheta * dTheta_dSigma);
  if (Tdim == 2) _dF_dSigma(4) = _dF_dSigma(5) = 0.;

  // Compute dP/dSigma, dP/dJ
  double dP_dJ = 0.;
  if (_yield_type == 1) {
    double et_value = 0.6;
    double xit = 0.1;
    double Rt_den =
        2. * (1 - et_value * et_value) * cos(_theta) +
        (2. * et_value - 1) *
            sqrt(4. * (1 - et_value * et_value) * cos(_theta) * cos(_theta) +
                 5. * et_value * et_value - 4. * et_value);
    double Rt_num = 4. * (1 - et_value * et_value) * cos(_theta) * cos(_theta) +
                    (2. * et_value - 1) * (2. * et_value - 1);
    double Rt = Rt_num / (3. * Rt_den);
    double dP_dRt = 1.5 * _rho * _rho * Rt /
                    sqrt(xit * xit * tension_cutoff_ * tension_cutoff_ +
                         1.5 * Rt * Rt * _rho * _rho);
    double dP_dRho = 1.5 * _rho * Rt * Rt /
                     sqrt(xit * xit * tension_cutoff_ * tension_cutoff_ +
                          1.5 * Rt * Rt * _rho * _rho);
    double dP_dEpsilon = 1. / sqrt(3.);
    double dRtden_dTheta =
        -2. * (1 - et_value * et_value) * sin(_theta) -
        (2. * et_value - 1) * 4. * (1 - et_value * et_value) * cos(_theta) *
            sin(_theta) /
            sqrt(4. * (1 - et_value * et_value) * cos(_theta) * cos(_theta) +
                 5. * et_value * et_value - 4. * et_value);
    double dRtnum_dTheta =
        -8. * (1 - et_value * et_value) * cos(_theta) * sin(_theta);
    double dRt_dTheta = (dRtnum_dTheta * Rt_den - dRtden_dTheta * Rt_num) /
                        (3. * Rt_den * Rt_den);
    // Compute the value of dP/dSigma and dP/dJ in tension yield
    _dP_dSigma = (dP_dEpsilon * dEpsilon_dSigma) + (dP_dRho * dRho_dSigma) +
                 (dP_dRt * dRt_dTheta * dTheta_dSigma);
    dP_dJ = dP_dRho * sqrt(2.);
  } else {
    double R_mc = (3. - sin(phi_)) / (6 * cos(phi_));
    double e_val = (3. - sin(phi_)) / (3. + sin(phi_));
    if ((e_val - 0.5) < 0.) e_val = 0.501;
    if ((e_val - 1.) > 0.) e_val = 1.0;
    double sqpart = (4. * (1 - e_val * e_val) * pow(cos(_theta), 2)) +
                    (5 * e_val * e_val) - (4. * e_val);
    if (sqpart < 0.) sqpart = 0.00001;
    double R_mw_den = (2. * (1 - e_val * e_val) * cos(_theta)) +
                      ((2. * e_val - 1) * sqrt(sqpart));
    if (fabs(R_mw_den) < 1.E-22) R_mw_den = 0.001;
    double R_mw_num = (4. * (1. - e_val * e_val) * pow(cos(_theta), 2)) +
                      pow((2. * e_val - 1.), 2);
    double R_mw = (R_mw_num / R_mw_den) * R_mc;
    double xi = 0.1;
    double omega =
        pow((xi * c_ * tan(psi_)), 2) + pow((R_mw * sqrt(3. / 2.) * _rho), 2);
    if (omega < 1.E-22) omega = 0.001;
    double L = R_mw_num;
    double M = R_mw_den;
    double dL_dTheta = -8. * (1. - e_val * e_val) * cos(_theta) * sin(_theta);
    double dM_dTheta = (-2. * (1. - e_val * e_val) * sin(_theta)) +
                       (0.5 * (2. * e_val - 1.) * dL_dTheta) / sqrt(sqpart);
    double dRmw_dTheta = ((M * dL_dTheta) - (L * dM_dTheta)) / (M * M);
    double dP_dEpsilon = tan(psi_) / sqrt(3.);
    double dP_dRho = 3. * _rho * R_mw * R_mw / (2. * sqrt(omega));
    double dP_dTheta =
        (3. * _rho * _rho * R_mw * R_mc * dRmw_dTheta) / (2. * sqrt(omega));
    // Compute the value of dP/dSigma and dP/dJ in shear yield
    _dP_dSigma = (dP_dEpsilon * dEpsilon_dSigma) + (dP_dRho * dRho_dSigma) +
                 (dP_dTheta * dTheta_dSigma);
    dP_dJ = dP_dRho * sqrt(2.);
  }

  // compute softening part
  double dPhi_dPstrain = 0.;
  double dC_dPstrain = 0.;
  if (_yield_type == 2 && _epds > peak_epds_ && _epds < crit_epds_) {
    dPhi_dPstrain = (residual_friction_angle_ - friction_angle_) * PI / 180. /
                    (crit_epds_ - peak_epds_);
    dC_dPstrain = (residual_cohesion_ - cohesion_) / (crit_epds_ - peak_epds_);
  }
  double epsilon = (1. / sqrt(3.)) * (stress(0) + stress(1) + stress(2));
  double dF_dPhi =
      sqrt(3. / 2.) * _rho *
          ((sin(phi_) * sin(_theta + PI / 3.) /
            (sqrt(3.) * cos(phi_) * cos(phi_))) +
           (cos(_theta + PI / 3.) / (3. * cos(phi_) * cos(phi_)))) +
      (epsilon / (sqrt(3.) * cos(phi_) * cos(phi_)));
  double dF_dC = -1.;
  _softening =
      (-1.) * ((dF_dPhi * dPhi_dPstrain) + (dF_dC * dC_dPstrain)) * dP_dJ;
  return true;
}

//! Compute stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::MohrCoulomb<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr, std::map<std::string, double>* state_vars) {
  // Friction and dilation in radians
  const double phi_max = friction_angle_ * PI / 180.;
  const double psi_max = dilation_angle_ * PI / 180.;
  const double c_max = cohesion_;
  const double phi_min = residual_friction_angle_ * PI / 180.;
  const double psi_min = residual_dilation_angle_ * PI / 180.;
  const double c_min = residual_cohesion_;

  // Current MC parameters using a linear softening rule
  // plastic deviatoric strain
  Eigen::Matrix<double, 6, 1> PDS;
  PDS(0) = (*state_vars).at("epds0");
  PDS(1) = (*state_vars).at("epds1");
  PDS(2) = (*state_vars).at("epds2");
  PDS(3) = (*state_vars).at("epds3");
  PDS(4) = (*state_vars).at("epds4");
  PDS(5) = (*state_vars).at("epds5");

  const double epds =
      (2. / 3.) *
      sqrt(3. / 2. * (PDS(0) * PDS(0) + PDS(1) * PDS(1) + PDS(2) * PDS(2)) +
           3. * (PDS(3) * PDS(3) + PDS(4) * PDS(4) + PDS(5) * PDS(5)));
  double j2, j3, rho, theta;
  Vector6d trial_stress = stress + (this->de_ * dstrain);
  this->compute_rho_theta(stress, j2, j3, rho, theta);
  double epsilon = (1. / sqrt(3.)) * (stress(0) + stress(1) + stress(2));
  Eigen::Matrix<double, 2, 1> yield_function, yield_function_trial;
  phi_ = phi_max;
  psi_ = psi_max;
  c_ = c_max;
  int yield_type = this->check_yield(yield_function, epsilon, rho, theta);
  if (yield_type != 1 && (epds - peak_epds_) > 0. && (crit_epds_ - epds) > 0.) {
    phi_ = phi_min + ((phi_max - phi_min) * (epds - crit_epds_) /
                      (peak_epds_ - crit_epds_));
    psi_ = psi_min + ((psi_max - psi_min) * (epds - crit_epds_) /
                      (peak_epds_ - crit_epds_));
    c_ = c_min +
         ((c_max - c_min) * (epds - crit_epds_) / (peak_epds_ - crit_epds_));
  } else if (yield_type != 1 && (epds - crit_epds_) >= 0.) {
    phi_ = phi_min;
    psi_ = psi_min;
    c_ = c_min;
  }
  // Yield function for the current stress state
  yield_type = this->check_yield(yield_function, epsilon, rho, theta);
  int yield_type_trial = 0;
  double lambda = 0.;
  double lambda_trial = 0.;
  Vector6d dF_dSigma = Vector6d::Zero();
  Vector6d dP_dSigma = Vector6d::Zero();
  Vector6d dF_dSigma_trial = Vector6d::Zero();
  Vector6d dP_dSigma_trial = Vector6d::Zero();
  double softening = 0.;
  double softening_trial = 0.;
  // Compute plastic multiplier from the current stress state
  this->compute_df_dp(yield_type, j2, j3, rho, theta, stress, dF_dSigma,
                      dP_dSigma, epds, softening, ptr);
  // Check the epds
  if (epds_last < peak_epds_ && epds > peak_epds_) softening = 0;
  lambda = dF_dSigma.dot(this->de_ * dstrain) /
           ((dF_dSigma.dot(this->de_ * dP_dSigma)) + softening);
  //}
  // Compute the trial stress
  j2 = j3 = rho = theta = 0.;
  this->compute_rho_theta(trial_stress, j2, j3, rho, theta);
  epsilon =
      (1. / sqrt(3.)) * (trial_stress(0) + trial_stress(1) + trial_stress(2));
  yield_type_trial =
      this->check_yield(yield_function_trial, epsilon, rho, theta);
  this->compute_df_dp(yield_type_trial, j2, j3, rho, theta, trial_stress,
                      dF_dSigma_trial, dP_dSigma_trial, epds, softening_trial,
                      ptr);
  // Check the epds
  if (epds_last < peak_epds_ && epds > peak_epds_) softening_trial = 0;
  if (yield_type_trial == 1)
    lambda_trial =
        yield_function_trial(0) /
        ((dF_dSigma_trial.transpose() * de_).dot(dP_dSigma_trial.transpose()) +
         softening_trial);
  else
    lambda_trial =
        yield_function_trial(1) /
        ((dF_dSigma_trial.transpose() * de_).dot(dP_dSigma_trial.transpose()) +
         softening_trial);
  // Compute the correction stress
  double p_multiplier = 0.;
  Vector6d dP_dSigma_Final = Vector6d::Zero();
  bool PDS_update = false;
  // check if it is a load process or an unload process
  double dyfun_t = yield_function_trial(0) - yield_function(0);
  double dyfun_s = yield_function_trial(1) - yield_function(1);
  if (yield_type != 0) {
    if (yield_type == 1 && dyfun_t > 0) {
      p_multiplier = lambda;
      dP_dSigma_Final = dP_dSigma;
      PDS_update = true;
    } else if (yield_type == 2 && dyfun_s > 0) {
      p_multiplier = lambda;
      dP_dSigma_Final = dP_dSigma;
      PDS_update = true;
    }
  } else if (yield_type_trial != 0) {
    p_multiplier = lambda_trial;
    dP_dSigma_Final = dP_dSigma_trial;
    PDS_update = true;
  }
  // update stress (plastic correction)
  Vector6d stress_update =
      trial_stress - (p_multiplier * this->de_ * dP_dSigma_Final);
  // compute plastic deviatoric strain
  Vector6d dstress = stress_update - stress;
  Vector6d dpstrain = dstrain - (this->de_.inverse()) * dstress;
  if (Tdim == 2) dpstrain(4) = dpstrain(5) = 0.;
  // Record the epds
  epds_last = epds;
  // Update PDS if it is a load process
  if (PDS_update) {
    Vector6d dPDS = dpstrain;
    double dpvstrain = p_multiplier * (dP_dSigma_Final(0) + dP_dSigma_Final(1) +
                                       dP_dSigma_Final(2));
    dPDS(0) -= (1. / 3. * dpvstrain);
    dPDS(1) -= (1. / 3. * dpvstrain);
    dPDS(2) -= (1. / 3. * dpvstrain);
    dPDS(3) = 0.5 * dpstrain(3);
    dPDS(4) = 0.5 * dpstrain(4);
    dPDS(5) = 0.5 * dpstrain(5);
    PDS += dPDS;
    (*state_vars).at("epds0") = PDS(0);
    (*state_vars).at("epds1") = PDS(1);
    (*state_vars).at("epds2") = PDS(2);
    (*state_vars).at("epds3") = PDS(3);
    (*state_vars).at("epds4") = PDS(4);
    (*state_vars).at("epds5") = PDS(5);
  }
  return stress_update;
}

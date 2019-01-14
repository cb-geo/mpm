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
    friction_angle_ =
      material_properties["friction_angle"].template get<double>();
    dilation_angle_ =
      material_properties["dilation_angle"].template get<double>();
    cohesion_ = material_properties["cohesion"].template get<double>();
    residual_friction_angle_ =
      material_properties["residual_friction_angle"].template get<double>();
    residual_dilation_angle_ =
      material_properties["residual_dilation_angle"].template get<double>();
    residual_cohesion_ =
      material_properties["residual_cohesion"].template get<double>();
    peak_EPDS_ =
      material_properties["peak_EPDS"].template get<double>();
    crit_EPDS_ =
      material_properties["crit_EPDS"].template get<double>();
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
bool mpm::MohrCoulomb<Tdim>::compute_rho_theta(const Vector6d stress,
					       double& _j2, double& _j3,
					       double& _rho, double& _theta) {
  const double ONETHIRDPI = 1.047197551;
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
  _j2 = (pow((stress(0) - stress(1)),2) + pow((stress(1) - stress(2)),2) + pow((stress(0) - stress(2)),2)) / 6.0 + pow(stress(3),2);
  if(Tdim == 3)
    _j2 += pow(stress(4),2) + pow(stress(5),2);
  // compute j3
  _j3 = (dev_stress(0) * dev_stress(1) * dev_stress(2))
    - (dev_stress(2) * pow(dev_stress(3),2));
  if(Tdim == 3)
    _j3 += ((2 * dev_stress(3)*dev_stress(4)*dev_stress(5)) - (dev_stress(0) * pow(dev_stress(4),2) - dev_stress(1)*pow(dev_stress(5),2)));

  // compute theta value
  double theta_val = 0.;
  if(fabs(_j2) > 0.0)
    theta_val = (3. * sqrt(3.) / 2.) * (_j3 / pow(_j2,1.5));


  if(theta_val > 0.99)
    theta_val = 1.0;
  if(theta_val < -0.99)
    theta_val= -1.0;
  _theta = (1./3.) * acos(theta_val);
  if(_theta > ONETHIRDPI) _theta = ONETHIRDPI;
  if(_theta < 0.0) _theta = 0.;
  
  _rho = sqrt(2 * _j2);
  return true;
}

//! Return dF/dSigma and dP/dSigma
template <unsigned Tdim>
bool mpm::MohrCoulomb<Tdim>::compute_df_dp(const double _j2, const double _j3,
					   const double _rho,
					   const double _theta,
					   const Vector6d stress,
					   Vector6d& _dF_dSigma,
					   Vector6d& _dP_dSigma,
					   double& _softening) {
  const double ONETHIRDPI = 1.047197551;
  // Deviatoric stress
  double mean_p = (stress(0) + stress(1) + stress(2)) / 3.0;
  if (mean_p >= 0.0)
    mean_p = 1.0;
  Vector6d dev_stress = Vector6d::Zero();
  dev_stress(0) = stress(0) - mean_p;
  dev_stress(1) = stress(1) - mean_p;
  dev_stress(2) = stress(2) - mean_p;
  dev_stress(3) = stress(3);
  if(Tdim == 3) {
    dev_stress(4) = stress(4);
    dev_stress(5) = stress(5);
  }
 
  // compute dF / dEpsilon,  dF / dRho, dF / dTheta
  double dF_dEpsilon = tan(phi_) / sqrt(3.);
  double dF_dRho = sqrt(3. / 2.) * ( (sin(_theta + ONETHIRDPI) / (sqrt(3.)*cos(phi_))) + (cos(_theta + ONETHIRDPI) * tan(phi_) / 3.) );
  double dF_dTheta = sqrt(3. / 2.) * _rho *( (cos(_theta + ONETHIRDPI) / (sqrt(3.)*cos(phi_))) - (sin(_theta + ONETHIRDPI) * tan(phi_) / 3.) );
  
  Vector6d dEpsilon_dSigma, dRho_dSigma, dTheta_dSigma;
  dEpsilon_dSigma = dRho_dSigma = dTheta_dSigma = Vector6d::Zero();
  // compute dEpsilon / dSigma
  dEpsilon_dSigma(0) = dEpsilon_dSigma(1) = dEpsilon_dSigma(2) = 1./sqrt(3.);
  // compute dRho / dSigma
  double multiplier = 1.;
  if (fabs(_rho) > 0.)
    multiplier = 1. / _rho;
  dRho_dSigma = multiplier * dev_stress;
  if(Tdim == 2)
    dRho_dSigma(4) = dRho_dSigma(5) = 0.;
  // compute dTheta / dSigma
  double r_val = 0.;
  if(fabs(_j2) > 1.E-22)
    r_val = (3. * sqrt(3.) / 2.) * (_j3 / pow(_j2,1.5));
  double devider = 1 - (r_val *r_val);
  if(devider <= 0.) devider = 0.001;
  double dTheta_dR = -1 / (3. * sqrt(devider));
  double dR_dJ2 = (-9*sqrt(3.) / 4.) * _j3;
  if(fabs(_j2) > 1.E-22)
    dR_dJ2 = dR_dJ2 / pow(_j2, 2.5);
  double dR_dJ3 = 1.5 * sqrt(3.);
  if(fabs(_j2) > 1.E-22)
    dR_dJ3 = dR_dJ3 / pow(_j2, 1.5);

  Vector6d dJ2_dSigma = dev_stress;
  Vector6d dJ3_dSigma = Vector6d::Zero();
  Eigen::Matrix<double,3,1> dev1, dev2, dev3;
  dev1(0) = dev_stress(0); dev1(1) = dev_stress(3); dev1(2) = dev_stress(5);
  dev2(0) = dev_stress(3); dev2(1) = dev_stress(1); dev2(2) = dev_stress(4);
  dev3(0) = dev_stress(5); dev3(1) = dev_stress(4); dev3(2) = dev_stress(2);
  dJ3_dSigma(0) = dev1.dot(dev1) - (2. / 3.) * _j2;
  dJ3_dSigma(1) = dev2.dot(dev2) - (2. / 3.) * _j2;
  dJ3_dSigma(2) = dev3.dot(dev3) - (2. / 3.) * _j2;
  dJ3_dSigma(3) = dev1.dot(dev2);
  if(Tdim == 3) {
    dJ3_dSigma(4) = dev2.dot(dev3);
    dJ3_dSigma(5) = dev1.dot(dev3);
  }
  dTheta_dSigma = dTheta_dR * ((dR_dJ2 * dJ2_dSigma)+(dR_dJ3 * dJ3_dSigma));
  if(Tdim == 2)
    dTheta_dSigma(4) = dTheta_dSigma(5) = 0.;

  _dF_dSigma = (dF_dEpsilon * dEpsilon_dSigma) + (dF_dRho * dRho_dSigma) + (dF_dTheta * dTheta_dSigma);
  if(Tdim == 2)
    _dF_dSigma(4) = _dF_dSigma(5) = 0.;

  // compute dP/dSigma
  double R_mc = (3. - sin(phi_)) / (6 * cos(phi_));
  double e_val = (3. - sin(phi_)) / (3. + sin(phi_));
  if((e_val - 0.5) < 0.) e_val = 0.501;
  if((e_val - 1.) > 0.) e_val = 1.0;
  double sqpart = (4.*(1 - e_val*e_val)*pow(cos(_theta),2))
    + (5*e_val*e_val) - (4.*e_val);
  if(sqpart < 0.) sqpart = 0.00001;
  double R_mw_den = (2.*(1-e_val*e_val)*cos(_theta))
    + ((2.*e_val - 1)*sqrt(sqpart));
  if(fabs(R_mw_den) < 1.E-22) R_mw_den = 0.001;
  double R_mw_num = (4. * (1. -e_val *e_val) * pow(cos(_theta),2))
    + pow((2. * e_val - 1.),2);
  double R_mw = (R_mw_num / R_mw_den) * R_mc;

  double xi = 0.1;
  double omega = pow((xi * c_ * tan(psi_)),2) + pow((R_mw * sqrt(3./2.)*_rho),2);
  if(omega < 1.E-22) omega = 0.001;

  double L = R_mw_num;
  double M = R_mw_den;
  double dL_dTheta = -8. * (1. - e_val * e_val)*cos(_theta)*sin(_theta);
  double dM_dTheta = (-2. * (1. - e_val * e_val)*sin(_theta))
    + (0.5 * (2. *e_val - 1.) *dL_dTheta)/sqrt(sqpart);
  double dRmw_dTheta = ((M * dL_dTheta)-(L * dM_dTheta)) / (M * M);

  double dP_dEpsilon = tan(psi_) / sqrt(3.);
  double dP_dRho = 3. * _rho * R_mw * R_mw / (2.*sqrt(omega));
  double dP_dTheta = (3.*_rho * _rho * R_mw * R_mc * dRmw_dTheta) / (2.*sqrt(omega));

  _dP_dSigma = (dP_dEpsilon * dEpsilon_dSigma)
    + (dP_dRho * dRho_dSigma) + (dP_dTheta * dTheta_dSigma);

  // compute softening part
  double dPhi_dPstrain = 0.;
  double dC_dPstrain = 0.;
  // if(epds_ > epds_peak_ && epds_ < epds_crit_) {
  //   _dPhidPstrain = (phi_resd_ - phi_) / (epds_crit_ - epds_peak_);
  //   _dCdPstrain = (c_resd_ - c_) / (epds_crit_ - epds_peak_);
  // }
  // double  epsilon = (1. / sqrt(3.)) * (stress(0) + stress(1) + stress(2));
  // double dF_dPhi = sqrt(3./2.) * _rho * ((sin(phi_) * sin(_theta + ONETHIRDPI) / (sqrt(3.)*cos(phi_) * cos(phi_))) + (cos(_theta +ONETHIRDPI)/(3.*cos(phi_)*cos(phi_))) ) + (epsilon / (sqrt(3.)*cos(phi_)*cos(phi_)));

  // double dF_dC = -1.;
  // _softening = (-1.) * ((dF_dPhi*dPhi_dPstrain) + (dF_dC*dC_dPstrain)) * dP_dRho;

  return true;
}

//! Compute stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::MohrCoulomb<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr) {

  const double PI = std::atan(1.0) * 4.;
  const double ONETHIRDPI = PI / 3.;

  // Friction and dilation in radians
  const double phi_max = friction_angle_ * PI / 180.;
  const double psi_max = dilation_angle_ * PI / 180.;
  const double c_max   = cohesion_;
  const double phi_min = residual_friction_angle_ * PI / 180.;
  const double psi_min = residual_dilation_angle_ * PI / 180.;
  const double c_min   = residual_cohesion_;

  // Current MC parameters using a linear softening rule
  double epds = 0;
  if((peak_EPDS_ - epds) >= 0.) {
    phi_ = phi_max; psi_ = psi_max; c_ = c_max;
  }
  else if ((epds - peak_EPDS_) > 0. && (crit_EPDS_ - epds) > 0.) {
    phi_ = phi_min
      + ((phi_max - phi_min) * (epds - crit_EPDS_) / (peak_EPDS_ - crit_EPDS_));
    psi_ = psi_min
      + ((psi_max - psi_min) * (epds - crit_EPDS_) / (peak_EPDS_ - crit_EPDS_));
    c_ = c_min
      + ((c_max - c_min) * (epds - crit_EPDS_) / (peak_EPDS_ - crit_EPDS_));
  }
  else if ((epds - crit_EPDS_) >= 0.) {
    phi_ = phi_max;
    psi_ = psi_max;
    c_ = c_max;
  }

  // Yield function for the current stress state
  double j2, j3, rho, theta;
  this->compute_rho_theta(stress, j2, j3, rho, theta);
  double epsilon = (1. / sqrt(3.)) * (stress(0) + stress(1) + stress(2));
  double yield_func = sqrt(3. / 2.) * rho * ( (sin(theta + ONETHIRDPI) / (sqrt(3.)*cos(phi_))) + (cos(theta + ONETHIRDPI) * tan(phi_) / 3.) ) + (epsilon / 3.)*tan(phi_) - c_;
  bool yield_state;
  if(yield_func > 1.E-22)
    yield_state = true;
  else
    yield_state = false;

  // compute plastic multiplier from the current stress state
  Vector6d dF_dSigma, dP_dSigma;
  double softening = 0;
  this->compute_df_dp(j2,j3,rho,theta,stress,dF_dSigma, dP_dSigma, softening);
  double lambda = dF_dSigma.dot(this->de_ * dstrain) / ((dF_dSigma.dot(this->de_ * dP_dSigma)) + softening);
  if(!yield_state)
    lambda = 0.;

  // compute the trial stress
  // sigma_trial = sigma + De * dstrain
  j2 = j3 = rho = theta = 0.;
  Vector6d trial_stress = stress + (this->de_ * dstrain);
  this->compute_rho_theta(trial_stress, j2, j3, rho, theta);
  epsilon = (1. / sqrt(3.)) * (trial_stress(0) + trial_stress(1) + trial_stress(2));
  double yield_func_trial = sqrt(3. / 2.) * rho * ( (sin(theta + ONETHIRDPI) / (sqrt(3.)*cos(phi_))) + (cos(theta + ONETHIRDPI) * tan(phi_) / 3.) ) + (epsilon / 3.)*tan(phi_) - c_;

  bool yield_state_trial;
  if(yield_func_trial > 1.E-22)
    yield_state_trial = true;
  else
    yield_state_trial = false;

  Vector6d dF_dSigma_trial, dP_dSigma_trial;
  double softening_trial = 0;
  this->compute_df_dp(j2,j3,rho,theta,trial_stress,dF_dSigma_trial, dP_dSigma_trial,softening_trial);
  double lambda_trial = yield_func_trial / ((dF_dSigma_trial.transpose() * de_).dot(dP_dSigma_trial.transpose()) + softening_trial);

  double p_multiplier;
  if(yield_state)
    p_multiplier = lambda;
  else if (!yield_state) {
    if(yield_state_trial)
      p_multiplier = lambda_trial;
    else if (!yield_state_trial)
      p_multiplier = 0.;
  }

  // update stress (plastic correction)
  Vector6d stress_update = trial_stress -  (p_multiplier * de_ * dP_dSigma);
  // compute plastic deviatoric strain
  Vector6d dstress = stress - stress_update;
  Vector6d dpstrain = dstrain - (de_.inverse()) * dstress;
  if(Tdim == 2)
    dpstrain(4) = dpstrain(5) = 0.;
  //PDS_ += dPstrain;

  // compute equivalent plastic deviatoric strain
  //double epds_inc = (2./3.)*sqrt(0.5*((dPstrain(0)-dPstrain(1))*(dPstrain(0)-dPstrain(1)) + (dPstrain(1)-dPstrain(2))*(dPstrain(1)-dPstrain(2)) + (dPstrain(0)-dPstrain(2))*(dPstrain(0)-dPstrain(2))) + 3.*((dPstrain(3)*dPstrain(3)) + (dPstrain(4)*dPstrain(4)) + (dPstrain(5)*dPstrain(5))));
  //epds_ += epds_inc;

  return stress_update;
}

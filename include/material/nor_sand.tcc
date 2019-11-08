//! Constructor with id and material properties
template <unsigned Tdim>
mpm::NorSand<Tdim>::NorSand(unsigned id, const Json& material_properties)
    : Material<Tdim>(id, material_properties) {
  try {
    // Density
    density_ = material_properties["density"].template get<double>();
    // Youngs modulus E
    youngs_modulus_ =
        material_properties["youngs_modulus"].template get<double>();
    // Shear modulus constant An
    shear_modulus_constant_ =
        material_properties["shear_modulus_constant"].template get<double>();
    // Shear modulus exponent Gn
    shear_modulus_exponent_ =
        material_properties["shear_modulus_exponent"].template get<double>();
    // Reference pressure pref
    reference_pressure_ =
        material_properties["reference_pressure"].template get<double>();
    // Poisson ratio
    poisson_ratio_ =
        material_properties["poisson_ratio"].template get<double>();
    // Critical state friction angle
    friction_cs_ =
        material_properties["friction_cs"].template get<double>() * M_PI / 180.;
    // Volumetric coupling (dilatancy) parameter N
    N_ = material_properties["N"].template get<double>();
    // Minimum void ratio
    e_min_ = material_properties["e_min"].template get<double>();
    // Maximum void ratio
    e_max_ = material_properties["e_max"].template get<double>();
    // Crushing pressure
    crushing_pressure_ =
        material_properties["crushing_pressure"].template get<double>();
    // Dilatancy coefficient chi
    chi_ = material_properties["chi"].template get<double>();
    // Hardening modulus
    hardening_modulus_ =
        material_properties["hardening_modulus"].template get<double>();
    // Initial void ratio
    void_ratio_initial_ =
        material_properties["void_ratio_initial"].template get<double>();
    // Initial image pressure
    p_image_initial_ =
        material_properties["p_image_initial"].template get<double>();
    // Flag for bonded model
    bond_model_ = material_properties["bond_model"].template get<bool>();
    // Initial p_cohesion
    p_cohesion_initial_ =
        material_properties["p_cohesion_initial"].template get<double>();
    // Degradation parameter m upon shearing
    m_ = material_properties["m"].template get<double>();

    double sin_friction_cs = sin(friction_cs_);
    Mtc_ = (6 * sin_friction_cs) / (3 - sin_friction_cs);
    Mte_ = (6 * sin_friction_cs) / (3 + sin_friction_cs);

    // Properties
    properties_ = material_properties;

  } catch (std::exception& except) {
    console_->error("Material parameter not set: {}\n", except.what());
  }
}

//! Initialise state variables
template <unsigned Tdim>
mpm::dense_map mpm::NorSand<Tdim>::initialise_state_variables() {
  mpm::dense_map state_vars = {
      // Mean stress
      {"p", 0.},
      // Deviatoric stress
      {"q", 0.},
      // Lode angle
      {"lode_angle", 0.},
      // J2
      {"j2", 0.},
      // J3
      {"j3", 0.},
      // M_theta
      {"M_theta", Mtc_},
      // Current void ratio
      {"void_ratio", void_ratio_initial_},
      // Void ratio image
      {"e_image",
       e_max_ - (e_max_ - e_min_) / log(crushing_pressure_ / p_image_initial_)},
      // Image pressure
      {"p_image", p_image_initial_},
      // State variable psi
      {"psi_image", void_ratio_initial_ -
                        (e_max_ - (e_max_ - e_min_) / log(crushing_pressure_ /
                                                          p_image_initial_))},
      // p_cohesion
      {"p_cohesion", p_cohesion_initial_},
      // p_cohesion degradation function
      {"zeta", 1.},
      // Equivalent plastic deviatoric strain
      {"epds", 0.},
      // Plastic strain components
      {"plastic_strain0", 0.},
      {"plastic_strain1", 0.},
      {"plastic_strain2", 0.},
      {"plastic_strain3", 0.},
      {"plastic_strain4", 0.},
      {"plastic_strain5", 0.}};

  return state_vars;
}

//! Compute elastic tensor
template <unsigned Tdim>
bool mpm::NorSand<Tdim>::compute_elastic_tensor() {
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
bool mpm::NorSand<Tdim>::compute_stress_invariants(const Vector6d& stress,
                                                   mpm::dense_map* state_vars) {

  // Note that in this subroutine, stress is compression positive

  // Compute mean stress p
  double mean_p = (stress(0) + stress(1) + stress(2)) / 3.;
  if (mean_p < 1.0E-15) mean_p = 1.0E-15;
  (*state_vars).at("p") = mean_p;

  // Compute the deviatoric stress
  Vector6d dev_stress = Vector6d::Zero();
  dev_stress(0) = stress(0) - mean_p;
  dev_stress(1) = stress(1) - mean_p;
  dev_stress(2) = stress(2) - mean_p;
  dev_stress(3) = stress(3);
  if (Tdim == 3) {
    dev_stress(4) = stress(4);
    dev_stress(5) = stress(5);
  }

  // Compute J2
  double j2 =
      (pow((stress(0) - stress(1)), 2) + pow((stress(1) - stress(2)), 2) +
       pow((stress(0) - stress(2)), 2)) /
          6.0 +
      pow(stress(3), 2);
  if (Tdim == 3) j2 += pow(stress(4), 2) + pow(stress(5), 2);
  if (fabs(j2) < 1.0E-15) j2 = 1.0E-15;
  (*state_vars).at("j2") = j2;

  // Compute J3
  double j3 = (dev_stress(0) * dev_stress(1) * dev_stress(2)) -
              (dev_stress(2) * pow(dev_stress(3), 2));
  if (Tdim == 3)
    j3 += ((2 * dev_stress(3) * dev_stress(4) * dev_stress(5)) -
           (dev_stress(0) * pow(dev_stress(4), 2)) -
           (dev_stress(1) * pow(dev_stress(5), 2)));
  (*state_vars).at("j3") = j3;

  // Compute q
  double deviatoric_q = sqrt(3 * j2);
  if (deviatoric_q < 1.0E-15) deviatoric_q = 1.0E-15;
  (*state_vars).at("q") = deviatoric_q;

  // Compute Lode angle value
  double lode_angle_val = (3. * sqrt(3.) / 2.) * (j3 / pow(j2, 1.5));
  if (lode_angle_val > 1.0) lode_angle_val = 1.0;
  if (lode_angle_val < -1.0) lode_angle_val = -1.0;

  // Compute Lode angle (sin convention)
  double lode_angle = (1. / 3.) * asin(lode_angle_val);
  if (lode_angle > M_PI / 6.) lode_angle = M_PI / 6.;
  if (lode_angle < -M_PI / 6.) lode_angle = -M_PI / 6.;
  (*state_vars)["lode_angle"] = lode_angle;

  // Compute M_theta (Jefferies and Shuttle, 2011)
  const double cos_lode_angle = cos(3. / 2. * lode_angle + M_PI / 4.);
  const double sqrt_3 = sqrt(3.);
  double M_theta = Mtc_ - pow(Mtc_, 2) / (3. + Mtc_) * cos_lode_angle;
  (*state_vars)["M_theta"] = M_theta;

  return true;
}

//! Compute state parameters
template <unsigned Tdim>
bool mpm::NorSand<Tdim>::compute_state_variables(const Vector6d& stress,
                                                 const Vector6d& dstrain,
                                                 mpm::dense_map* state_vars) {

  // Get state variables
  const double mean_p = (*state_vars).at("p");
  const double deviatoric_q = (*state_vars).at("q");
  const double M_theta = (*state_vars).at("M_theta");
  const double p_cohesion = (*state_vars).at("p_cohesion");

  // Compute and update pressure image
  double p_image;
  if (bond_model_) {
    p_image =
        (mean_p + p_cohesion) *
            pow(((1 - N_ / M_theta * deviatoric_q / (mean_p + p_cohesion)) /
                 (1 - N_)),
                ((N_ - 1) / N_)) -
        p_cohesion;
  } else {
    p_image =
        mean_p * pow(((1 - N_ / M_theta * deviatoric_q / mean_p) / (1 - N_)),
                     ((N_ - 1) / N_));
  }

  (*state_vars).at("p_image") = p_image;

  // Compute and update void ratio image
  double e_image;
  if (bond_model_) {
    e_image = e_max_ - (e_max_ - e_min_) /
                           log(crushing_pressure_ / (p_image + p_cohesion));
  } else {
    e_image = e_max_ - (e_max_ - e_min_) / log(crushing_pressure_ / p_image);
  }
  (*state_vars).at("e_image") = e_image;

  // Update void ratio
  // Note that dstrain is in tension positive - depsv = de / (1 + e)
  double dvolumetric_strain = dstrain(0) + dstrain(1) + dstrain(2);
  double void_ratio = (*state_vars).at("void_ratio") -
                      (1 + (*state_vars).at("void_ratio")) * dvolumetric_strain;
  (*state_vars).at("void_ratio") = void_ratio;

  // Compute and update psi
  double psi_image = void_ratio - e_image;
  (*state_vars).at("psi_image") = psi_image;

  return true;
}

//! Compute elastic tensor
template <unsigned Tdim>
bool mpm::NorSand<Tdim>::compute_p_cohesion(mpm::dense_map* state_vars) {

  // Compute current zeta
  double zeta = exp(-m_ * (*state_vars).at("epds"));

  if (zeta > 1.)
    (*state_vars).at("zeta") = 1.;
  else if (zeta < 1.0E-15)
    (*state_vars).at("zeta") = 0.;
  else
    (*state_vars).at("zeta") = zeta;

  // Update p_cohesion
  double p_cohesion = p_cohesion_initial_ * (*state_vars).at("zeta");
  (*state_vars).at("p_cohesion") = p_cohesion;

  return true;
}

//! Compute yield function and yield state
template <unsigned Tdim>
typename mpm::NorSand<Tdim>::FailureState
    mpm::NorSand<Tdim>::compute_yield_state(double* yield_function,
                                            const Vector6d& stress,
                                            const mpm::dense_map* state_vars) {

  // Get stress invariants
  // Compute mean stress p
  double mean_p = (stress(0) + stress(1) + stress(2)) / 3.;
  if (mean_p < 1.0E-15) mean_p = 1.0E-15;

  // Compute J2
  double j2 =
      (pow((stress(0) - stress(1)), 2) + pow((stress(1) - stress(2)), 2) +
       pow((stress(0) - stress(2)), 2)) /
          6.0 +
      pow(stress(3), 2);
  if (Tdim == 3) j2 += pow(stress(4), 2) + pow(stress(5), 2);

  // Compute q
  double deviatoric_q = sqrt(3 * j2);
  if (deviatoric_q < 1.0E-15) deviatoric_q = 1.0E-15;

  // Get image pressure and M_theta
  const double p_image = (*state_vars).at("p_image");
  const double M_theta = (*state_vars).at("M_theta");
  const double p_cohesion = (*state_vars).at("p_cohesion");

  // Initialise yield status (0: elastic, 1: yield)
  auto yield_type = FailureState::Elastic;

  // Compute yield functions
  if (bond_model_) {
    (*yield_function) =
        deviatoric_q / (mean_p + p_cohesion) -
        M_theta / N_ *
            (1 +
             (N_ - 1) * pow(((mean_p + p_cohesion) / (p_image + p_cohesion)),
                            (N_ / (1 - N_))));
  } else {
    (*yield_function) =
        deviatoric_q / mean_p -
        M_theta / N_ *
            (1 + (N_ - 1) * pow((mean_p / p_image), (N_ / (1 - N_))));
  }

  // Yield criterion
  if ((*yield_function) > 1.E-15) yield_type = FailureState::Yield;

  return yield_type;
}

//! Compute plastic tensor
template <unsigned Tdim>
bool mpm::NorSand<Tdim>::compute_plastic_tensor(const Vector6d& stress,
                                                mpm::dense_map* state_vars) {

  // Note that in this subroutine, stress is compression positive
  // Get state variables
  const double mean_p = (*state_vars).at("p");
  const double deviatoric_q = (*state_vars).at("q");
  const double M_theta = (*state_vars).at("M_theta");
  const double p_image = (*state_vars).at("p_image");
  const double psi_image = (*state_vars).at("psi_image");
  const double p_cohesion = (*state_vars).at("p_cohesion");

  // Estimate dilatancy at peak
  const double D_min = chi_ * psi_image;

  // Estimate maximum image pressure
  double p_image_max;
  if (bond_model_) {
    p_image_max = (mean_p + p_cohesion) *
                  pow((1 + D_min * N_ / M_theta), ((N_ - 1) / N_));
  } else {
    p_image_max = mean_p * pow((1 + D_min * N_ / M_theta), ((N_ - 1) / N_));
  }
  // Compute derivatives
  // Compute dF / dp
  double dF_dp = 0;

  if (bond_model_) {
    dF_dp = -1. * M_theta / N_ *
            (1 - pow(((mean_p + p_cohesion) / (p_image + p_cohesion)),
                     (N_ / (1 - N_))));
  } else {
    dF_dp = -1. * M_theta / N_ * (1 - pow((mean_p / p_image), (N_ / (1 - N_))));
  }

  // Compute dp / dsigma
  Vector6d dp_dsigma = Vector6d::Zero();
  dp_dsigma(0) = 1. / 3.;
  dp_dsigma(1) = 1. / 3.;
  dp_dsigma(2) = 1. / 3.;

  // Compute dF / dq
  const double dF_dq = 1.;

  // Compute the deviatoric stress
  Vector6d dev_stress = Vector6d::Zero();
  dev_stress(0) = stress(0) - mean_p;
  dev_stress(1) = stress(1) - mean_p;
  dev_stress(2) = stress(2) - mean_p;
  dev_stress(3) = stress(3);
  if (Tdim == 3) {
    dev_stress(4) = stress(4);
    dev_stress(5) = stress(5);
  }

  // Compute dq / dsigma
  Vector6d dq_dsigma = Vector6d::Zero();
  dq_dsigma(0) = 3. / 2. / deviatoric_q * dev_stress(0);
  dq_dsigma(1) = 3. / 2. / deviatoric_q * dev_stress(1);
  dq_dsigma(2) = 3. / 2. / deviatoric_q * dev_stress(2);
  dq_dsigma(3) = 3. / deviatoric_q * dev_stress(3);
  if (Tdim == 3) {
    dq_dsigma(4) = 3. / deviatoric_q * dev_stress(4);
    dq_dsigma(5) = 3. / deviatoric_q * dev_stress(5);
  }

  const double lode_angle = (*state_vars)["lode_angle"];
  const double j2 = (*state_vars)["j2"];
  const double j3 = (*state_vars)["j3"];
  const double sin_lode_angle = sin(3. / 2. * lode_angle + M_PI / 4.);
  const double cos_lode_angle = cos(3. / 2. * lode_angle + M_PI / 4.);

  // Compute dF / dM
  double dF_dM;
  if (bond_model_) {
    dF_dM =
        -1. / N_ *
        (1 + (N_ - 1) * pow(((mean_p + p_cohesion) / (p_image + p_cohesion)),
                            (N_ / (1 - N_))));
  } else {
    dF_dM =
        -1. / N_ * (1 + (N_ - 1) * pow((mean_p / p_image), (N_ / (1 - N_))));
  }

  // Compute dM / dtehta
  const double dM_dtheta =
      3. / 2. * pow(Mtc_, 3) / (3. + Mtc_) * sin_lode_angle;

  // Compute dj2 / dsigma
  Vector6d dj2_dsigma = dev_stress;

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
  Vector6d dj3_dsigma = Vector6d::Zero();
  dj3_dsigma(0) = dev1.dot(dev1) - (2. / 3.) * j2;
  dj3_dsigma(1) = dev2.dot(dev2) - (2. / 3.) * j2;
  dj3_dsigma(2) = dev3.dot(dev3) - (2. / 3.) * j2;
  dj3_dsigma(3) = dev1.dot(dev2);
  if (Tdim == 3) {
    dj3_dsigma(4) = dev2.dot(dev3);
    dj3_dsigma(5) = dev1.dot(dev3);
  }

  // Compute dtheta / dsigma
  Vector6d dtheta_dsigma = sqrt(3.) / 2. / cos(3 * lode_angle) / pow(j2, 1.5) *
                           (dj3_dsigma - 3. / 2. * j3 / j2 * dj2_dsigma);
  if (Tdim == 2) {
    dtheta_dsigma(4) = 0.;
    dtheta_dsigma(5) = 0.;
  }

  Vector6d dF_dsigma = (dF_dp * dp_dsigma) + (dF_dq * dq_dsigma) +
                       (dF_dM * dM_dtheta * dtheta_dsigma);

  // Compute hardening terms
  double dF_dpi = 0;
  if (bond_model_) {
    dF_dpi =
        -1. * M_theta *
        pow(((mean_p + p_cohesion) / (p_image + p_cohesion)), (1 / (1 - N_)));
  } else {
    dF_dpi = -1. * M_theta * pow((mean_p / p_image), (1 / (1 - N_)));
  }

  double dpi_depsd;
  if (bond_model_) {
    dpi_depsd = hardening_modulus_ * (p_image_max - p_image - p_cohesion);
  } else {
    dpi_depsd = hardening_modulus_ * (p_image_max - p_image);
  }

  const double dF_dsigma_v = (dF_dsigma(0) + dF_dsigma(1) + dF_dsigma(2)) / 3;
  double dF_dsigma_deviatoric =
      sqrt(2. / 3.) *
      sqrt(pow(dF_dsigma(0) - dF_dsigma_v, 2) +
           pow(dF_dsigma(1) - dF_dsigma_v, 2) +
           pow(dF_dsigma(2) - dF_dsigma_v, 2) + 2 * pow(dF_dsigma(3), 2) +
           2 * pow(dF_dsigma(4), 2) + 2 * pow(dF_dsigma(5), 2));

  const double dF_dpcohesion =
      -1. * M_theta / N_ *
      (1 +
       pow(((mean_p + p_cohesion) / (p_image + p_cohesion)), (1 / (1 - N_))) -
       pow(((mean_p + p_cohesion) / (p_image + p_cohesion)), (N_ / (1 - N_))));

  const double dpcohesion_depsd =
      -p_cohesion_initial_ * m_ * exp(-m_ * (*state_vars).at("epds"));

  double hardening_term;
  if (bond_model_) {
    hardening_term = dF_dpi * dpi_depsd * dF_dsigma_deviatoric +
                     dF_dpcohesion * dpcohesion_depsd * dF_dsigma_deviatoric;
  } else {
    hardening_term = dF_dpi * dpi_depsd * dF_dsigma_deviatoric;
  }

  // Construct Dp matrix
  this->dp_ = (de_ * dF_dsigma * dF_dsigma.transpose() * de_) /
              (dF_dsigma.transpose() * de_ * dF_dsigma - hardening_term);

  return true;
}

//! Compute stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::NorSand<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr, mpm::dense_map* state_vars) {

  // Note: compression positive in all derivations
  Vector6d stress_neg = -1 * stress;
  Vector6d dstrain_neg = -1 * dstrain;

  // Compute current stress invariants (need to remove this, only needed in the
  // first step)
  this->compute_stress_invariants(stress_neg, state_vars);

  // Elastic step
  // --------------------------------------------------------------------------------------
  // Shear modulus
  // shear_modulus_ = shear_modulus_constant_ *
  //                  pow(((*state_vars).at("p") / reference_pressure_),
  //                      shear_modulus_exponent_);
  shear_modulus_ = youngs_modulus_ / (2.0 * (1. + poisson_ratio_));
  // Bulk modulus
  bulk_modulus_ = shear_modulus_ * (2.0 * (1 + poisson_ratio_)) /
                  (3.0 * (1. - 2. * poisson_ratio_));
  // Set elastic tensor
  this->compute_elastic_tensor();

  // Trial stress - elastic
  Vector6d trial_stress = stress_neg + (this->de_ * dstrain_neg);

  // Initialise value for yield function
  double yield_function;
  auto yield_type =
      this->compute_yield_state(&yield_function, trial_stress, state_vars);

  // Return the updated stress in elastic state
  if (yield_type == FailureState::Elastic) return -trial_stress;
  // --------------------------------------------------------------------------------------

  // Set plastic tensor
  this->compute_plastic_tensor(stress_neg, state_vars);

  // Plastic step
  // Compute D matrix used in stress update
  Matrix6x6 D_matrix = this->de_ - this->dp_;

  // Update stress
  Vector6d updated_stress = stress_neg + D_matrix * dstrain_neg;

  // Update state variables
  this->compute_state_variables(-1 * updated_stress, dstrain_neg, state_vars);

  // Update stress invariants
  this->compute_stress_invariants(updated_stress, state_vars);

  // Compute incremental plastic strain
  Vector6d dstress = updated_stress - stress;
  Vector6d dpstrain = dstrain - (this->de_.inverse() * dstress);
  if (Tdim == 2) dpstrain(4) = dpstrain(5) = 0.;

  // Update plastic strain
  (*state_vars).at("plastic_strain0") += dpstrain(0);
  (*state_vars).at("plastic_strain1") += dpstrain(1);
  (*state_vars).at("plastic_strain2") += dpstrain(2);
  (*state_vars).at("plastic_strain3") += dpstrain(3);
  (*state_vars).at("plastic_strain4") += dpstrain(4);
  (*state_vars).at("plastic_strain5") += dpstrain(5);

  // Update equivalent plastic deviatoric strain
  (*state_vars).at("epds") =
      sqrt(2. / 9. *
               (pow(((*state_vars).at("plastic_strain0") -
                     (*state_vars).at("plastic_strain1")),
                    2.) +
                pow(((*state_vars).at("plastic_strain1") -
                     (*state_vars).at("plastic_strain2")),
                    2.) +
                pow(((*state_vars).at("plastic_strain2") -
                     (*state_vars).at("plastic_strain0")),
                    2.)) +
           1. / 3. *
               (pow(((*state_vars).at("plastic_strain3")), 2.) +
                pow(((*state_vars).at("plastic_strain4")), 2.) +
                pow(((*state_vars).at("plastic_strain5")), 2.)));

  // Update p_cohesion
  this->compute_p_cohesion(state_vars);

  return (-updated_stress);
}
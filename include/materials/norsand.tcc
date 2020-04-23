//! Constructor with id and material properties
template <unsigned Tdim>
mpm::NorSand<Tdim>::NorSand(unsigned id, const Json& material_properties)
    : Material<Tdim>(id, material_properties) {
  try {
    // Density
    density_ = material_properties.at("density").template get<double>();
    // Poisson ratio
    poisson_ratio_ =
        material_properties.at("poisson_ratio").template get<double>();
    // Reference pressure pref
    reference_pressure_ =
        material_properties.at("reference_pressure").template get<double>();
    // Critical state friction angle
    friction_cs_ =
        material_properties.at("friction_cs").template get<double>() * M_PI /
        180.;
    // Volumetric coupling (dilatancy) parameter N
    N_ = material_properties.at("N").template get<double>();
    // Lambda volumetric
    lambda_ = material_properties.at("lambda").template get<double>();
    // Kappa swelling volumetric
    kappa_ = material_properties.at("kappa").template get<double>();
    // Gamma void ratio at reference pressure
    gamma_ = material_properties.at("gamma").template get<double>();
    // Dilatancy coefficient chi
    chi_ = material_properties.at("chi").template get<double>();
    // Hardening modulus
    hardening_modulus_ =
        material_properties.at("hardening_modulus").template get<double>();
    // Initial void ratio
    void_ratio_initial_ =
        material_properties.at("void_ratio_initial").template get<double>();
    // Initial image pressure
    p_image_initial_ =
        material_properties.at("p_image_initial").template get<double>();
    // Flag for bonded model
    bond_model_ = material_properties.at("bond_model").template get<bool>();
    // Initial p_cohesion
    p_cohesion_initial_ =
        material_properties.at("p_cohesion_initial").template get<double>();
    // Initial p_dilation
    p_dilation_initial_ =
        material_properties.at("p_dilation_initial").template get<double>();
    // Cohesion degradation parameter m upon shearing
    m_cohesion_ = material_properties.at("m_cohesion").template get<double>();
    // Dilation degradation parameter m upon shearing
    m_dilation_ = material_properties.at("m_dilation").template get<double>();
    // Parameter for shear modulus
    m_modulus_ = material_properties.at("m_modulus").template get<double>();
    // Default tolerance
    tolerance_ = material_properties.at("tolerance").template get<double>();

    const double sin_friction_cs = sin(friction_cs_);
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
      // M_theta
      {"M_theta", Mtc_},
      // Current void ratio
      {"void_ratio", void_ratio_initial_},
      // Void ratio image
      {"e_image",
       gamma_ - lambda_ * log(p_image_initial_ / reference_pressure_)},
      // Image pressure
      {"p_image", p_image_initial_},
      // p_cohesion
      {"p_cohesion", p_cohesion_initial_},
      // p_dilation
      {"p_dilation", p_dilation_initial_},
      // Equivalent plastic deviatoric strain
      {"pdstrain", 0.},
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
Eigen::Matrix<double, 4, 1> mpm::NorSand<Tdim>::compute_stress_invariants(
    const Vector6d& stress, mpm::dense_map* state_vars) {

  // Note that in this subroutine, stress is compression positive

  // Compute mean stress p
  double mean_p = -1. * mpm::materials::p(-stress);
  mean_p = check_low(mean_p);

  // Compute q
  double deviatoric_q = mpm::materials::q(-stress);
  deviatoric_q = check_low(deviatoric_q);

  // Compute Lode angle (cos convetion)
  const double lode_angle = mpm::materials::lode_angle(-stress, tolerance_);

  // Compute M_theta (Jefferies and Shuttle, 2011)
  const double cos_lode_angle = cos(3. / 2. * lode_angle);
  double M_theta = Mtc_ - std::pow(Mtc_, 2) / (3. + Mtc_) * cos_lode_angle;

  // Store to return
  Eigen::Matrix<double, 4, 1> invariants;
  invariants << mean_p, deviatoric_q, lode_angle, M_theta;

  return invariants;
}

//! Compute state parameters
template <unsigned Tdim>
void mpm::NorSand<Tdim>::compute_state_variables(
    const Vector6d& stress, const Vector6d& dstrain, mpm::dense_map* state_vars,
    mpm::norsand::FailureState yield_type) {

  // Get invariants
  auto invariants = this->compute_stress_invariants(stress, state_vars);

  const double mean_p = invariants(0);
  const double deviatoric_q = invariants(1);

  // Get state variables (note that M_theta used is at current stress)
  const double M_theta = (*state_vars).at("M_theta");
  const double p_cohesion = (*state_vars).at("p_cohesion");
  const double p_dilation = (*state_vars).at("p_dilation");
  double p_image;
  double e_image;

  if (yield_type == mpm::norsand::FailureState::Elastic) {
    // Keep the same pressure image and void ratio image at critical state
    p_image = (*state_vars).at("p_image");
    e_image = (*state_vars).at("e_image");
  } else {
    // Compute and update pressure image
    p_image =
        (mean_p + p_cohesion) *
            std::pow(
                ((1 - N_ / M_theta * deviatoric_q / (mean_p + p_cohesion)) /
                 (1 - N_)),
                ((N_ - 1) / N_)) -
        p_cohesion - p_dilation;
    (*state_vars).at("p_image") = p_image;

    // Compute and update void ratio image
    // e_image = e_max_ - (e_max_ - e_min_) / log(crushing_pressure_ / p_image);
    e_image = gamma_ - lambda_ * log(p_image / reference_pressure_);
    e_image = check_low(e_image);

    (*state_vars).at("e_image") = e_image;
  }

  // Update M_theta at the updated stress state
  (*state_vars).at("M_theta") = invariants(3);

  // Update void ratio
  // Note that dstrain is in tension positive - depsv = de / (1 + e_initial)
  double dvolumetric_strain = dstrain(0) + dstrain(1) + dstrain(2);
  double void_ratio = (*state_vars).at("void_ratio") -
                      (1 + void_ratio_initial_) * dvolumetric_strain;
  void_ratio = check_low(void_ratio);
  (*state_vars).at("void_ratio") = void_ratio;
}

//! Compute elastic tensor
template <unsigned Tdim>
void mpm::NorSand<Tdim>::compute_p_bond(mpm::dense_map* state_vars) {

  // Compute current zeta cohesion
  double zeta_cohesion = exp(-m_cohesion_ * (*state_vars).at("pdstrain"));
  zeta_cohesion = check_one(zeta_cohesion);
  zeta_cohesion = check_low(zeta_cohesion);

  // Update p_cohesion
  double p_cohesion = p_cohesion_initial_ * zeta_cohesion;
  (*state_vars).at("p_cohesion") = p_cohesion;

  // Compute current zeta dilation
  double zeta_dilation = exp(-m_dilation_ * (*state_vars).at("pdstrain"));
  zeta_dilation = check_one(zeta_dilation);
  zeta_dilation = check_low(zeta_dilation);

  // Update p_dilation
  double p_dilation = p_dilation_initial_ * zeta_dilation;
  (*state_vars).at("p_dilation") = p_dilation;
}

//! Compute yield function and yield state
template <unsigned Tdim>
typename mpm::norsand::FailureState mpm::NorSand<Tdim>::compute_yield_state(
    double* yield_function, const Vector6d& stress,
    mpm::dense_map* state_vars) {

  // Get stress invariants
  auto invariants = this->compute_stress_invariants(stress, state_vars);

  double mean_p = invariants(0);
  double deviatoric_q = invariants(1);

  // Get state variables
  const double p_image = (*state_vars).at("p_image");
  const double M_theta = (*state_vars).at("M_theta");
  const double p_cohesion = (*state_vars).at("p_cohesion");
  const double p_dilation = (*state_vars).at("p_dilation");

  // Initialise yield status (Elastic, Yield)
  auto yield_type = mpm::norsand::FailureState::Elastic;

  // Compute yield functions
  (*yield_function) =
      deviatoric_q / (mean_p + p_cohesion) -
      M_theta / N_ *
          (1 + (N_ - 1) * std::pow(((mean_p + p_cohesion) /
                                    (p_image + p_cohesion + p_dilation)),
                                   (N_ / (1 - N_))));

  // Yield criterion
  if ((*yield_function) > tolerance_)
    yield_type = mpm::norsand::FailureState::Yield;

  return yield_type;
}

//! Compute plastic tensor
template <unsigned Tdim>
void mpm::NorSand<Tdim>::compute_plastic_tensor(const Vector6d& stress,
                                                mpm::dense_map* state_vars) {

  // Note that in this subroutine, stress is compression positive

  // Get stress invariants
  auto invariants = this->compute_stress_invariants(stress, state_vars);

  const double mean_p = invariants(0);
  const double deviatoric_q = invariants(1);
  const double lode_angle = invariants(2);

  // Get state variables
  const double M_theta = (*state_vars).at("M_theta");
  const double p_image = (*state_vars).at("p_image");
  const double e_image = (*state_vars).at("e_image");
  const double void_ratio = (*state_vars).at("void_ratio");
  const double p_cohesion = (*state_vars).at("p_cohesion");
  const double p_dilation = (*state_vars).at("p_dilation");

  // Estimate dilatancy at peak
  const double D_min = chi_ * (void_ratio - e_image);

  // Estimate maximum image pressure
  const double p_image_max =
      (mean_p + p_cohesion) *
      std::pow((1 + D_min * N_ / M_theta), ((N_ - 1) / N_));

  // Compute derivatives
  // Compute dF / dp
  const double dF_dp =
      -1. * M_theta / N_ *
      (1 - std::pow(((mean_p + p_cohesion) / (p_image + p_cohesion)),
                    (N_ / (1 - N_))));

  // Compute dp / dsigma
  const Vector6d dp_dsigma = mpm::materials::dp_dsigma(-stress);

  // Compute dF / dq
  const double dF_dq = 1.;

  // Compute dq / dsigma
  const Vector6d dq_dsigma = mpm::materials::dq_dsigma(-stress);

  // Compute dF / dM
  const double dF_dM =
      -1.0 / N_ * (mean_p + p_cohesion) *
      (1 + (N_ - 1) * std::pow(((mean_p + p_cohesion) /
                                (p_image + p_cohesion + p_dilation)),
                               (N_ / (1 - N_))));

  // Use current lode angle to compute dtheta
  const double sin_lode_angle = sin(3. / 2. * lode_angle);

  // Compute dM / dtehta
  const double dM_dtheta =
      3. / 2. * std::pow(Mtc_, 2) / (3. + Mtc_) * sin_lode_angle;

  // Compute dtheta / dsigma
  const Vector6d dtheta_dsigma = mpm::materials::dtheta_dsigma(-stress);

  // dF_dsigma is in compression negative
  const Vector6d dF_dsigma = (dF_dp * dp_dsigma) + (-1. * dF_dq * dq_dsigma) +
                             (-1. * dF_dM * dM_dtheta * dtheta_dsigma);

  // Derivatives in respect to p_image
  const double dF_dpi =
      -1. * M_theta *
      std::pow(((mean_p + p_cohesion) / (p_image + p_cohesion + p_dilation)),
               (1 / (1 - N_)));

  const double dpi_depsd = hardening_modulus_ * (p_image_max - p_image);

  const double dF_dsigma_v = (dF_dsigma(0) + dF_dsigma(1) + dF_dsigma(2)) / 3;
  const double dF_dsigma_deviatoric =
      std::sqrt(2. / 3.) *
      std::sqrt(std::pow(dF_dsigma(0) - dF_dsigma_v, 2) +
                std::pow(dF_dsigma(1) - dF_dsigma_v, 2) +
                std::pow(dF_dsigma(2) - dF_dsigma_v, 2) +
                2 * std::pow(dF_dsigma(3), 2) + 2 * std::pow(dF_dsigma(4), 2) +
                2 * std::pow(dF_dsigma(5), 2));

  // Compute hardering term
  double hardening_term;
  if (bond_model_) {
    // Derivatives in respect to p_cohesion
    const double dF_dpcohesion =
        M_theta / N_ *
        (1 +
         (N_ - 1) * std::pow(((mean_p + p_cohesion) /
                              (p_image + p_cohesion + p_dilation)),
                             (N_ / (1 - N_))) -
         N_ * (p_image + p_dilation - mean_p) /
             (p_image + p_cohesion + p_dilation) *
             std::pow(
                 ((mean_p + p_cohesion) / (p_image + p_cohesion + p_dilation)),
                 (N_ / (1 - N_))));

    const double dpcohesion_depsd =
        -p_cohesion_initial_ * m_cohesion_ *
        exp(-m_cohesion_ * (*state_vars).at("pdstrain"));

    // Derivatives in respect to p_dilation
    const double dF_dpdilation =
        -1. * M_theta *
        std::pow(((mean_p + p_cohesion) / (p_image + p_cohesion + p_dilation)),
                 (1 / (1 - N_)));

    const double dpdilation_depsd =
        -p_dilation_initial_ * m_dilation_ *
        exp(-m_dilation_ * (*state_vars).at("pdstrain"));

    hardening_term = dF_dpi * dpi_depsd * dF_dsigma_deviatoric +
                     dF_dpcohesion * dpcohesion_depsd * dF_dsigma_deviatoric +
                     dF_dpdilation * dpdilation_depsd * dF_dsigma_deviatoric;
  } else {
    hardening_term = dF_dpi * dpi_depsd * dF_dsigma_deviatoric;
  }

  // Construct Dp matrix
  this->dp_ = (de_ * dF_dsigma * dF_dsigma.transpose() * de_) /
              (dF_dsigma.transpose() * de_ * dF_dsigma - hardening_term);
}

//! Compute stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::NorSand<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr, mpm::dense_map* state_vars) {

  // Zero parameters for non-bond model
  if (!bond_model_) {
    (*state_vars).at("p_cohesion") = 0.0;
    (*state_vars).at("p_dilation") = 0.0;
    m_cohesion_ = 0.0;
    m_dilation_ = 0.0;
    m_modulus_ = 0.0;
  }

  // Note: compression positive in all derivations
  Vector6d stress_neg = -1 * stress;
  Vector6d dstrain_neg = -1 * dstrain;

  // Get stress invariants
  auto invariants = this->compute_stress_invariants(stress_neg, state_vars);
  double mean_p = invariants(0);

  // Elastic step
  // Bulk modulus computation
  bulk_modulus_ = (1. + (*state_vars).at("void_ratio")) / kappa_ * mean_p +
                  m_modulus_ * ((*state_vars).at("p_cohesion") +
                                (*state_vars).at("p_dilation"));
  // Shear modulus computation
  shear_modulus_ = 3. * bulk_modulus_ * (1. - 2. * poisson_ratio_) /
                   (2.0 * (1. + poisson_ratio_));

  // Set elastic tensor
  this->compute_elastic_tensor();

  // Trial stress - elastic
  Vector6d trial_stress = stress_neg + (this->de_ * dstrain_neg);

  // Initialise value for yield function
  double yield_function;
  auto yield_type =
      this->compute_yield_state(&yield_function, trial_stress, state_vars);

  // Return the updated stress in elastic state
  if (yield_type == mpm::norsand::FailureState::Elastic) {

    // Update state variables
    this->compute_state_variables(trial_stress, dstrain_neg, state_vars,
                                  yield_type);

    // Update p_cohesion
    this->compute_p_bond(state_vars);

    // Return elastic stress
    return (-trial_stress);
  }

  // Set plastic tensor
  this->compute_plastic_tensor(stress_neg, state_vars);

  // Plastic step
  // Compute D matrix used in stress update
  Matrix6x6 D_matrix = this->de_ - this->dp_;

  // Update stress
  Vector6d updated_stress = stress_neg + D_matrix * dstrain_neg;

  // Update state variables
  this->compute_state_variables(updated_stress, dstrain_neg, state_vars,
                                yield_type);

  // Compute incremental plastic strain, still in tension positive
  Vector6d dstress = updated_stress - stress_neg;
  Vector6d dpstrain = dstrain_neg - (this->de_.inverse() * dstress);
  if (Tdim == 2) dpstrain(4) = dpstrain(5) = 0.;

  // Update plastic strain
  (*state_vars).at("plastic_strain0") += dpstrain(0);
  (*state_vars).at("plastic_strain1") += dpstrain(1);
  (*state_vars).at("plastic_strain2") += dpstrain(2);
  (*state_vars).at("plastic_strain3") += dpstrain(3);
  (*state_vars).at("plastic_strain4") += dpstrain(4);
  (*state_vars).at("plastic_strain5") += dpstrain(5);

  // Update equivalent plastic deviatoric strain
  (*state_vars).at("pdstrain") =
      std::sqrt(2. / 9. *
                    (std::pow(((*state_vars).at("plastic_strain0") -
                               (*state_vars).at("plastic_strain1")),
                              2.) +
                     std::pow(((*state_vars).at("plastic_strain1") -
                               (*state_vars).at("plastic_strain2")),
                              2.) +
                     std::pow(((*state_vars).at("plastic_strain2") -
                               (*state_vars).at("plastic_strain0")),
                              2.)) +
                1. / 3. *
                    (std::pow(((*state_vars).at("plastic_strain3")), 2.) +
                     std::pow(((*state_vars).at("plastic_strain4")), 2.) +
                     std::pow(((*state_vars).at("plastic_strain5")), 2.)));

  // Update p_cohesion
  this->compute_p_bond(state_vars);

  // Return updated stress
  return (-updated_stress);
}
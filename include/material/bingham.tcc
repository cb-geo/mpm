//! Read material properties
template <unsigned Tdim>
void mpm::Bingham<Tdim>::properties(const Json& material_properties) {
  try {
    density_ = material_properties["density"].template get<double>();
    youngs_modulus_ =
        material_properties["youngs_modulus"].template get<double>();
    poisson_ratio_ =
        material_properties["poisson_ratio"].template get<double>();
    tau0_ = material_properties["tau0"].template get<double>();
    mu_ = material_properties["mu"].template get<double>();
    critical_shear_rate_ =
        material_properties["critical_shear_rate"].template get<double>();
    properties_ = material_properties;
    status_ = true;
  } catch (std::exception& except) {
    std::cerr << "Material parameter not set: " << except.what() << '\n';
  }
}

//! Elastic tensor is not defined in Bingham model, throws an error
template <unsigned Tdim>
Eigen::Matrix<double, 6, 6> mpm::Bingham<Tdim>::elastic_tensor() {

  Eigen::Matrix<double, 6, 6> de;
  de.setZero();

  throw std::runtime_error("Elastic tensor is not used for this material");

  return de;
}

//! Compute stress without a particle handle is undefined in the Bingham model,
//! throws an error
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::Bingham<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain) {

  Vector6d stress_results;
  stress_results.setZero();

  throw std::runtime_error("Stress computation for this material is not valid");

  return stress_results;
}

//! Compute stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::Bingham<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr) {

  const unsigned phase = 0;
  const auto strain_rate = ptr->strain_rate(phase);

  // Bulk modulus
  const double K = youngs_modulus_ / (3.0 * (1. - 2. * poisson_ratio_));

  // Get volumetric change and update pressure
  // pressure_new = pressure_old + dpressure
  // dpressure = K * strain_volumetric
  const double dpressure = K * (dstrain(0) + dstrain(1) + dstrain(2));
  const double pressure_old = (stress(0) + stress(1) + stress(2)) / 3.0;
  const double pressure_new = pressure_old + dpressure;

  // Determine accuracy of minimum critical shear rate
  const double shear_rate_threshold = 1.0E-15;
  if (critical_shear_rate_ < shear_rate_threshold)
    critical_shear_rate_ = shear_rate_threshold;

  // Checking yielding from strain rate vs critical yielding shear rate
  // rate of shear = sqrt(2 * strain_rate * strain_rate)
  // yielding is defined: rate of shear > critical_shear_rate_^2
  // modulus maps shear rate to shear stress
  const double shear_rate = 2 * strain_rate.dot(strain_rate);
  double modulus;
  if (shear_rate > critical_shear_rate_ * critical_shear_rate_)
    modulus = 2 * ((tau0_ / (std::sqrt(shear_rate))) + mu_);
  else
    modulus = 0.;

  // Compute shear change to volumetric
  // tau deviatoric part of cauchy stress tensor
  Eigen::Matrix<double, 6, 1> tau;
  tau.setZero();
  tau = modulus * strain_rate;

  // Use von Mises criterion
  // second invariant of tau > 2 tau0^2
  double invariant2 = tau.dot(tau);
  if (invariant2 < 2 * (tau0_ * tau0_)) tau.setZero();

  // Update volumetric and deviatoric stress
  Eigen::Matrix<double, 6, 1> stress_results;
  stress_results.setZero();
  // Get dirac delta function in Voigt notation
  const auto dirac_delta = this->dirac_delta();

  stress_results = pressure_new * dirac_delta + tau;

  return stress_results;
}

//! Dirac delta 2D
template <>
Eigen::Matrix<double, 6, 1> mpm::Bingham<2>::dirac_delta() { 

  Eigen::Matrix<double, 6, 1> dirac_delta;  
  dirac_delta << 1, 1, 0, 0, 0, 0;

  return dirac_delta;
}

//! Dirac delta 3D
template <>
Eigen::Matrix<double, 6, 1> mpm::Bingham<3>::dirac_delta() { 

  Eigen::Matrix<double, 6, 1> dirac_delta;  
  dirac_delta << 1, 1, 1, 0, 0, 0;

  return dirac_delta;
}
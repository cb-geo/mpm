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
  auto strain_rate = ptr->strain_rate(phase);

  // Get defintion of D for Bingham
  strain_rate.tail(3) *= 0.5;

  // Determine accuracy of minimum critical shear rate
  const double shear_rate_threshold = 1.0E-15;
  if (critical_shear_rate_ < shear_rate_threshold)
    critical_shear_rate_ = shear_rate_threshold;

  // Checking yielding from strain rate vs critical yielding shear rate
  // rate of shear = sqrt(2 * strain_rate * strain_rate)
  // yielding is defined: rate of shear > critical_shear_rate_^2
  // modulus maps shear rate to shear stress
  const double shear_rate = 2 * strain_rate.dot(strain_rate);
  double modulus = 0;
  if (shear_rate > critical_shear_rate_ * critical_shear_rate_)
    modulus = 2 * ((tau0_ / (std::sqrt(shear_rate))) + mu_);

  // Compute shear change to volumetric
  // tau deviatoric part of cauchy stress tensor
  Eigen::Matrix<double, 6, 1> tau;
  tau = modulus * strain_rate;

  // Use von Mises criterion
  // second invariant of tau > 2 tau0^2
  double invariant2 = 0.5 * tau.dot(tau);
  if (invariant2 < (tau0_ * tau0_)) tau.setZero();

  // Get pressure
  const double pressure = ptr->pressure(phase);

  // Get dirac delta function in Voigt notation
  const auto dirac_delta = this->dirac_delta();

  // Update volumetric and deviatoric stress
  Eigen::Matrix<double, 6, 1> updated_stress;
  updated_stress = -pressure * dirac_delta + tau;

  return updated_stress;
}

//! Dirac delta 2D
template <>
inline Eigen::Matrix<double, 6, 1> mpm::Bingham<2>::dirac_delta() const {

  return (Eigen::Matrix<double, 6, 1>() << 1.f, 1.f, 0.f, 0.f, 0.f, 0.f)
      .finished();
}

//! Dirac delta 3D
template <>
inline Eigen::Matrix<double, 6, 1> mpm::Bingham<3>::dirac_delta() const {

  return (Eigen::Matrix<double, 6, 1>() << 1.f, 1.f, 1.f, 0.f, 0.f, 0.f)
      .finished();
}
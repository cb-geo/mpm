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
    strain_cutoff_ =
        material_properties["strain_cutoff"].template get<double>();
    properties_ = material_properties;
    status_ = true;
  } catch (std::exception& except) {
    std::cerr << "Material parameter not set: " << except.what() << '\n';
  }
}

//! Not used in this model, thus return error
template <unsigned Tdim>
Eigen::Matrix<double, 6, 6> mpm::Bingham<Tdim>::elastic_tensor() {

  Eigen::Matrix<double, 6, 6> de;
  de.setZero();

  try {
    throw std::runtime_error("Elastic tensor is not used for this material");
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
  }

  return de;
}

//! Not used in this model, thus return error
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::Bingham<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain) {

  Vector6d stress_results;
  stress_results.setZero();

  try {
    throw std::runtime_error(
        "Stress computation for this material is not valid");
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
  }

  return stress_results;
}

//! Compute stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::Bingham<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr) {

  unsigned phase = 0;
  auto strain_rate = ptr->strain_rate(phase);

  // Bulk modulus
  const double K = youngs_modulus_ / (3.0 * (1. - 2. * poisson_ratio_));

  // Make minimum of strain_cutoff accuracy
  const double strain_threshold = 1.0E-15;
  if (strain_cutoff_ < strain_threshold) strain_cutoff_ = strain_threshold;

  // Get volumetric change and update pressure
  const double pressure_old = (stress(0) + stress(1) + stress(2)) / 3.0;
  const double dpressure = K * (dstrain(0) + dstrain(1) + dstrain(2));
  const double pressure_new = pressure_old + dpressure;

  const double invariant2 = 0.5 * strain_rate.dot(strain_rate);

  // Compute deviatoric change
  double factor;
  if (std::sqrt(invariant2) > strain_cutoff_)
    factor = ((tau0_ / (std::sqrt(invariant2))) + 2 * mu_);
  else
    factor = 0.;

  // Compute deviatoric strain
  Eigen::Vector3d tau;
  tau.setZero();
  tau(0) = factor * strain_rate(0);
  tau(1) = factor * strain_rate(1);
  tau(2) = factor * strain_rate(2);

  // double sum_squared_tau = 0.5 * (tau(0) * tau(0) + tau(1) * tau(1) + tau(2)
  // * tau(2));
  double sum_squared_tau = 0.5 * tau.dot(tau);
  if (sum_squared_tau < (tau0_ * tau0_)) tau.setZero();

  // Update stress
  Eigen::Matrix<double, 6, 1> stress_results;
  stress_results.setZero();
  stress_results(0) = tau(0) + pressure_new;
  stress_results(1) = tau(1) + pressure_new;
  stress_results(2) = tau(2) + pressure_new;
  stress_results(3) = factor * strain_rate(3);
  stress_results(4) = factor * strain_rate(4);
  stress_results(5) = factor * strain_rate(5);

  return stress_results;
}

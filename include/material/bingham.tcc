//! Read material properties
template <unsigned Tdim>
void mpm::Bingham<Tdim>::properties(const Json& material_properties) {
  try {
    density_ = material_properties["density"].template get<double>();
    youngs_modulus_ =
        material_properties["youngs_modulus"].template get<double>();
    poisson_ratio_ =
        material_properties["poisson_ratio"].template get<double>();
    tau0_ =
        material_properties["tau0"].template get<double>();
    mu_ =
        material_properties["mu"].template get<double>();
    strain_cutoff_ =
        material_properties["strain_cutoff"].template get<double>();        
    properties_ = material_properties;
    status_ = true;
  } catch (std::exception& except) {
    std::cerr << "Material parameter not set: " << except.what() << '\n';
  }
}

//! Compute stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::Bingham<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr) {

  Eigen::Matrix<double, 6, 1> stress_results;

  unsigned phase = 0;
  auto strain_rate = ptr->strain_rate(phase);

  // Bulk and shear modulus
  const double K = youngs_modulus_ / (3.0 * (1. - 2. * poisson_ratio_));
  const double G = youngs_modulus_ / (2.0 * (1. + poisson_ratio_));

  // Make minimum of strain_cutoff accuracy 
  if (strain_cutoff_ < 1.E-15) strain_cutoff_ = 1.0E-15;

  // Get volumetric change
  double p_old = (stress(0) + stress(1) + stress(2)) / 3.0;
  double dp = K * (dstrain(0) + dstrain(1) + dstrain(2));
  double p_new = p_old + dp;

  // double invariant2 =
  //     0.5 * ((strain_rate[0, 0) * strain_rate[0, 0)) + (strain_rate[0, 1) * strain_rate[0, 1)) +
  //            (strain_rate[0, 2) * strain_rate[0, 2)) + (strain_rate[1, 0) * strain_rate[1, 0)) +
  //            (strain_rate[1, 1) * strain_rate[1, 1)) + (strain_rate[1, 2) * strain_rate[1, 2)) +
  //            (strain_rate[2, 0) * strain_rate[2, 0)) + (strain_rate[2, 1) * strain_rate[2, 1)) +
  //            (strain_rate[2, 2) * strain_rate[2, 2)));

  double invariant2 = 0.5 * strain_rate.dot(strain_rate);

  // Compute deviatoric change  
  double factor;
  if (std::sqrt(invariant2) > strain_cutoff_)
    factor = ((tau0_ / (std::sqrt(invariant2))) + 2 * mu_);
  else
    factor = 0.;

  Eigen::Vector3d tau;
  tau(0) = factor * strain_rate(0);
  tau(1) = factor * strain_rate(1);
  tau(2) = factor * strain_rate(2);

  // double sum_squared_tau = 0.5 * (tau(0) * tau(0) + tau(1) * tau(1) + tau(2) * tau(2));
  double sum_squared_tau = 0.5 * tau.dot(tau);
  if (sum_squared_tau < (tau0_ * tau0_)) {
    tau(0) = 0;
    tau(1) = 0;
    tau(2) = 0;
  }

  // Update stress
  stress_results(0) = tau(0) + p_new;
  stress_results(1) = tau(1) + p_new;
  stress_results(2) = tau(2) + p_new;
  stress_results(3) = factor * strain_rate(3);
  stress_results(4) = factor * strain_rate(4);
  stress_results(5) = factor * strain_rate(5);

  return stress_results;

}

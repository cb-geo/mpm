//! Read material properties
void mpm::Bingham::properties(const Json& material_properties) {
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
void mpm::Bingham::compute_stress(Vector6d& stress,
                                  const Vector6d& dstrain,
                                  const Eigen::Matrix<double, 3, 3>& strain_rate) {

  // Bulk and shear modulus
  const double K = youngs_modulus_ / (3.0 * (1. - 2. * poisson_ratio_));
  const double G = youngs_modulus_ / (2.0 * (1. + poisson_ratio_));

  // Make minimum of strain_cutoff accuracy 
  if (strain_cutoff_ < 1.E-15) strain_cutoff_ = 1.0E-15;

  // Get volumetric change
  double pOld = (stress(0) + stress(1) + stress(2)) / 3.0;
  double dP = K * (dstrain(0) + dstrain(1) + dstrain(2));
  double pNew = pOld + dP;

  // double I2 =
  //     0.5 * ((strain_rate[0, 0] * strain_rate[0, 0]) + (strain_rate[0, 1] * strain_rate[0, 1]) +
  //            (strain_rate[0, 2] * strain_rate[0, 2]) + (strain_rate[1, 0] * strain_rate[1, 0]) +
  //            (strain_rate[1, 1] * strain_rate[1, 1]) + (strain_rate[1, 2] * strain_rate[1, 2]) +
  //            (strain_rate[2, 0] * strain_rate[2, 0]) + (strain_rate[2, 1] * strain_rate[2, 1]) +
  //            (strain_rate[2, 2] * strain_rate[2, 2]));

  double I2 = 0.5 * (strain_rate.cwiseProduct(strain_rate)).sum();

  // Compute deviatoric change  
  double factor;
  if (std::sqrt(I2) > strain_cutoff_)
    factor = ((tau0_ / (std::sqrt(I2))) + 2 * mu_);
  else
    factor = 0.;

  Eigen::Vector3d tau;
  tau[0] = factor * strain_rate[0, 0];
  tau[1] = factor * strain_rate[1, 1];
  tau[2] = factor * strain_rate[2, 2];

  // double traceTau2 = 0.5 * (tau[0] * tau[0] + tau[1] * tau[1] + tau[2] * tau[2]);
  double tracetau2 = 0.5 * tau.dot(tau);
  if (tracetau2 < (tau0_ * tau0_)) {
    tau(0) = 0;
    tau(1) = 0;
    tau(2) = 0;
  }

  // Update stress
  stress[0] = tau[0] + pNew;
  stress[1] = tau[1] + pNew;
  stress[2] = tau[2] + pNew;
  stress[3] = factor * strain_rate[0, 1];
  stress[4] = factor * strain_rate[1, 2];
  stress[5] = factor * strain_rate[0, 2];

}

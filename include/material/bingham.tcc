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
    console_->error("Material parameter not set: {}\n", except.what());
  }
}

//! Elastic tensor is not defined in Bingham model, throws an error
template <unsigned Tdim>
Eigen::Matrix<double, 6, 6> mpm::Bingham<Tdim>::elastic_tensor() {

  throw std::runtime_error("Elastic tensor is not used for this material");

  return Eigen::Matrix<double, 6, 6>::Zero();
}

//! Compute stress without a particle handle is undefined in the Bingham model,
//! throws an error
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::Bingham<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain) {

  throw std::runtime_error(
      "Stress computation for this material requires a particle handle");

  return Vector6d::Zero();
}

//! Compute stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::Bingham<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr) {

  const unsigned phase = 0;

  // Bulk modulus
  const double bulk_modulus =
      youngs_modulus_ / (3.0 * (1. - 2. * poisson_ratio_));

  // Get volumetric strain from particle
  const double volumetric_strain = ptr->volumetric_strain_centroid(phase);

  // Compute thermodynamic pressure
  // thermodynamic_pressure = -bulk_modulus * volstrain
  // Expansion causes a decrease in thermodynamic pressure
  const double thermodynamic_pressure = -bulk_modulus * volumetric_strain;

  auto strain_rate = ptr->strain_rate(phase);

  // Convert strain rate to rate of deformation tensor
  strain_rate.tail(3) *= 0.5;

  // Set threshold for minimum critical shear rate
  const double shear_rate_threshold = 1.0E-15;
  if (critical_shear_rate_ < shear_rate_threshold)
    critical_shear_rate_ = shear_rate_threshold;

  // Rate of shear = sqrt(2 * D_ij * D_ij)
  // Since D (D_ij) is in Voigt notation (D_i), and the definition above is in
  // matrix, the last 3 components have to be doubled D_ij * D_ij = D_0^2 +
  // D_1^2 + D_2^2 + 2*D_3^2 + 2*D_4^2 + 2*D_5^2 Yielding is defined: rate of
  // shear > critical_shear_rate_^2 Checking yielding from strain rate vs
  // critical yielding shear rate
  double shear_rate =
      std::sqrt(2 * (strain_rate.dot(strain_rate) +
                     strain_rate.tail(3).dot(strain_rate.tail(3))));

  // Apparent_viscosity maps shear rate to shear stress
  // Check if shear rate is 0
  double apparent_viscosity = 0;
  if (shear_rate * shear_rate > critical_shear_rate_ * critical_shear_rate_) {
    if (shear_rate < shear_rate_threshold) shear_rate = shear_rate_threshold;
    apparent_viscosity = 2 * ((tau0_ / shear_rate) + mu_);
  }

  // Compute shear change to volumetric
  // tau deviatoric part of cauchy stress tensor
  Eigen::Matrix<double, 6, 1> tau = apparent_viscosity * strain_rate;

  // von Mises criterion
  // trace of second invariant J2 of deviatoric stress in matrix form
  // Since tau is in Voigt notation, only the first three numbers matter
  // yield condition trace of the invariant > tau0^2
  const double trace_invariant2 = 0.5 * (tau.head(3)).dot(tau.head(3));
  if (trace_invariant2 < (tau0_ * tau0_)) tau.setZero();

  // Update volumetric and deviatoric stress
  // stress = -thermodynamic_pressure I + tau, where I is identity matrix or
  // direc_delta in Voigt notation
  const Eigen::Matrix<double, 6, 1> updated_stress =
      -thermodynamic_pressure * this->dirac_delta() + tau;

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

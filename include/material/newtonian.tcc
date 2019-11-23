//! Constructor with material properties
template <unsigned Tdim>
mpm::Newtonian<Tdim>::Newtonian(unsigned id, const Json& material_properties)
    : Material<Tdim>(id, material_properties) {
  try {
    density_ = material_properties.at("density").template get<double>();
    bulk_modulus_ =
        material_properties.at("bulk_modulus").template get<double>();
    mu_ = material_properties.at("mu").template get<double>();

    properties_ = material_properties;
  } catch (Json::exception& except) {
    console_->error("Material parameter not set: {} {}\n", except.what(),
                    except.id);
  }
}

//! Compute pressure
template <unsigned Tdim>
double mpm::Newtonian<Tdim>::thermodynamic_pressure(
    double volumetric_strain) const {
  // Bulk modulus
  return (-bulk_modulus_ * volumetric_strain);
}

//! Compute stress in 2D
template <>
Eigen::Matrix<double, 6, 1> mpm::Newtonian<2>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain, const ParticleBase<2>* ptr,
    mpm::dense_map* state_vars) {

  const unsigned phase = 0;

  const double pressure = ptr->pressure();

  const double volumetric_strain = dstrain(0) + dstrain(1);

  // Update volumetric and deviatoric stress
  Eigen::Matrix<double, 6, 1> pstress;
  pstress(0) =
      -pressure + (2. * mu_ * dstrain(0)) - (2. * mu_ * volumetric_strain / 3.);
  pstress(1) =
      -pressure + (2. * mu_ * dstrain(1)) - (2. * mu_ * volumetric_strain / 3.);
  pstress(2) = -pressure - (2. * mu_ * volumetric_strain / 3.);
  pstress(3) = mu_ * dstrain(2);
  pstress(4) = 0.0;
  pstress(5) = 0.0;

  return pstress;
}

//! Compute stress in 3D
template <>
Eigen::Matrix<double, 6, 1> mpm::Newtonian<3>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain, const ParticleBase<3>* ptr,
    mpm::dense_map* state_vars) {

  const unsigned phase = 0;

  const double pressure = ptr->pressure();

  const double volumetric_strain = dstrain(0) + dstrain(1) + dstrain(2);

  // Update volumetric and deviatoric stress
  Eigen::Matrix<double, 6, 1> pstress;
  pstress(0) =
      -pressure + (2. * mu_ * dstrain(0)) - (2. * mu_ * volumetric_strain / 3.);
  pstress(1) =
      -pressure + (2. * mu_ * dstrain(1)) - (2. * mu_ * volumetric_strain / 3.);
  pstress(2) =
      -pressure + (2. * mu_ * dstrain(2)) - (2. * mu_ * volumetric_strain / 3.);
  pstress(3) = mu_ * dstrain(3);
  pstress(4) = mu_ * dstrain(4);
  pstress(5) = mu_ * dstrain(5);

  return pstress;
}

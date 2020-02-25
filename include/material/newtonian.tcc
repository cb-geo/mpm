//! Constructor with material properties
template <unsigned Tdim>
mpm::Newtonian<Tdim>::Newtonian(unsigned id, const Json& material_properties)
    : Material<Tdim>(id, material_properties) {
  try {
    density_ = material_properties.at("density").template get<double>();
    bulk_modulus_ =
        material_properties.at("bulk_modulus").template get<double>();
    dynamic_viscosity_ =
        material_properties.at("dynamic_viscosity").template get<double>();

    // Special material properties
    if (material_properties.contains("incompressible")) {
      bool incompressible =
          material_properties.at("incompressible").template get<bool>();
      if (incompressible) compressibility_multiplier_ = 0.0;
    }

    properties_ = material_properties;
  } catch (Json::exception& except) {
    console_->error("Material parameter not set: {} {}\n", except.what(),
                    except.id);
  }
}

//! Initialise history variables
template <unsigned Tdim>
mpm::dense_map mpm::Newtonian<Tdim>::initialise_state_variables() {
  mpm::dense_map state_vars = {{"pressure", 0.0}};
  return state_vars;
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

  // Get strain rate
  const auto& strain_rate = ptr->strain_rate();
  const double volumetric_strain_rate = strain_rate(0) + strain_rate(1);

  // Update pressure
  (*state_vars).at("pressure") +=
      (compressibility_multiplier_ *
       this->thermodynamic_pressure(ptr->dvolumetric_strain()));

  // Volumetric stress component
  const double volumetric_component =
      compressibility_multiplier_ *
      (-(*state_vars).at("pressure") -
       (2. * dynamic_viscosity_ * volumetric_strain_rate / 3.));

  // Update stress component
  Eigen::Matrix<double, 6, 1> pstress;
  pstress(0) = volumetric_component + 2. * dynamic_viscosity_ * strain_rate(0);
  pstress(1) = volumetric_component + 2. * dynamic_viscosity_ * strain_rate(1);
  pstress(2) = volumetric_component;
  pstress(3) = dynamic_viscosity_ * strain_rate(3);
  pstress(4) = 0.;
  pstress(5) = 0.;

  return pstress;
}

//! Compute stress in 3D
template <>
Eigen::Matrix<double, 6, 1> mpm::Newtonian<3>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain, const ParticleBase<3>* ptr,
    mpm::dense_map* state_vars) {

  // Get strain rate
  const auto& strain_rate = ptr->strain_rate();
  const double volumetric_strain_rate =
      strain_rate(0) + strain_rate(1) + strain_rate(2);

  // Update pressure
  (*state_vars).at("pressure") +=
      (compressibility_multiplier_ *
       this->thermodynamic_pressure(ptr->dvolumetric_strain()));

  // Volumetric stress component
  const double volumetric_component =
      compressibility_multiplier_ *
      (-(*state_vars).at("pressure") -
       (2. * dynamic_viscosity_ * volumetric_strain_rate / 3.));

  // Update stress component
  Eigen::Matrix<double, 6, 1> pstress;
  pstress(0) = volumetric_component + 2. * dynamic_viscosity_ * strain_rate(0);
  pstress(1) = volumetric_component + 2. * dynamic_viscosity_ * strain_rate(1);
  pstress(2) = volumetric_component + 2. * dynamic_viscosity_ * strain_rate(2);
  pstress(3) = dynamic_viscosity_ * strain_rate(3);
  pstress(4) = dynamic_viscosity_ * strain_rate(4);
  pstress(5) = dynamic_viscosity_ * strain_rate(5);

  return pstress;
}

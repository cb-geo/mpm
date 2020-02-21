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
    incompressible_ =
        material_properties.at("incompressible").template get<bool>();

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

  // Update volumetric and deviatoric stress
  Eigen::Matrix<double, 6, 1> pstress;
  pstress(0) = 2. * dynamic_viscosity_ * strain_rate(0);
  pstress(1) = 2. * dynamic_viscosity_ * strain_rate(1);
  pstress(2) = 0.;
  pstress(3) = dynamic_viscosity_ * strain_rate(3);
  pstress(4) = 0.;
  pstress(5) = 0.;

  // If weakly compressible
  if (!incompressible_) {
    // Update pressure
    const double volumetric_strain_rate = strain_rate(0) + strain_rate(1);

    (*state_vars).at("pressure") +=
        this->thermodynamic_pressure(ptr->dvolumetric_strain());
    double pressure = (*state_vars).at("pressure");

    // Update volumetric and deviatoric stress
    pstress(0) +=
        -pressure - (2. * dynamic_viscosity_ * volumetric_strain_rate / 3.);
    pstress(1) +=
        -pressure - (2. * dynamic_viscosity_ * volumetric_strain_rate / 3.);
    pstress(2) +=
        -pressure - (2. * dynamic_viscosity_ * volumetric_strain_rate / 3.);
  }

  return pstress;
}

//! Compute stress in 3D
template <>
Eigen::Matrix<double, 6, 1> mpm::Newtonian<3>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain, const ParticleBase<3>* ptr,
    mpm::dense_map* state_vars) {

  // Get strain rate
  const auto& strain_rate = ptr->strain_rate();

  // Update volumetric and deviatoric stress
  Eigen::Matrix<double, 6, 1> pstress;
  pstress(0) = 2. * dynamic_viscosity_ * strain_rate(0);
  pstress(1) = 2. * dynamic_viscosity_ * strain_rate(1);
  pstress(2) = 2. * dynamic_viscosity_ * strain_rate(2);
  pstress(3) = dynamic_viscosity_ * strain_rate(3);
  pstress(4) = dynamic_viscosity_ * strain_rate(4);
  pstress(5) = dynamic_viscosity_ * strain_rate(5);

  // If weakly compressible
  if (!incompressible_) {
    // Update pressure
    const double volumetric_strain_rate =
        strain_rate(0) + strain_rate(1) + strain_rate(2);

    (*state_vars).at("pressure") +=
        this->thermodynamic_pressure(ptr->dvolumetric_strain());
    double pressure = (*state_vars).at("pressure");

    // Update volumetric and deviatoric stress
    pstress(0) +=
        -pressure - (2. * dynamic_viscosity_ * volumetric_strain_rate / 3.);
    pstress(1) +=
        -pressure - (2. * dynamic_viscosity_ * volumetric_strain_rate / 3.);
    pstress(2) +=
        -pressure - (2. * dynamic_viscosity_ * volumetric_strain_rate / 3.);
  }

  return pstress;
}

//! Read material properties
template <unsigned Tdim>
mpm::Rigid<Tdim>::Rigid(unsigned id, const Json& material_properties)
    : Material<Tdim>(id, material_properties) {
  try {
    density_ = material_properties.at("density").template get<double>();
    properties_ = material_properties;

  } catch (Json::exception& except) {
    console_->error("Material parameter not set: {} {}\n", except.what(),
                    except.id);
  }
}

//! Initialise state variables
template <unsigned Tdim>
std::vector<std::string> mpm::Rigid<Tdim>::state_variables() const {
  const std::vector<std::string> state_vars = {};
  return state_vars;
}

//! Compute stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::Rigid<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr, mpm::dense_map* state_vars) {
  const Vector6d dstress = Vector6d::Zero();
  return dstress;
}

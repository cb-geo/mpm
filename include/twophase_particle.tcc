//! Construct a two phase particle with id and coordinates
template <unsigned Tdim>
mpm::TwoPhaseParticle<Tdim>::TwoPhaseParticle(Index id, const VectorDim& coord)
    : mpm::Particl<Tdim>(id, coord) {
  this->initialise_liquid_phase();

  // Set material pointer to null
  liquid_material_ = nullptr;
  // Logger
  std::string logger =
      "twophaseparticle" + std::to_string(Tdim) + "d::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

// Initialise liquid phase particle properties
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::initialise_liquid_phase() {
  liquid_mass_ = 0.;
  liquid_displacement_.setZero();
  liquid_velocity_.setZero();
  liquid_strain_rate_.setZero();
  liquid_strain_.setZero();
  set_traction_ = false;
  liquid_traction_.setZero();
  pore_pressure_ = 0.;
  liquid_saturation_ = 1.;

  this->liquid_properties_["velocities"] = [&]() { return liquid_velocity(); };
  this->liquid_properties_["pressure"] = [&]() {
    Eigen::VectorXd vec_pressure(1);
    vec_pressure << this->pore_pressure();
    return vec_pressure;
  };
}

// Assign a liquid material to particle
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::assign_liquid_material(
    const std::shared_ptr<Material<Tdim>>& material) {
  bool status = false;
  try {
    // Check if material is valid and properties are set
    if (material != nullptr) {
      liquid_material_ = material;
      liquid_material_id_ = material_->id();
      status = true;
    } else {
      throw std::runtime_error("Material is undefined!");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return status;
}

// Assign degree of saturation to the liquid phase
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::assign_saturation_degree() {
  bool status = true;
  try {
    if (material_ != nullptr) {
      liquid_saturation_ = material_
                      ->template property<double>(std::string("saturation"));
      if (liquid_saturation < 0. || liquid_saturation_ > 1.)
        throw std::runtime_error(
            "Particle saturation degree is negative or larger than one");
    } else {
      throw std::runtime_error("Liquid material is invalid");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

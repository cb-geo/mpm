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

// Assign traction to the liquid phase
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::assign_liquid_traction(unsigned direction,
                                                         double traction) {
  bool status = false;
  try {
    if (direction >= Tdim ||
        this->volume_ == std::numeric_limits<double>::max()) {
      throw std::runtime_error(
          "Particle liquid traction property: volume / direction is invalid");
    }
    // Assign liquid traction
    liquid_traction_(direction) =
        traction * this->volume_ / this->size_(direction);
    status = true;
    this->set_liquid_traction_ = true;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Compute mass of particle
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::compute_liquid_mass() {
  bool status = true;
  try {
    // Check if particle volume is set and liquid material ptr is valid
    if (volume_ != std::numeric_limits<double>::max() &&
        liquid_material_ != nullptr) {
      // Mass = volume of particle * bulk_density
      this->liquid_mass_density_ =
          liquid_saturation_ * porosity_ *
          liquid_material_->template property<double>(std::string("density"));
      this->liquid_mass_ = volume_ * liquid_mass_density_;
    } else {
      throw std::runtime_error(
          "Particle volume or density is invalid! cannot compute mass for the "
          "liquid particle");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Map particle mass and momentum to nodes
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::map_liquid_mass_momentum_to_nodes() {
  bool status = true;
  try {
    // Check if particle mass is set and positive
    if (liquid_mass_ != std::numeric_limits<double>::max()) {
      // Map particle mass and momentum to nodes
      this->cell_->map_mass_momentum_to_nodes(
          this->shapefn_, mpm::ParticlePhase::Liquid, liquid_mass_,
          liquid_velocity_);
    } else {
      throw std::runtime_error("Particle mass is not set or negative");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Compute pore pressure
bool mpm::TwoPhaseParticle<Tdim>::compute_pore_pressure(double dt) {
  bool status = true;
  try {
    if (liquid_material_ != nullptr && cell_ != nullptr) {
      // get the bulk modulus of liquid
      double K = liquid_material_->template property<double>(
          std::string("bulk_modulus"));
      // get liquid phase strain rate at cell centre
      Eigen::VectorXd liquid_strain_rate_centroid =
          cell_->compute_strain_rate_centroid(mpm::ParticlePhase::Liquid);
      // update pressure
      this->pore_pressure_ +=
          -dt * (K / porosity_) *
          ((1 - porosity_) * strain_rate_.head(Tdim).sum() +
           porosity_ * liquid_strain_rate_centroid.head(Tdim).sum());
    } else
      throw std::runtime_error("Liquid material is not set");
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Map particle pore liquid pressure to nodes
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::map_pore_pressure_to_nodes() {
  bool status = true;
  try {
    // Check if particle mass is set
    if (liquid_mass_ != std::numeric_limits<double>::max()) {
      // Map particle liquid mass and pore pressure to nodes
      this->cell_->map_pressure_to_nodes(this->shapefn_,
                                         mpm::ParticlePhase::Liquid,
                                         liquid_mass_, pore_pressure_);
    } else {
      throw std::runtime_error("Particleliquid mass has not been computed");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Compute pore liquid pressure smoothing based on nodal pressure
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::compute_pore_pressure_smoothing() {
  bool status = true;
  try {
    // Check if particle has a valid cell ptr
    if (cell_ != nullptr)
      // Update pore liquid pressure to interpolated nodal pressure
      this->pore_pressure_ = cell_->interpolate_nodal_pressure(
          shapefn_, mpm::ParticlePhase::Liquid);
    else
      throw std::runtime_error(
          "Cell is not initialised! "
          "cannot compute pore liquid pressure smoothing of the particle");

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Map liquid phase body force
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_liquid_body_force(
    const VectorDim& pgravity) {
  // Compute nodal liquid body forces
  cell_->compute_nodal_body_force(this->shapefn_, mpm::ParticlePhase::Liquid,
                                  this->liquid_mass_, pgravity);
}

//! Map liquid phase traction force
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_liquid_traction_force() {
  if (this->set_liquid_traction_)
    // Map particle liquid phase traction forces to nodes
    cell_->compute_nodal_traction_force(this->shapefn_,
                                        mpm::ParticlePhase::Liquid,
                                        -1. * porosity_ * this->traction_);
}


//! Map liquid phase internal force
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::map_liquid_internal_force() {
  bool status = true;
  try {
    // initialise a vector of pore pressure
    Eigen::Matrix<double, 6, 1> pressure;
    pressure.setZero();
    pressure(0) = pressure(1) = pressure(2) = pore_pressure_;
    // Compute nodal liquid phase  internal forces
    // porosity * pressure * volume
    cell_->compute_nodal_internal_force(
        this->bmatrix_, mpm::ParticlePhase::Liquid, this->volume_,
        porosity_ * this->stress_);
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

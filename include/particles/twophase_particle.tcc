//! Construct a twophase particle with id and coordinates
template <unsigned Tdim>
mpm::TwoPhaseParticle<Tdim>::TwoPhaseParticle(Index id, const VectorDim& coord)
    : mpm::Particle<Tdim>(id, coord) {
  // Initialise variables for solid phase
  mpm::Particle<Tdim>::initialise();
  // Initialise variables for liquid phase
  this->initialise_liquid_phase();
  // Clear cell ptr
  cell_ = nullptr;
  // Nodes
  nodes_.clear();
  // Set material pointer to null
  material_ = nullptr;
  liquid_material_ = nullptr;
  // Logger
  std::string logger =
      "twophaseparticle" + std::to_string(Tdim) + "d::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

//! Construct a twophase particle with id, coordinates and status
template <unsigned Tdim>
mpm::TwoPhaseParticle<Tdim>::TwoPhaseParticle(Index id, const VectorDim& coord,
                                              bool status)
    : mpm::Particle<Tdim>(id, coord, status) {
  // Initialise variables for solid phase
  mpm::Particle<Tdim>::initialise();
  // Initialise variables for liquid phase
  this->initialise_liquid_phase();
  // Clear cell ptr
  cell_ = nullptr;
  // Nodes
  nodes_.clear();
  // Set material pointer to null
  material_ = nullptr;
  liquid_material_ = nullptr;
  // Logger
  std::string logger =
      "twophaseparticle" + std::to_string(Tdim) + "d::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

//! Initialise particle data from HDF5
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::initialise_particle(
    const HDF5Particle& particle) {
  // Derive from particle
  mpm::Particle<Tdim>::initialise_particle(particle);

  // TODO:HDF5
  // // Liquid mass
  // this->liquid_mass_ = particle.liquid_mass;
  // // Liquid mass Density
  // this->liquid_mass_density_ = particle.mass / particle.volume;
  // // Pore pressure
  // this->pore_pressure_ = particle.pore_pressure;
  // // Liquid velocity
  // Eigen::Vector3d liquid_velocity;
  // liquid_velocity << particle.liquid_velocity_x, particle.liquid_velocity_y,
  //     particle.liquid_velocity_z;
  // // Initialise velocity
  // for (unsigned i = 0; i < Tdim; ++i)
  //   this->liquid_velocity_(i) = liquid_velocity(i);
  // // Liquid strain
  // this->liquid_strain_[0] = particle.liquid_strain_xx;
  // this->liquid_strain_[1] = particle.liquid_strain_yy;
  // this->liquid_strain_[2] = particle.liquid_strain_zz;
  // this->liquid_strain_[3] = particle.liquid_gamma_xy;
  // this->liquid_strain_[4] = particle.liquid_gamma_yz;
  // this->liquid_strain_[5] = particle.liquid_gamma_xz;
  // // Liquid material id
  // this->liquid_material_id_ = particle.liquid_material_id;

  return true;
}

//! Initialise particle data from HDF5
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::initialise_particle(
    const HDF5Particle& particle,
    const std::shared_ptr<mpm::Material<Tdim>>& material) {
  bool status = this->initialise_particle(particle);
  if (material != nullptr) {
    if (this->material_id_ == material->id() ||
        this->material_id_ == std::numeric_limits<unsigned>::max()) {
      material_ = material;
      // Reinitialize state variables
      auto mat_state_vars = material_->initialise_state_variables();
      if (mat_state_vars.size() == particle.nstate_vars) {
        unsigned i = 0;
        auto state_variables = material_->state_variables();
        for (const auto& state_var : state_variables) {
          this->state_variables_.at(state_var) = particle.svars[i];
          ++i;
        }
      }
    } else {
      status = false;
      throw std::runtime_error("Material is invalid to assign to particle!");
    }
  }
  return status;
}

//! Return particle data in HDF5 format
template <unsigned Tdim>
// cppcheck-suppress *
mpm::HDF5Particle mpm::TwoPhaseParticle<Tdim>::hdf5() const {
  // Derive from particle
  auto particle_data = mpm::Particle<Tdim>::hdf5();
  // TODO:HDF5
  // // Particle liquid velocity
  // Eigen::Vector3d liquid_velocity;
  // liquid_velocity.setZero();
  // for (unsigned j = 0; j < Tdim; ++j)
  //   liquid_velocity[j] = this->liquid_velocity_[j];
  // // Particle liquid strain
  // Eigen::Matrix<double, 6, 1> liquid_strain = this->liquid_strain_;
  // // Particle liquid mass
  // particle_data.liquid_mass = this->liquid_mass_;
  // // Particle pore pressure
  // particle_data.pore_pressure = this->pore_pressure_;
  // // Particle liquid velocity
  // particle_data.liquid_velocity_x = liquid_velocity[0];
  // particle_data.liquid_velocity_y = liquid_velocity[1];
  // particle_data.liquid_velocity_z = liquid_velocity[2];
  // // Particle liquid strain
  // particle_data.liquid_strain_xx = liquid_strain[0];
  // particle_data.liquid_strain_yy = liquid_strain[1];
  // particle_data.liquid_strain_zz = liquid_strain[2];
  // particle_data.liquid_gamma_xy = liquid_strain[3];
  // particle_data.liquid_gamma_yz = liquid_strain[4];
  // particle_data.liquid_gamma_xz = liquid_strain[5];
  // // Particle liquid material id
  // particle_data.liquid_material_id = this->liquid_material_id_;

  return particle_data;
}

// Initialise liquid phase particle properties
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::initialise_liquid_phase() {
  liquid_mass_ = 0.;
  pore_pressure_ = 0.;
  porosity_ = 0.;
  liquid_velocity_.setZero();
  liquid_strain_rate_.setZero();
  liquid_strain_.setZero();
  liquid_traction_.setZero();
  permeability_c1_.setZero();
  set_liquid_traction_ = false;

  // Initialize vector data properties
  this->properties_["liquid_strains"] = [&]() { return liquid_strain(); };
  this->properties_["liquid_velocities"] = [&]() { return liquid_velocity(); };
  this->properties_["pore_pressure"] = [&]() {
    Eigen::VectorXd pore_pressure = Eigen::VectorXd::Zero(3);
    // Total pore pressure
    pore_pressure[0] = this->pore_pressure();
    // Excessive pore pressure
    pore_pressure[1] = this->excessive_pore_pressure();
    // Free surface
    pore_pressure[2] = this->free_surface();

    return pore_pressure;
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
      liquid_material_id_ = liquid_material_->id();
      status = true;
    } else {
      throw std::runtime_error("Liquid material is undefined!");
    }

    // Assign porosity
    this->assign_porosity();
    // Assign permeability
    this->assign_permeability();

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return status;
}

//! Assign particle permeability
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::assign_permeability() {
  bool status = true;
  try {
    // Check if material ptr is valid
    if (material_ != nullptr) {
      // Porosity parameter k_p
      const double k_p =
          std::pow(this->porosity_, 3) / std::pow((1. - this->porosity_), 2);
      // Different dimensions permeability
      switch (Tdim) {
        case (1): {
          permeability_c1_(0) =
              material_->template property<double>("k_x") / k_p;
          break;
        }
        case (2): {
          permeability_c1_(0) =
              material_->template property<double>("k_x") / k_p;
          permeability_c1_(1) =
              material_->template property<double>("k_y") / k_p;
          break;
        }
        default: {
          permeability_c1_(0) =
              material_->template property<double>("k_x") / k_p;
          permeability_c1_(1) =
              material_->template property<double>("k_y") / k_p;
          permeability_c1_(2) =
              material_->template property<double>("k_z") / k_p;
          break;
        }
      }
    } else {
      throw std::runtime_error(
          "Material is invalid, could not assign permeability");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Assign particle porosity
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::assign_porosity() {
  bool status = true;
  try {
    // Check if material ptr is valid
    if (material_ != nullptr) {
      porosity_ = material_->template property<double>(std::string("porosity"));
      // Check if the porosity value is valid
      if (porosity_ < 0. || porosity_ > 1.)
        throw std::runtime_error(
            "Particle porosity is negative or larger than one");
    } else {
      throw std::runtime_error(
          "Material is invalid, could not assign porosity");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Assign velocity to the particle liquid phase
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::assign_liquid_velocity(
    const Eigen::Matrix<double, Tdim, 1>& velocity) {
  bool status = false;
  try {
    // Assign velocity
    liquid_velocity_ = velocity;
    status = true;
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

// Compute mass of particle (both solid and liquid)
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::compute_mass() noexcept {
  // Compute mass of particle (solid phase)
  mpm::Particle<Tdim>::compute_mass();
  // Compute mass of particle (liquid phase)
  this->compute_liquid_mass();
}

// Compute liquid mass of particle
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::compute_liquid_mass() noexcept {
  // Check if particle volume is set and liquid material ptr is valid
  assert(volume_ != std::numeric_limits<double>::max() &&
         liquid_material_ != nullptr);

  // Mass = volume of particle * bulk_density
  this->liquid_mass_density_ =
      porosity_ *
      liquid_material_->template property<double>(std::string("density"));
  this->liquid_mass_ = volume_ * liquid_mass_density_;
}

//! Map particle mass and momentum to nodes
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_mass_momentum_to_nodes() noexcept {
  // Map particle mass and momentum to nodes (solid phase)
  mpm::Particle<Tdim>::map_mass_momentum_to_nodes();
  // Map particle mass and momentum to nodes (liquid phase)
  this->map_liquid_mass_momentum_to_nodes();
}

//! Map liquid mass and momentum to nodes
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_liquid_mass_momentum_to_nodes() noexcept {
  // Check if liquid mass is set and positive
  assert(liquid_mass_ != std::numeric_limits<double>::max());

  // Map liquid mass and momentum to nodes
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    nodes_[i]->update_mass(true, mpm::ParticlePhase::Liquid,
                           liquid_mass_ * shapefn_[i]);
    nodes_[i]->update_momentum(true, mpm::ParticlePhase::Liquid,
                               liquid_mass_ * shapefn_[i] * liquid_velocity_);
  }
}

//! Compute pore pressure (compressible fluid)
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::compute_pore_pressure(double dt) noexcept {
  // Check if liquid material and cell pointer are set and positive
  assert(liquid_material_ != nullptr && cell_ != nullptr);
  // Apply free surface
  if (this->free_surface()) this->pore_pressure_ = 0.0;
  // Compute pore pressure
  else {
    // get the bulk modulus of liquid
    double K = liquid_material_->template property<double>(
        std::string("bulk_modulus"));
    // Compute strain rate of liquid phase at centroid
    auto liquid_strain_rate_centroid = mpm::Particle<Tdim>::compute_strain_rate(
        dn_dx_centroid_, mpm::ParticlePhase::Liquid);
    // update pressure
    this->pore_pressure_ +=
        -dt * (K / porosity_) *
        ((1 - porosity_) * strain_rate_.head(Tdim).sum() +
         porosity_ * liquid_strain_rate_centroid.head(Tdim).sum());
  }
}

//! Map particle pore liquid pressure to nodes
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_pore_pressure_to_nodes() noexcept {
  // Check if particle mass is set
  assert(liquid_mass_ != std::numeric_limits<double>::max());
  // Map particle liquid mass and pore pressure to nodes
  for (unsigned i = 0; i < nodes_.size(); ++i)
    nodes_[i]->update_mass_pressure(
        mpm::ParticlePhase::Liquid,
        shapefn_[i] * liquid_mass_ * pore_pressure_);
}

// Compute pore liquid pressure smoothing based on nodal pressure
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::compute_pore_pressure_smoothing() noexcept {
  // Check if particle has a valid cell ptr
  assert(cell_ != nullptr);

  bool status = false;
  // Check if particle has a valid cell ptr
  if (cell_ != nullptr) {
    // Apply free surface
    if (this->free_surface()) this->pore_pressure_ = 0.0;
    // Compute particle pore pressure from nodal values
    else {
      double pore_pressure = 0;
      for (unsigned i = 0; i < nodes_.size(); ++i)
        pore_pressure +=
            shapefn_(i) * nodes_[i]->pressure(mpm::ParticlePhase::Liquid);
      // Update pore liquid pressure to interpolated nodal pressure
      this->pore_pressure_ = pore_pressure;
    }
  }
  return status;
}

//! Map body force for both mixture and liquid
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_body_force(
    const VectorDim& pgravity) noexcept {
  //! Map body force for mixture
  this->map_mixture_body_force(mpm::ParticlePhase::Mixture, pgravity);
  //! Map body force for liquid
  this->map_liquid_body_force(pgravity);
}

//! Map liquid phase body force
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_liquid_body_force(
    const VectorDim& pgravity) noexcept {
  // Compute nodal liquid body forces
  for (unsigned i = 0; i < nodes_.size(); ++i)
    nodes_[i]->update_external_force(
        true, mpm::ParticlePhase::Liquid,
        (pgravity * this->liquid_mass_ * shapefn_(i)));
}

//! Map mixture body force
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_mixture_body_force(
    unsigned mixture, const VectorDim& pgravity) noexcept {
  // Compute nodal mixture body forces
  for (unsigned i = 0; i < nodes_.size(); ++i)
    nodes_[i]->update_external_force(
        true, mixture,
        (pgravity * (this->liquid_mass_ + this->mass_) * shapefn_(i)));
}

//! Map liquid traction force
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_liquid_traction_force() noexcept {
  if (this->set_liquid_traction_) {
    // Map particle liquid traction forces to nodes
    for (unsigned i = 0; i < nodes_.size(); ++i)
      nodes_[i]->update_external_force(
          true, mpm::ParticlePhase::Liquid,
          (-1. * shapefn_[i] * porosity_ * this->liquid_traction_));
  }
}

//! Map both mixture and liquid internal force
template <unsigned Tdim>
inline void mpm::TwoPhaseParticle<Tdim>::map_internal_force() noexcept {
  //! Map mixture internal force
  mpm::TwoPhaseParticle<Tdim>::map_mixture_internal_force(
      mpm::ParticlePhase::Mixture);
  //! Map liquid internal force
  mpm::TwoPhaseParticle<Tdim>::map_liquid_internal_force();
}

//! Map liquid phase internal force
template <>
void mpm::TwoPhaseParticle<2>::map_liquid_internal_force() noexcept {
  // initialise a vector of pore pressure
  Eigen::Matrix<double, 6, 1> pressure;
  pressure.setZero();
  pressure(0) = -this->pore_pressure_;
  pressure(1) = -this->pore_pressure_;

  // Compute nodal internal forces
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    // Compute force: -pstress * volume
    Eigen::Matrix<double, 2, 1> force;
    force[0] = dn_dx_(i, 0) * pressure[0] * this->porosity_;
    force[1] = dn_dx_(i, 1) * pressure[1] * this->porosity_;

    force *= -1. * this->volume_;

    nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Liquid, force);
  }
}

template <>
void mpm::TwoPhaseParticle<3>::map_liquid_internal_force() noexcept {
  // initialise a vector of pore pressure
  Eigen::Matrix<double, 6, 1> pressure;
  pressure.setZero();
  pressure(0) = -this->pore_pressure_;
  pressure(1) = -this->pore_pressure_;
  pressure(2) = -this->pore_pressure_;

  // Compute nodal internal forces
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    // Compute force: -pstress * volume
    Eigen::Matrix<double, 3, 1> force;
    force[0] = dn_dx_(i, 0) * pressure[0] * this->porosity_;
    force[1] = dn_dx_(i, 1) * pressure[1] * this->porosity_;
    force[2] = dn_dx_(i, 2) * pressure[2] * this->porosity_;

    force *= -1. * this->volume_;

    nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Liquid, force);
  }
}

//! Map mixture internal force
template <>
void mpm::TwoPhaseParticle<2>::map_mixture_internal_force(
    unsigned mixture) noexcept {
  // initialise a vector of pore pressure
  Eigen::Matrix<double, 6, 1> total_stress = this->stress_;
  total_stress(0) -= this->pore_pressure_;
  total_stress(1) -= this->pore_pressure_;

  // Compute nodal internal forces
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    // Compute force: -pstress * volume
    Eigen::Matrix<double, 2, 1> force;
    force[0] = dn_dx_(i, 0) * total_stress[0] + dn_dx_(i, 1) * total_stress[3];
    force[1] = dn_dx_(i, 1) * total_stress[1] + dn_dx_(i, 0) * total_stress[3];

    force *= -1. * this->volume_;

    nodes_[i]->update_internal_force(true, mixture, force);
  }
}

template <>
void mpm::TwoPhaseParticle<3>::map_mixture_internal_force(
    unsigned mixture) noexcept {
  // initialise a vector of pore pressure
  Eigen::Matrix<double, 6, 1> total_stress = this->stress_;
  total_stress(0) -= this->pore_pressure_;
  total_stress(1) -= this->pore_pressure_;
  total_stress(2) -= this->pore_pressure_;

  // Compute nodal internal forces
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    // Compute force: -pstress * volume
    Eigen::Matrix<double, 3, 1> force;
    force[0] = dn_dx_(i, 0) * total_stress[0] + dn_dx_(i, 1) * total_stress[3] +
               dn_dx_(i, 2) * total_stress[5];

    force[1] = dn_dx_(i, 1) * total_stress[1] + dn_dx_(i, 0) * total_stress[3] +
               dn_dx_(i, 2) * total_stress[4];

    force[2] = dn_dx_(i, 2) * total_stress[2] + dn_dx_(i, 1) * total_stress[4] +
               dn_dx_(i, 0) * total_stress[5];

    force *= -1. * this->volume_;

    nodes_[i]->update_internal_force(true, mixture, force);
  }
}

//! Map drag force coefficient
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::map_drag_force_coefficient() {
  bool status = true;
  try {
    // Update permeability
    auto permeability = this->update_permeability();
    // Initialise drag force coefficient
    VectorDim drag_force_coefficient;
    drag_force_coefficient.setZero();

    // Check if permeability coefficient is valid
    for (unsigned i = 0; i < Tdim; ++i) {
      drag_force_coefficient(i) =
          porosity_ * porosity_ * 9.81 *
          liquid_material_->template property<double>(std::string("density")) /
          permeability(i);
    }

    // Map drag forces from particle to nodes
    for (unsigned j = 0; j < nodes_.size(); ++j)
      nodes_[j]->update_drag_force_coefficient(
          true, drag_force_coefficient * this->volume_ * shapefn_(j));

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Update particle permeability
template <unsigned Tdim>
Eigen::Matrix<double, Tdim, 1>
    mpm::TwoPhaseParticle<Tdim>::update_permeability() {
  // Initialise permeability
  Eigen::Matrix<double, Tdim, 1> permeability;
  permeability.setZero();
  try {
    // Porosity parameter
    const double k_p =
        std::pow(this->porosity_, 3) / std::pow((1. - this->porosity_), 2);
    // Update permeability by KC equation
    permeability = k_p * permeability_c1_;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return permeability;
}

// Update particle porosity
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::update_porosity(double dt) {
  bool status = true;
  try {
    // Update particle porosity
    this->porosity_ =
        1 - (1 - this->porosity_) / (1 + dt * strain_rate_.head(Tdim).sum());
    // Check if the value is valid
    if (this->porosity_ < 0) this->porosity_ = 1E-5;
    if (this->porosity_ < 1) this->porosity_ = 1 - 1E-5;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Compute updated position and velocities
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::compute_updated_position(
    double dt, bool velocity_update) noexcept {
  // Compute updated position
  mpm::Particle<Tdim>::compute_updated_position(dt, velocity_update);
  // Compute updated liquid velocities
  this->compute_updated_liquid_velocity(dt, velocity_update);
}

// Compute updated velocity of the liquid phase based on nodal velocity
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::compute_updated_liquid_velocity(
    double dt, bool velocity_update) noexcept {
  // Check if particle has a valid cell ptr
  assert(cell_ != nullptr);

  if (!velocity_update) {
    // Get interpolated nodal acceleration
    Eigen::Matrix<double, Tdim, 1> acceleration;
    acceleration.setZero();

    for (unsigned i = 0; i < nodes_.size(); ++i)
      acceleration +=
          shapefn_(i) * nodes_[i]->acceleration(mpm::ParticlePhase::Liquid);

    // Update particle velocity from interpolated nodal acceleration
    this->liquid_velocity_ += acceleration * dt;
  } else {
    // Get interpolated nodal velocity
    Eigen::Matrix<double, Tdim, 1> velocity;
    velocity.setZero();

    for (unsigned i = 0; i < nodes_.size(); ++i)
      velocity += shapefn_(i) * nodes_[i]->velocity(mpm::ParticlePhase::Liquid);

    // Update particle velocity to interpolated nodal velocity
    this->liquid_velocity_ = velocity;
  }

  // Apply particle velocity constraints
  this->apply_particle_liquid_velocity_constraints();
}

//! Assign particle liquid phase velocity constraint
//! Constrain directions can take values between 0 and Dim
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::assign_particle_liquid_velocity_constraint(
    unsigned dir, double velocity) {
  bool status = true;
  try {
    //! Constrain directions can take values between 0 and Dim
    if (dir < Tdim)
      this->liquid_velocity_constraints_.insert(
          std::make_pair<unsigned, double>(static_cast<unsigned>(dir),
                                           static_cast<double>(velocity)));
    else
      throw std::runtime_error(
          "Particle liquid velocity constraint direction is out of bounds");

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Apply particle velocity constraints
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::apply_particle_liquid_velocity_constraints() {
  // Set particle velocity constraint
  for (const auto& constraint : this->liquid_velocity_constraints_) {
    // Direction value in the constraint (0, Dim)
    const unsigned dir = constraint.first;
    // Direction: dir % Tdim (modulus)
    const auto direction = static_cast<unsigned>(dir % Tdim);
    this->liquid_velocity_(direction) = constraint.second;
  }
}

//! Assign particle pressure constraints
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::assign_particle_pore_pressure_constraint(
    double pressure) {
  bool status = true;
  try {
    this->pore_pressure_constraint_ = pressure;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Initial pore pressure by water table
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::initialise_pore_pressure_watertable(
    const unsigned dir_v, const unsigned dir_h,
    std::map<double, double>& refernece_points) {
  bool status = true;
  try {
    // Initialise left boundary position (coordinate) and h0
    double left_boundary = std::numeric_limits<double>::lowest();
    double h0_left = 0.;
    // Initialise right boundary position (coordinate) and h0
    double right_boundary = std::numeric_limits<double>::max();
    double h0_right = 0.;
    // Position and h0 of particle (coordinate)
    const double position = this->coordinates_(dir_h);
    // Iterate over each refernece_points
    for (const auto& refernece_point : refernece_points) {
      // Find boundary
      if (refernece_point.first > left_boundary &&
          refernece_point.first <= position) {
        // Left boundary position and h0
        left_boundary = refernece_point.first;
        h0_left = refernece_point.second;
      } else if (refernece_point.first > position &&
                 refernece_point.first <= right_boundary) {
        // Right boundary position and h0
        right_boundary = refernece_point.first;
        h0_right = refernece_point.second;
      }
    }
    // Check if the boundaries are assigned
    if (left_boundary != std::numeric_limits<double>::lowest()) {
      // Particle with left and right boundary
      if (right_boundary != std::numeric_limits<double>::max()) {
        this->pore_pressure_ =
            ((h0_right - h0_left) / (right_boundary - left_boundary) *
                 (position - left_boundary) +
             h0_left - this->coordinates_(dir_v)) *
            1000 * 9.81;
      } else
        // Particle with only left boundary
        this->pore_pressure_ =
            (h0_left - this->coordinates_(dir_v)) * 1000 * 9.81;
    }
    // Particle with only right boundary
    else if (right_boundary != std::numeric_limits<double>::max())
      this->pore_pressure_ =
          (h0_right - this->coordinates_(dir_v)) * 1000 * 9.81;

    else
      throw std::runtime_error(
          "Particle pore pressure can not be initialised by water table");
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}
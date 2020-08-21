//! Construct a two phase particle with id and coordinates
template <unsigned Tdim>
mpm::TwoPhaseParticle<Tdim>::TwoPhaseParticle(Index id, const VectorDim& coord)
    : mpm::Particle<Tdim>(id, coord) {
  this->initialise();
  // Clear cell ptr
  cell_ = nullptr;
  // Nodes
  nodes_.clear();
  // Set material containers
  this->initialise_material(2);
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
  this->initialise();
  cell_ = nullptr;
  nodes_.clear();
  // Set material containers
  this->initialise_material(2);
  //! Logger
  std::string logger =
      "twophaseparticle" + std::to_string(Tdim) + "d::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

//! Return particle data as HDF5 pointer
template <unsigned Tdim>
// cppcheck-suppress *
std::shared_ptr<void> mpm::TwoPhaseParticle<Tdim>::hdf5_ptr() {
  // Initialise particle_data
  auto particle_data = std::make_shared<mpm::HDF5ParticleTwoPhase>();

  Eigen::Vector3d coordinates;
  coordinates.setZero();
  for (unsigned j = 0; j < Tdim; ++j) coordinates[j] = this->coordinates_[j];

  Eigen::Vector3d displacement;
  displacement.setZero();
  for (unsigned j = 0; j < Tdim; ++j) displacement[j] = this->displacement_[j];

  Eigen::Vector3d velocity;
  velocity.setZero();
  for (unsigned j = 0; j < Tdim; ++j) velocity[j] = this->velocity_[j];

  // Particle local size
  Eigen::Vector3d nsize;
  nsize.setZero();
  Eigen::VectorXd size = this->natural_size();
  for (unsigned j = 0; j < Tdim; ++j) nsize[j] = size[j];

  Eigen::Matrix<double, 6, 1> stress = this->stress_;

  Eigen::Matrix<double, 6, 1> strain = this->strain_;

  particle_data->id = this->id();
  particle_data->mass = this->mass();
  particle_data->volume = this->volume();
  particle_data->pressure =
      (state_variables_[mpm::ParticlePhase::Solid].find("pressure") !=
       state_variables_[mpm::ParticlePhase::Solid].end())
          ? state_variables_[mpm::ParticlePhase::Solid].at("pressure")
          : 0.;

  particle_data->coord_x = coordinates[0];
  particle_data->coord_y = coordinates[1];
  particle_data->coord_z = coordinates[2];

  particle_data->displacement_x = displacement[0];
  particle_data->displacement_y = displacement[1];
  particle_data->displacement_z = displacement[2];

  particle_data->nsize_x = nsize[0];
  particle_data->nsize_y = nsize[1];
  particle_data->nsize_z = nsize[2];

  particle_data->velocity_x = velocity[0];
  particle_data->velocity_y = velocity[1];
  particle_data->velocity_z = velocity[2];

  particle_data->stress_xx = stress[0];
  particle_data->stress_yy = stress[1];
  particle_data->stress_zz = stress[2];
  particle_data->tau_xy = stress[3];
  particle_data->tau_yz = stress[4];
  particle_data->tau_xz = stress[5];

  particle_data->strain_xx = strain[0];
  particle_data->strain_yy = strain[1];
  particle_data->strain_zz = strain[2];
  particle_data->gamma_xy = strain[3];
  particle_data->gamma_yz = strain[4];
  particle_data->gamma_xz = strain[5];

  particle_data->epsilon_v = this->volumetric_strain_centroid_;

  particle_data->status = this->status();

  particle_data->cell_id = this->cell_id();

  particle_data->material_id = this->material_id();

  // Write state variables
  if (this->material() != nullptr) {
    particle_data->nstate_vars =
        state_variables_[mpm::ParticlePhase::Solid].size();
    if (state_variables_[mpm::ParticlePhase::Solid].size() > 20)
      throw std::runtime_error("# of state variables cannot be more than 20");
    unsigned i = 0;
    auto state_variables = (this->material())->state_variables();
    for (const auto& state_var : state_variables) {
      particle_data->svars[i] =
          state_variables_[mpm::ParticlePhase::Solid].at(state_var);
      ++i;
    }
  }

  // Particle liquid mass
  particle_data->liquid_mass = this->liquid_mass_;

  // Particle liquid velocity
  Eigen::Vector3d liquid_velocity;
  liquid_velocity.setZero();
  for (unsigned j = 0; j < Tdim; ++j)
    liquid_velocity[j] = this->liquid_velocity_[j];

  particle_data->liquid_velocity_x = liquid_velocity[0];
  particle_data->liquid_velocity_y = liquid_velocity[1];
  particle_data->liquid_velocity_z = liquid_velocity[2];

  // Particle porosity and saturation
  particle_data->porosity = this->porosity_;
  particle_data->liquid_saturation = this->liquid_saturation_;

  // Particle liquid material id
  particle_data->liquid_material_id =
      this->material_id(mpm::ParticlePhase::Liquid);

  // Write state variables
  if (this->material(mpm::ParticlePhase::Liquid) != nullptr) {
    particle_data->nliquid_state_vars =
        state_variables_[mpm::ParticlePhase::Liquid].size();
    if (state_variables_[mpm::ParticlePhase::Liquid].size() > 5)
      throw std::runtime_error("# of state variables cannot be more than 5");
    unsigned i = 0;
    auto state_variables =
        (this->material(mpm::ParticlePhase::Liquid))->state_variables();
    for (const auto& state_var : state_variables) {
      particle_data->liquid_svars[i] =
          state_variables_[mpm::ParticlePhase::Liquid].at(state_var);
      ++i;
    }
  }

  return particle_data;
}

//! Initialise particle data from HDF5
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::initialise_particle(HDF5Particle& particle) {
  // Initialise solid phase
  bool status = mpm::Particle<Tdim>::initialise_particle(particle);
  auto twophase_particle = reinterpret_cast<HDF5ParticleTwoPhase*>(&particle);

  // Liquid mass
  this->liquid_mass_ = twophase_particle->liquid_mass;
  // Liquid mass Density
  this->liquid_mass_density_ = twophase_particle->liquid_mass / particle.volume;

  // Liquid velocity
  Eigen::Vector3d liquid_velocity;
  liquid_velocity << twophase_particle->liquid_velocity_x,
      twophase_particle->liquid_velocity_y,
      twophase_particle->liquid_velocity_z;
  // Initialise velocity
  for (unsigned i = 0; i < Tdim; ++i)
    this->liquid_velocity_(i) = liquid_velocity(i);

  // Particle porosity and saturation
  this->porosity_ = twophase_particle->porosity;
  this->liquid_saturation_ = twophase_particle->liquid_saturation;
  this->assign_permeability();

  // Liquid material id
  this->material_id_[mpm::ParticlePhase::Liquid] =
      twophase_particle->liquid_material_id;

  return status;
}

//! Initialise particle data from HDF5
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::initialise_particle(
    HDF5Particle& particle,
    const std::vector<std::shared_ptr<mpm::Material<Tdim>>>& materials) {
  auto twophase_particle = reinterpret_cast<HDF5ParticleTwoPhase*>(&particle);
  bool status = this->initialise_particle(*twophase_particle);

  assert(materials.size() == 2);

  // Solid Phase
  const auto& solid_material = materials.at(mpm::ParticlePhase::Solid);
  if (solid_material != nullptr) {
    if (this->material_id(mpm::ParticlePhase::Solid) == solid_material->id() ||
        this->material_id(mpm::ParticlePhase::Solid) ==
            std::numeric_limits<unsigned>::max()) {
      bool assign_mat =
          this->assign_material(solid_material, mpm::ParticlePhase::Solid);
      if (!assign_mat) throw std::runtime_error("Material assignment failed");
      // Reinitialize state variables
      auto mat_state_vars = (this->material(mpm::ParticlePhase::Solid))
                                ->initialise_state_variables();
      if (mat_state_vars.size() == twophase_particle->nstate_vars) {
        unsigned i = 0;
        auto state_variables =
            (this->material(mpm::ParticlePhase::Solid))->state_variables();
        for (const auto& state_var : state_variables) {
          this->state_variables_[mpm::ParticlePhase::Solid].at(state_var) =
              twophase_particle->svars[i];
          ++i;
        }
      }
    } else {
      status = false;
      throw std::runtime_error("Material is invalid to assign to particle!");
    }
  }

  // Fluid Phase
  const auto& liquid_material = materials.at(mpm::ParticlePhase::Liquid);
  if (liquid_material != nullptr) {
    if (this->material_id(mpm::ParticlePhase::Liquid) ==
            liquid_material->id() ||
        this->material_id(mpm::ParticlePhase::Liquid) ==
            std::numeric_limits<unsigned>::max()) {
      bool assign_mat =
          this->assign_material(liquid_material, mpm::ParticlePhase::Liquid);
      if (!assign_mat) throw std::runtime_error("Material assignment failed");
      // Reinitialize state variables
      auto mat_state_vars = (this->material(mpm::ParticlePhase::Liquid))
                                ->initialise_state_variables();
      if (mat_state_vars.size() == twophase_particle->nliquid_state_vars) {
        unsigned i = 0;
        auto state_variables =
            (this->material(mpm::ParticlePhase::Liquid))->state_variables();
        for (const auto& state_var : state_variables) {
          this->state_variables_[mpm::ParticlePhase::Liquid].at(state_var) =
              twophase_particle->liquid_svars[i];
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

// Initialise liquid phase particle properties
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::initialise() {
  mpm::Particle<Tdim>::initialise();
  liquid_mass_ = 0.;
  liquid_velocity_.setZero();
  set_liquid_traction_ = false;
  permeability_.setZero();
  liquid_traction_.setZero();
  liquid_saturation_ = 1.;

  // Initialize vector data properties
  this->vector_properties_["liquid_velocities"] = [&]() {
    return liquid_velocity();
  };
}

// Assign degree of saturation to the liquid phase
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::assign_saturation_degree() {
  bool status = true;
  try {
    if (this->material(mpm::ParticlePhase::Liquid) != nullptr) {
      liquid_saturation_ =
          this->material(mpm::ParticlePhase::Liquid)
              ->template property<double>(std::string("saturation"));
      if (liquid_saturation_ < 0. || liquid_saturation_ > 1.)
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

// Compute mass of particle (both solid and fluid)
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::compute_mass() noexcept {
  // Check if particle volume is set and material ptr is valid
  assert(volume_ != std::numeric_limits<double>::max() &&
         this->material(mpm::ParticlePhase::Solid) != nullptr &&
         this->material(mpm::ParticlePhase::Liquid) != nullptr);
  // Mass = volume of particle * mass_density
  // Solid mass
  this->mass_density_ =
      (this->material(mpm::ParticlePhase::Solid))
          ->template property<double>(std::string("density")) *
      (1 - this->porosity_);
  this->mass_ = volume_ * mass_density_;

  // Liquid mass
  this->liquid_mass_density_ =
      liquid_saturation_ * porosity_ *
      (this->material(mpm::ParticlePhase::Liquid))
          ->template property<double>(std::string("density"));
  this->liquid_mass_ = volume_ * liquid_mass_density_;
}

//! Map particle mass and momentum to nodes
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_mass_momentum_to_nodes() noexcept {
  mpm::Particle<Tdim>::map_mass_momentum_to_nodes();
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

//! Compute pore pressure
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::compute_pore_pressure(double dt) noexcept {
  // Check if liquid material and cell pointer are set and positive
  assert(this->material(mpm::ParticlePhase::Liquid) != nullptr &&
         cell_ != nullptr);

  // get the bulk modulus of liquid
  double K = this->material(mpm::ParticlePhase::Liquid)
                 ->template property<double>(std::string("bulk_modulus"));

  // Compute at centroid
  // get liquid phase strain rate at cell centre
  auto liquid_strain_rate_centroid =
      this->compute_strain_rate(dn_dx_centroid_, mpm::ParticlePhase::Liquid);

  // update pressure
  this->state_variables_[mpm::ParticlePhase::Liquid].at("pressure") +=
      -dt * (K / porosity_) *
      ((1 - porosity_) * strain_rate_.head(Tdim).sum() +
       porosity_ * liquid_strain_rate_centroid.head(Tdim).sum());

  // Apply free surface
  if (this->free_surface())
    this->assign_state_variable("pressure", 0., mpm::ParticlePhase::Liquid);
}

//! Map body force for both mixture and liquid
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_body_force(
    const VectorDim& pgravity) noexcept {
  this->map_mixture_body_force(mpm::ParticlePhase::Mixture, pgravity);
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

//! Map traction force
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_traction_force() noexcept {
  if (this->set_traction_) this->map_mixture_traction_force();
  if (this->set_liquid_traction_) this->map_liquid_traction_force();
}

//! Map mixture traction force
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_mixture_traction_force() noexcept {
  // Map particle traction forces to nodes
  for (unsigned i = 0; i < nodes_.size(); ++i)
    nodes_[i]->update_external_force(true, mpm::ParticlePhase::Mixture,
                                     (shapefn_[i] * traction_));
}

//! Map liquid traction force
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_liquid_traction_force() noexcept {
  // Map particle liquid traction forces to nodes
  for (unsigned i = 0; i < nodes_.size(); ++i)
    nodes_[i]->update_external_force(true, mpm::ParticlePhase::Liquid,
                                     (shapefn_[i] * liquid_traction_));
}

//! Map both mixture and liquid internal force
template <unsigned Tdim>
inline void mpm::TwoPhaseParticle<Tdim>::map_internal_force() noexcept {
  mpm::TwoPhaseParticle<Tdim>::map_mixture_internal_force(
      mpm::ParticlePhase::Mixture);
  mpm::TwoPhaseParticle<Tdim>::map_liquid_internal_force();
}

//! Map liquid phase internal force
template <>
inline void mpm::TwoPhaseParticle<1>::map_liquid_internal_force() noexcept {
  // initialise a vector of pore pressure
  Eigen::Matrix<double, 6, 1> pressure;
  pressure.setZero();
  pressure(0) = -this->state_variable("pressure", mpm::ParticlePhase::Liquid);

  // Compute nodal internal forces
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    // Compute force: -pstress * volume
    Eigen::Matrix<double, 1, 1> force;
    force[0] = dn_dx_(i, 0) * pressure[0] * this->porosity_;

    force *= -1. * this->volume_;

    nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Liquid, force);
  }
}

//! Map liquid phase internal force
template <>
inline void mpm::TwoPhaseParticle<2>::map_liquid_internal_force() noexcept {
  // initialise a vector of pore pressure
  Eigen::Matrix<double, 6, 1> pressure;
  pressure.setZero();
  pressure(0) = -this->state_variable("pressure", mpm::ParticlePhase::Liquid);
  pressure(1) = -this->state_variable("pressure", mpm::ParticlePhase::Liquid);

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
inline void mpm::TwoPhaseParticle<3>::map_liquid_internal_force() noexcept {
  // initialise a vector of pore pressure
  Eigen::Matrix<double, 6, 1> pressure;
  pressure.setZero();
  pressure(0) = -this->state_variable("pressure", mpm::ParticlePhase::Liquid);
  pressure(1) = -this->state_variable("pressure", mpm::ParticlePhase::Liquid);
  pressure(2) = -this->state_variable("pressure", mpm::ParticlePhase::Liquid);

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
inline void mpm::TwoPhaseParticle<1>::map_mixture_internal_force(
    unsigned mixture) noexcept {
  // initialise a vector of pore pressure
  Eigen::Matrix<double, 6, 1> total_stress = this->stress_;
  total_stress(0) -=
      this->state_variable("pressure", mpm::ParticlePhase::Liquid);

  // Compute nodal internal forces
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    // Compute force: -pstress * volume
    Eigen::Matrix<double, 1, 1> force;
    force[0] = dn_dx_(i, 0) * total_stress[0] + dn_dx_(i, 1) * total_stress[3];

    force *= -1. * this->volume_;

    nodes_[i]->update_internal_force(true, mixture, force);
  }
}

//! Map mixture internal force
template <>
inline void mpm::TwoPhaseParticle<2>::map_mixture_internal_force(
    unsigned mixture) noexcept {
  // initialise a vector of pore pressure
  Eigen::Matrix<double, 6, 1> total_stress = this->stress_;
  total_stress(0) -=
      this->state_variable("pressure", mpm::ParticlePhase::Liquid);
  total_stress(1) -=
      this->state_variable("pressure", mpm::ParticlePhase::Liquid);

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
inline void mpm::TwoPhaseParticle<3>::map_mixture_internal_force(
    unsigned mixture) noexcept {
  // initialise a vector of pore pressure
  Eigen::Matrix<double, 6, 1> total_stress = this->stress_;
  total_stress(0) -=
      this->state_variable("pressure", mpm::ParticlePhase::Liquid);
  total_stress(1) -=
      this->state_variable("pressure", mpm::ParticlePhase::Liquid);
  total_stress(2) -=
      this->state_variable("pressure", mpm::ParticlePhase::Liquid);

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

// Compute updated position of the particle and kinematics of both solid and
// liquid phase
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::compute_updated_position(
    double dt, bool velocity_update) noexcept {
  mpm::Particle<Tdim>::compute_updated_position(dt, velocity_update);
  this->compute_updated_liquid_velocity(dt, velocity_update);
}

//! Map particle pressure to nodes
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::map_pressure_to_nodes(
    unsigned phase) noexcept {
  // Mass is initialized
  assert(liquid_mass_ != std::numeric_limits<double>::max());

  bool status = false;
  // If phase is Solid, use the default map_pressure_to_nodes
  if (phase == mpm::ParticlePhase::Solid)
    status = mpm::Particle<Tdim>::map_pressure_to_nodes(phase);
  else {
    // Check if particle liquid mass is set and state variable pressure is found
    if (liquid_mass_ != std::numeric_limits<double>::max() &&
        (state_variables_[phase].find("pressure") !=
         state_variables_[phase].end())) {
      // Map particle pressure to nodes
      for (unsigned i = 0; i < nodes_.size(); ++i)
        nodes_[i]->update_mass_pressure(
            phase,
            shapefn_[i] * liquid_mass_ * state_variables_[phase]["pressure"]);

      status = true;
    }
  }
  return status;
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
}

//! Apply particle velocity constraints
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::apply_particle_velocity_constraints(
    unsigned dir, double velocity) {
  // Set particle velocity constraint for solid phase
  if (dir < Tdim)
    mpm::Particle<Tdim>::apply_particle_velocity_constraints(dir, velocity);

  // Set particle velocity constraint for liquid phase
  else
    this->liquid_velocity_(static_cast<unsigned>(dir % Tdim)) = velocity;
}

// Assign porosity to the particle
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::assign_porosity() {
  bool status = true;
  try {
    if (this->material(mpm::ParticlePhase::Solid) != nullptr) {
      this->porosity_ =
          this->material(mpm::ParticlePhase::Solid)
              ->template property<double>(std::string("porosity"));
      if (porosity_ < 0. || porosity_ > 1.)
        throw std::runtime_error(
            "Particle porosity is negative or larger than one");
    } else {
      throw std::runtime_error("Material is invalid");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assign particle permeability
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::assign_permeability() {
  bool status = true;
  try {
    // Check if material ptr is valid
    if (this->material(mpm::ParticlePhase::Solid) != nullptr) {
      // Get initial porosity
      double porosity =
          this->material(mpm::ParticlePhase::Solid)
              ->template property<double>(std::string("porosity"));
      if (porosity < 0. || porosity > 1.)
        throw std::runtime_error(
            "Particle porosity is negative or larger than one, can not assign "
            "permeability");
      // Porosity parameter
      const double k_p = std::pow(porosity, 3) / std::pow((1. - porosity), 2);
      // Assign permeability
      switch (Tdim) {
        case (3):
          // Check if the permeability is valid
          if (this->material(mpm::ParticlePhase::Solid)
                  ->template property<double>("k_z") < 0)
            throw std::runtime_error("Material's permeability is invalid");
          permeability_(2) = this->material(mpm::ParticlePhase::Solid)
                                 ->template property<double>("k_z") /
                             k_p;
        case (2):
          // Check if the permeability is valid
          if (this->material(mpm::ParticlePhase::Solid)
                  ->template property<double>("k_y") < 0)
            throw std::runtime_error("Material's permeability is invalid");
          permeability_(1) = this->material(mpm::ParticlePhase::Solid)
                                 ->template property<double>("k_y") /
                             k_p;
        case (1):
          // Check if the permeability is valid
          if (this->material(mpm::ParticlePhase::Solid)
                  ->template property<double>("k_x") < 0)
            throw std::runtime_error("Material's permeability is invalid");
          // Assign permeability
          permeability_(0) = this->material(mpm::ParticlePhase::Solid)
                                 ->template property<double>("k_x") /
                             k_p;
      }
    } else {
      throw std::runtime_error("Material is invalid");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Map drag force
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::map_drag_force_coefficient() {
  bool status = true;
  try {
    // Porosity parameter
    const double k_p =
        std::pow(this->porosity_, 3) / std::pow((1. - this->porosity_), 2);
    // Initialise drag force coefficient
    VectorDim drag_force_coefficient;
    drag_force_coefficient.setZero();

    // Check if permeability coefficient is valid
    for (unsigned i = 0; i < Tdim; ++i) {
      if (k_p > 0.)
        drag_force_coefficient(i) =
            porosity_ * porosity_ * 9.81 *
            this->material(mpm::ParticlePhase::Liquid)
                ->template property<double>(std::string("density")) /
            (k_p * permeability_(i));
      else
        throw std::runtime_error("Porosity coefficient is invalid");
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

//! Initial pore pressure
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::initialise_pore_pressure_watertable(
    const unsigned dir_v, const unsigned dir_h, VectorDim& gravity,
    std::map<double, double>& reference_points) {
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
    // Iterate over each reference_points
    for (const auto& reference_point : reference_points) {
      // Find boundary
      if (reference_point.first > left_boundary &&
          reference_point.first <= position) {
        // Left boundary position and h0
        left_boundary = reference_point.first;
        h0_left = reference_point.second;
      } else if (reference_point.first > position &&
                 reference_point.first <= right_boundary) {
        // Right boundary position and h0
        right_boundary = reference_point.first;
        h0_right = reference_point.second;
      }
    }

    // Initialise pore pressure
    double pore_pressure = 0;

    if (left_boundary != std::numeric_limits<double>::lowest()) {
      // Particle with left and right boundary
      if (right_boundary != std::numeric_limits<double>::max()) {
        pore_pressure =
            ((h0_right - h0_left) / (right_boundary - left_boundary) *
                 (position - left_boundary) +
             h0_left - this->coordinates_(dir_v)) *
            (this->material(mpm::ParticlePhase::Liquid))
                ->template property<double>(std::string("density")) *
            (-gravity(dir_v));
      } else
        // Particle with only left boundary
        pore_pressure =
            (h0_left - this->coordinates_(dir_v)) *
            (this->material(mpm::ParticlePhase::Liquid))
                ->template property<double>(std::string("density")) *
            (-gravity(dir_v));

    }
    // Particle with only right boundary
    else if (right_boundary != std::numeric_limits<double>::max())
      pore_pressure = (h0_right - this->coordinates_(dir_v)) *
                      (this->material(mpm::ParticlePhase::Liquid))
                          ->template property<double>(std::string("density")) *
                      (-gravity(dir_v));
    else
      throw std::runtime_error(
          "Particle pore pressure can not be initialised by water table");

    // Check negative pore pressure
    if (pore_pressure < 0) pore_pressure = 0;

    // Assign pore pressure
    this->assign_state_variable("pressure", pore_pressure,
                                mpm::ParticlePhase::Liquid);
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Update material point porosity
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::update_porosity(double dt) {
  // Update particle porosity
  const double porosity =
      1 - (1 - this->porosity_) / (1 + dt * strain_rate_.head(Tdim).sum());
  // Check if the value is valid
  if (porosity < 0.)
    this->porosity_ = 1.E-5;
  else if (porosity > 1.)
    this->porosity_ = 1 - 1.E-5;
  else
    this->porosity_ = porosity;
}

// Assign traction to the particle
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::assign_traction(unsigned direction,
                                                  double traction) {
  bool status = false;
  try {
    if (direction >= Tdim * 2 ||
        this->volume_ == std::numeric_limits<double>::max()) {
      throw std::runtime_error(
          "Particle traction property: volume / direction is invalid");
    }
    // Assign mixture traction
    if (direction < Tdim) {
      this->set_traction_ = true;
      traction_(direction) = traction * this->volume_ / this->size_(direction);
    }
    // Assign liquid traction
    else {
      this->set_liquid_traction_ = true;
      liquid_traction_(direction - Tdim) =
          traction * this->volume_ / this->size_(direction - Tdim);
    }
    status = true;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Compute size of serialized particle data
template <unsigned Tdim>
int mpm::TwoPhaseParticle<Tdim>::compute_pack_size() const {
  int total_size = mpm::Particle<Tdim>::compute_pack_size();
  int partial_size;
#ifdef USE_MPI
  // material id for liquid phase
  MPI_Pack_size(1, MPI_UNSIGNED, MPI_COMM_WORLD, &partial_size);
  total_size += partial_size;

  // liquid mass
  MPI_Pack_size(1, MPI_DOUBLE, MPI_COMM_WORLD, &partial_size);
  total_size += partial_size;

  // liquid velocity
  MPI_Pack_size(Tdim, MPI_DOUBLE, MPI_COMM_WORLD, &partial_size);
  total_size += partial_size;

  // porosity, liquid saturation
  MPI_Pack_size(2 * 1, MPI_DOUBLE, MPI_COMM_WORLD, &partial_size);
  total_size += partial_size;

  // nliquid state variables
  unsigned nliquid_state_vars =
      state_variables_[mpm::ParticlePhase::Liquid].size();
  MPI_Pack_size(1, MPI_UNSIGNED, MPI_COMM_WORLD, &partial_size);
  total_size += partial_size;

  // liquid state variables
  MPI_Pack_size(nliquid_state_vars, MPI_DOUBLE, MPI_COMM_WORLD, &partial_size);
  total_size += partial_size;
#endif
  return total_size;
}

//! Serialize particle data
template <unsigned Tdim>
std::vector<uint8_t> mpm::TwoPhaseParticle<Tdim>::serialize() {
  // Compute pack size
  if (pack_size_ == 0) pack_size_ = this->compute_pack_size();
  // Initialize data buffer
  std::vector<uint8_t> data;
  data.resize(pack_size_);
  uint8_t* data_ptr = &data[0];
  int position = 0;

#ifdef USE_MPI
  // Type
  int type = ParticleType.at(this->type());
  MPI_Pack(&type, 1, MPI_INT, data_ptr, data.size(), &position, MPI_COMM_WORLD);

  // Material ID
  unsigned nmaterials = material_id_.size();
  MPI_Pack(&nmaterials, 1, MPI_UNSIGNED, data_ptr, data.size(), &position,
           MPI_COMM_WORLD);
  MPI_Pack(&material_id_[mpm::ParticlePhase::Solid], 1, MPI_UNSIGNED, data_ptr,
           data.size(), &position, MPI_COMM_WORLD);
  MPI_Pack(&material_id_[mpm::ParticlePhase::Liquid], 1, MPI_UNSIGNED, data_ptr,
           data.size(), &position, MPI_COMM_WORLD);

  // ID
  MPI_Pack(&this->id_, 1, MPI_UNSIGNED_LONG_LONG, data_ptr, data.size(),
           &position, MPI_COMM_WORLD);

  // Solid Phase
  // Mass
  MPI_Pack(&mass_, 1, MPI_DOUBLE, data_ptr, data.size(), &position,
           MPI_COMM_WORLD);
  // Volume
  MPI_Pack(&volume_, 1, MPI_DOUBLE, data_ptr, data.size(), &position,
           MPI_COMM_WORLD);
  // Pressure
  double pressure =
      (state_variables_[mpm::ParticlePhase::Solid].find("pressure") !=
       state_variables_[mpm::ParticlePhase::Solid].end())
          ? state_variables_[mpm::ParticlePhase::Solid].at("pressure")
          : 0.;
  MPI_Pack(&pressure, 1, MPI_DOUBLE, data_ptr, data.size(), &position,
           MPI_COMM_WORLD);

  // Coordinates
  MPI_Pack(coordinates_.data(), Tdim, MPI_DOUBLE, data_ptr, data.size(),
           &position, MPI_COMM_WORLD);
  // Displacement
  MPI_Pack(this->displacement_.data(), Tdim, MPI_DOUBLE, data_ptr, data.size(),
           &position, MPI_COMM_WORLD);
  // Natural size
  MPI_Pack(this->natural_size_.data(), Tdim, MPI_DOUBLE, data_ptr, data.size(),
           &position, MPI_COMM_WORLD);
  // Velocity
  MPI_Pack(this->velocity_.data(), Tdim, MPI_DOUBLE, data_ptr, data.size(),
           &position, MPI_COMM_WORLD);
  // Stress
  MPI_Pack(stress_.data(), 6, MPI_DOUBLE, data_ptr, data.size(), &position,
           MPI_COMM_WORLD);
  // Strain
  MPI_Pack(this->strain_.data(), 6, MPI_DOUBLE, data_ptr, data.size(),
           &position, MPI_COMM_WORLD);

  // epsv
  MPI_Pack(&this->volumetric_strain_centroid_, 1, MPI_DOUBLE, data_ptr,
           data.size(), &position, MPI_COMM_WORLD);

  // Cell id
  MPI_Pack(&this->cell_id_, 1, MPI_UNSIGNED_LONG_LONG, data_ptr, data.size(),
           &position, MPI_COMM_WORLD);

  // Status
  MPI_Pack(&this->status_, 1, MPI_C_BOOL, data_ptr, data.size(), &position,
           MPI_COMM_WORLD);

  // nstate variables
  unsigned nstate_vars = state_variables_[mpm::ParticlePhase::Solid].size();
  MPI_Pack(&nstate_vars, 1, MPI_UNSIGNED, data_ptr, data.size(), &position,
           MPI_COMM_WORLD);

  // state variables
  if (this->material(mpm::ParticlePhase::Solid) != nullptr) {
    std::vector<double> svars;
    auto state_variables =
        (this->material(mpm::ParticlePhase::Solid))->state_variables();
    for (const auto& state_var : state_variables)
      svars.emplace_back(
          state_variables_[mpm::ParticlePhase::Solid].at(state_var));

    // Write state vars
    MPI_Pack(&svars[0], nstate_vars, MPI_DOUBLE, data_ptr, data.size(),
             &position, MPI_COMM_WORLD);
  }

  // Liquid Phase
  // Mass
  MPI_Pack(&liquid_mass_, 1, MPI_DOUBLE, data_ptr, data.size(), &position,
           MPI_COMM_WORLD);
  // Velocity
  MPI_Pack(liquid_velocity_.data(), Tdim, MPI_DOUBLE, data_ptr, data.size(),
           &position, MPI_COMM_WORLD);
  // Porosity
  MPI_Pack(&porosity_, 1, MPI_DOUBLE, data_ptr, data.size(), &position,
           MPI_COMM_WORLD);
  // Liquid Saturation
  MPI_Pack(&liquid_saturation_, 1, MPI_DOUBLE, data_ptr, data.size(), &position,
           MPI_COMM_WORLD);

  // nstate variables
  unsigned nliquid_state_vars =
      state_variables_[mpm::ParticlePhase::Liquid].size();
  MPI_Pack(&nliquid_state_vars, 1, MPI_UNSIGNED, data_ptr, data.size(),
           &position, MPI_COMM_WORLD);

  // state variables
  if (this->material(mpm::ParticlePhase::Liquid) != nullptr) {
    std::vector<double> svars;
    auto state_variables =
        (this->material(mpm::ParticlePhase::Liquid))->state_variables();
    for (const auto& state_var : state_variables)
      svars.emplace_back(
          state_variables_[mpm::ParticlePhase::Liquid].at(state_var));

    // Write state vars
    MPI_Pack(&svars[0], nliquid_state_vars, MPI_DOUBLE, data_ptr, data.size(),
             &position, MPI_COMM_WORLD);
  }
#endif
  return data;
}

//! Deserialize particle data
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::deserialize(
    const std::vector<uint8_t>& data,
    std::vector<std::shared_ptr<mpm::Material<Tdim>>>& materials) {
  uint8_t* data_ptr = const_cast<uint8_t*>(&data[0]);
  int position = 0;

#ifdef USE_MPI
  // Type
  int type;
  MPI_Unpack(data_ptr, data.size(), &position, &type, 1, MPI_INT,
             MPI_COMM_WORLD);
  if (type != ParticleType.at(this->type()))
    throw std::runtime_error("Deserialize particle(): particle type mismatch");

  // nmaterials
  int nmaterials = 0;
  MPI_Unpack(data_ptr, data.size(), &position, &nmaterials, 1, MPI_UNSIGNED,
             MPI_COMM_WORLD);
  if (nmaterials != materials.size())
    throw std::runtime_error(
        "Deserialize particle(): nmaterials mismatch with the input materials "
        "size");

  // Material ID
  MPI_Unpack(data_ptr, data.size(), &position,
             &material_id_[mpm::ParticlePhase::Solid], 1, MPI_UNSIGNED,
             MPI_COMM_WORLD);
  MPI_Unpack(data_ptr, data.size(), &position,
             &material_id_[mpm::ParticlePhase::Liquid], 1, MPI_UNSIGNED,
             MPI_COMM_WORLD);

  // ID
  MPI_Unpack(data_ptr, data.size(), &position, &this->id_, 1,
             MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);

  // Solid Phase
  // mass
  MPI_Unpack(data_ptr, data.size(), &position, &mass_, 1, MPI_DOUBLE,
             MPI_COMM_WORLD);
  // volume
  MPI_Unpack(data_ptr, data.size(), &position, &volume_, 1, MPI_DOUBLE,
             MPI_COMM_WORLD);
  // pressure
  double pressure;
  MPI_Unpack(data_ptr, data.size(), &position, &pressure, 1, MPI_DOUBLE,
             MPI_COMM_WORLD);
  this->assign_pressure(pressure);

  // Coordinates
  MPI_Unpack(data_ptr, data.size(), &position, coordinates_.data(), Tdim,
             MPI_DOUBLE, MPI_COMM_WORLD);
  // Displacement
  MPI_Unpack(data_ptr, data.size(), &position, this->displacement_.data(), Tdim,
             MPI_DOUBLE, MPI_COMM_WORLD);
  // Natural size
  MPI_Unpack(data_ptr, data.size(), &position, this->natural_size_.data(), Tdim,
             MPI_DOUBLE, MPI_COMM_WORLD);
  // Velocity
  MPI_Unpack(data_ptr, data.size(), &position, this->velocity_.data(), Tdim,
             MPI_DOUBLE, MPI_COMM_WORLD);
  // Stress
  MPI_Unpack(data_ptr, data.size(), &position, stress_.data(), 6, MPI_DOUBLE,
             MPI_COMM_WORLD);
  // Strain
  MPI_Unpack(data_ptr, data.size(), &position, this->strain_.data(), 6,
             MPI_DOUBLE, MPI_COMM_WORLD);

  // epsv
  MPI_Unpack(data_ptr, data.size(), &position,
             &this->volumetric_strain_centroid_, 1, MPI_DOUBLE, MPI_COMM_WORLD);
  // cell id
  MPI_Unpack(data_ptr, data.size(), &position, &this->cell_id_, 1,
             MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);
  // status
  MPI_Unpack(data_ptr, data.size(), &position, &this->status_, 1, MPI_C_BOOL,
             MPI_COMM_WORLD);

  // Assign materials
  if (material_id_[mpm::ParticlePhase::Solid] ==
      materials.at(mpm::ParticlePhase::Solid)->id()) {
    bool assign_mat =
        this->assign_material(materials.at(mpm::ParticlePhase::Solid));
    if (!assign_mat) throw std::runtime_error("Material assignment failed");
  }

  // nstate vars
  unsigned nstate_vars;
  MPI_Unpack(data_ptr, data.size(), &position, &nstate_vars, 1, MPI_UNSIGNED,
             MPI_COMM_WORLD);

  if (nstate_vars > 0) {
    std::vector<double> svars;
    svars.reserve(nstate_vars);
    MPI_Unpack(data_ptr, data.size(), &position, &svars, nstate_vars,
               MPI_DOUBLE, MPI_COMM_WORLD);

    // Reinitialize state variables
    auto mat_state_vars = (this->material(mpm::ParticlePhase::Solid))
                              ->initialise_state_variables();
    if (mat_state_vars.size() == nstate_vars) {
      unsigned i = 0;
      auto state_variables =
          (this->material(mpm::ParticlePhase::Solid))->state_variables();
      for (const auto& state_var : state_variables) {
        this->state_variables_[mpm::ParticlePhase::Solid].at(state_var) =
            svars[i];
        ++i;
      }
    } else
      throw std::runtime_error(
          "Deserialize particle(): state_vars size mismatch");
  }

  // Liquid Phase
  // liquid mass
  MPI_Unpack(data_ptr, data.size(), &position, &liquid_mass_, 1, MPI_DOUBLE,
             MPI_COMM_WORLD);
  // Liquid velocity
  MPI_Unpack(data_ptr, data.size(), &position, liquid_velocity_.data(), Tdim,
             MPI_DOUBLE, MPI_COMM_WORLD);
  // porosity
  MPI_Unpack(data_ptr, data.size(), &position, &porosity_, 1, MPI_DOUBLE,
             MPI_COMM_WORLD);
  // liquid Saturation
  MPI_Unpack(data_ptr, data.size(), &position, &liquid_saturation_, 1,
             MPI_DOUBLE, MPI_COMM_WORLD);

  // Assign permeability
  this->assign_permeability();

  // Assign liquid materials
  if (material_id_[mpm::ParticlePhase::Liquid] ==
      materials.at(mpm::ParticlePhase::Liquid)->id()) {
    bool assign_mat =
        this->assign_material(materials.at(mpm::ParticlePhase::Liquid));
    if (!assign_mat) throw std::runtime_error("Material assignment failed");
  }

  // nliquid state vars
  unsigned nliquid_state_vars;
  MPI_Unpack(data_ptr, data.size(), &position, &nliquid_state_vars, 1,
             MPI_UNSIGNED, MPI_COMM_WORLD);

  if (nliquid_state_vars > 0) {
    std::vector<double> svars;
    svars.reserve(nliquid_state_vars);
    MPI_Unpack(data_ptr, data.size(), &position, &svars, nliquid_state_vars,
               MPI_DOUBLE, MPI_COMM_WORLD);

    // Reinitialize state variables
    auto mat_state_vars = (this->material(mpm::ParticlePhase::Liquid))
                              ->initialise_state_variables();
    if (mat_state_vars.size() == nliquid_state_vars) {
      unsigned i = 0;
      auto state_variables =
          (this->material(mpm::ParticlePhase::Liquid))->state_variables();
      for (const auto& state_var : state_variables) {
        this->state_variables_[mpm::ParticlePhase::Liquid].at(state_var) =
            svars[i];
        ++i;
      }
    } else
      throw std::runtime_error(
          "Deserialize particle(): liquid_state_vars size mismatch");
  }
#endif
}
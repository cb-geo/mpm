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

//! Return particle data in HDF5 format
template <unsigned Tdim>
// cppcheck-suppress *
mpm::HDF5Particle mpm::TwoPhaseParticle<Tdim>::hdf5() const {
  // Derive from particle
  auto particle_data = mpm::Particle<Tdim>::hdf5();
  // // Particle liquid mass
  // particle_data.liquid_mass = this->liquid_mass_;

  // // Particle liquid velocity
  // particle_data.liquid_velocity_x = liquid_velocity[0];
  // particle_data.liquid_velocity_y = liquid_velocity[1];
  // particle_data.liquid_velocity_z = liquid_velocity[2];

  // // Particle liquid material id
  // particle_data.liquid_material_id =
  // this->material_id_(mpm::ParticlePhase::Liquid);

  if (this->material(mpm::ParticlePhase::Liquid) != nullptr) {
    if ((state_variables_[mpm::ParticlePhase::Solid].size() +
         state_variables_[mpm::ParticlePhase::Liquid].size()) > 20)
      throw std::runtime_error("# of state variables cannot be more than 20");
    // Assign number of state variables
    particle_data.nstate_vars =
        state_variables_[mpm::ParticlePhase::Solid].size() +
        state_variables_[mpm::ParticlePhase::Liquid].size();
    // First id
    unsigned i = state_variables_[mpm::ParticlePhase::Solid].size();
    // Liquid state variables
    auto state_variables =
        (this->material(mpm::ParticlePhase::Liquid))->state_variables();
    for (const auto& state_var : state_variables) {
      particle_data.svars[i] =
          state_variables_[mpm::ParticlePhase::Liquid].at(state_var);
      ++i;
    }
  }

  return particle_data;
}

//! Initialise particle data from HDF5
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::initialise_particle(
    const HDF5Particle& particle) {
  // Derive from particle
  mpm::Particle<Tdim>::initialise_particle(particle);

  // // Liquid mass
  // this->liquid_mass_ = particle.liquid_mass;
  // // Liquid mass Density
  // this->liquid_mass_density_ = particle.liquid_mass / particle.volume;

  // // Liquid velocity
  // Eigen::Vector3d liquid_velocity;
  // liquid_velocity << particle.liquid_velocity_x, particle.liquid_velocity_y,
  //     particle.liquid_velocity_z;
  // // Initialise velocity
  // for (unsigned i = 0; i < Tdim; ++i)
  //   this->liquid_velocity_(i) = liquid_velocity(i);

  // // Liquid material id
  // this->liquid_material_id_ = particle.liquid_material_id;

  return true;
}

//! Initialise particle data from HDF5
//! TODO
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::initialise_particle(
    const HDF5Particle& particle,
    const std::shared_ptr<mpm::Material<Tdim>>& material) {
  bool status = mpm::Particle<Tdim>::initialise_particle(particle, material);
  return status;
}

// Initialise liquid phase particle properties
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::initialise() {
  mpm::Particle<Tdim>::initialise();
  liquid_mass_ = 0.;
  liquid_velocity_.setZero();
  set_liquid_traction_ = false;
  liquid_traction_.setZero();
  liquid_saturation_ = 1.;

  // Initialize vector data properties
  this->properties_["liquid_velocities"] = [&]() { return liquid_velocity(); };
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
         this->material() != nullptr &&
         this->material(mpm::ParticlePhase::Liquid) != nullptr);
  // Mass = volume of particle * mass_density
  // Solid mass
  this->mass_density_ =
      (this->material())->template property<double>(std::string("density")) *
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
    if (this->material() != nullptr) {
      this->porosity_ =
          this->material()->template property<double>(std::string("porosity"));
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
    if (this->material() != nullptr) {
      // Porosity parameter
      const double k_p =
          std::pow(this->porosity_, 3) / std::pow((1. - this->porosity_), 2);
      switch (Tdim) {
        case (1): {
          // Check if the permeability is valid
          if (this->material()->template property<double>("k_x") < 0)
            throw std::runtime_error("Material's permeability is invalid");
          // Assign permeability
          permeability_(0) =
              this->material()->template property<double>("k_x") / k_p;
          break;
        }
        case (2): {
          // Check if the permeability is valid
          if (this->material()->template property<double>("k_x") < 0 ||
              this->material()->template property<double>("k_y") < 0)
            throw std::runtime_error("Material's permeability is invalid");
          // Assign permeability
          permeability_(0) =
              this->material()->template property<double>("k_x") / k_p;
          permeability_(1) =
              this->material()->template property<double>("k_y") / k_p;
          break;
        }
        default: {
          // Check if the permeability is valid
          if (this->material()->template property<double>("k_x") < 0 ||
              this->material()->template property<double>("k_y") < 0 ||
              this->material()->template property<double>("k_z") < 0)
            throw std::runtime_error("Material's permeability is invalid");
          // Assign permeability
          permeability_(0) =
              this->material()->template property<double>("k_x") / k_p;
          permeability_(1) =
              this->material()->template property<double>("k_y") / k_p;
          permeability_(2) =
              this->material()->template property<double>("k_z") / k_p;
          break;
        }
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
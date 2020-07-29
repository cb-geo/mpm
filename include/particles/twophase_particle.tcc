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
    : mpm::ParticleBase<Tdim>(id, coord, status) {
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
  // // Particle liquid velocity
  // Eigen::Vector3d liquid_velocity;
  // liquid_velocity.setZero();
  // for (unsigned j = 0; j < Tdim; ++j)
  //   liquid_velocity[j] = this->liquid_velocity_[j];
  // // Particle liquid strain
  // Eigen::Matrix<double, 6, 1> liquid_strain = this->liquid_strain_;
  // // Particle liquid mass
  // particle_data.liquid_mass = this->liquid_mass_;
  // Particle pore pressure
  particle_data.pressure = this->pore_pressure_;
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

//! Initialise particle data from HDF5
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::initialise_particle(
    const HDF5Particle& particle) {
  // Derive from particle
  mpm::Particle<Tdim>::initialise_particle(particle);

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
  liquid_strain_rate_.setZero();
  liquid_strain_.setZero();
  set_liquid_traction_ = false;
  set_mixture_traction_ = false;
  liquid_traction_.setZero();
  pore_pressure_ = 0.;
  liquid_saturation_ = 1.;

  // Initialize vector data properties
  this->properties_["liquid_strains"] = [&]() { return liquid_strain(); };
  this->properties_["liquid_velocities"] = [&]() { return liquid_velocity(); };
  this->properties_["pore_pressure"] = [&]() {
    Eigen::VectorXd vec_pressure = Eigen::VectorXd::Zero(3);
    vec_pressure[0] = this->pore_pressure();
    // FIXME: This is to check free surface particles
    // TODO: To be removed somewhere
    // vec_pressure[1] = this->free_surface();
    return vec_pressure;
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

// Assign traction
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::assign_traction(unsigned direction,
                                                  double traction) {
  bool status = true;
  this->assign_mixture_traction(direction, traction);
  return status;
}

// Assign traction to the mixture
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::assign_mixture_traction(unsigned direction,
                                                          double traction) {
  bool status = false;
  try {
    if (direction >= Tdim ||
        this->volume_ == std::numeric_limits<double>::max()) {
      throw std::runtime_error(
          "Particle mixture traction property: volume / direction is invalid");
    }
    // Assign mixture traction
    mixture_traction_(direction) =
        traction * this->volume_ / this->size_(direction);
    status = true;
    this->set_mixture_traction_ = true;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Compute mass of particle (both solid and fluid)
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::compute_mass() noexcept {
  mpm::Particle<Tdim>::compute_mass();
  this->compute_liquid_mass();
}

// Compute fluid mass of particle
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::compute_liquid_mass() noexcept {
  // Check if particle volume is set and liquid material ptr is valid
  assert(volume_ != std::numeric_limits<double>::max() &&
         this->material(mpm::ParticlePhase::Liquid) != nullptr);

  // Mass = volume of particle * bulk_density
  this->liquid_mass_density_ =
      liquid_saturation_ * porosity_ *
      this->material(mpm::ParticlePhase::Liquid)
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
  this->pore_pressure_ +=
      -dt * (K / porosity_) *
      ((1 - porosity_) * strain_rate_.head(Tdim).sum() +
       porosity_ * liquid_strain_rate_centroid.head(Tdim).sum());

  // Apply free surface
  // if (this->free_surface()) this->pore_pressure_ = 0.0;
}

// Compute pore liquid pressure smoothing based on nodal pressure
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::compute_pore_pressure_smoothing() noexcept {
  // Check if particle has a valid cell ptr
  assert(cell_ != nullptr);

  bool status = true;
  double pressure = 0;
  for (unsigned i = 0; i < nodes_.size(); ++i)
    pressure += shapefn_(i) * nodes_[i]->pressure(mpm::ParticlePhase::Liquid);

  // Update pore liquid pressure to interpolated nodal pressure
  this->pore_pressure_ = pressure;
  // Apply free surface
  // if (this->free_surface()) this->pore_pressure_ = 0.0;
  return status;
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
  this->map_mixture_traction_force(mpm::ParticlePhase::Mixture);
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

//! Map mixture traction force
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_mixture_traction_force(
    unsigned mixture) noexcept {
  if (this->set_mixture_traction_) {
    // Map particle mixture traction forces to nodes
    for (unsigned i = 0; i < nodes_.size(); ++i)
      nodes_[i]->update_external_force(true, mixture,
                                       (shapefn_[i] * this->mixture_traction_));
  }
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
  assert(mass_ != std::numeric_limits<double>::max() &&
         liquid_mass_ != std::numeric_limits<double>::max());

  bool status = false;
  if (phase == mpm::ParticlePhase::Solid) {
    // Check if particle mass is set and state variable pressure is found
    if (mass_ != std::numeric_limits<double>::max() &&
        (state_variables_[phase].find("pressure") !=
         state_variables_[phase].end())) {
      // Map particle pressure to nodes
      for (unsigned i = 0; i < nodes_.size(); ++i)
        nodes_[i]->update_mass_pressure(
            phase, shapefn_[i] * mass_ * state_variables_[phase]["pressure"]);

      status = true;
    }
  } else {
    // Map particle pore pressure to nodes
    for (unsigned i = 0; i < nodes_.size(); ++i)
      nodes_[i]->update_mass_pressure(
          phase, shapefn_[i] * liquid_mass_ * this->pore_pressure_);

    status = true;
  }
  return status;
}

// Compute pressure smoothing of the particle based on nodal pressure
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::compute_pressure_smoothing(
    unsigned phase) noexcept {
  // Assert
  assert(cell_ != nullptr);

  bool status = false;
  if (phase == mpm::ParticlePhase::Solid) {
    // Check if particle has a valid cell ptr
    if (cell_ != nullptr && (state_variables_[phase].find("pressure") !=
                             state_variables_[phase].end())) {

      double pressure = 0.;
      // Update particle pressure to interpolated nodal pressure
      for (unsigned i = 0; i < this->nodes_.size(); ++i)
        pressure += shapefn_[i] * nodes_[i]->pressure(phase);

      state_variables_[phase]["pressure"] = pressure;
      status = true;
    }
  } else {
    double pressure = 0.;
    // Update particle pressure to interpolated nodal pressure
    for (unsigned i = 0; i < this->nodes_.size(); ++i)
      pressure += shapefn_[i] * nodes_[i]->pressure(phase);

    this->pore_pressure_ = pressure;
    status = true;
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
          permeability_(0) = this->material()->template property<double>("k_x");
          c1_(0) = permeability_(0) / k_p;
          break;
        }
        case (2): {
          permeability_(0) = this->material()->template property<double>("k_x");
          permeability_(1) = this->material()->template property<double>("k_y");
          c1_(0) = permeability_(0) / k_p;
          c1_(1) = permeability_(1) / k_p;
          break;
        }
        default: {
          permeability_(0) = this->material()->template property<double>("k_x");
          permeability_(1) = this->material()->template property<double>("k_y");
          permeability_(2) = this->material()->template property<double>("k_z");
          c1_(0) = permeability_(0) / k_p;
          c1_(1) = permeability_(1) / k_p;
          c1_(2) = permeability_(2) / k_p;
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

//! Update particle permeability
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::update_permeability() {
  bool status = true;
  try {
    // Porosity parameter
    const double k_p =
        std::pow(this->porosity_, 3) / std::pow((1. - this->porosity_), 2);

    // Update permeability by KC equation
    permeability_ = k_p * c1_;

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
    // Update permeability
    this->update_permeability();
    // Initialise drag force coefficient
    VectorDim drag_force_coefficient;
    drag_force_coefficient.setZero();

    // Check if permeability coefficient is valid
    for (unsigned i = 0; i < Tdim; ++i) {
      if (permeability_(i) > 0.)
        drag_force_coefficient(i) =
            porosity_ * porosity_ * 9.81 *
            this->material(mpm::ParticlePhase::Liquid)
                ->template property<double>(std::string("density")) /
            permeability_(i);
      else
        throw std::runtime_error("Permeability coefficient is invalid");
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

// Update material point porosity
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::update_porosity(double dt) {
  bool status = true;
  try {
    // Update particle porosity
    const double porosity =
        1 - (1 - this->porosity_) / (1 + dt * strain_rate_.head(Tdim).sum());
    // Check if the value is valid
    if (porosity > 0 && porosity < 1) this->porosity_ = porosity;
    // Invalid value
    else
      throw std::runtime_error(
          "Invalid porosity, less than zero or greater than one");

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}
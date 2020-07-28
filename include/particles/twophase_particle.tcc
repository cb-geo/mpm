//! Construct a twophase particle with id and coordinates
template <unsigned Tdim>
mpm::TwoPhaseParticle<Tdim>::TwoPhaseParticle(Index id, const VectorDim& coord)
    : mpm::Particle<Tdim>(id, coord) {
  // Initialise variables for solid phase and liquid phase
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
  // Initialise variables for solid phase and liquid phase
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

//! Initialise particle data from HDF5
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::initialise_particle(
    const HDF5Particle& particle) {
  // Derive from particle
  mpm::Particle<Tdim>::initialise_particle(particle);

  // // Liquid mass
  // scalar_properties_.at(mpm::properties::Scalar::LiquidMass) =
  //     particle.liquid_mass;
  // // Liquid mass Density
  // scalar_properties_.at(mpm::properties::Scalar::LiquidMassDensity) =
  //     particle.mass / particle.volume;
  // // Pore pressure
  // scalar_properties_.at(mpm::properties::Scalar::PorePressure) =
  //     particle.pore_pressure;
  // // Liquid velocity
  // Eigen::Vector3d liquid_velocity;
  // liquid_velocity << particle.liquid_velocity_x, particle.liquid_velocity_y,
  //     particle.liquid_velocity_z;
  // // Initialise velocity
  // for (unsigned i = 0; i < Tdim; ++i)
  //   vector_properties_.at(mpm::properties::Vector::LiquidVelocity)(i) =
  //       liquid_velocity(i);
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
    if (this->material_id() == material->id() ||
        this->material_id() == std::numeric_limits<unsigned>::max()) {
      bool assign_mat = mpm::Particle<Tdim>::assign_material(material);
      if (!assign_mat) throw std::runtime_error("Material assignment failed");
      // Reinitialize state variables
      auto mat_state_vars = (this->material())->initialise_state_variables();
      if (mat_state_vars.size() == particle.nstate_vars) {
        unsigned i = 0;
        auto state_variables = (this->material())->state_variables();
        for (const auto& state_var : state_variables) {
          this->state_variables_[mpm::ParticlePhase::Solid].at(state_var) =
              particle.svars[i];
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
  // particle_data.pore_pressure =
  // scalar_properties_.at(mpm::properties::Scalar::PorePressure);
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
void mpm::TwoPhaseParticle<Tdim>::initialise() {
  mpm::Particle<Tdim>::initialise();
  liquid_strain_rate_.setZero();
  liquid_strain_.setZero();
  liquid_traction_.setZero();
  set_liquid_traction_ = false;

  // Initialize scalar properties
  scalar_properties_.emplace(
      std::make_pair(mpm::properties::Scalar::LiquidMass, double(0.)));
  scalar_properties_.emplace(
      std::make_pair(mpm::properties::Scalar::Porosity, double(0.)));
  scalar_properties_.emplace(
      std::make_pair(mpm::properties::Scalar::PorePressure, double(0.)));
  scalar_properties_.emplace(
      std::make_pair(mpm::properties::Scalar::LiquidMassDensity, double(0.)));

  // Initialize vector properties
  vector_properties_.emplace(
      std::make_pair(mpm::properties::Vector::Permeability, VectorDim::Zero()));
  vector_properties_.emplace(std::make_pair(
      mpm::properties::Vector::LiquidVelocity, VectorDim::Zero()));
  vector_properties_.emplace(
      std::make_pair(mpm::properties::Vector::DragForce, VectorDim::Zero()));

  // Initialize vector data properties
  this->properties_["liquid_strains"] = [&]() { return liquid_strain(); };
  this->properties_["liquid_velocities"] = [&]() { return liquid_velocity(); };
  this->properties_["pore_pressure"] = [&]() {
    Eigen::VectorXd pore_pressure = Eigen::VectorXd::Zero(3);
    // Total pore pressure
    pore_pressure[0] = this->pore_pressure();
    // Excessive pore pressure
    // pore_pressure[1] = this->excessive_pore_pressure();
    // Free surface
    // pore_pressure[2] = this->free_surface();

    return pore_pressure;
  };
}

//! Map both mixture and liquid internal force
template <>
inline void mpm::TwoPhaseParticle<1>::map_internal_force() noexcept {
  // Map mixture internal force
  // Initialise a vector of pore pressure
  Eigen::Matrix<double, 6, 1> total_stress = this->stress_;
  total_stress(0) -=
      scalar_properties_.at(mpm::properties::Scalar::PorePressure);
  // Compute nodal internal forces
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    // Compute force: -pstress * volume
    Eigen::Matrix<double, 1, 1> force;
    force[0] = dn_dx_(i, 0) * total_stress[0] + dn_dx_(i, 1) * total_stress[3];

    force *= -1. * mpm::Particle<1>::volume();

    nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Liquid, force);
  }

  // Map liquid internal force
  // Initialise a vector of pore pressure
  Eigen::Matrix<double, 6, 1> pressure;
  pressure.setZero();
  pressure(0) = -scalar_properties_.at(mpm::properties::Scalar::PorePressure);
  pressure(1) = -scalar_properties_.at(mpm::properties::Scalar::PorePressure);
  // Compute nodal internal forces
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    // Compute force: -pstress * volume
    Eigen::Matrix<double, 1, 1> force;
    force[0] = dn_dx_(i, 0) * pressure[0] * this->porosity();

    force *= -1. * mpm::Particle<1>::volume();

    nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Liquid, force);
  }
}

//! Map both mixture and liquid internal force
template <>
inline void mpm::TwoPhaseParticle<2>::map_internal_force() noexcept {
  // Map mixture internal force
  // Initialise a vector of pore pressure
  Eigen::Matrix<double, 6, 1> total_stress = this->stress_;
  total_stress(0) -=
      scalar_properties_.at(mpm::properties::Scalar::PorePressure);
  total_stress(1) -=
      scalar_properties_.at(mpm::properties::Scalar::PorePressure);
  // Compute nodal internal forces
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    // Compute force: -pstress * volume
    Eigen::Matrix<double, 2, 1> force;
    force[0] = dn_dx_(i, 0) * total_stress[0] + dn_dx_(i, 1) * total_stress[3];
    force[1] = dn_dx_(i, 1) * total_stress[1] + dn_dx_(i, 0) * total_stress[3];

    force *= -1. * mpm::Particle<2>::volume();

    nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Liquid, force);
  }

  // Map liquid internal force
  // Initialise a vector of pore pressure
  Eigen::Matrix<double, 6, 1> pressure;
  pressure.setZero();
  pressure(0) = -scalar_properties_.at(mpm::properties::Scalar::PorePressure);
  pressure(1) = -scalar_properties_.at(mpm::properties::Scalar::PorePressure);
  // Compute nodal internal forces
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    // Compute force: -pstress * volume
    Eigen::Matrix<double, 2, 1> force;
    force[0] = dn_dx_(i, 0) * pressure[0] * this->porosity();
    force[1] = dn_dx_(i, 1) * pressure[1] * this->porosity();

    force *= -1. * mpm::Particle<2>::volume();

    nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Liquid, force);
  }
}

//! Map both mixture and liquid internal force
template <>
inline void mpm::TwoPhaseParticle<3>::map_internal_force() noexcept {
  // Map mixture internal force
  // Initialise a vector of pore pressure
  Eigen::Matrix<double, 6, 1> total_stress = this->stress_;
  total_stress(0) -=
      scalar_properties_.at(mpm::properties::Scalar::PorePressure);
  total_stress(1) -=
      scalar_properties_.at(mpm::properties::Scalar::PorePressure);
  total_stress(2) -=
      scalar_properties_.at(mpm::properties::Scalar::PorePressure);
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

    force *= -1. * mpm::Particle<3>::volume();

    nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Liquid, force);
  }

  // Map liquid internal force
  // Initialise a vector of pore pressure
  Eigen::Matrix<double, 6, 1> pressure;
  pressure.setZero();
  pressure(0) = -scalar_properties_.at(mpm::properties::Scalar::PorePressure);
  pressure(1) = -scalar_properties_.at(mpm::properties::Scalar::PorePressure);
  pressure(2) = -scalar_properties_.at(mpm::properties::Scalar::PorePressure);
  // Compute nodal internal forces
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    // Compute force: -pstress * volume
    Eigen::Matrix<double, 3, 1> force;
    force[0] = dn_dx_(i, 0) * pressure[0] * this->porosity();
    force[1] = dn_dx_(i, 1) * pressure[1] * this->porosity();
    force[2] = dn_dx_(i, 2) * pressure[2] * this->porosity();

    force *= -1. * mpm::Particle<3>::volume();

    nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Liquid, force);
  }
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
    vector_properties_.at(mpm::properties::Vector::LiquidVelocity) +=
        acceleration * dt;
  } else {
    // Get interpolated nodal velocity
    Eigen::Matrix<double, Tdim, 1> velocity;
    velocity.setZero();

    for (unsigned i = 0; i < nodes_.size(); ++i)
      velocity += shapefn_(i) * nodes_[i]->velocity(mpm::ParticlePhase::Liquid);

    // Update particle velocity to interpolated nodal velocity
    vector_properties_.at(mpm::properties::Vector::LiquidVelocity) = velocity;
  }

  // Apply particle velocity constraints
  this->apply_particle_liquid_velocity_constraints();
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
    vector_properties_.at(mpm::properties::Vector::LiquidVelocity)(direction) =
        constraint.second;
  }
}

//! Map liquid traction force
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::map_liquid_traction_force() noexcept {
  if (this->set_liquid_traction_) {
    // Map particle liquid traction forces to nodes
    for (unsigned i = 0; i < nodes_.size(); ++i)
      nodes_[i]->update_external_force(
          true, mpm::ParticlePhase::Liquid,
          (-1. * shapefn_[i] * this->porosity() * this->liquid_traction_));
  }
}

//----------------------------------------------------------------------------
//! TODO
// Assign traction to the liquid phase
template <unsigned Tdim>
bool mpm::TwoPhaseParticle<Tdim>::assign_liquid_traction(unsigned direction,
                                                         double traction) {
  bool status = false;
  try {
    if (direction >= Tdim ||
        mpm::Particle<Tdim>::volume() == std::numeric_limits<double>::max()) {
      throw std::runtime_error(
          "Particle liquid traction property: volume / direction is invalid");
    }
    // Assign liquid traction
    liquid_traction_(direction) =
        traction * mpm::Particle<Tdim>::volume() / this->size_(direction);
    status = true;
    this->set_liquid_traction_ = true;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Compute pore pressure (compressible fluid)
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::compute_pore_pressure(double dt) noexcept {
  // Check if liquid material and cell pointer are set and positive
  assert(material_.at(mpm::ParticlePhase::Liquid) != nullptr &&
         cell_ != nullptr);
  // Apply free surface
  if (this->free_surface())
    scalar_properties_.at(mpm::properties::Scalar::PorePressure) = 0.0;
  // Compute pore pressure
  else {
    // get the bulk modulus of liquid
    double K = material_.at(mpm::ParticlePhase::Liquid)
                   ->template property<double>(std::string("bulk_modulus"));
    // Compute strain rate of liquid phase at centroid
    auto liquid_strain_rate_centroid = mpm::Particle<Tdim>::compute_strain_rate(
        dn_dx_centroid_, mpm::ParticlePhase::Liquid);
    // update pressure
    scalar_properties_.at(mpm::properties::Scalar::PorePressure) +=
        -dt * (K / this->scalar_property(mpm::properties::Scalar::Porosity)) *
        ((1 - this->scalar_property(mpm::properties::Scalar::Porosity)) *
             strain_rate_.head(Tdim).sum() +
         this->scalar_property(mpm::properties::Scalar::Porosity) *
             liquid_strain_rate_centroid.head(Tdim).sum());
  }
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
    if (this->free_surface())
      scalar_properties_.at(mpm::properties::Scalar::PorePressure) = 0.0;
    // Compute particle pore pressure from nodal values
    else {
      double pore_pressure = 0;
      for (unsigned i = 0; i < nodes_.size(); ++i)
        pore_pressure +=
            shapefn_(i) * nodes_[i]->pressure(mpm::ParticlePhase::Liquid);
      // Update pore liquid pressure to interpolated nodal pressure
      scalar_properties_.at(mpm::properties::Scalar::PorePressure) =
          pore_pressure;
    }
  }
  return status;
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

//! Assign particle pressure constraints
//! TODO
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
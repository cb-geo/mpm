//! Construct a particle with id and coordinates
template <unsigned Tdim, unsigned Tnphases>
mpm::Particle<Tdim, Tnphases>::Particle(Index id, const VectorDim& coord)
    : mpm::ParticleBase<Tdim>(id, coord) {
  this->initialise();
  cell_ = nullptr;

  //! Set material pointer to null
  material_.clear();
  for (unsigned i = 0; i < Tnphases; ++i) material_[i] = nullptr;

  //! Logger
  std::string logger =
      "particle" + std::to_string(Tdim) + "d::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

//! Construct a particle with id, coordinates and status
template <unsigned Tdim, unsigned Tnphases>
mpm::Particle<Tdim, Tnphases>::Particle(Index id, const VectorDim& coord,
                                        bool status)
    : mpm::ParticleBase<Tdim>(id, coord, status) {
  this->initialise();
  cell_ = nullptr;
  material_.clear();
  //! Logger
  std::string logger =
      "particle" + std::to_string(Tdim) + "d::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

//! Initialise particle data from HDF5
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::initialise_particle(
    const HDF5Particle& particle) {

  // TODO: Set phase
  const unsigned phase = 0;

  // Assign id
  this->id_ = particle.id;
  // Mass
  this->mass_(phase) = particle.mass;
  // Volume
  this->assign_volume(phase, particle.volume);
  // Set local size of particle
  Eigen::Vector3d psize;
  psize << particle.nsize_x, particle.nsize_y, particle.nsize_z;
  // Initialise particle size
  for (unsigned i = 0; i < Tdim; ++i) this->natural_size_(i) = psize(i);

  // Coordinates
  Eigen::Vector3d coordinates;
  coordinates << particle.coord_x, particle.coord_y, particle.coord_z;
  // Initialise coordinates
  for (unsigned i = 0; i < Tdim; ++i)
    this->coordinates_(i, phase) = coordinates(i);

  // Displacement
  Eigen::Vector3d displacement;
  displacement << particle.displacement_x, particle.displacement_y,
      particle.displacement_z;
  // Initialise displacement
  for (unsigned i = 0; i < Tdim; ++i)
    this->displacement_(i, phase) = displacement(i);

  // Velocity
  Eigen::Vector3d velocity;
  velocity << particle.velocity_x, particle.velocity_y, particle.velocity_z;
  // Initialise velocity
  for (unsigned i = 0; i < Tdim; ++i) this->velocity_(i, phase) = velocity(i);

  // Stress
  this->stress_.col(phase) << particle.stress_xx, particle.stress_yy,
      particle.stress_zz, particle.tau_xy, particle.tau_yz, particle.tau_xz;

  // Strain
  this->strain_.col(phase) << particle.strain_xx, particle.strain_yy,
      particle.strain_zz, particle.gamma_xy, particle.gamma_yz,
      particle.gamma_xz;

  // Volumetric strain
  this->volumetric_strain_centroid_(phase) = particle.epsilon_v;

  // Status
  this->status_ = particle.status;

  // Cell id
  this->cell_id_ = particle.cell_id;
  this->cell_ = nullptr;
  return true;
}

//! Return particle data in HDF5 format
template <unsigned Tdim, unsigned Tnphases>
mpm::HDF5Particle mpm::Particle<Tdim, Tnphases>::hdf5(unsigned phase) const {

  mpm::HDF5Particle particle_data;

  Eigen::Vector3d coordinates;
  coordinates.setZero();
  for (unsigned j = 0; j < Tdim; ++j) coordinates[j] = this->coordinates()[j];

  Eigen::Vector3d displacement;
  displacement.setZero();
  for (unsigned j = 0; j < Tdim; ++j)
    displacement[j] = this->displacement(phase)[j];

  Eigen::Vector3d velocity;
  velocity.setZero();
  for (unsigned j = 0; j < Tdim; ++j) velocity[j] = this->velocity(phase)[j];

  // Particle local size
  Eigen::Vector3d nsize;
  nsize.setZero();
  Eigen::VectorXd size = this->natural_size();
  for (unsigned j = 0; j < Tdim; ++j) nsize[j] = size[j];

  Eigen::Matrix<double, 6, 1> stress = this->stress(phase);

  Eigen::Matrix<double, 6, 1> strain = this->strain(phase);

  particle_data.id = this->id();
  particle_data.mass = this->mass(phase);
  particle_data.volume = this->volume(phase);
  particle_data.pressure = this->pressure(phase);

  particle_data.coord_x = coordinates[0];
  particle_data.coord_y = coordinates[1];
  particle_data.coord_z = coordinates[2];

  particle_data.displacement_x = displacement[0];
  particle_data.displacement_y = displacement[1];
  particle_data.displacement_z = displacement[2];

  particle_data.nsize_x = nsize[0];
  particle_data.nsize_y = nsize[1];
  particle_data.nsize_z = nsize[2];

  particle_data.velocity_x = velocity[0];
  particle_data.velocity_y = velocity[1];
  particle_data.velocity_z = velocity[2];

  particle_data.stress_xx = stress[0];
  particle_data.stress_yy = stress[1];
  particle_data.stress_zz = stress[2];
  particle_data.tau_xy = stress[3];
  particle_data.tau_yz = stress[4];
  particle_data.tau_xz = stress[5];

  particle_data.strain_xx = strain[0];
  particle_data.strain_yy = strain[1];
  particle_data.strain_zz = strain[2];
  particle_data.gamma_xy = strain[3];
  particle_data.gamma_yz = strain[4];
  particle_data.gamma_xz = strain[5];

  particle_data.epsilon_v = this->volumetric_strain_centroid(phase);

  particle_data.status = this->status();

  particle_data.cell_id = this->cell_id();

  return particle_data;
}

// Initialise particle properties
template <unsigned Tdim, unsigned Tnphases>
void mpm::Particle<Tdim, Tnphases>::initialise() {
  displacement_.setZero();
  dstrain_.setZero();
  mass_.setZero();
  natural_size_.setZero();
  pressure_.setZero();
  set_traction_ = false;
  size_.setZero();
  strain_rate_.setZero();
  strain_.setZero();
  stress_.setZero();
  traction_.setZero();
  velocity_.setZero();
  volume_.fill(std::numeric_limits<double>::max());
  volumetric_strain_centroid_.setZero();

  // Initialize vector data properties
  this->properties_["stresses"] = [&](unsigned phase) { return stress(phase); };
  this->properties_["strains"] = [&](unsigned phase) { return strain(phase); };
  this->properties_["velocities"] = [&](unsigned phase) {
    return velocity(phase);
  };
  this->properties_["displacements"] = [&](unsigned phase) {
    return displacement(phase);
  };
}

// Assign a cell to particle
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::assign_cell(
    const std::shared_ptr<Cell<Tdim>>& cellptr) {
  bool status = true;
  try {
    Eigen::Matrix<double, Tdim, 1> xi;
    // Assign cell to the new cell ptr, if point can be found in new cell
    if (cellptr->is_point_in_cell(this->coordinates_, &xi)) {
      // if a cell already exists remove particle from that cell
      if (cell_ != nullptr) cell_->remove_particle_id(this->id_);

      cell_ = cellptr;
      cell_id_ = cellptr->id();
      // Compute reference location of particle
      bool xi_status = this->compute_reference_location();
      if (!xi_status) return false;
      status = cell_->add_particle_id(this->id());
    } else {
      throw std::runtime_error("Point cannot be found in cell!");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Assign a cell to particle
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::assign_cell_xi(
    const std::shared_ptr<Cell<Tdim>>& cellptr,
    const Eigen::Matrix<double, Tdim, 1>& xi) {
  bool status = true;
  try {
    // Assign cell to the new cell ptr, if point can be found in new cell
    if (cellptr != nullptr) {
      // if a cell already exists remove particle from that cell
      if (cell_ != nullptr) cell_->remove_particle_id(this->id_);

      cell_ = cellptr;
      cell_id_ = cellptr->id();
      // Assign the reference location of particle
      bool xi_nan = false;

      // Check if point is within the cell
      for (unsigned i = 0; i < xi.size(); ++i)
        if (xi(i) < -1. || xi(i) > 1. || std::isnan(xi(i))) xi_nan = true;

      if (xi_nan == false)
        this->xi_ = xi;
      else
        return false;

      status = cell_->add_particle_id(this->id());
    } else {
      throw std::runtime_error("Point cannot be found in cell!");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Assign a cell id to particle
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::assign_cell_id(mpm::Index id) {
  bool status = false;
  try {
    // if a cell ptr is null
    if (cell_ == nullptr && id != std::numeric_limits<Index>::max()) {
      cell_id_ = id;
      status = true;
    } else {
      throw std::runtime_error("Invalid cell id or cell is already assigned!");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Remove cell for the particle
template <unsigned Tdim, unsigned Tnphases>
void mpm::Particle<Tdim, Tnphases>::remove_cell() {
  // if a cell is not nullptr
  if (cell_ != nullptr) cell_->remove_particle_id(this->id_);
  cell_id_ = std::numeric_limits<Index>::max();
}

// Assign a material to particle
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::assign_material(
    unsigned phase, const std::shared_ptr<Material<Tdim>>& material) {
  bool status = false;
  try {
    // Check if material is valid and properties are set
    if (material != nullptr) {
      material_.at(phase) = material;
      state_variables_ = material_.at(phase)->initialise_state_variables();
      status = true;
    } else {
      throw std::runtime_error("Material is undefined!");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return status;
}

// Compute reference location cell to particle
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::compute_reference_location() {
  bool status = true;
  try {
    // Check if particle has a valid cell ptr
    if (cell_ != nullptr) {
      // Compute local coordinates
      Eigen::Matrix<double, Tdim, 1> xi;
      // Check if the point is in cell
      if (cell_->is_point_in_cell(this->coordinates_, &xi)) {
        this->xi_ = xi;
        status = true;
      } else
        status = false;
    } else {
      throw std::runtime_error(
          "Cell is not initialised! "
          "cannot compute local reference coordinates of the particle");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Compute shape functions and gradients
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::compute_shapefn() {
  bool status = true;
  try {
    // Check if particle has a valid cell ptr
    if (cell_ != nullptr) {
      // Get element ptr of a cell
      const auto element = cell_->element_ptr();

      // Zero matrix
      Eigen::Matrix<double, Tdim, 1> zero;
      zero.setZero();

      // Compute shape function of the particle
      shapefn_ = element->shapefn(this->xi_, this->natural_size_, zero);

      // Compute bmatrix of the particle for reference cell
      bmatrix_ = element->bmatrix(this->xi_, cell_->nodal_coordinates(),
                                  this->natural_size_, zero);
    } else {
      throw std::runtime_error(
          "Cell is not initialised! "
          "cannot compute shapefns for the particle");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Assign volume to the particle
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::assign_volume(unsigned phase,
                                                  double volume) {
  bool status = true;
  try {
    if (volume <= 0.)
      throw std::runtime_error("Particle volume cannot be negative");

    this->volume_(phase) = volume;
    // Compute size of particle in each direction
    const double length =
        std::pow(this->volume_(phase), static_cast<double>(1. / Tdim));
    // Set particle size as length on each side
    this->size_.fill(length);

    if (cell_ != nullptr) {
      // Get element ptr of a cell
      const auto element = cell_->element_ptr();

      // Set local particle size based on length of element in natural
      // coordinates (cpGIMP Bardenhagen 2008 (pp485))
      this->natural_size_.fill(element->unit_element_length() /
                               cell_->nparticles());
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Compute volume of the particle
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::compute_volume(unsigned phase) {
  bool status = true;
  try {
    // Check if particle has a valid cell ptr
    if (cell_ != nullptr) {
      // Volume of the cell / # of particles
      this->assign_volume(phase, cell_->volume() / cell_->nparticles());
    } else {
      throw std::runtime_error(
          "Cell is not initialised! "
          "cannot compute volume for the particle");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Update volume based on the central strain rate
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::update_volume_strainrate(unsigned phase,
                                                             double dt) {
  bool status = true;
  try {
    // Check if particle has a valid cell ptr and a valid volume
    if (cell_ != nullptr &&
        volume_(phase) != std::numeric_limits<double>::max()) {
      // Compute at centroid
      // Strain rate for reduced integration
      Eigen::VectorXd strain_rate_centroid =
          cell_->compute_strain_rate_centroid(phase);
      this->volume_(phase) *= (1. + dt * strain_rate_centroid.head(Tdim).sum());
    } else {
      throw std::runtime_error(
          "Cell or volume is not initialised! cannot update particle volume");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Compute mass of particle
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::compute_mass(unsigned phase) {
  bool status = true;
  try {
    // Check if particle volume is set and material ptr is valid
    if (volume_(phase) != std::numeric_limits<double>::max() &&
        material_.at(phase) != nullptr) {
      // Mass = volume of particle * mass_density
      mass_density_(phase) = material_.at(phase)->template property<double>(
          std::string("density"));
      this->mass_(phase) = volume_(phase) * mass_density_(phase);
    } else {
      throw std::runtime_error(
          "Cell or material is invalid! cannot compute mass for the particle");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Map particle mass and momentum to nodes
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::map_mass_momentum_to_nodes(unsigned phase) {
  bool status = true;
  try {
    // Check if particle mass is set
    if (mass_(phase) != std::numeric_limits<double>::max()) {
      // Map particle mass and momentum to nodes
      this->cell_->map_mass_momentum_to_nodes(
          this->shapefn_, phase, mass_(phase), velocity_.col(phase));
    } else {
      throw std::runtime_error("Particle mass has not been computed");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Compute strain of the particle
template <unsigned Tdim, unsigned Tnphases>
void mpm::Particle<Tdim, Tnphases>::compute_strain(unsigned phase, double dt) {
  // Strain rate
  const auto strain_rate = cell_->compute_strain_rate(bmatrix_, phase);
  // particle_strain_rate
  Eigen::Matrix<double, 6, 1> particle_strain_rate;
  particle_strain_rate.setZero();
  // Set dimension of strain rate
  switch (Tdim) {
    case (1): {
      particle_strain_rate(0) = strain_rate(0);
      break;
    }
    case (2): {
      particle_strain_rate(0) = strain_rate(0);
      particle_strain_rate(1) = strain_rate(1);
      particle_strain_rate(3) = strain_rate(2);
      break;
    }
    default: {
      particle_strain_rate = strain_rate;
      break;
    }
  }

  // Check to see if value is below threshold
  for (unsigned i = 0; i < particle_strain_rate.size(); ++i)
    if (std::fabs(particle_strain_rate(i)) < 1.E-15)
      particle_strain_rate(i) = 0.;

  // Assign strain rate
  strain_rate_.col(phase) = particle_strain_rate;
  // Update dstrain
  dstrain_.col(phase) = particle_strain_rate * dt;
  // Update strain
  strain_.col(phase) += particle_strain_rate * dt;

  // Compute at centroid
  // Strain rate for reduced integration
  Eigen::VectorXd strain_rate_centroid =
      cell_->compute_strain_rate_centroid(phase);

  // Check to see if value is below threshold
  for (unsigned i = 0; i < strain_rate_centroid.size(); ++i)
    if (std::fabs(strain_rate_centroid(i)) < 1.E-15)
      strain_rate_centroid(i) = 0.;

  // Assign volumetric strain at centroid
  const double dvolumetric_strain = dt * strain_rate_centroid.head(Tdim).sum();
  volumetric_strain_centroid_(phase) += dvolumetric_strain;

  // Update thermodynamic pressure
  this->update_pressure(phase, dvolumetric_strain);
}

// Compute stress
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::compute_stress(unsigned phase) {
  bool status = true;
  try {
    // Check if material ptr is valid
    if (material_.at(phase) != nullptr) {
      Eigen::Matrix<double, 6, 1> dstrain = this->dstrain_.col(phase);
      // Calculate stress
      this->stress_.col(phase) = material_.at(phase)->compute_stress(
          this->stress_.col(phase), dstrain, this, &state_variables_);
    } else {
      throw std::runtime_error("Material is invalid");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Map body force
//! \param[in] phase Index corresponding to the phase
//! \param[in] pgravity Gravity of a particle
template <unsigned Tdim, unsigned Tnphases>
void mpm::Particle<Tdim, Tnphases>::map_body_force(unsigned phase,
                                                   const VectorDim& pgravity) {
  // Compute nodal body forces
  cell_->compute_nodal_body_force(this->shapefn_, phase, this->mass_(phase),
                                  pgravity);
}

//! Map internal force
//! \param[in] phase Index corresponding to the phase
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::map_internal_force(unsigned phase) {
  bool status = true;
  try {
    // Check if  material ptr is valid
    if (material_.at(phase) != nullptr) {
      // Compute nodal internal forces
      // -pstress * volume
      cell_->compute_nodal_internal_force(this->bmatrix_, phase,
                                          this->volume_(phase),
                                          -1. * this->stress_.col(phase));
    } else {
      throw std::runtime_error("Material is invalid");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Assign velocity to the particle
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::assign_velocity(
    unsigned phase, const Eigen::Matrix<double, Tdim, 1>& velocity) {
  bool status = false;
  try {
    if (phase >= Tnphases)
      throw std::runtime_error("Particle velocity: Invalid phase");

    // Assign velocity
    velocity_.col(phase) = velocity;
    status = true;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Assign traction to the particle
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::assign_traction(unsigned phase,
                                                    unsigned direction,
                                                    double traction) {
  bool status = false;
  try {
    if (phase >= Tnphases || direction >= Tdim ||
        this->volume_(phase) == std::numeric_limits<double>::max()) {
      throw std::runtime_error(
          "Particle traction property: volume / direction / phase is invalid");
    }
    // Assign traction
    traction_(direction, phase) =
        traction * this->volume_(phase) / this->size_(direction);
    status = true;
    this->set_traction_ = true;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Map traction force
//! \param[in] phase Index corresponding to the phase
template <unsigned Tdim, unsigned Tnphases>
void mpm::Particle<Tdim, Tnphases>::map_traction_force(unsigned phase) {
  if (this->set_traction_)
    // Map particle traction forces to nodes
    cell_->compute_nodal_traction_force(this->shapefn_, phase,
                                        this->traction_.col(phase));
}

// Compute updated position of the particle
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::compute_updated_position(unsigned phase,
                                                             double dt) {
  bool status = true;
  try {
    // Check if particle has a valid cell ptr
    if (cell_ != nullptr) {
      // Get interpolated nodal acceleration
      const Eigen::Matrix<double, Tdim, 1> nodal_acceleration =
          cell_->interpolate_nodal_acceleration(this->shapefn_, phase);

      // Update particle velocity from interpolated nodal acceleration
      this->velocity_.col(phase) += nodal_acceleration * dt;

      // Apply particle velocity constraints
      this->apply_particle_velocity_constraints();

      // Get interpolated nodal velocity
      const Eigen::Matrix<double, Tdim, 1> nodal_velocity =
          cell_->interpolate_nodal_velocity(this->shapefn_, phase);

      // New position  current position + velocity * dt
      this->coordinates_ += nodal_velocity * dt;
      // Update displacement (displacement is initialized from zero)
      this->displacement_ += nodal_velocity * dt;
    } else {
      throw std::runtime_error(
          "Cell is not initialised! "
          "cannot compute updated coordinates of the particle");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Compute updated position of the particle based on nodal velocity
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::compute_updated_position_velocity(
    unsigned phase, double dt) {
  bool status = true;
  try {
    // Check if particle has a valid cell ptr
    if (cell_ != nullptr) {
      // Get interpolated nodal velocity
      const Eigen::Matrix<double, Tdim, 1> nodal_velocity =
          cell_->interpolate_nodal_velocity(this->shapefn_, phase);

      // Update particle velocity to interpolated nodal velocity
      this->velocity_.col(phase) = nodal_velocity;

      // Apply particle velocity constraints
      this->apply_particle_velocity_constraints();

      // New position current position + velocity * dt
      this->coordinates_ += nodal_velocity * dt;
      // Update displacement (displacement is initialized from zero)
      this->displacement_ += nodal_velocity * dt;
    } else {
      throw std::runtime_error(
          "Cell is not initialised! "
          "cannot compute updated coordinates of the particle");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Update pressure
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::update_pressure(unsigned phase,
                                                    double dvolumetric_strain) {
  bool status = true;
  try {
    // Check if material ptr is valid
    if (material_.at(phase) != nullptr) {
      // Update pressure
      this->pressure_(phase) +=
          material_.at(phase)->thermodynamic_pressure(dvolumetric_strain);
    } else {
      throw std::runtime_error("Material is invalid");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Map particle pressure to nodes
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::map_pressure_to_nodes(unsigned phase) {
  bool status = true;
  try {
    // Check if particle mass is set
    if (mass_(phase) != std::numeric_limits<double>::max()) {
      // Map particle mass and momentum to nodes
      this->cell_->map_pressure_to_nodes(this->shapefn_, phase, mass_(phase),
                                         pressure_(phase));
    } else {
      throw std::runtime_error("Particle mass has not been computed");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Compute pressure smoothing of the particle based on nodal pressure
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::compute_pressure_smoothing(unsigned phase) {
  bool status = true;
  try {
    // Check if particle has a valid cell ptr
    if (cell_ != nullptr)
      // Update particle pressure to interpolated nodal pressure
      this->pressure_(phase) =
          cell_->interpolate_nodal_pressure(this->shapefn_, phase);
    else
      throw std::runtime_error(
          "Cell is not initialised! "
          "cannot compute pressure smoothing of the particle");

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assign particle velocity constraint
//! Constrain directions can take values between 0 and Dim * Nphases
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::assign_particle_velocity_constraint(
    unsigned dir, double velocity) {
  bool status = true;
  try {
    //! Constrain directions can take values between 0 and Dim * Nphases
    if (dir < (Tdim * Tnphases))
      this->particle_velocity_constraints_.insert(
          std::make_pair<unsigned, double>(static_cast<unsigned>(dir),
                                           static_cast<double>(velocity)));
    else
      throw std::runtime_error(
          "Particle velocity constraint direction is out of bounds");

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Apply particle velocity constraints
template <unsigned Tdim, unsigned Tnphases>
void mpm::Particle<Tdim, Tnphases>::apply_particle_velocity_constraints() {
  // Set particle velocity constraint
  for (const auto& constraint : this->particle_velocity_constraints_) {
    // Direction value in the constraint (0, Dim * Nphases)
    const unsigned dir = constraint.first;
    // Direction: dir % Tdim (modulus)
    const auto direction = static_cast<unsigned>(dir % Tdim);
    // Phase: Integer value of division (dir / Tdim)
    const auto phase = static_cast<unsigned>(dir / Tdim);
    this->velocity_(direction, phase) = constraint.second;
  }
}

//! Return particle vector data
template <unsigned Tdim, unsigned Tnphases>
Eigen::VectorXd mpm::Particle<Tdim, Tnphases>::vector_data(
    unsigned phase, const std::string& property) {
  return this->properties_.at(property)(phase);
}

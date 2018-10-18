//! Construct a particle with id and coordinates
template <unsigned Tdim, unsigned Tnphases>
mpm::Particle<Tdim, Tnphases>::Particle(Index id, const VectorDim& coord)
    : mpm::ParticleBase<Tdim>(id, coord) {
  this->initialise();
  cell_ = nullptr;
  material_ = nullptr;
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
  material_ = nullptr;
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

  // Coordinates
  Eigen::Vector3d coordinates;
  coordinates << particle.coord_x, particle.coord_y, particle.coord_z;
  // Initialise coordinates
  for (unsigned i = 0; i < Tdim; ++i)
    this->coordinates_(i, phase) = coordinates(i);

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
  return true;
}

// Initialise particle properties
template <unsigned Tdim, unsigned Tnphases>
void mpm::Particle<Tdim, Tnphases>::initialise() {
  dstrain_.setZero();
  mass_.setZero();
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
}

// Assign a cell to particle
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::assign_cell(
    const std::shared_ptr<Cell<Tdim>>& cellptr) {
  bool status = true;
  try {
    // Assign cell to the new cell ptr, if point can be found in new cell
    if (cellptr->is_point_in_cell(this->coordinates_)) {
      // if a cell already exists remove particle from that cell
      if (cell_ != nullptr) cell_->remove_particle_id(this->id_);

      cell_ = cellptr;
      cell_id_ = cellptr->id();
      // Calculate the reference location of particle
      this->compute_reference_location();
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
    const std::shared_ptr<Material<Tdim>>& material) {
  bool status = false;
  try {
    // Check if material is valid and properties are set
    if (material != nullptr && material->status()) {
      material_ = material;
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
      //#ifdef _MPM_ISOPARAMETRIC_
      if (cell_->is_point_in_cell(this->coordinates_))
        // Get reference location of a particle with isoparametric
        // transformation
        this->xi_ = cell_->transform_real_to_unit_cell(this->coordinates_);
      else
        return false;
      //#else
      // Get reference location of a particle on cartesian grid
      // this->xi_ = cell_->local_coordinates_point(this->coordinates_);
      //#endif
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
      // Compute local coordinates
      this->compute_reference_location();

      // Get element ptr of a cell
      const auto element = cell_->element_ptr();

      // Compute shape function of the particle
      shapefn_ = element->shapefn(this->xi_, this->size_,
                                  Eigen::Matrix<double, Tdim, 1>::Zero());
      // Compute bmatrix of the particle for reference cell
      bmatrix_ =
          element->bmatrix(this->xi_, cell_->nodal_coordinates(), this->size_,
                           Eigen::Matrix<double, Tdim, 1>::Zero());
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
void mpm::Particle<Tdim, Tnphases>::assign_volume(unsigned phase,
                                                  double volume) {
  this->volume_(phase) = volume;
  // Compute size of particle in each direction
  const double length =
      std::pow(this->volume_(phase), static_cast<double>(1. / Tdim));
  // Set particle size as length on each side
  this->size_ = Eigen::Matrix<double, Tdim, 1>::Constant(length);
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
        material_ != nullptr) {
      // Mass = volume of particle * density
      this->mass_(phase) = volume_(phase) * material_->property("density");
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
      throw std::runtime_error("Particle mass has not be computed");
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
  Eigen::VectorXd strain_rate = cell_->compute_strain_rate(bmatrix_, phase);
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
    if (material_ != nullptr) {
      Eigen::Matrix<double, 6, 1> dstrain = this->dstrain_.col(phase);
      // Calculate stress
      this->stress_.col(phase) =
          material_->compute_stress(this->stress_.col(phase), dstrain, this);
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
    if (material_ != nullptr) {
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
    if (phase < 0 || phase >= Tnphases || direction < 0 || direction >= Tdim ||
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
      Eigen::Matrix<double, Tdim, 1> nodal_acceleration =
          cell_->interpolate_nodal_acceleration(this->shapefn_, phase);

      // Update particle velocity from interpolated nodal acceleration
      this->velocity_.col(phase) += nodal_acceleration * dt;

      // Get interpolated nodal velocity
      Eigen::Matrix<double, Tdim, 1> nodal_velocity =
          cell_->interpolate_nodal_velocity(this->shapefn_, phase);

      // New position  current position + velocity * dt
      this->coordinates_ += nodal_velocity * dt;
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
      Eigen::Matrix<double, Tdim, 1> velocity =
          cell_->interpolate_nodal_velocity(this->shapefn_, phase);

      // Update particle velocity to interpolated nodal velocity
      this->velocity_.col(phase) = velocity;

      // New position current position + velocity * dt
      this->coordinates_ += this->velocity_.col(phase) * dt;
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
    if (material_ != nullptr) {
      // Update pressure
      this->pressure_(phase) +=
          material_->thermodynamic_pressure(dvolumetric_strain);
    } else {
      throw std::runtime_error("Material is invalid");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

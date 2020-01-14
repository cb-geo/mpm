//! Construct a particle with id and coordinates
template <unsigned Tdim>
mpm::Particle<Tdim>::Particle(Index id, const VectorDim& coord)
    : mpm::ParticleBase<Tdim>(id, coord) {
  this->initialise();
  // Clear cell ptr
  cell_ = nullptr;
  // Nodes
  nodes_.clear();
  // Set material pointer to null
  material_ = nullptr;
  // Logger
  std::string logger =
      "particle" + std::to_string(Tdim) + "d::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

//! Construct a particle with id, coordinates and status
template <unsigned Tdim>
mpm::Particle<Tdim>::Particle(Index id, const VectorDim& coord, bool status)
    : mpm::ParticleBase<Tdim>(id, coord, status) {
  this->initialise();
  cell_ = nullptr;
  nodes_.clear();
  material_ = nullptr;
  //! Logger
  std::string logger =
      "particle" + std::to_string(Tdim) + "d::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

//! Initialise particle data from HDF5
template <unsigned Tdim>
bool mpm::Particle<Tdim>::initialise_particle(const HDF5Particle& particle) {

  // Assign id
  this->id_ = particle.id;
  // Mass
  this->mass_ = particle.mass;
  // Volume
  this->assign_volume(particle.volume);
  // Mass Density
  this->mass_density_ = particle.mass / particle.volume;
  // Set local size of particle
  Eigen::Vector3d psize;
  psize << particle.nsize_x, particle.nsize_y, particle.nsize_z;
  // Initialise particle size
  for (unsigned i = 0; i < Tdim; ++i) this->natural_size_(i) = psize(i);

  // Coordinates
  Eigen::Vector3d coordinates;
  coordinates << particle.coord_x, particle.coord_y, particle.coord_z;
  // Initialise coordinates
  for (unsigned i = 0; i < Tdim; ++i) this->coordinates_(i) = coordinates(i);

  // Displacement
  Eigen::Vector3d displacement;
  displacement << particle.displacement_x, particle.displacement_y,
      particle.displacement_z;
  // Initialise displacement
  for (unsigned i = 0; i < Tdim; ++i) this->displacement_(i) = displacement(i);

  // Velocity
  Eigen::Vector3d velocity;
  velocity << particle.velocity_x, particle.velocity_y, particle.velocity_z;
  // Initialise velocity
  for (unsigned i = 0; i < Tdim; ++i) this->velocity_(i) = velocity(i);

  // Stress
  this->stress_[0] = particle.stress_xx;
  this->stress_[1] = particle.stress_yy;
  this->stress_[2] = particle.stress_zz;
  this->stress_[3] = particle.tau_xy;
  this->stress_[4] = particle.tau_yz;
  this->stress_[5] = particle.tau_xz;

  // Strain
  this->strain_[0] = particle.strain_xx;
  this->strain_[1] = particle.strain_yy;
  this->strain_[2] = particle.strain_zz;
  this->strain_[3] = particle.gamma_xy;
  this->strain_[4] = particle.gamma_yz;
  this->strain_[5] = particle.gamma_xz;

  // Volumetric strain
  this->volumetric_strain_centroid_ = particle.epsilon_v;

  // Status
  this->status_ = particle.status;

  // Cell id
  this->cell_id_ = particle.cell_id;
  this->cell_ = nullptr;

  // Clear nodes
  this->nodes_.clear();

  // Material id
  this->material_id_ = particle.material_id;

  return true;
}

//! Initialise particle data from HDF5
template <unsigned Tdim>
bool mpm::Particle<Tdim>::initialise_particle(
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
        for (const auto& mat_state_var : mat_state_vars) {
          this->state_variables_[mat_state_var.first] = particle.svars[i];
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
mpm::HDF5Particle mpm::Particle<Tdim>::hdf5() const {

  mpm::HDF5Particle particle_data;

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

  particle_data.id = this->id();
  particle_data.mass = this->mass();
  particle_data.volume = this->volume();
  particle_data.pressure =
      (state_variables_.find("pressure") != state_variables_.end())
          ? state_variables_.at("pressure")
          : 0.;

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

  particle_data.epsilon_v = this->volumetric_strain_centroid_;

  particle_data.status = this->status();

  particle_data.cell_id = this->cell_id();

  particle_data.material_id = this->material_id();

  // Write state variables
  if (material_ != nullptr) {
    particle_data.nstate_vars = state_variables_.size();
    if (state_variables_.size() > 20)
      throw std::runtime_error("# of state variables cannot be more than 20");
    unsigned i = 0;
    for (const auto& state_var : this->state_variables_) {
      particle_data.svars[i] = state_var.second;
      ++i;
    }
  }

  return particle_data;
}

// Initialise particle properties
template <unsigned Tdim>
void mpm::Particle<Tdim>::initialise() {
  displacement_.setZero();
  dstrain_.setZero();
  mass_ = 0.;
  natural_size_.setZero();
  set_traction_ = false;
  size_.setZero();
  strain_rate_.setZero();
  strain_.setZero();
  stress_.setZero();
  traction_.setZero();
  velocity_.setZero();
  volume_ = std::numeric_limits<double>::max();
  volumetric_strain_centroid_ = 0.;

  // Initialize vector data properties
  this->properties_["stresses"] = [&]() { return stress(); };
  this->properties_["strains"] = [&]() { return strain(); };
  this->properties_["velocities"] = [&]() { return velocity(); };
  this->properties_["displacements"] = [&]() { return displacement(); };
}

// Assign a cell to particle
template <unsigned Tdim>
bool mpm::Particle<Tdim>::assign_cell(
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
      // dn_dx centroid
      dn_dx_centroid_ = cell_->dn_dx_centroid();
      // Copy nodal pointer to cell
      nodes_.clear();
      nodes_ = cell_->nodes();

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
template <unsigned Tdim>
bool mpm::Particle<Tdim>::assign_cell_xi(
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
      // dn_dx centroid
      dn_dx_centroid_ = cell_->dn_dx_centroid();
      // Copy nodal pointer to cell
      nodes_.clear();
      nodes_ = cell_->nodes();

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
template <unsigned Tdim>
bool mpm::Particle<Tdim>::assign_cell_id(mpm::Index id) {
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
template <unsigned Tdim>
void mpm::Particle<Tdim>::remove_cell() {
  // if a cell is not nullptr
  if (cell_ != nullptr) cell_->remove_particle_id(this->id_);
  cell_id_ = std::numeric_limits<Index>::max();
  // Clear all the nodes
  nodes_.clear();
}

// Assign a material to particle
template <unsigned Tdim>
bool mpm::Particle<Tdim>::assign_material(
    const std::shared_ptr<Material<Tdim>>& material) {
  bool status = false;
  try {
    // Check if material is valid and properties are set
    if (material != nullptr) {
      material_ = material;
      material_id_ = material_->id();
      state_variables_ = material_->initialise_state_variables();
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
template <unsigned Tdim>
bool mpm::Particle<Tdim>::compute_reference_location() {
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
template <unsigned Tdim>
void mpm::Particle<Tdim>::compute_shapefn() noexcept {
  // Check if particle has a valid cell ptr
  assert(cell_ != nullptr);
  // Get element ptr of a cell
  const auto element = cell_->element_ptr();

  // Zero matrix
  Eigen::Matrix<double, Tdim, 1> zero = Eigen::Matrix<double, Tdim, 1>::Zero();

  // Compute shape function of the particle
  shapefn_ = element->shapefn(this->xi_, this->natural_size_, zero);

  // Compute dN/dx
  dn_dx_ = element->dn_dx(this->xi_, cell_->nodal_coordinates(),
                          this->natural_size_, zero);
}

// Assign volume to the particle
template <unsigned Tdim>
bool mpm::Particle<Tdim>::assign_volume(double volume) {
  bool status = true;
  try {
    if (volume <= 0.)
      throw std::runtime_error("Particle volume cannot be negative");

    this->volume_ = volume;
    // Compute size of particle in each direction
    const double length =
        std::pow(this->volume_, static_cast<double>(1. / Tdim));
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
template <unsigned Tdim>
bool mpm::Particle<Tdim>::compute_volume() {
  bool status = true;
  try {
    // Check if particle has a valid cell ptr
    if (cell_ != nullptr) {
      // Volume of the cell / # of particles
      this->assign_volume(cell_->volume() / cell_->nparticles());
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
template <unsigned Tdim>
bool mpm::Particle<Tdim>::update_volume() {
  bool status = true;
  try {
    // Check if particle has a valid cell ptr and a valid volume
    if (cell_ != nullptr && volume_ != std::numeric_limits<double>::max()) {
      // Compute at centroid
      // Strain rate for reduced integration
      this->volume_ *= (1. + dvolumetric_strain_);
      this->mass_density_ = this->mass_density_ / (1. + dvolumetric_strain_);
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
template <unsigned Tdim>
bool mpm::Particle<Tdim>::compute_mass() {
  bool status = true;
  try {
    // Check if particle volume is set and material ptr is valid
    if (volume_ != std::numeric_limits<double>::max() && material_ != nullptr) {
      // Mass = volume of particle * mass_density
      this->mass_density_ =
          material_->template property<double>(std::string("density"));
      this->mass_ = volume_ * mass_density_;
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
template <unsigned Tdim>
bool mpm::Particle<Tdim>::map_mass_momentum_to_nodes() {
  bool status = true;
  try {
    // Check if particle mass is set
    if (mass_ != std::numeric_limits<double>::max()) {
      // Map mass and momentum to nodes
      for (unsigned i = 0; i < nodes_.size(); ++i) {
        nodes_[i]->update_mass(true, mpm::ParticlePhase::Solid,
                               mass_ * shapefn_[i]);
        nodes_[i]->update_momentum(true, mpm::ParticlePhase::Solid,
                                   mass_ * shapefn_[i] * velocity_);
      }
    } else {
      throw std::runtime_error("Particle mass has not been computed");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Compute strain rate of the particle
template <>
inline Eigen::Matrix<double, 6, 1> mpm::Particle<1>::compute_strain_rate(
    const Eigen::MatrixXd& dn_dx, unsigned phase) {
  // Define strain rate
  Eigen::Matrix<double, 6, 1> strain_rate = Eigen::Matrix<double, 6, 1>::Zero();

  for (unsigned i = 0; i < this->nodes_.size(); ++i) {
    Eigen::Matrix<double, 1, 1> vel = nodes_[i]->velocity(phase);
    strain_rate[0] += dn_dx(i, 0) * vel[0];
  }

  if (std::fabs(strain_rate(0)) < 1.E-15) strain_rate[0] = 0.;
  return strain_rate;
}

// Compute strain rate of the particle
template <>
inline Eigen::Matrix<double, 6, 1> mpm::Particle<2>::compute_strain_rate(
    const Eigen::MatrixXd& dn_dx, unsigned phase) {
  // Define strain rate
  Eigen::Matrix<double, 6, 1> strain_rate = Eigen::Matrix<double, 6, 1>::Zero();

  for (unsigned i = 0; i < this->nodes_.size(); ++i) {
    Eigen::Matrix<double, 2, 1> vel = nodes_[i]->velocity(phase);
    strain_rate[0] += dn_dx(i, 0) * vel[0];
    strain_rate[1] += dn_dx(i, 1) * vel[1];
    strain_rate[3] += dn_dx(i, 1) * vel[0] + dn_dx(i, 0) * vel[1];
  }

  if (std::fabs(strain_rate[0]) < 1.E-15) strain_rate[0] = 0.;
  if (std::fabs(strain_rate[1]) < 1.E-15) strain_rate[1] = 0.;
  if (std::fabs(strain_rate[3]) < 1.E-15) strain_rate[3] = 0.;
  return strain_rate;
}

// Compute strain rate of the particle
template <>
inline Eigen::Matrix<double, 6, 1> mpm::Particle<3>::compute_strain_rate(
    const Eigen::MatrixXd& dn_dx, unsigned phase) {
  // Define strain rate
  Eigen::Matrix<double, 6, 1> strain_rate = Eigen::Matrix<double, 6, 1>::Zero();

  for (unsigned i = 0; i < this->nodes_.size(); ++i) {
    Eigen::Matrix<double, 3, 1> vel = nodes_[i]->velocity(phase);
    strain_rate[0] += dn_dx(i, 0) * vel[0];
    strain_rate[1] += dn_dx(i, 1) * vel[1];
    strain_rate[2] += dn_dx(i, 2) * vel[2];
    strain_rate[3] += dn_dx(i, 1) * vel[0] + dn_dx(i, 0) * vel[1];
    strain_rate[4] += dn_dx(i, 2) * vel[1] + dn_dx(i, 1) * vel[2];
    strain_rate[5] += dn_dx(i, 2) * vel[0] + dn_dx(i, 0) * vel[2];
  }

  for (unsigned i = 0; i < strain_rate.size(); ++i)
    if (std::fabs(strain_rate[i]) < 1.E-15) strain_rate[i] = 0.;
  return strain_rate;
}

// Compute strain of the particle
template <unsigned Tdim>
void mpm::Particle<Tdim>::compute_strain(double dt) {
  // Assign strain rate
  strain_rate_ = this->compute_strain_rate(dn_dx_, mpm::ParticlePhase::Solid);
  // Update dstrain
  dstrain_ = strain_rate_ * dt;
  // Update strain
  strain_ += dstrain_;

  // Compute at centroid
  // Strain rate for reduced integration
  const Eigen::Matrix<double, 6, 1> strain_rate_centroid =
      this->compute_strain_rate(dn_dx_centroid_, mpm::ParticlePhase::Solid);

  // Assign volumetric strain at centroid
  dvolumetric_strain_ = dt * strain_rate_centroid.head(Tdim).sum();
  volumetric_strain_centroid_ += dvolumetric_strain_;
}

// Compute stress
template <unsigned Tdim>
bool mpm::Particle<Tdim>::compute_stress() {
  bool status = true;
  try {
    // Check if material ptr is valid
    if (material_ != nullptr) {
      Eigen::Matrix<double, 6, 1> dstrain = this->dstrain_;
      // Calculate stress
      this->stress_ = material_->compute_stress(this->stress_, dstrain, this,
                                                &state_variables_);
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
template <unsigned Tdim>
void mpm::Particle<Tdim>::map_body_force(const VectorDim& pgravity) {
  // Compute nodal body forces
  for (unsigned i = 0; i < nodes_.size(); ++i)
    nodes_[i]->update_external_force(true, mpm::ParticlePhase::Solid,
                                     (pgravity * mass_ * shapefn_(i)));
}

//! Map internal force
template <>
inline void mpm::Particle<1>::map_internal_force() {
  // Compute nodal internal forces
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    // Compute force: -pstress * volume
    Eigen::Matrix<double, 1, 1> force;
    force[0] = -1. * dn_dx_(i, 0) * volume_ * stress_[0];

    nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Solid, force);
  }
}

//! Map internal force
template <>
inline void mpm::Particle<2>::map_internal_force() {
  // Compute nodal internal forces
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    // Compute force: -pstress * volume
    Eigen::Matrix<double, 2, 1> force;
    force[0] = dn_dx_(i, 0) * stress_[0] + dn_dx_(i, 1) * stress_[3];
    force[1] = dn_dx_(i, 1) * stress_[1] + dn_dx_(i, 0) * stress_[3];

    force *= -1. * this->volume_;

    nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Solid, force);
  }
}

//! Map internal force
template <>
inline void mpm::Particle<3>::map_internal_force() {
  // Compute nodal internal forces
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    // Compute force: -pstress * volume
    Eigen::Matrix<double, 3, 1> force;
    force[0] = dn_dx_(i, 0) * stress_[0] + dn_dx_(i, 1) * stress_[3] +
               dn_dx_(i, 2) * stress_[5];

    force[1] = dn_dx_(i, 1) * stress_[1] + dn_dx_(i, 0) * stress_[3] +
               dn_dx_(i, 2) * stress_[4];

    force[2] = dn_dx_(i, 2) * stress_[2] + dn_dx_(i, 1) * stress_[4] +
               dn_dx_(i, 0) * stress_[5];

    force *= -1. * this->volume_;

    nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Solid, force);
  }
}

// Assign velocity to the particle
template <unsigned Tdim>
bool mpm::Particle<Tdim>::assign_velocity(
    const Eigen::Matrix<double, Tdim, 1>& velocity) {
  // Assign velocity
  velocity_ = velocity;
  return true;
}

// Assign traction to the particle
template <unsigned Tdim>
bool mpm::Particle<Tdim>::assign_traction(unsigned direction, double traction) {
  bool status = false;
  try {
    if (direction >= Tdim ||
        this->volume_ == std::numeric_limits<double>::max()) {
      throw std::runtime_error(
          "Particle traction property: volume / direction is invalid");
    }
    // Assign traction
    traction_(direction) = traction * this->volume_ / this->size_(direction);
    status = true;
    this->set_traction_ = true;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Map traction force
template <unsigned Tdim>
void mpm::Particle<Tdim>::map_traction_force() {
  if (this->set_traction_) {
    // Map particle traction forces to nodes
    for (unsigned i = 0; i < nodes_.size(); ++i)
      nodes_[i]->update_external_force(true, mpm::ParticlePhase::Solid,
                                       (shapefn_[i] * traction_));
  }
}

// Compute updated position of the particle
template <unsigned Tdim>
bool mpm::Particle<Tdim>::compute_updated_position(double dt,
                                                   bool velocity_update) {
  bool status = true;
  try {
    // Check if particle has a valid cell ptr
    if (cell_ != nullptr) {
      // Get interpolated nodal velocity
      Eigen::Matrix<double, Tdim, 1> nodal_velocity;
      nodal_velocity.setZero();

      for (unsigned i = 0; i < nodes_.size(); ++i)
        nodal_velocity +=
            shapefn_[i] * nodes_[i]->velocity(mpm::ParticlePhase::Solid);

      // Acceleration update
      if (!velocity_update) {
        // Get interpolated nodal acceleration
        Eigen::Matrix<double, Tdim, 1> nodal_acceleration;
        nodal_acceleration.setZero();
        for (unsigned i = 0; i < nodes_.size(); ++i)
          nodal_acceleration +=
              shapefn_[i] * nodes_[i]->acceleration(mpm::ParticlePhase::Solid);

        // Update particle velocity from interpolated nodal acceleration
        this->velocity_ += nodal_acceleration * dt;
      }
      // Update particle velocity using interpolated nodal velocity
      else
        this->velocity_ = nodal_velocity;

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

//! Map particle pressure to nodes
template <unsigned Tdim>
bool mpm::Particle<Tdim>::map_pressure_to_nodes() {
  bool status = true;
  try {
    // Check if particle mass is set and state variable pressure is found
    if (mass_ != std::numeric_limits<double>::max() &&
        (state_variables_.find("pressure") != state_variables_.end())) {
      // Map particle pressure to nodes
      for (unsigned i = 0; i < nodes_.size(); ++i)
        nodes_[i]->update_mass_pressure(
            mpm::ParticlePhase::Solid,
            shapefn_[i] * mass_ * state_variables_.at("pressure"));

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
template <unsigned Tdim>
bool mpm::Particle<Tdim>::compute_pressure_smoothing() {
  bool status = true;
  try {
    // Check if particle has a valid cell ptr
    if (cell_ != nullptr &&
        (state_variables_.find("pressure") != state_variables_.end())) {
      double pressure = 0.;
      // Update particle pressure to interpolated nodal pressure
      for (unsigned i = 0; i < this->nodes_.size(); ++i)
        pressure +=
            shapefn_[i] * nodes_[i]->pressure(mpm::ParticlePhase::Solid);

      state_variables_.at("pressure") = pressure;
    } else {
      throw std::runtime_error(
          "Cell is not initialised! "
          "cannot compute pressure smoothing of the particle");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Apply particle velocity constraints
template <unsigned Tdim>
void mpm::Particle<Tdim>::apply_particle_velocity_constraints(unsigned dir,
                                                              double velocity) {
  // Set particle velocity constraint
  this->velocity_(dir) = velocity;
}

//! Return particle vector data
template <unsigned Tdim>
Eigen::VectorXd mpm::Particle<Tdim>::vector_data(const std::string& property) {
  return this->properties_.at(property)();
}

//! Assign material id of this particle to nodes
template <unsigned Tdim>
void mpm::Particle<Tdim>::append_material_id_to_nodes() const {
  for (unsigned i = 0; i < nodes_.size(); ++i)
    nodes_[i]->append_material_id(material_id_);
}

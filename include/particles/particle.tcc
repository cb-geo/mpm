//! Construct a particle with id and coordinates
template <unsigned Tdim>
mpm::Particle<Tdim>::Particle(Index id, const VectorDim& coord)
    : mpm::ParticleBase<Tdim>(id, coord) {
  this->initialise();
  // Clear cell ptr
  cell_ = nullptr;
  // Nodes
  nodes_.clear();
  // Set material containers
  this->initialise_material(1);
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
  // Set material containers
  this->initialise_material(1);
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
  this->volume_ = particle.volume;
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
  this->material_id_[mpm::ParticlePhase::Solid] = particle.material_id;

  return true;
}

//! Initialise particle data from HDF5
template <unsigned Tdim>
bool mpm::Particle<Tdim>::initialise_particle(
    const HDF5Particle& particle,
    const std::shared_ptr<mpm::Material<Tdim>>& material) {
  bool status = this->initialise_particle(particle);
  if (material != nullptr) {
    if (this->material_id() == material->id() ||
        this->material_id() == std::numeric_limits<unsigned>::max()) {
      bool assign_mat = this->assign_material(material);
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
      (state_variables_[mpm::ParticlePhase::Solid].find("pressure") !=
       state_variables_[mpm::ParticlePhase::Solid].end())
          ? state_variables_[mpm::ParticlePhase::Solid].at("pressure")
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
  if (this->material() != nullptr) {
    particle_data.nstate_vars =
        state_variables_[mpm::ParticlePhase::Solid].size();
    if (state_variables_[mpm::ParticlePhase::Solid].size() > 20)
      throw std::runtime_error("# of state variables cannot be more than 20");
    unsigned i = 0;
    auto state_variables = (this->material())->state_variables();
    for (const auto& state_var : state_variables) {
      particle_data.svars[i] =
          state_variables_[mpm::ParticlePhase::Solid].at(state_var);
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

//! Initialise particle material container
template <unsigned Tdim>
void mpm::Particle<Tdim>::initialise_material(unsigned phase_size) {
  material_.resize(phase_size);
  material_id_.resize(phase_size);
  state_variables_.resize(phase_size);
  std::fill(material_.begin(), material_.end(), nullptr);
  std::fill(material_id_.begin(), material_id_.end(),
            std::numeric_limits<unsigned>::max());
  std::fill(state_variables_.begin(), state_variables_.end(), mpm::dense_map());
}

//! Assign material history variables
template <unsigned Tdim>
bool mpm::Particle<Tdim>::assign_material_state_vars(
    const mpm::dense_map& state_vars,
    const std::shared_ptr<mpm::Material<Tdim>>& material, unsigned phase) {
  bool status = false;
  if (material != nullptr && this->material(phase) != nullptr &&
      this->material_id(phase) == material->id()) {
    // Clone state variables
    auto mat_state_vars = (this->material(phase))->initialise_state_variables();
    if (state_variables_[phase].size() == state_vars.size() &&
        mat_state_vars.size() == state_vars.size()) {
      this->state_variables_[phase] = state_vars;
      status = true;
    }
  }
  return status;
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
    const std::shared_ptr<Material<Tdim>>& material, unsigned phase) {
  bool status = false;
  try {
    // Check if material is valid and properties are set
    if (material != nullptr) {
      material_.at(phase) = material;
      material_id_.at(phase) = material_[phase]->id();
      state_variables_.at(phase) =
          material_[phase]->initialise_state_variables();
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
bool mpm::Particle<Tdim>::compute_reference_location() noexcept {
  // Set status of compute reference location
  bool status = false;
  // Compute local coordinates
  Eigen::Matrix<double, Tdim, 1> xi;
  // Check if the point is in cell
  if (cell_ != nullptr && cell_->is_point_in_cell(this->coordinates_, &xi)) {
    this->xi_ = xi;
    status = true;
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

      // Set local particle length based on length of element in natural
      // coordinates. Length/(npartices^(1/Dimension))
      this->natural_size_.fill(
          element->unit_element_length() /
          std::pow(cell_->nparticles(), static_cast<double>(1. / Tdim)));
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Compute volume of the particle
template <unsigned Tdim>
void mpm::Particle<Tdim>::compute_volume() noexcept {
  // Check if particle has a valid cell ptr
  assert(cell_ != nullptr);
  // Volume of the cell / # of particles
  this->assign_volume(cell_->volume() / cell_->nparticles());
}

// Update volume based on the central strain rate
template <unsigned Tdim>
void mpm::Particle<Tdim>::update_volume() noexcept {
  // Check if particle has a valid cell ptr and a valid volume
  assert(cell_ != nullptr && volume_ != std::numeric_limits<double>::max());
  // Compute at centroid
  // Strain rate for reduced integration
  this->volume_ *= (1. + dvolumetric_strain_);
  this->mass_density_ = this->mass_density_ / (1. + dvolumetric_strain_);
}

// Compute mass of particle
template <unsigned Tdim>
void mpm::Particle<Tdim>::compute_mass() noexcept {
  // Check if particle volume is set and material ptr is valid
  assert(volume_ != std::numeric_limits<double>::max() &&
         this->material() != nullptr);
  // Mass = volume of particle * mass_density
  this->mass_density_ =
      (this->material())->template property<double>(std::string("density"));
  this->mass_ = volume_ * mass_density_;
}

//! Map particle mass and momentum to nodes
template <unsigned Tdim>
void mpm::Particle<Tdim>::map_mass_momentum_to_nodes() noexcept {
  // Check if particle mass is set
  assert(mass_ != std::numeric_limits<double>::max());

  // Map mass and momentum to nodes
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    nodes_[i]->update_mass(true, mpm::ParticlePhase::Solid,
                           mass_ * shapefn_[i]);
    nodes_[i]->update_momentum(true, mpm::ParticlePhase::Solid,
                               mass_ * shapefn_[i] * velocity_);
  }
}

//! Map multimaterial properties to nodes
template <unsigned Tdim>
void mpm::Particle<Tdim>::map_multimaterial_mass_momentum_to_nodes() noexcept {
  // Check if particle mass is set
  assert(mass_ != std::numeric_limits<double>::max());

  // Unit 1x1 Eigen matrix to be used with scalar quantities
  Eigen::Matrix<double, 1, 1> nodal_mass;

  // Map mass and momentum to nodal property taking into account the material id
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    nodal_mass(0, 0) = mass_ * shapefn_[i];
    nodes_[i]->update_property(true, "masses", nodal_mass, this->material_id(),
                               1);
    nodes_[i]->update_property(true, "momenta", velocity_ * nodal_mass,
                               this->material_id(), Tdim);
  }
}

//! Map multimaterial displacements to nodes
template <unsigned Tdim>
void mpm::Particle<Tdim>::map_multimaterial_displacements_to_nodes() noexcept {
  // Check if particle mass is set
  assert(mass_ != std::numeric_limits<double>::max());

  // Map displacements to nodal property and divide it by the respective
  // nodal-material mass
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    const auto& displacement = mass_ * shapefn_[i] * displacement_;
    nodes_[i]->update_property(true, "displacements", displacement,
                               this->material_id(), Tdim);
  }
}

//! Map multimaterial domain gradients to nodes
template <unsigned Tdim>
void mpm::Particle<
    Tdim>::map_multimaterial_domain_gradients_to_nodes() noexcept {
  // Check if particle volume is set
  assert(volume_ != std::numeric_limits<double>::max());

  // Map domain gradients to nodal property. The domain gradients is defined as
  // the gradient of the particle volume
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    Eigen::Matrix<double, Tdim, 1> gradient;
    for (unsigned j = 0; j < Tdim; ++j) gradient[j] = volume_ * dn_dx_(i, j);
    nodes_[i]->update_property(true, "domain_gradients", gradient,
                               this->material_id(), Tdim);
  }
}

// Compute strain rate of the particle
template <>
inline Eigen::Matrix<double, 6, 1> mpm::Particle<1>::compute_strain_rate(
    const Eigen::MatrixXd& dn_dx, unsigned phase) noexcept {
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
    const Eigen::MatrixXd& dn_dx, unsigned phase) noexcept {
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
    const Eigen::MatrixXd& dn_dx, unsigned phase) noexcept {
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
void mpm::Particle<Tdim>::compute_strain(double dt) noexcept {
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
void mpm::Particle<Tdim>::compute_stress() noexcept {
  // Check if material ptr is valid
  assert(this->material() != nullptr);
  // Calculate stress
  this->stress_ =
      (this->material())
          ->compute_stress(stress_, dstrain_, this,
                           &state_variables_[mpm::ParticlePhase::Solid]);
}

//! Map body force
template <unsigned Tdim>
void mpm::Particle<Tdim>::map_body_force(const VectorDim& pgravity) noexcept {
  // Compute nodal body forces
  for (unsigned i = 0; i < nodes_.size(); ++i)
    nodes_[i]->update_external_force(true, mpm::ParticlePhase::Solid,
                                     (pgravity * mass_ * shapefn_(i)));
}

//! Map internal force
template <>
inline void mpm::Particle<1>::map_internal_force() noexcept {
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
inline void mpm::Particle<2>::map_internal_force() noexcept {
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
inline void mpm::Particle<3>::map_internal_force() noexcept {
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
void mpm::Particle<Tdim>::map_traction_force() noexcept {
  if (this->set_traction_) {
    // Map particle traction forces to nodes
    for (unsigned i = 0; i < nodes_.size(); ++i)
      nodes_[i]->update_external_force(true, mpm::ParticlePhase::Solid,
                                       (shapefn_[i] * traction_));
  }
}

// Compute updated position of the particle
template <unsigned Tdim>
void mpm::Particle<Tdim>::compute_updated_position(
    double dt, bool velocity_update) noexcept {
  // Check if particle has a valid cell ptr
  assert(cell_ != nullptr);
  // Get interpolated nodal velocity
  Eigen::Matrix<double, Tdim, 1> nodal_velocity =
      Eigen::Matrix<double, Tdim, 1>::Zero();

  for (unsigned i = 0; i < nodes_.size(); ++i)
    nodal_velocity +=
        shapefn_[i] * nodes_[i]->velocity(mpm::ParticlePhase::Solid);

  // Acceleration update
  if (!velocity_update) {
    // Get interpolated nodal acceleration
    Eigen::Matrix<double, Tdim, 1> nodal_acceleration =
        Eigen::Matrix<double, Tdim, 1>::Zero();
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
}

//! Map particle pressure to nodes
template <unsigned Tdim>
bool mpm::Particle<Tdim>::map_pressure_to_nodes(unsigned phase) noexcept {
  // Mass is initialized
  assert(mass_ != std::numeric_limits<double>::max());

  bool status = false;
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
  return status;
}

// Compute pressure smoothing of the particle based on nodal pressure
template <unsigned Tdim>
bool mpm::Particle<Tdim>::compute_pressure_smoothing(unsigned phase) noexcept {
  // Assert
  assert(cell_ != nullptr);

  bool status = false;
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
  return status;
}

//! Apply particle velocity constraints
template <unsigned Tdim>
void mpm::Particle<Tdim>::apply_particle_velocity_constraints(unsigned dir,
                                                              double velocity) {
  // Set particle velocity constraint
  this->velocity_(dir) = velocity;
}

//! Return particle tensor data
template <unsigned Tdim>
Eigen::VectorXd mpm::Particle<Tdim>::tensor_data(const std::string& property) {
  return this->properties_.at(property)();
}

//! Assign material id of this particle to nodes
template <unsigned Tdim>
void mpm::Particle<Tdim>::append_material_id_to_nodes() const {
  for (unsigned i = 0; i < nodes_.size(); ++i)
    nodes_[i]->append_material_id(this->material_id());
}

//! Assign neighbour particles
template <unsigned Tdim>
void mpm::Particle<Tdim>::assign_neighbours(
    const std::vector<mpm::Index>& neighbours) {
  neighbours_ = neighbours;
  neighbours_.erase(std::remove(neighbours_.begin(), neighbours_.end(), id_),
                    neighbours_.end());
}

//! Compute size of serialized particle data
template <unsigned Tdim>
int mpm::Particle<Tdim>::compute_pack_size() const {
  int total_size = 0;
  int partial_size;

  // ID
  MPI_Pack_size(1, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD, &partial_size);
  total_size += partial_size;
  // mass, volume, pressure
  MPI_Pack_size(3 * 1, MPI_DOUBLE, MPI_COMM_WORLD, &partial_size);
  total_size += partial_size;

  // Coordinates, displacement, natural size, velocity
  MPI_Pack_size(4 * Tdim, MPI_DOUBLE, MPI_COMM_WORLD, &partial_size);
  total_size += partial_size;
  // Stress & strain
  MPI_Pack_size(6 * 2, MPI_DOUBLE, MPI_COMM_WORLD, &partial_size);
  total_size += partial_size;

  // epsv
  MPI_Pack_size(1, MPI_DOUBLE, MPI_COMM_WORLD, &partial_size);
  total_size += partial_size;

  // Cell id
  MPI_Pack_size(1, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD, &partial_size);
  total_size += partial_size;

  // Status
  MPI_Pack_size(1, MPI_C_BOOL, MPI_COMM_WORLD, &partial_size);
  total_size += partial_size;

  // Material id
  MPI_Pack_size(1, MPI_UNSIGNED, MPI_COMM_WORLD, &partial_size);
  total_size += partial_size;

  // nstate variables
  unsigned nstate_vars = state_variables_[mpm::ParticlePhase::Solid].size();
  MPI_Pack_size(1, MPI_UNSIGNED, MPI_COMM_WORLD, &partial_size);
  total_size += partial_size;

  // state variables
  MPI_Pack_size(nstate_vars, MPI_DOUBLE, MPI_COMM_WORLD, &partial_size);
  total_size += partial_size;

  return total_size;
}

//! Serialize particle data
template <unsigned Tdim>
std::vector<uint8_t> mpm::Particle<Tdim>::serialize() {
  // Compute pack size
  if (pack_size_ == 0) pack_size_ = compute_pack_size();
  // Initialize data buffer
  // std::string data(pack_size_, ' ');
  std::vector<uint8_t> data;
  data.resize(pack_size_);
  uint8_t* data_ptr = &data[0];
  int position = 0;

  // ID
  MPI_Pack(&id_, 1, MPI_UNSIGNED_LONG_LONG, data_ptr, data.size(), &position,
           MPI_COMM_WORLD);
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
  MPI_Pack(displacement_.data(), Tdim, MPI_DOUBLE, data_ptr, data.size(),
           &position, MPI_COMM_WORLD);
  // Natural size
  MPI_Pack(natural_size_.data(), Tdim, MPI_DOUBLE, data_ptr, data.size(),
           &position, MPI_COMM_WORLD);
  // Velocity
  MPI_Pack(velocity_.data(), Tdim, MPI_DOUBLE, data_ptr, data.size(), &position,
           MPI_COMM_WORLD);
  // Stress
  MPI_Pack(stress_.data(), 6, MPI_DOUBLE, data_ptr, data.size(), &position,
           MPI_COMM_WORLD);
  // Strain
  MPI_Pack(strain_.data(), 6, MPI_DOUBLE, data_ptr, data.size(), &position,
           MPI_COMM_WORLD);

  // epsv
  MPI_Pack(&volumetric_strain_centroid_, 1, MPI_DOUBLE, data_ptr, data.size(),
           &position, MPI_COMM_WORLD);

  // Cell id
  MPI_Pack(&cell_id_, 1, MPI_UNSIGNED_LONG_LONG, data_ptr, data.size(),
           &position, MPI_COMM_WORLD);

  // Status
  MPI_Pack(&status_, 1, MPI_C_BOOL, data_ptr, data.size(), &position,
           MPI_COMM_WORLD);

  // Material id
  MPI_Pack(&material_id_[0], 1, MPI_UNSIGNED, data_ptr, data.size(), &position,
           MPI_COMM_WORLD);

  // nstate variables
  unsigned nstate_vars = state_variables_[mpm::ParticlePhase::Solid].size();
  MPI_Pack(&nstate_vars, 1, MPI_UNSIGNED, data_ptr, data.size(), &position,
           MPI_COMM_WORLD);

  // state variables
  if (this->material() != nullptr) {
    std::vector<double> svars;
    unsigned i = 0;
    auto state_variables = (this->material())->state_variables();
    for (const auto& state_var : state_variables)
      svars.emplace_back(
          state_variables_[mpm::ParticlePhase::Solid].at(state_var));

    // Write state vars
    MPI_Pack(&svars[0], nstate_vars, MPI_DOUBLE, data_ptr, data.size(),
             &position, MPI_COMM_WORLD);
  }

  return data;
}

//! Deserialize particle data
template <unsigned Tdim>
void mpm::Particle<Tdim>::deserialize(const std::vector<uint8_t>& data) {
  uint8_t* data_ptr = const_cast<uint8_t*>(&data[0]);
  int position = 0;

  // ID
  MPI_Unpack(data_ptr, data.size(), &position, &id_, 1, MPI_UNSIGNED_LONG_LONG,
             MPI_COMM_WORLD);
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

  // Coordinates
  MPI_Unpack(data_ptr, data.size(), &position, coordinates_.data(), Tdim,
             MPI_DOUBLE, MPI_COMM_WORLD);
  // Displacement
  MPI_Unpack(data_ptr, data.size(), &position, displacement_.data(), Tdim,
             MPI_DOUBLE, MPI_COMM_WORLD);
  // Natural size
  MPI_Unpack(data_ptr, data.size(), &position, natural_size_.data(), Tdim,
             MPI_DOUBLE, MPI_COMM_WORLD);
  // Velocity
  MPI_Unpack(data_ptr, data.size(), &position, velocity_.data(), Tdim,
             MPI_DOUBLE, MPI_COMM_WORLD);
  // Stress
  MPI_Unpack(data_ptr, data.size(), &position, stress_.data(), 6, MPI_DOUBLE,
             MPI_COMM_WORLD);
  // Strain
  MPI_Unpack(data_ptr, data.size(), &position, strain_.data(), 6, MPI_DOUBLE,
             MPI_COMM_WORLD);

  // epsv
  MPI_Unpack(data_ptr, data.size(), &position, &volumetric_strain_centroid_, 1,
             MPI_DOUBLE, MPI_COMM_WORLD);
  // cell id
  MPI_Unpack(data_ptr, data.size(), &position, &cell_id_, 1,
             MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);
  // status
  MPI_Unpack(data_ptr, data.size(), &position, &status_, 1, MPI_C_BOOL,
             MPI_COMM_WORLD);
  // material id
  MPI_Unpack(data_ptr, data.size(), &position, &material_id_[0], 1,
             MPI_UNSIGNED, MPI_COMM_WORLD);

  // nstate vars
  unsigned nstate_vars;
  MPI_Unpack(data_ptr, data.size(), &position, &nstate_vars, 1, MPI_UNSIGNED,
             MPI_COMM_WORLD);
}

//! Constructor
template <unsigned Tdim>
mpm::MPMExplicitUSL<Tdim>::MPMExplicitUSL(std::unique_ptr<IO>&& io)
    : mpm::MPMExplicit<Tdim>(std::move(io)) {
  //! Logger
  console_ = spdlog::get("MPMExplicitUSL");
}

//! MPM Explicit USL solver
template <unsigned Tdim>
bool mpm::MPMExplicitUSL<Tdim>::solve() {
  bool status = true;

  // Phase
  const unsigned phase = 0;
  // Initialise material
  bool mat_status = this->initialise_materials();
  if (!mat_status) status = false;

  // Initialise mesh and materials
  bool mesh_status = this->initialise_mesh_particles();
  if (!mesh_status) status = false;

  // Assign material to particles
  // Get mesh properties
  auto mesh_props = io_->json_object("mesh");
  // Material id
  const auto material_id = mesh_props["material_id"].template get<unsigned>();

  // Get material from list of materials
  auto material = materials_.at(material_id);

  // Iterate over each particle to assign material
  meshes_.at(0)->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::assign_material,
                std::placeholders::_1, material));

  // Compute volume
  meshes_.at(0)->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::compute_volume, std::placeholders::_1, phase));

  // Compute mass
  meshes_.at(0)->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::compute_mass, std::placeholders::_1, phase));

  // Test if checkpoint resume is needed
  bool resume = false;
  if (analysis_.find("resume") != analysis_.end())
    resume = analysis_["resume"]["resume"].template get<bool>();
  if (resume) this->checkpoint_resume();

  for (; step_ < nsteps_; ++step_) {
    console_->info("Step: {} of {}.\n", step_, nsteps_);
    // Initialise nodes
    meshes_.at(0)->iterate_over_nodes(
        std::bind(&mpm::NodeBase<Tdim>::initialise, std::placeholders::_1));

    meshes_.at(0)->iterate_over_cells(
        std::bind(&mpm::Cell<Tdim>::activate_nodes, std::placeholders::_1));

    // Iterate over each particle to compute shapefn
    meshes_.at(0)->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::compute_shapefn, std::placeholders::_1));

    // Assign mass and momentum to nodes
    meshes_.at(0)->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_mass_momentum_to_nodes,
                  std::placeholders::_1, phase));

    // Compute nodal velocity
    meshes_.at(0)->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_velocity,
                  std::placeholders::_1),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Iterate over each particle to compute nodal body force
    meshes_.at(0)->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_body_force,
                  std::placeholders::_1, phase, this->gravity_));

    // Iterate over each particle to compute nodal internal force
    meshes_.at(0)->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_internal_force,
                  std::placeholders::_1, phase));

    // Iterate over active nodes to compute acceleratation and velocity
    meshes_.at(0)->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_acceleration_velocity,
                  std::placeholders::_1, phase, this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Iterate over each particle to compute updated position
    meshes_.at(0)->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::compute_updated_position,
                  std::placeholders::_1, phase, this->dt_));

    // Iterate over each particle to calculate strain
    meshes_.at(0)->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::compute_strain,
                  std::placeholders::_1, phase, dt_));

    // Iterate over each particle to compute stress
    meshes_.at(0)->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::compute_stress,
                  std::placeholders::_1, phase));

    // Locate particles
    auto unlocatable_particles = meshes_.at(0)->locate_particles_mesh();

    if (!unlocatable_particles.empty())
      throw std::runtime_error("Particle outside the mesh domain");

    if (step_ % output_steps_ == 0) {
      // HDF5 outputs
      this->write_hdf5(step_, this->nsteps_);
#ifdef USE_VTK
      // VTK outputs
      this->write_vtk(this->step_, this->nsteps_);
#endif
    }
  }
  return status;
}

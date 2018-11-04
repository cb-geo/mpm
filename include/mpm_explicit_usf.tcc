//! Constructor
template <unsigned Tdim>
mpm::MPMExplicitUSF<Tdim>::MPMExplicitUSF(std::unique_ptr<IO>&& io)
    : mpm::MPMExplicit<Tdim>(std::move(io)) {
  //! Logger
  console_ = spdlog::get("MPMExplicitUSF");
}

//! MPM Explicit USF solver
template <unsigned Tdim>
bool mpm::MPMExplicitUSF<Tdim>::solve() {
  bool status = true;

  // Initialise MPI rank and size
  int mpi_rank = 0;
  int mpi_size = 1;

#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  // Get number of MPI ranks
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif

  // Phase
  const unsigned phase = 0;

  // Test if checkpoint resume is needed
  bool resume = false;
  if (analysis_.find("resume") != analysis_.end())
    resume = analysis_["resume"]["resume"].template get<bool>();

  // Initialise material
  bool mat_status = this->initialise_materials();
  if (!mat_status) status = false;

  // Initialise mesh
  bool mesh_status = this->initialise_mesh();
  if (!mesh_status) status = false;

  // Initialise particles
  bool particle_status = this->initialise_particles();
  if (!particle_status) status = false;

  // Assign material to particles
  // Get mesh properties
  auto mesh_props = io_->json_object("mesh");
  // Material id
  const auto material_id = mesh_props["material_id"].template get<unsigned>();

  // Get material from list of materials
  auto material = materials_.at(material_id);

  // Iterate over each particle to assign material
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::assign_material,
                std::placeholders::_1, material));

  // Check point resume
  if (resume) this->checkpoint_resume();

  // Compute mass
  mesh_->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::compute_mass, std::placeholders::_1, phase));

  auto solver_begin = std::chrono::steady_clock::now();
  // Main loop
  for (; step_ < nsteps_; ++step_) {

    if (mpi_rank == 0) console_->info("Step: {} of {}.\n", step_, nsteps_);

    // Create a TBB task group
    tbb::task_group task_group;

    // Spawn a task for initialising nodes and cells
    task_group.run([&] {
      // Initialise nodes
      mesh_->iterate_over_nodes(
          std::bind(&mpm::NodeBase<Tdim>::initialise, std::placeholders::_1));

      mesh_->iterate_over_cells(
          std::bind(&mpm::Cell<Tdim>::activate_nodes, std::placeholders::_1));

      mesh_->find_active_nodes();
    });

    // Spawn a task for particles
    task_group.run([&] {
      // Iterate over each particle to compute shapefn
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::compute_shapefn, std::placeholders::_1));
    });

    task_group.wait();

    // Assign mass and momentum to nodes
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_mass_momentum_to_nodes,
                  std::placeholders::_1, phase));

#ifdef USE_MPI
    // Run if there is more than a single MPI task
    if (mpi_size > 1) {
      // MPI all reduce nodal mass
      mesh_->allreduce_nodal_scalar_property(
          std::bind(&mpm::NodeBase<Tdim>::mass, std::placeholders::_1, phase),
          std::bind(&mpm::NodeBase<Tdim>::update_mass, std::placeholders::_1,
                    false, phase, std::placeholders::_2));
      // MPI all reduce nodal momentum
      mesh_->allreduce_nodal_vector_property(
          std::bind(&mpm::NodeBase<Tdim>::momentum, std::placeholders::_1,
                    phase),
          std::bind(&mpm::NodeBase<Tdim>::update_momentum,
                    std::placeholders::_1, false, phase,
                    std::placeholders::_2));
    }
#endif

    // Compute nodal velocity
    mesh_->iterate_over_active_nodes(std::bind(
        &mpm::NodeBase<Tdim>::compute_velocity, std::placeholders::_1));

    // Iterate over each particle to calculate strain
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::compute_strain,
                  std::placeholders::_1, phase, dt_));

    // Iterate over each particle to compute stress
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::compute_stress,
                  std::placeholders::_1, phase));

    // Spawn a task for external force
    task_group.run([&] {
      // Iterate over each particle to compute nodal body force
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::map_body_force,
                    std::placeholders::_1, phase, this->gravity_));

      // Iterate over each particle to map traction force to nodes
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::map_traction_force,
                    std::placeholders::_1, phase));
    });

    // Spawn a task for internal force
    task_group.run([&] {
      // Iterate over each particle to compute nodal internal force
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::map_internal_force,
                    std::placeholders::_1, phase));
    });
    task_group.wait();

#ifdef USE_MPI
    // Run if there is more than a single MPI task
    if (mpi_size > 1) {
      // MPI all reduce external force
      mesh_->allreduce_nodal_vector_property(
          std::bind(&mpm::NodeBase<Tdim>::external_force, std::placeholders::_1,
                    phase),
          std::bind(&mpm::NodeBase<Tdim>::update_external_force,
                    std::placeholders::_1, false, phase,
                    std::placeholders::_2));
      // MPI all reduce internal force
      mesh_->allreduce_nodal_vector_property(
          std::bind(&mpm::NodeBase<Tdim>::internal_force, std::placeholders::_1,
                    phase),
          std::bind(&mpm::NodeBase<Tdim>::update_internal_force,
                    std::placeholders::_1, false, phase,
                    std::placeholders::_2));
    }
#endif

    // Iterate over active nodes to compute acceleratation and velocity
    mesh_->iterate_over_active_nodes(
        std::bind(&mpm::NodeBase<Tdim>::compute_acceleration_velocity,
                  std::placeholders::_1, phase, this->dt_));

    // Iterate over each particle to compute updated position
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::compute_updated_position,
                  std::placeholders::_1, phase, this->dt_));

    // Iterate over each particle to update particle volume
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::update_volume_strainrate,
                  std::placeholders::_1, phase, this->dt_));

    // Locate particles
    auto unlocatable_particles = mesh_->locate_particles_mesh();

    if (!unlocatable_particles.empty())
      throw std::runtime_error("Particle outside the mesh domain");

    if (step_ % output_steps_ == 0) {
      // HDF5 outputs
      this->write_hdf5(this->step_, this->nsteps_);
#ifdef USE_VTK
      // VTK outputs
      this->write_vtk(this->step_, this->nsteps_);
#endif
    }
  }
  auto solver_end = std::chrono::steady_clock::now();
  console_->info("Rank {}, USF solver duration: {} ms", mpi_rank,
                 std::chrono::duration_cast<std::chrono::milliseconds>(
                     solver_end - solver_begin)
                     .count());

  return status;
}

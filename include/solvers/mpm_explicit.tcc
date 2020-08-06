//! Constructor
template <unsigned Tdim>
mpm::MPMExplicit<Tdim>::MPMExplicit(const std::shared_ptr<IO>& io)
    : mpm::MPMBase<Tdim>(io) {
  //! Logger
  console_ = spdlog::get("MPMExplicit");
  //! Stress update
  if (this->stress_update_ == "usl")
    stress_update_scheme_ =
        std::make_shared<mpm::StressUpdateUSL<Tdim>>(mesh_, dt_);
  else
    stress_update_scheme_ =
        std::make_shared<mpm::StressUpdateUSF<Tdim>>(mesh_, dt_);

  //! Interface scheme
  if (this->interface_)
    interface_scheme_ = std::make_shared<mpm::InterfaceContact<Tdim>>(mesh_);
  else
    interface_scheme_ = std::make_shared<mpm::Interface<Tdim>>(mesh_);
}

//! MPM Explicit compute stress strain
template <unsigned Tdim>
void mpm::MPMExplicit<Tdim>::compute_stress_strain(unsigned phase) {
  // Iterate over each particle to calculate strain
  mesh_->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::compute_strain, std::placeholders::_1, dt_));

  // Iterate over each particle to update particle volume
  mesh_->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::update_volume, std::placeholders::_1));

  // Pressure smoothing
  if (pressure_smoothing_) this->pressure_smoothing(phase);

  // Iterate over each particle to compute stress
  mesh_->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::compute_stress, std::placeholders::_1));
}

//! MPM Explicit solver
template <unsigned Tdim>
bool mpm::MPMExplicit<Tdim>::solve() {
  bool status = true;

  console_->info("MPM analysis type {}", io_->analysis_type());

  // Initialise MPI rank and size
  int mpi_rank = 0;
  int mpi_size = 1;

#ifdef USE_MPI
  // Get MPI rank
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

  // Pressure smoothing
  pressure_smoothing_ = io_->analysis_bool("pressure_smoothing");

  // Interface
  if (analysis_.find("interface") != analysis_.end())
    interface_ = analysis_.at("interface").template get<bool>();

  // Initialise material
  this->initialise_materials();

  // Initialise mesh
  this->initialise_mesh();

  // Initialise particles
  this->initialise_particles();

  // Initialise loading conditions
  this->initialise_loads();

  // Create nodal properties
  if (interface_) mesh_->create_nodal_properties();

  // Compute mass
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::compute_mass, std::placeholders::_1));

  // Check point resume
  if (resume) this->checkpoint_resume();

  // Domain decompose
  bool initial_step = (resume == true) ? false : true;
  this->mpi_domain_decompose(initial_step);

  auto solver_begin = std::chrono::steady_clock::now();
  // Main loop
  for (; step_ < nsteps_; ++step_) {

    if (mpi_rank == 0) console_->info("Step: {} of {}.\n", step_, nsteps_);

#ifdef USE_MPI
#ifdef USE_GRAPH_PARTITIONING
    // Run load balancer at a specified frequency
    if (step_ % nload_balance_steps_ == 0 && step_ != 0)
      this->mpi_domain_decompose(false);
#endif
#endif

    // Inject particles
    mesh_->inject_particles(step_ * dt_);

    // Initialise nodes, cells and shape functions
    stress_update_scheme_->initialise();

    // Initialise nodal properties and append material ids to node
    interface_scheme_->initialise();

    // Mass momentum and compute velocity at nodes
    stress_update_scheme_->momentum_nodes(phase);

    // Map material properties to nodes
    interface_scheme_->compute_contact_forces();

    // Update stress first
    stress_update_scheme_->precompute_stress_strain(phase, pressure_smoothing_);

    // Compute forces
    stress_update_scheme_->compute_forces(gravity_, phase, step_,
                                          set_node_concentrated_force_);

    // Particle kinematics
    stress_update_scheme_->compute_particle_kinematics(
        velocity_update_, phase, "Cundall", damping_factor_);

    // Update Stress Last
    stress_update_scheme_->postcompute_stress_strain(phase,
                                                     pressure_smoothing_);

    // Locate particles
    stress_update_scheme_->locate_particles(this->locate_particles_);

#ifdef USE_MPI
#ifdef USE_GRAPH_PARTITIONING
    mesh_->transfer_halo_particles();
#endif
#endif

    if (step_ % output_steps_ == 0) {
      // HDF5 outputs
      this->write_hdf5(this->step_, this->nsteps_);
#ifdef USE_VTK
      // VTK outputs
      this->write_vtk(this->step_, this->nsteps_);
#endif
#ifdef USE_PARTIO
      // Partio outputs
      this->write_partio(this->step_, this->nsteps_);
#endif
    }
  }
  auto solver_end = std::chrono::steady_clock::now();
  console_->info("Rank {}, Explicit {} solver duration: {} ms", mpi_rank,
                 stress_update_scheme_->scheme(),
                 std::chrono::duration_cast<std::chrono::milliseconds>(
                     solver_end - solver_begin)
                     .count());

  return status;
}

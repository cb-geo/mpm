//! Constructor
template <unsigned Tdim>
mpm::MPMExplicit<Tdim>::MPMExplicit(const std::shared_ptr<IO>& io)
    : mpm::MPMBase<Tdim>(io) {
  //! Logger
  console_ = spdlog::get("MPMExplicit");
  //! Stress update
  if (this->stress_update_ == "usl")
    mpm_scheme_ = std::make_shared<mpm::MPMSchemeUSL<Tdim>>(mesh_, dt_);
  else
    mpm_scheme_ = std::make_shared<mpm::MPMSchemeUSF<Tdim>>(mesh_, dt_);
}

//! MPM Explicit compute stress strain
template <unsigned Tdim>
void mpm::MPMExplicit<Tdim>::compute_stress_strain(unsigned phase) {
  // Iterate over each particle to calculate strain
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::compute_strain, std::placeholders::_1,
                dt_, interface_));

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
    interface_ = analysis_["interface"]["interface"].template get<bool>();

  // Interface scheme
  if (interface_) {
    double friction = analysis_["interface"]["friction"].template get<double>();
    std::string normal_type =
        analysis_["interface"]["normal_type"].template get<std::string>();
    contact_ = std::make_shared<mpm::ContactFriction<Tdim>>(mesh_, friction,
                                                            normal_type);
  } else {
    contact_ = std::make_shared<mpm::Contact<Tdim>>(mesh_);
  }

  // Initialise material
  this->initialise_materials();

  // Initialise mesh
  this->initialise_mesh();

  // Initialise particles
  if (!resume) this->initialise_particles();

  // Create nodal properties
  if (interface_) mesh_->create_nodal_properties();

  // Compute mass
  if (!resume)
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::compute_mass, std::placeholders::_1));

  // Check point resume
  if (resume) {
    this->checkpoint_resume();
    mesh_->resume_domain_cell_ranks();
#ifdef USE_MPI
#ifdef USE_GRAPH_PARTITIONING
    MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
  } else {
    // Domain decompose
    bool initial_step = (resume == true) ? false : true;
    this->mpi_domain_decompose(initial_step);
  }

  //! Particle entity sets and velocity constraints
  if (resume) {
    this->particle_entity_sets(false);
    this->particle_velocity_constraints();
  }

  // Initialise loading conditions
  this->initialise_loads();

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
    mpm_scheme_->initialise();

    // Initialise nodal properties and append material ids to node
    contact_->initialise();

    // Mass momentum and compute velocity at nodes
    mpm_scheme_->compute_nodal_kinematics(phase);

    // Compute mass, momentum and velocity at the contact nodes
    contact_->compute_nodal_kinematics();

    // Update stress first
    mpm_scheme_->precompute_stress_strain(phase, pressure_smoothing_,
                                          interface_);

    // Compute forces
    mpm_scheme_->compute_forces(gravity_, phase, step_,
                                set_node_concentrated_force_);

    // Compute internal and external forces at contact nodes
    contact_->compute_contact_forces(gravity_, phase, (step_ * dt_),
                                     set_node_concentrated_force_);

    // Particle kinematics
    mpm_scheme_->compute_particle_kinematics(phase, "Cundall", damping_factor_);

    // Nodal kinematics at contact nodes
    contact_->compute_contact_kinematics(dt_);

    // Compute particle updated position
    if (interface_) {
      contact_->update_particles_contact(dt_, velocity_update_);
    } else {
      mpm_scheme_->update_particles(velocity_update_);
    }

    // Update Stress Last
    mpm_scheme_->postcompute_stress_strain(phase, pressure_smoothing_,
                                           interface_);

    // Locate particles
    mpm_scheme_->locate_particles(this->locate_particles_);

#ifdef USE_MPI
#ifdef USE_GRAPH_PARTITIONING
    mesh_->transfer_halo_particles();
    MPI_Barrier(MPI_COMM_WORLD);
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
                 mpm_scheme_->scheme(),
                 std::chrono::duration_cast<std::chrono::milliseconds>(
                     solver_end - solver_begin)
                     .count());

  return status;
}

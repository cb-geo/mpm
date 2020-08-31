//! Constructor
template <unsigned Tdim>
mpm::XMPMExplicit<Tdim>::XMPMExplicit(const std::shared_ptr<IO>& io)
    : mpm::MPMBase<Tdim>(io) {
  //! Logger
  console_ = spdlog::get("XMPMExplicit");
  //! Stress update
  if (this->stress_update_ == "usl")
    mpm_scheme_ = std::make_shared<mpm::MPMSchemeUSL<Tdim>>(mesh_, dt_);
  else
    mpm_scheme_ = std::make_shared<mpm::MPMSchemeUSF<Tdim>>(mesh_, dt_);

  //! Interface scheme
  if (this->interface_)
    contact_ = std::make_shared<mpm::ContactFriction<Tdim>>(mesh_);
  else
    contact_ = std::make_shared<mpm::Contact<Tdim>>(mesh_);
}

//! MPM Explicit compute stress strain
template <unsigned Tdim>
void mpm::XMPMExplicit<Tdim>::compute_stress_strain(unsigned phase) {
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
bool mpm::XMPMExplicit<Tdim>::solve() {
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
  interface_ = io_->analysis_bool("interface");

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

  // Initialise discontinuity
  this->initialise_discontinuities();

  // Initialise the levelset values for particles
  if (discontinuity_) mesh_->initialise_levelset_discontinuity();

  // Create nodal properties for discontinuity
  if (discontinuity_) mesh_->create_nodal_properties_discontinuity();
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
    mpm_scheme_->initialise();

    // Initialise nodal properties and append material ids to node
    contact_->initialise();

    if (discontinuity_) {
      // Initialise nodal properties
      mesh_->initialise_nodal_properties();

      // locate points of discontinuity
      mesh_->locate_discontinuity();

      // Iterate over each points to compute shapefn
      mesh_->compute_shapefn_discontinuity();

      // obtain the normal direction of each enrich nodes
      mesh_->compute_normal_vector_discontinuity();
    }

    // Mass momentum and compute velocity at nodes
    mpm_scheme_->compute_nodal_kinematics(phase);

    // Map material properties to nodes
    contact_->compute_contact_forces();

    // Update stress first
    mpm_scheme_->precompute_stress_strain(phase, pressure_smoothing_);

    // Compute forces
    mpm_scheme_->compute_forces(gravity_, phase, step_,
                                set_node_concentrated_force_);

    // intergrate momentum Iterate over
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::intergrate_momentum_discontinuity,
                  std::placeholders::_1, phase, this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Update the discontinuity position
    if (discontinuity_)
      mesh_->compute_updated_position_discontinuity(this->dt_);

    // // Particle kinematics
    // mpm_scheme_->compute_particle_kinematics(velocity_update_, phase, "Cundall",
    //                                          damping_factor_);
    // Iterate over each particle to compute updated position
    mesh_->iterate_over_particles(
    std::bind(&mpm::ParticleBase<Tdim>::compute_updated_position,
              std::placeholders::_1, dt_, velocity_update_));

    // Update Stress Last
    mpm_scheme_->postcompute_stress_strain(phase, pressure_smoothing_);

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

// Initialise discontinuities
template <unsigned Tdim>
void mpm::XMPMExplicit<Tdim>::initialise_discontinuities() {
  try {
    // Get discontinuities data
    auto json_discontinuities = io_->json_object("discontinuity");
    if (!json_discontinuities.empty()) {
      discontinuity_ = true;
      for (const auto discontinuity_props : json_discontinuities) {
        // Get discontinuity type
        const std::string discontunity_type =
            discontinuity_props["type"].template get<std::string>();

        // Get discontinuity id
        auto discontinuity_id =
            discontinuity_props["id"].template get<unsigned>();

        // Create a new discontinuity surface from JSON object
        auto discontinuity =
            Factory<mpm::DiscontinuityBase<Tdim>, unsigned,
                    const Json&>::instance()
                ->create(discontunity_type, std::move(discontinuity_id),
                         discontinuity_props);

        // Get discontinuity  input type
        auto io_type =
            discontinuity_props["io_type"].template get<std::string>();

        // discontinuity file
        std::string discontinuity_file = io_->file_name(
            discontinuity_props["file"].template get<std::string>());
        // Create a mesh reader
        auto discontunity_io =
            Factory<mpm::IOMesh<Tdim>>::instance()->create(io_type);

        // Create points and cells from file
        discontinuity->initialize(
            discontunity_io->read_mesh_nodes(discontinuity_file),
            discontunity_io->read_mesh_cells(discontinuity_file));
        // Add discontinuity to list
        auto result = discontinuities_.insert(
            std::make_pair(discontinuity_id, discontinuity));

        // If insert discontinuity failed
        if (!result.second) {
          throw std::runtime_error(
              "New discontinuity cannot be added, insertion failed");
        }
      }
      // Copy discontinuities to mesh
      mesh_->initialise_discontinuities(this->discontinuities_);
    }
  } catch (std::exception& exception) {
    console_->warn("{} #{}: No discontinuity is defined", __FILE__, __LINE__,
                   exception.what());
  }
}

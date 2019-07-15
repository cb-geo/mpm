//! Constructor
template <unsigned Tdim>
mpm::MPMExplicitTwoPhase<Tdim>::MPMExplicitTwoPhase(std::unique_ptr<IO>&& io)
    : mpm::MPMBase<Tdim>(std::move(io)) {
  //! Logger
  console_ = spdlog::get("MPMExplicitTwoPhase");
}

//! MPM Explicit solver
template <unsigned Tdim>
bool mpm::MPMExplicitTwoPhase<Tdim>::solve() {
  bool status = true;

  console_->error("Analysis{} {}", io_->analysis_type());

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
  const unsigned solid_skeleton = 0;
  const unsigned pore_fluid = 1;

  // Test if checkpoint resume is needed
  bool resume = false;
  if (analysis_.find("resume") != analysis_.end())
    resume = analysis_["resume"]["resume"].template get<bool>();

  // Pressure smoothing
  if (analysis_.find("pressure_smoothing") != analysis_.end())
    pressure_smoothing_ = analysis_["pressure_smoothing"].template get<bool>();

  // Initialise material
  bool mat_status = this->initialise_materials();
  if (!mat_status) status = false;

  // Initialise mesh
  bool mesh_status = this->initialise_mesh();
  if (!mesh_status) status = false;

  // Initialise particles
  bool particle_status = this->initialise_particles();
  if (!particle_status) status = false;

  // Assign material to particles for each phases
  // Get particle properties
  auto particle_props = io_->json_object("particle");

  // Get material ids from input file
  if (!particle_props.at("material_id").is_array() ||
      particle_props.at("material_id").size() != 2)
    throw std::runtime_error("Unable to assign material for each phase");

  // Initialise material ids for each phase
  const auto solid_skeleton_mid =
      particle_props["material_id"][0].template get<unsigned>();
  const auto pore_fluid_mid =
      particle_props["material_id"][1].template get<unsigned>();

  // Get material from list of materials for each phase
  auto solid_skeleton_material = materials_.at(solid_skeleton_mid);
  auto pore_fluid_material = materials_.at(pore_fluid_mid);

  // Iterate over each particle to assign material to each phase
  mesh_->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::assign_material, std::placeholders::_1,
      solid_skeleton, solid_skeleton_material));
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::assign_material,
                std::placeholders::_1, pore_fluid, pore_fluid_material));

  // Assign material to particle sets
  if (particle_props["particle_sets"].size() != 0) {
    // Assign material to particles in the specific sets
    bool set_material_status = this->apply_properties_to_particles_sets();
  }

  // Read and assign porosity (volume fraction and phase volume)
  const double porosity = particle_props["porosity"].template get<double>();
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::assign_porosity,
                std::placeholders::_1, porosity));

  // Compute mass for each phase
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::compute_mass, std::placeholders::_1,
                solid_skeleton));
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::compute_mass, std::placeholders::_1,
                pore_fluid));

  // Check point resume
  if (resume) this->checkpoint_resume();

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

      // mesh_->find_active_nodes();
    });

    // Spawn a task for particles
    task_group.run([&] {
      // Iterate over each particle to compute shapefn
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::compute_shapefn, std::placeholders::_1));
    });

    task_group.wait();

    // Assign mass and momentum to nodes
    // Solid phase
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_mass_momentum_to_nodes,
                  std::placeholders::_1, solid_skeleton));
    // Fluid phase
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_mass_momentum_to_nodes,
                  std::placeholders::_1, pore_fluid));

#ifdef USE_MPI
    // Run if there is more than a single MPI task
    if (mpi_size > 1) {
      // MPI all reduce nodal mass
      // Solid phase
      mesh_->allreduce_nodal_scalar_property(
          std::bind(&mpm::NodeBase<Tdim>::mass, std::placeholders::_1,
                    solid_skeleton),
          std::bind(&mpm::NodeBase<Tdim>::update_mass, std::placeholders::_1,
                    false, solid_skeleton, std::placeholders::_2));
      // Fluid phase
      mesh_->allreduce_nodal_scalar_property(
          std::bind(&mpm::NodeBase<Tdim>::mass, std::placeholders::_1,
                    pore_fluid),
          std::bind(&mpm::NodeBase<Tdim>::update_mass, std::placeholders::_1,
                    false, pore_fluid, std::placeholders::_2));

      // MPI all reduce nodal momentum
      // Solid phase
      mesh_->allreduce_nodal_vector_property(
          std::bind(&mpm::NodeBase<Tdim>::momentum, std::placeholders::_1,
                    solid_skeleton),
          std::bind(&mpm::NodeBase<Tdim>::update_momentum,
                    std::placeholders::_1, false, solid_skeleton,
                    std::placeholders::_2));
      // Fluid phase
      mesh_->allreduce_nodal_vector_property(
          std::bind(&mpm::NodeBase<Tdim>::momentum, std::placeholders::_1,
                    pore_fluid),
          std::bind(&mpm::NodeBase<Tdim>::update_momentum,
                    std::placeholders::_1, false, pore_fluid,
                    std::placeholders::_2));
    }
#endif

    // Compute nodal velocity
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_velocity,
                  std::placeholders::_1),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Update stress first
    if (stress_update_ == mpm::StressUpdate::usf) {
      // Iterate over each particle to calculate strain of solid_skeleton
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_strain,
                    std::placeholders::_1, solid_skeleton, dt_));

      // Iterate over each particle to calculate strain of solid_skeleton
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_strain,
                    std::placeholders::_1, pore_fluid, dt_));

      // Iterate over each particle to update particle volume
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::update_volume,
                    std::placeholders::_1, solid_skeleton, dt_));

      // Iterate over each particle to update porosity
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::update_porosity,
                    std::placeholders::_1, solid_skeleton, dt_));

      // Pressure smoothing
      if (pressure_smoothing_) {
        // Assign pressure to nodes
        mesh_->iterate_over_particles(
            std::bind(&mpm::ParticleBase<Tdim>::map_pressure_to_nodes,
                      std::placeholders::_1, solid_skeleton));

#ifdef USE_MPI
        // Run if there is more than a single MPI task
        if (mpi_size > 1) {
          // MPI all reduce nodal pressure
          mesh_->allreduce_nodal_scalar_property(
              std::bind(&mpm::NodeBase<Tdim>::pressure, std::placeholders::_1,
                        solid_skeleton),
              std::bind(&mpm::NodeBase<Tdim>::assign_pressure,
                        std::placeholders::_1, solid_skeleton,
                        std::placeholders::_2));
        }
#endif

        // Smooth pressure over particles
        mesh_->iterate_over_particles(
            std::bind(&mpm::ParticleBase<Tdim>::compute_pressure_smoothing,
                      std::placeholders::_1, solid_skeleton));
      }

      // Iterate over each particle to compute stress of solid skeleton
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_stress,
                    std::placeholders::_1, solid_skeleton));

      // Iterate over each particle to compute stress of pore fluid
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_pore_pressure,
                    std::placeholders::_1, solid_skeleton, pore_fluid));
    }

    // Spawn a task for external force
    task_group.run([&] {
      // Iterate over each particle to compute nodal body force of soild
      // skeleton
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::map_body_force,
                    std::placeholders::_1, solid_skeleton, this->gravity_));
      // Iterate over each particle to compute nodal body force of pore fluid
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::map_body_force,
                    std::placeholders::_1, pore_fluid, this->gravity_));

      // Iterate over each particle to map traction force of solid skeleton to
      // nodes
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::map_traction_force,
                    std::placeholders::_1, solid_skeleton));

      // Iterate over each particle to map traction force of pore fluid to nodes
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::map_traction_force,
                    std::placeholders::_1, pore_fluid));

      //! Apply nodal tractions
      if (nodal_tractions_) this->apply_nodal_tractions();
    });

    // Spawn a task for internal force
    task_group.run([&] {
      // Iterate over each particle to compute mixture nodal internal
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::map_mixture_internal_force,
                    std::placeholders::_1));

      // Iterate over each particle to compute nodal internal force of pore
      // fluid
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::map_internal_force,
                    std::placeholders::_1, pore_fluid));

      // Iterate over each particle to compute nodal drag force
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::map_drag_force_coefficient,
                    std::placeholders::_1, pore_fluid));
    });
    task_group.wait();

#ifdef USE_MPI
    // Run if there is more than a single MPI task
    if (mpi_size > 1) {
      // MPI all reduce external force of solid skeleton
      mesh_->allreduce_nodal_vector_property(
          std::bind(&mpm::NodeBase<Tdim>::external_force, std::placeholders::_1,
                    solid_skeleton),
          std::bind(&mpm::NodeBase<Tdim>::update_external_force,
                    std::placeholders::_1, false, solid_skeleton,
                    std::placeholders::_2));
      // MPI all reduce external force of pore fluid
      mesh_->allreduce_nodal_vector_property(
          std::bind(&mpm::NodeBase<Tdim>::external_force, std::placeholders::_1,
                    pore_fluid),
          std::bind(&mpm::NodeBase<Tdim>::update_external_force,
                    std::placeholders::_1, false, pore_fluid,
                    std::placeholders::_2));

      // MPI all reduce internal force of solid skeleton
      mesh_->allreduce_nodal_vector_property(
          std::bind(&mpm::NodeBase<Tdim>::mixture_internal_force,
                    std::placeholders::_1),
          std::bind(&mpm::NodeBase<Tdim>::update_mixture_internal_force,
                    std::placeholders::_1, false, std::placeholders::_2));
      // MPI all reduce internal force of pore fluid
      mesh_->allreduce_nodal_vector_property(
          std::bind(&mpm::NodeBase<Tdim>::internal_force, std::placeholders::_1,
                    pore_fluid),
          std::bind(&mpm::NodeBase<Tdim>::update_internal_force,
                    std::placeholders::_1, false, pore_fluid,
                    std::placeholders::_2));

      // MPI all reduce drag force
      mesh_->allreduce_nodal_vector_property(
          std::bind(&mpm::NodeBase<Tdim>::drag_force_coefficient,
                    std::placeholders::_1),
          std::bind(&mpm::NodeBase<Tdim>::update_drag_force_coefficient,
                    std::placeholders::_1, false, std::placeholders::_2));
    }
#endif

    // Iterate over active nodes to compute acceleratation and velocity
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_acceleration_velocity_two_phase,
                  std::placeholders::_1, solid_skeleton, pore_fluid, this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Use nodal velocity to update position
    if (velocity_update_)
      // Iterate over each particle to compute updated position and velocity
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::compute_updated_position_velocity_two_phase,
          std::placeholders::_1, solid_skeleton, pore_fluid, this->dt_));

    else
      // Iterate over each particle to compute updated position and velocity
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::compute_updated_position_two_phase,
          std::placeholders::_1, solid_skeleton, pore_fluid, this->dt_));

    // Update Stress Last
    if (stress_update_ == mpm::StressUpdate::usl) {
      // Iterate over each particle to calculate strain of solid_skeleton
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_strain,
                    std::placeholders::_1, solid_skeleton, dt_));

      // Iterate over each particle to calculate strain of solid_skeleton
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_strain,
                    std::placeholders::_1, pore_fluid, dt_));

      // Iterate over each particle to update particle volume of solid_skeleton
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::update_volume,
                    std::placeholders::_1, solid_skeleton, dt_));

      // Iterate over each particle to update porosity
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::update_porosity,
                    std::placeholders::_1, solid_skeleton, dt_));

      // Pressure smoothing
      if (pressure_smoothing_) {
        // Assign pressure to nodes
        mesh_->iterate_over_particles(
            std::bind(&mpm::ParticleBase<Tdim>::map_pressure_to_nodes,
                      std::placeholders::_1, solid_skeleton));

#ifdef USE_MPI
        // Run if there is more than a single MPI task
        if (mpi_size > 1) {
          // MPI all reduce nodal pressure
          mesh_->allreduce_nodal_scalar_property(
              std::bind(&mpm::NodeBase<Tdim>::pressure, std::placeholders::_1,
                        solid_skeleton),
              std::bind(&mpm::NodeBase<Tdim>::assign_pressure,
                        std::placeholders::_1, solid_skeleton,
                        std::placeholders::_2));
        }
#endif

        // Smooth pressure over particles
        mesh_->iterate_over_particles(
            std::bind(&mpm::ParticleBase<Tdim>::compute_pressure_smoothing,
                      std::placeholders::_1, solid_skeleton));
      }

      // Iterate over each particle to compute stress
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_stress,
                    std::placeholders::_1, solid_skeleton));

      // Iterate over each particle to compute stress of pore fluid
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_pore_pressure,
                    std::placeholders::_1, solid_skeleton, pore_fluid));
    }

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
  console_->info("Rank {}, Explicit {} solver duration: {} ms", mpi_rank,
                 (stress_update_ == mpm::StressUpdate::usl ? "USL" : "USF"),
                 std::chrono::duration_cast<std::chrono::milliseconds>(
                     solver_end - solver_begin)
                     .count());

  return status;
}

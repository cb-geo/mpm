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

  // Two phases (soil skeleton and pore liquid)
  const unsigned soil_skeleton = 0;
  const unsigned pore_liquid = 1;
  const unsigned mixture = 0;

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
  const auto pore_liquid_mid =
      particle_props["material_id"][1].template get<unsigned>();

  // Get material from list of materials for each phase
  auto solid_skeleton_material = materials_.at(solid_skeleton_mid);
  auto pore_liquid_material = materials_.at(pore_liquid_mid);

  // Iterate over each particle to assign material to each phase
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::assign_material,
                std::placeholders::_1, solid_skeleton_material));
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::assign_liquid_material,
                std::placeholders::_1, pore_liquid_material));

  // Assign material to particle sets
  if (particle_props["particle_sets"].size() != 0) {
    // Assign material to particles in the specific sets
    bool set_material_status = this->apply_properties_to_particles_sets();
  }

  // Assign porosity
  mesh_->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::assign_porosity, std::placeholders::_1));

  // Compute mass for each phase
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::compute_mass, std::placeholders::_1));
  mesh_->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::compute_liquid_mass, std::placeholders::_1));

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
    // Soil skeleton
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_mass_momentum_to_nodes,
                  std::placeholders::_1));
    // Liquid phase
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_liquid_mass_momentum_to_nodes,
                  std::placeholders::_1));

#ifdef USE_MPI
    // Run if there is more than a single MPI task
    if (mpi_size > 1) {
      // MPI all reduce nodal mass
      // Soil skeleton
      mesh_->allreduce_nodal_scalar_property(
          std::bind(&mpm::NodeBase<Tdim>::mass, std::placeholders::_1,
                    soil_skeleton),
          std::bind(&mpm::NodeBase<Tdim>::update_mass, std::placeholders::_1,
                    false, soil_skeleton, std::placeholders::_2));
      // Liquid phase
      mesh_->allreduce_nodal_scalar_property(
          std::bind(&mpm::NodeBase<Tdim>::mass, std::placeholders::_1,
                    pore_liquid),
          std::bind(&mpm::NodeBase<Tdim>::update_mass, std::placeholders::_1,
                    false, pore_liquid, std::placeholders::_2));

      // MPI all reduce nodal momentum
      // Soil skeleton
      mesh_->allreduce_nodal_vector_property(
          std::bind(&mpm::NodeBase<Tdim>::momentum, std::placeholders::_1,
                    soil_skeleton),
          std::bind(&mpm::NodeBase<Tdim>::update_momentum,
                    std::placeholders::_1, false, soil_skeleton,
                    std::placeholders::_2));
      // Fluid phase
      mesh_->allreduce_nodal_vector_property(
          std::bind(&mpm::NodeBase<Tdim>::momentum, std::placeholders::_1,
                    pore_liquid),
          std::bind(&mpm::NodeBase<Tdim>::update_momentum,
                    std::placeholders::_1, false, pore_liquid,
                    std::placeholders::_2));
    }
#endif

    // Compute nodal velocity at the begining of time step
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_velocity,
                  std::placeholders::_1),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Update stress first
    if (this->stress_update_ == mpm::StressUpdate::USF) {
      // Iterate over each particle to calculate strain of soil_skeleton
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_strain,
                    std::placeholders::_1, dt_));

      // Iterate over each particle to update particle volume
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::update_volume_strainrate,
                    std::placeholders::_1, dt_));

      // Iterate over each particle to update porosity
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::update_porosity,
                    std::placeholders::_1, dt_));

      // Iterate over each particle to compute stress of soil skeleton
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::compute_stress, std::placeholders::_1));

      // Iterate over each particle to compute pore pressure
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_pore_pressure,
                    std::placeholders::_1, dt_));

      // Pressure smoothing for liquid phase
      if (pressure_smoothing_) {
        // Assign pressure to nodes
        mesh_->iterate_over_particles(
            std::bind(&mpm::ParticleBase<Tdim>::map_pore_pressure_to_nodes,
                      std::placeholders::_1));

#ifdef USE_MPI
        // Run if there is more than a single MPI task
        if (mpi_size > 1) {
          // MPI all reduce nodal pressure
          mesh_->allreduce_nodal_scalar_property(
              std::bind(&mpm::NodeBase<Tdim>::pressure, std::placeholders::_1,
                        pore_liquid),
              std::bind(&mpm::NodeBase<Tdim>::assign_pressure,
                        std::placeholders::_1, pore_liquid,
                        std::placeholders::_2));
        }
#endif

        // Smooth pressure over particles
        mesh_->iterate_over_particles(
            std::bind(&mpm::ParticleBase<Tdim>::compute_pore_pressure_smoothing,
                      std::placeholders::_1));
      }
    }

    // Spawn a task for external force
    task_group.run([&] {
      // Iterate over particles to compute nodal body force of soil skeleton
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::map_mixture_body_force,
                    std::placeholders::_1, mixture, this->gravity_));
      // Iterate over particles to compute nodal body force of pore liquid
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::map_liquid_body_force,
                    std::placeholders::_1, this->gravity_));

      // Iterate over particles to compute nodal traction force of mixture
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::map_mixture_traction_force,
                    std::placeholders::_1, mixture));

      // Iterate over particles to compute nodal traction force of liquid phase
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::map_liquid_traction_force,
                    std::placeholders::_1));
    });

    // Spawn a task for internal force
    task_group.run([&] {
      // Iterate over particles to compute nodal mixture internal force
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::map_mixture_internal_force,
                    std::placeholders::_1, mixture));

      // Iterate over particles to compute nodal liquid phase internal force
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::map_liquid_internal_force,
                    std::placeholders::_1));

      // Iterate over particles to compute nodal drag force coefficient
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::map_drag_force_coefficient,
                    std::placeholders::_1));
    });
    task_group.wait();

#ifdef USE_MPI
    // Run if there is more than a single MPI task
    if (mpi_size > 1) {
      // MPI all reduce external force of mixture
      mesh_->allreduce_nodal_vector_property(
          std::bind(&mpm::NodeBase<Tdim>::external_force, std::placeholders::_1,
                    mixture),
          std::bind(&mpm::NodeBase<Tdim>::update_external_force,
                    std::placeholders::_1, false, mixture,
                    std::placeholders::_2));
      // MPI all reduce external force of pore fluid
      mesh_->allreduce_nodal_vector_property(
          std::bind(&mpm::NodeBase<Tdim>::external_force, std::placeholders::_1,
                    pore_liquid),
          std::bind(&mpm::NodeBase<Tdim>::update_external_force,
                    std::placeholders::_1, false, pore_liquid,
                    std::placeholders::_2));

      // MPI all reduce internal force of mixture
      mesh_->allreduce_nodal_vector_property(
          std::bind(&mpm::NodeBase<Tdim>::internal_force, std::placeholders::_1,
                    mixture),
          std::bind(&mpm::NodeBase<Tdim>::update_internal_force,
                    std::placeholders::_1, false, mixture,
                    std::placeholders::_2));
      // MPI all reduce internal force of pore liquid
      mesh_->allreduce_nodal_vector_property(
          std::bind(&mpm::NodeBase<Tdim>::internal_force, std::placeholders::_1,
                    pore_liquid),
          std::bind(&mpm::NodeBase<Tdim>::update_internal_force,
                    std::placeholders::_1, false, pore_liquid,
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
                  std::placeholders::_1, soil_skeleton, pore_liquid, mixture,
                  this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // TODO
    /*
    // Use nodal velocity to update particle velocity
    if (velocity_update_) {
      // Soil skeleton
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_updated_position_velocity,
                    std::placeholders::_1, this->dt_));
      // Pore liquid
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_updated_liquid_velocity,
                    std::placeholders::_1, this->dt_));
    }
    // Use nodal acceleration to update particle velocity
    else {
      // Solid skeleton
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_updated_position,
                    std::placeholders::_1, this->dt_));
      // Pore fluid
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_updated_liquid_kinematics,
                    std::placeholders::_1, this->dt_));
    }
    */

    // Update Stress Last
    if (this->stress_update_ == mpm::StressUpdate::USL) {
      // Iterate over each particle to calculate strain of soil_skeleton
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_strain,
                    std::placeholders::_1, dt_));

      // Iterate over each particle to update particle volume
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::update_volume_strainrate,
                    std::placeholders::_1, dt_));

      // Iterate over each particle to update porosity
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::update_porosity,
                    std::placeholders::_1, dt_));

      // Iterate over each particle to compute stress of soil skeleton
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::compute_stress, std::placeholders::_1));

      // Iterate over each particle to compute pore pressure
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_pore_pressure,
                    std::placeholders::_1, dt_));

      // Pressure smoothing for liquid phase
      if (pressure_smoothing_) {
        // Assign pressure to nodes
        mesh_->iterate_over_particles(
            std::bind(&mpm::ParticleBase<Tdim>::map_pore_pressure_to_nodes,
                      std::placeholders::_1));

#ifdef USE_MPI
        // Run if there is more than a single MPI task
        if (mpi_size > 1) {
          // MPI all reduce nodal pressure
          mesh_->allreduce_nodal_scalar_property(
              std::bind(&mpm::NodeBase<Tdim>::pressure, std::placeholders::_1,
                        pore_liquid),
              std::bind(&mpm::NodeBase<Tdim>::assign_pressure,
                        std::placeholders::_1, pore_liquid,
                        std::placeholders::_2));
        }
#endif

        // Smooth pressure over particles
        mesh_->iterate_over_particles(
            std::bind(&mpm::ParticleBase<Tdim>::compute_pore_pressure_smoothing,
                      std::placeholders::_1));
      }
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
  console_->info(
      "Rank {}, Explicit {} solver duration: {} ms", mpi_rank,
      (this->stress_update_ == mpm::StressUpdate::USL ? "USL" : "USF"),
      std::chrono::duration_cast<std::chrono::milliseconds>(solver_end -
                                                            solver_begin)
          .count());

  return status;
}

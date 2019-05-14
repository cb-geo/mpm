//! Constructor
template <unsigned Tdim>
mpm::MPMExplicit<Tdim>::MPMExplicit(std::unique_ptr<IO>&& io)
    : mpm::MPMBase<Tdim>(std::move(io)) {
  //! Logger
  console_ = spdlog::get("MPMExplicit");
}

//! MPM Explicit solver
template <unsigned Tdim>
bool mpm::MPMExplicit<Tdim>::solve() {
  bool status = true;

  // Get analysis type USL/USF
  if (io_->analysis_type() == "MPMExplicitUSL2D" ||
      io_->analysis_type() == "MPMExplicitUSL3D")
    this->usl_ = true;

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
  const unsigned phase = 0;

  // Test if checkpoint resume is needed
  bool resume = false;
  if (analysis_.find("resume") != analysis_.end())
    resume = analysis_["resume"]["resume"].template get<bool>();

  // Test if contact computation is needed
  bool contact = false;
  if (analysis_.find("contact") != analysis_.end())
    contact = analysis_["contact"]["contact"].template get<bool>();

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

  // Assign material to particles
  // Get particle properties
  auto particle_props = io_->json_object("particle");
  // Material id
  const auto material_id =
      particle_props["material_id"].template get<unsigned>();

  // Get material from list of materials
  auto material = materials_.at(material_id);

  // Iterate over each particle to assign material
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::assign_material,
                std::placeholders::_1, material));

  // Assign material to particle sets
  if (particle_props["particle_sets"].size() != 0) {
    // Assign material to particles in the specific sets
    bool set_material_status = this->apply_properties_to_particles_sets();
  }

  // Compute mass
  mesh_->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::compute_mass, std::placeholders::_1, phase));

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

      // Initialise nodes for each subdomain for contact computation
      if (contact) {
        auto subdomains = io_->analysis()["contact"]["subdomains"];
        for (const auto& sub : subdomains) {
          unsigned mid = sub["material_id"];
          mesh_->iterate_over_nodes(
              std::bind(&mpm::NodeBase<Tdim>::initialise_subdomain,
                        std::placeholders::_1, mid));
        }
      }

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
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_mass_momentum_to_nodes,
                  std::placeholders::_1, phase));

    // Assign mass and momentum for each subdomain to nodes
    if (contact)
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::map_mass_momentum_to_nodes_subdomain,
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
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_velocity,
                  std::placeholders::_1),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Compute nodal velocity of each subdomain for contact computation
    if (contact) {
      // Map coordinates to nodes from particles
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::map_coordinates_to_nodes,
                    std::placeholders::_1, phase));

      // Map subdomain coordinates to nodes from particles
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::map_coordinates_to_nodes_subdomain,
          std::placeholders::_1, phase));

      // Compute normal vector
      auto subdomains = io_->analysis()["contact"]["subdomains"];
      for (const auto& sub : subdomains) {
        unsigned mid = sub["material_id"];
        // Get default normal vector in "SN" method
        Eigen::Matrix<double, Tdim, 1> normal_vector;
        normal_vector.setZero();
        if (sub["normal_vector_type"] == "SN") {
          for (unsigned i = 0; i < Tdim; ++i) {
            normal_vector[i] = sub["normal_vector"].at(i);
          }
        }
        // Compute normalised normal vector
        mesh_->iterate_over_nodes_predicate(
            std::bind(&mpm::NodeBase<Tdim>::compute_normal_vector,
                      std::placeholders::_1, sub["normal_vector_type"],
                      normal_vector, mid),
            std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));
      }

      // Compute contact nodes ("Correct momentume" or "Implement contact
      // force")
      for (const auto& sub : subdomains) {
        unsigned mid = sub["material_id"];
        mesh_->iterate_over_nodes_predicate(
            std::bind(&mpm::NodeBase<Tdim>::compute_contact_interface,
                      std::placeholders::_1,
                      io_->analysis()["contact"]["contact_force"],
                      sub["friction_type"], sub["friction_coefficient"],
                      io_->analysis()["contact"]["separation_cut_off"],
                      sub["dc_n"], sub["dc_t"], sub["material_id"]),
            std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));
      }

      // Compute nodal velocity
      for (const auto& sub : subdomains) {
        unsigned mid = sub["material_id"];
        mesh_->iterate_over_nodes_predicate(
            std::bind(&mpm::NodeBase<Tdim>::compute_velocity_subdomain,
                      std::placeholders::_1, mid),
            std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));
      }
    }

    // Update stress first
    if (!usl_) {
      // Iterate over each particle to calculate strain
      if (!contact)
        mesh_->iterate_over_particles(
            std::bind(&mpm::ParticleBase<Tdim>::compute_strain,
                      std::placeholders::_1, phase, dt_));

      // Iterate over each particle to calculate strain of each subdomain for
      // contact computation
      if (contact)
        mesh_->iterate_over_particles(
            std::bind(&mpm::ParticleBase<Tdim>::compute_strain_subdomain,
                      std::placeholders::_1, phase, dt_));

      // Iterate over each particle to update particle volume
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::update_volume_strainrate,
                    std::placeholders::_1, phase, this->dt_));

      // Pressure smoothing
      if (pressure_smoothing_) {
        // Assign pressure to nodes
        mesh_->iterate_over_particles(
            std::bind(&mpm::ParticleBase<Tdim>::map_pressure_to_nodes,
                      std::placeholders::_1, phase));

#ifdef USE_MPI
        // Run if there is more than a single MPI task
        if (mpi_size > 1) {
          // MPI all reduce nodal pressure
          mesh_->allreduce_nodal_scalar_property(
              std::bind(&mpm::NodeBase<Tdim>::pressure, std::placeholders::_1,
                        phase),
              std::bind(&mpm::NodeBase<Tdim>::update_pressure,
                        std::placeholders::_1, false, phase,
                        std::placeholders::_2));
        }
#endif

        // Smooth pressure over particles
        mesh_->iterate_over_particles(
            std::bind(&mpm::ParticleBase<Tdim>::compute_pressure_smoothing,
                      std::placeholders::_1, phase));
      }

      // Iterate over each particle to compute stress
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_stress,
                    std::placeholders::_1, phase));
    }

    // Spawn a task for external force
    task_group.run([&] {
      if (!contact) {
        // Iterate over each particle to compute nodal body force
        mesh_->iterate_over_particles(
            std::bind(&mpm::ParticleBase<Tdim>::map_body_force,
                      std::placeholders::_1, phase, this->gravity_));

        // Iterate over each particle to map traction force to nodes
        mesh_->iterate_over_particles(
            std::bind(&mpm::ParticleBase<Tdim>::map_traction_force,
                      std::placeholders::_1, phase));
      }

      if (contact) {
        // Iterate over each particle to compute nodal body force of each
        // subdomian for contact computation
        mesh_->iterate_over_particles(
            std::bind(&mpm::ParticleBase<Tdim>::map_body_force_subdomain,
                      std::placeholders::_1, phase, this->gravity_));
        // Iterate over each particle to map traction force of each subdomian
        // for contact computation
        mesh_->iterate_over_particles(
            std::bind(&mpm::ParticleBase<Tdim>::map_traction_force_subdomain,
                      std::placeholders::_1, phase));
      }

      //! Apply nodal tractions
      if (!contact && nodal_tractions_) this->apply_nodal_tractions();

      //! Apply nodal tractions of each subdomain for contact computation
      if (contact && nodal_tractions_) {
        auto subdomains = io_->analysis()["contact"]["subdomains"];
        for (const auto& sub : subdomains) {
          unsigned mid = sub["material_id"];
          // Get mesh properties
          auto mesh_props = io_->json_object("mesh");
          // Get Mesh reader from JSON object
          const std::string reader =
              mesh_props["mesh_reader"].template get<std::string>();
          auto node_reader =
              Factory<mpm::ReadMesh<Tdim>>::instance()->create(reader);
          bool nodal_tractions_subdomain =
              mesh_->assign_nodal_tractions_subdomain(
                  node_reader->read_particles_tractions(
                      io_->file_name("nodal_tractions")),
                  mid);
        }
      }
    });

    // Spawn a task for internal force
    task_group.run([&] {
      // Iterate over each particle to compute nodal internal force
      if (!contact)
        mesh_->iterate_over_particles(
            std::bind(&mpm::ParticleBase<Tdim>::map_internal_force,
                      std::placeholders::_1, phase));

      // Iterate over each particle to compute nodal internal force of each
      // subdomain for contact computation
      if (contact)
        mesh_->iterate_over_particles(
            std::bind(&mpm::ParticleBase<Tdim>::map_internal_force_subdomain,
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
    if (!contact)
      mesh_->iterate_over_nodes_predicate(
          std::bind(&mpm::NodeBase<Tdim>::compute_acceleration_velocity,
                    std::placeholders::_1, phase, this->dt_),
          std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Iterate over active nodes to compute acceleratation and velocity of each
    // subdomain for contact computation
    if (contact) {
      auto subdomains = io_->analysis()["contact"]["subdomains"];
      for (const auto& sub : subdomains) {
        unsigned mid = sub["material_id"];
        mesh_->iterate_over_nodes_predicate(
            std::bind(
                &mpm::NodeBase<Tdim>::compute_acceleration_velocity_subdomain,
                std::placeholders::_1, phase, this->dt_, mid),
            std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));
      }
    }

    // Use nodal velocity to update position
    if (velocity_update_) {
      // Iterate over each particle to compute updated position
      if (!contact)
        mesh_->iterate_over_particles(std::bind(
            &mpm::ParticleBase<Tdim>::compute_updated_position_velocity,
            std::placeholders::_1, phase, this->dt_));
      // Iterate over each particle to compute updated position of each
      // subdomain for contact computation
      if (contact)
        mesh_->iterate_over_particles(
            std::bind(&mpm::ParticleBase<
                          Tdim>::compute_updated_position_velocity_subdomain,
                      std::placeholders::_1, phase, this->dt_));
    } else {
      // Iterate over each particle to compute updated position
      if (!contact)
        mesh_->iterate_over_particles(
            std::bind(&mpm::ParticleBase<Tdim>::compute_updated_position,
                      std::placeholders::_1, phase, this->dt_));

      // Iterate over each particle to compute updated position of each
      // subdomain for contact computation
      if (contact)
        mesh_->iterate_over_particles(std::bind(
            &mpm::ParticleBase<Tdim>::compute_updated_position_subdomain,
            std::placeholders::_1, phase, this->dt_));
    }

    // Update Stress Last
    if (usl_ == true) {
      // Iterate over each particle to calculate strain
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_strain,
                    std::placeholders::_1, phase, dt_));

      // Iterate over each particle to update particle volume
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::update_volume_strainrate,
                    std::placeholders::_1, phase, this->dt_));

      // Pressure smoothing
      if (pressure_smoothing_) {
        // Assign pressure to nodes
        mesh_->iterate_over_particles(
            std::bind(&mpm::ParticleBase<Tdim>::map_pressure_to_nodes,
                      std::placeholders::_1, phase));

#ifdef USE_MPI
        // Run if there is more than a single MPI task
        if (mpi_size > 1) {
          // MPI all reduce nodal pressure
          mesh_->allreduce_nodal_scalar_property(
              std::bind(&mpm::NodeBase<Tdim>::pressure, std::placeholders::_1,
                        phase),
              std::bind(&mpm::NodeBase<Tdim>::update_pressure,
                        std::placeholders::_1, false, phase,
                        std::placeholders::_2));
        }
#endif

        // Smooth pressure over particles
        mesh_->iterate_over_particles(
            std::bind(&mpm::ParticleBase<Tdim>::compute_pressure_smoothing,
                      std::placeholders::_1, phase));
      }

      // Iterate over each particle to compute stress
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::compute_stress,
                    std::placeholders::_1, phase));
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
                 (this->usl_ ? "USL" : "USF"),
                 std::chrono::duration_cast<std::chrono::milliseconds>(
                     solver_end - solver_begin)
                     .count());

  return status;
}

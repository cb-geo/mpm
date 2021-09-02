//! Constructor
template <unsigned Tdim>
mpm::MPMImplicitLinear<Tdim>::MPMImplicitLinear(const std::shared_ptr<IO>& io)
    : mpm::MPMBase<Tdim>(io) {
  //! Logger
  console_ = spdlog::get("MPMImplicitLinear");

  // Check if stress update is not newmark
  if (stress_update_ != "newmark") {
    console_->warn(
        "The stress_update_ scheme chosen is not available and automatically "
        "set to default: \'newmark\'. Only \'newmark\' scheme is currently "
        "supported for implicit solver.");
    stress_update_ = "newmark";
  }

  // Initialise scheme
  if (stress_update_ == "newmark") {
    mpm_scheme_ = std::make_shared<mpm::MPMSchemeNewmark<Tdim>>(mesh_, dt_);
    // Read parameters of Newmark scheme
    if (analysis_.contains("newmark") && analysis_["newmark"].contains("beta"))
      newmark_beta_ = analysis_["newmark"].at("beta").template get<double>();
    if (analysis_.contains("newmark") && analysis_["newmark"].contains("gamma"))
      newmark_gamma_ = analysis_["newmark"].at("gamma").template get<double>();
  }
}

//! MPM Implicit Linear solver
template <unsigned Tdim>
bool mpm::MPMImplicitLinear<Tdim>::solve() {
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
  const unsigned phase = mpm::ParticlePhase::SinglePhase;

  // Test if checkpoint resume is needed
  bool resume = false;
  if (analysis_.find("resume") != analysis_.end())
    resume = analysis_["resume"]["resume"].template get<bool>();

  // Pressure smoothing
  pressure_smoothing_ = io_->analysis_bool("pressure_smoothing");

  // Initialise material
  this->initialise_materials();

  // Initialise mesh
  this->initialise_mesh();

  // Check point resume
  if (resume) {
    bool check_resume = this->checkpoint_resume();
    if (!check_resume) resume = false;
  }

  // Resume or Initialise
  if (resume) {
    mesh_->resume_domain_cell_ranks();
#ifdef USE_MPI
#ifdef USE_GRAPH_PARTITIONING
    MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif

    //! Particle entity sets and velocity constraints
    this->particle_entity_sets(false);
    this->particle_velocity_constraints();
  } else {
    // Initialise particles
    this->initialise_particles();

    // Compute mass
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::compute_mass, std::placeholders::_1));

    // Domain decompose
    bool initial_step = (resume == true) ? false : true;
    this->mpi_domain_decompose(initial_step);
  }

  // Initialise loading conditions
  this->initialise_loads();

  // Initialise matrix
  bool matrix_status = this->initialise_matrix();
  if (!matrix_status) {
    status = false;
    throw std::runtime_error("Initialisation of matrix failed");
  }

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

    // Mass momentum inertia and compute velocity and acceleration at nodes
    mpm_scheme_->compute_nodal_kinematics(phase);

    // Predict nodal velocity and acceleration -- Predictor step of Newmark
    // scheme
    mpm_scheme_->update_nodal_kinematics_newmark(phase, newmark_beta_,
                                                 newmark_gamma_);

    // Compute local residual force
    mpm_scheme_->compute_forces(gravity_, phase, step_,
                                set_node_concentrated_force_);

    // Reinitialise system matrix to construct equillibrium equation
    bool matrix_reinitialization_status = this->reinitialise_matrix();
    if (!matrix_reinitialization_status) {
      status = false;
      throw std::runtime_error("Reinitialisation of matrix failed");
    }

    // Compute equilibrium equation
    this->compute_equilibrium_equation();

    // Assign displacement increment to nodes
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::update_displacement_increment,
                  std::placeholders::_1, assembler_->displacement_increment(),
                  phase, assembler_->active_dof()),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Update nodal velocity and acceleration -- Corrector step of Newmark
    // scheme
    mpm_scheme_->update_nodal_kinematics_newmark(phase, newmark_beta_,
                                                 newmark_gamma_);

    // Update stress and strain
    mpm_scheme_->postcompute_stress_strain(phase, pressure_smoothing_);

    // Particle kinematics
    mpm_scheme_->compute_particle_kinematics(velocity_update_, phase, "Cundall",
                                             damping_factor_);

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
  console_->info("Rank {}, Implicit Linear {} solver duration: {} ms", mpi_rank,
                 mpm_scheme_->scheme(),
                 std::chrono::duration_cast<std::chrono::milliseconds>(
                     solver_end - solver_begin)
                     .count());

  return status;
}

// Initialise matrix
template <unsigned Tdim>
bool mpm::MPMImplicitLinear<Tdim>::initialise_matrix() {
  bool status = true;
  try {
    // Get matrix assembler type
    std::string assembler_type = analysis_["linear_solver"]["assembler_type"]
                                     .template get<std::string>();
    // Create matrix assembler
    assembler_ =
        Factory<mpm::AssemblerBase<Tdim>, unsigned>::instance()->create(
            assembler_type, std::move(node_neighbourhood_));

    // Solver settings
    if (analysis_["linear_solver"].contains("solver_settings") &&
        analysis_["linear_solver"].at("solver_settings").is_array() &&
        analysis_["linear_solver"].at("solver_settings").size() > 0) {
      mpm::MPMBase<Tdim>::initialise_linear_solver(
          analysis_["linear_solver"]["solver_settings"], linear_solver_);
    }
    // Default solver settings
    else {
      std::string solver_type = "IterativeEigen";
      unsigned max_iter = 1000;
      double tolerance = 1.E-7;

      // In case the default settings are specified in json
      if (analysis_["linear_solver"].contains("solver_type")) {
        solver_type = analysis_["linear_solver"]["solver_type"]
                          .template get<std::string>();
      }
      // Max iteration steps
      if (analysis_["linear_solver"].contains("max_iter")) {
        max_iter =
            analysis_["linear_solver"]["max_iter"].template get<unsigned>();
      }
      // Tolerance
      if (analysis_["linear_solver"].contains("tolerance")) {
        tolerance =
            analysis_["linear_solver"]["tolerance"].template get<double>();
      }

      // NOTE: Only KrylovPETSC solver is supported for MPI
#ifdef USE_MPI
      // Get number of MPI ranks
      int mpi_size = 1;
      MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

      if (solver_type != "KrylovPETSC" && mpi_size > 1) {
        console_->warn(
            "The linear solver in MPI setting is automatically set to default: "
            "\'KrylovPETSC\'. Only \'KrylovPETSC\' solver is supported for "
            "MPI.");
        solver_type = "KrylovPETSC";
      }
#endif

      // Create matrix solver
      auto lin_solver =
          Factory<mpm::SolverBase<Eigen::SparseMatrix<double>>, unsigned,
                  double>::instance()
              ->create(solver_type, std::move(max_iter), std::move(tolerance));
      // Add solver set to map
      linear_solver_.insert(
          std::pair<
              std::string,
              std::shared_ptr<mpm::SolverBase<Eigen::SparseMatrix<double>>>>(
              "displacement", lin_solver));
    }

    // Assign mesh pointer to assembler
    assembler_->assign_mesh_pointer(mesh_);

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Reinitialise and resize matrices at the beginning of every time step
template <unsigned Tdim>
bool mpm::MPMImplicitLinear<Tdim>::reinitialise_matrix() {

  bool status = true;
  try {
    // Assigning matrix id (in each MPI rank)
    const auto nactive_node = mesh_->assign_active_nodes_id();

    // Assigning matrix id globally (required for rank-to-global mapping)
    unsigned nglobal_active_node = nactive_node;
#ifdef USE_MPI
    nglobal_active_node = mesh_->assign_global_active_nodes_id();
#endif

    // Assign global node indice
    assembler_->assign_global_node_indices(nactive_node, nglobal_active_node);

    // Assign displacement constraints
    assembler_->assign_displacement_constraints(this->step_ * this->dt_);

    // Initialise element matrix
    mesh_->iterate_over_cells(
        std::bind(&mpm::Cell<Tdim>::initialise_element_stiffness_matrix,
                  std::placeholders::_1));

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Compute equilibrium equation
template <unsigned Tdim>
bool mpm::MPMImplicitLinear<Tdim>::compute_equilibrium_equation() {
  bool status = true;
  try {
    // Compute local cell stiffness matrices
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::map_material_stiffness_matrix_to_cell,
        std::placeholders::_1));
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_mass_matrix_to_cell,
                  std::placeholders::_1, newmark_beta_, dt_));

    // Assemble global stiffness matrix
    assembler_->assemble_stiffness_matrix();

    // Assemble global residual force RHS vector
    assembler_->assemble_residual_force_right();

    // Apply displacement constraints
    assembler_->apply_displacement_constraints();

#ifdef USE_MPI
    // Assign global active dof to solver
    linear_solver_["displacement"]->assign_global_active_dof(
        3 * assembler_->global_active_dof());

    // Prepare rank global mapper
    auto predictor_rgm = assembler_->rank_global_mapper();
    auto predictor_rgm_y = assembler_->rank_global_mapper();
    auto predictor_rgm_z = assembler_->rank_global_mapper();
    std::for_each(
        predictor_rgm_y.begin(), predictor_rgm_y.end(),
        [size = assembler_->global_active_dof()](int& rgm) { rgm += size; });
    std::for_each(predictor_rgm_z.begin(), predictor_rgm_z.end(),
                  [size = assembler_->global_active_dof()](int& rgm) {
                    rgm += 2 * size;
                  });
    predictor_rgm.insert(predictor_rgm.end(), predictor_rgm_y.begin(),
                         predictor_rgm_y.end());
    predictor_rgm.insert(predictor_rgm.end(), predictor_rgm_z.begin(),
                         predictor_rgm_z.end());

    // Assign rank global mapper to solver
    linear_solver_["displacement"]->assign_rank_global_mapper(predictor_rgm);
#endif

    // Solve matrix equation and assign solution to assembler
    assembler_->assign_displacement_increment(
        linear_solver_["displacement"]->solve(
            assembler_->stiffness_matrix(),
            assembler_->residual_force_rhs_vector()));

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}
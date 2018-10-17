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

#ifdef USE_MPI
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  // Get number of MPI ranks
  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif

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

  // Compute mass
  meshes_.at(0)->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::compute_mass, std::placeholders::_1, phase));

  // Test if checkpoint resume is needed
  bool resume = false;
  if (analysis_.find("resume") != analysis_.end())
    resume = analysis_["resume"]["resume"].template get<bool>();
  if (resume) this->checkpoint_resume();

  // Main loop
  for (; step_ < nsteps_; ++step_) {

#ifdef USE_MPI
    if (mpi_rank == 0) console_->info("Step: {} of {}.\n", step_, nsteps_);
#else
    console_->info("Step: {} of {}.\n", step_, nsteps_);
#endif

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

#ifdef USE_MPI
    // Run if there is more than a single MPI task
    if (mpi_size > 1) {
      // MPI all reduce nodal mass
      meshes_.at(0)->allreduce_nodal_scalar_property(
          std::bind(&mpm::NodeBase<Tdim>::mass, std::placeholders::_1, phase),
          std::bind(&mpm::NodeBase<Tdim>::update_mass, std::placeholders::_1,
                    false, phase, std::placeholders::_2));
      // MPI all reduce nodal momentum
      meshes_.at(0)->allreduce_nodal_vector_property(
          std::bind(&mpm::NodeBase<Tdim>::momentum, std::placeholders::_1,
                    phase),
          std::bind(&mpm::NodeBase<Tdim>::update_momentum,
                    std::placeholders::_1, false, phase,
                    std::placeholders::_2));
    }
#endif

    // Compute nodal velocity
    meshes_.at(0)->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_velocity,
                  std::placeholders::_1),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Iterate over each particle to calculate strain
    meshes_.at(0)->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::compute_strain,
                  std::placeholders::_1, phase, dt_));

    // Iterate over each particle to compute stress
    meshes_.at(0)->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::compute_stress,
                  std::placeholders::_1, phase));

    // Iterate over each particle to compute nodal body force
    meshes_.at(0)->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_body_force,
                  std::placeholders::_1, phase, this->gravity_));

    // Iterate over each particle to map traction force to nodes
    meshes_.at(0)->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_traction_force,
                  std::placeholders::_1, phase));

    // Iterate over each particle to compute nodal internal force
    meshes_.at(0)->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_internal_force,
                  std::placeholders::_1, phase));

#ifdef USE_MPI
    // Run if there is more than a single MPI task
    if (mpi_size > 1) {
      // MPI all reduce external force
      meshes_.at(0)->allreduce_nodal_vector_property(
          std::bind(&mpm::NodeBase<Tdim>::external_force, std::placeholders::_1,
                    phase),
          std::bind(&mpm::NodeBase<Tdim>::update_external_force,
                    std::placeholders::_1, false, phase,
                    std::placeholders::_2));
      // MPI all reduce internal force
      meshes_.at(0)->allreduce_nodal_vector_property(
          std::bind(&mpm::NodeBase<Tdim>::internal_force, std::placeholders::_1,
                    phase),
          std::bind(&mpm::NodeBase<Tdim>::update_internal_force,
                    std::placeholders::_1, false, phase,
                    std::placeholders::_2));
    }
#endif

    // Iterate over active nodes to compute acceleratation and velocity
    meshes_.at(0)->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_acceleration_velocity,
                  std::placeholders::_1, phase, this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Iterate over each particle to compute updated position
    meshes_.at(0)->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::compute_updated_position,
                  std::placeholders::_1, phase, this->dt_));

    // Iterate over each particle to update particle volume
    meshes_.at(0)->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::update_volume_strainrate,
                  std::placeholders::_1, phase, this->dt_));

    // Locate particles
    auto unlocatable_particles = meshes_.at(0)->locate_particles_mesh();

    if (!unlocatable_particles.empty())
      throw std::runtime_error("Particle outside the mesh domain");

    if (step_ % output_steps_ == 0) {
#ifdef USE_MPI
      // HDF5 outputs
      this->write_hdf5(this->step_, this->nsteps_);
#ifdef USE_VTK
      // VTK outputs
      this->write_vtk(this->step_, this->nsteps_);
#endif  // VTK

#else  // MPI Else
       // HDF5 outputs
      this->write_hdf5(this->step_, this->nsteps_);
#ifdef USE_VTK
      // VTK outputs
      this->write_vtk(this->step_, this->nsteps_);
#endif  // VTK

#endif  // MPI
    }
  }
  return status;
}

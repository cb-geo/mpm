//! Constructor
template <unsigned Tdim>
mpm::MPMSchemeMUSL<Tdim>::MPMSchemeMUSL(
    const std::shared_ptr<mpm::Mesh<Tdim>>& mesh, double dt)
    : mpm::MPMScheme<Tdim>(mesh, dt) {}

//! Precompute stresses and strains
template <unsigned Tdim>
inline void mpm::MPMSchemeMUSL<Tdim>::precompute_stress_strain(
    unsigned phase, bool pressure_smoothing) {}

//! Postcompute stresses and strains
template <unsigned Tdim>
inline void mpm::MPMSchemeMUSL<Tdim>::postcompute_stress_strain(
    unsigned phase, bool pressure_smoothing) {
  mpm::MPMScheme<Tdim>::compute_stress_strain(phase, pressure_smoothing);
}

// Compute particle kinematics
template <unsigned Tdim>
inline void mpm::MPMSchemeMUSL<Tdim>::compute_particle_kinematics(
    bool velocity_update, unsigned phase, const std::string& damping_type,
    double damping_factor) {

  // Check if damping has been specified and accordingly Iterate over
  // active nodes to compute acceleratation and velocity
  if (damping_type == "Cundall")
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_acceleration_velocity_cundall,
                  std::placeholders::_1, phase, dt_, damping_factor),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));
  else
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_acceleration_velocity,
                  std::placeholders::_1, phase, dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

  // Iterate over each particle to compute updated velocity
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::compute_updated_velocity_musl,
                std::placeholders::_1, dt_, velocity_update));

  // Apply particle velocity constraints
  mesh_->apply_particle_velocity_constraints();
}

//! Recompute nodal kinematics and compute particle position (only for musl)
template <unsigned Tdim>
inline void mpm::MPMSchemeMUSL<Tdim>::compute_particle_updated_position(
    bool velocity_update, unsigned phase) {

#pragma omp parallel sections
  {
    // Spawn a task for initialising nodal momentum
#pragma omp section
    {
      // Initialise nodal momentum
      mesh_->iterate_over_nodes_predicate(
          std::bind(&mpm::NodeBase<Tdim>::initialise_momentum, std::placeholders::_1),
          std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));
    }
  }

  // Assign momentum to nodes
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::map_momentum_to_nodes,
                std::placeholders::_1));

#ifdef USE_MPI
  // Run if there is more than a single MPI task
  if (mpi_size_ > 1) {
    // MPI all reduce nodal mass
    mesh_->template nodal_halo_exchange<double, 1>(
        std::bind(&mpm::NodeBase<Tdim>::mass, std::placeholders::_1, phase),
        std::bind(&mpm::NodeBase<Tdim>::update_mass, std::placeholders::_1,
                  false, phase, std::placeholders::_2));
    // MPI all reduce nodal momentum
    mesh_->template nodal_halo_exchange<Eigen::Matrix<double, Tdim, 1>, Tdim>(
        std::bind(&mpm::NodeBase<Tdim>::momentum, std::placeholders::_1, phase),
        std::bind(&mpm::NodeBase<Tdim>::update_momentum, std::placeholders::_1,
                  false, phase, std::placeholders::_2));
  }
#endif

  // Compute nodal velocity
  mesh_->iterate_over_nodes_predicate(
      std::bind(&mpm::NodeBase<Tdim>::compute_velocity, std::placeholders::_1),
      std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

  // Iterate over each particle to compute updated position
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::compute_updated_position_musl,
                std::placeholders::_1, dt_, velocity_update));
}

//! Stress update scheme
template <unsigned Tdim>
inline std::string mpm::MPMSchemeMUSL<Tdim>::scheme() const {
  return "MUSL";
}

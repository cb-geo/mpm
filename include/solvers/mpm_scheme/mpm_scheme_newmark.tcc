//! Constructor
template <unsigned Tdim>
mpm::MPMSchemeNewmark<Tdim>::MPMSchemeNewmark(
    const std::shared_ptr<mpm::Mesh<Tdim>>& mesh, double dt)
    : mpm::MPMScheme<Tdim>(mesh, dt) {}

//! Initialize nodes, cells and shape functions
template <unsigned Tdim>
inline void mpm::MPMSchemeNewmark<Tdim>::initialise() {
#pragma omp parallel sections
  {
    // Spawn a task for initialising nodes and cells
#pragma omp section
    {
      // Initialise nodes
      mesh_->iterate_over_nodes(
          std::bind(&mpm::NodeBase<Tdim>::initialise, std::placeholders::_1));

      mesh_->iterate_over_cells(
          std::bind(&mpm::Cell<Tdim>::activate_nodes, std::placeholders::_1));
    }
    // Spawn a task for particles
#pragma omp section
    {
      // Iterate over each particle to compute shapefn
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::compute_shapefn, std::placeholders::_1));

      // Iterate over each particle to initialise displacement
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::initialise_displacement,
                    std::placeholders::_1));
    }
  }  // Wait to complete
}

//! Compute nodal kinematics - map mass, momentum and inertia to nodes
template <unsigned Tdim>
inline void mpm::MPMSchemeNewmark<Tdim>::compute_nodal_kinematics(
    unsigned phase) {
  // Assign mass, momentum and inertia to nodes
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::map_mass_momentum_inertia_to_nodes,
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
    // MPI all reduce nodal inertia
    mesh_->template nodal_halo_exchange<Eigen::Matrix<double, Tdim, 1>, Tdim>(
        std::bind(&mpm::NodeBase<Tdim>::inertia, std::placeholders::_1, phase),
        std::bind(&mpm::NodeBase<Tdim>::update_inertia, std::placeholders::_1,
                  false, phase, std::placeholders::_2));
  }
#endif

  // Compute nodal velocity and acceleration
  mesh_->iterate_over_nodes_predicate(
      std::bind(&mpm::NodeBase<Tdim>::compute_velocity_acceleration,
                std::placeholders::_1),
      std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));
}

//! Update nodal kinematics by Newmark scheme
template <unsigned Tdim>
inline void mpm::MPMSchemeNewmark<Tdim>::update_nodal_kinematics_newmark(
    unsigned phase, double newmark_beta, double newmark_gamma) {

  // Update nodal velocity and acceleration
  mesh_->iterate_over_nodes_predicate(
      std::bind(&mpm::NodeBase<Tdim>::update_velocity_acceleration_newmark,
                std::placeholders::_1, phase, newmark_beta, newmark_gamma, dt_),
      std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));
}

//! Compute stress and strain by Newmark scheme
template <unsigned Tdim>
inline void mpm::MPMSchemeNewmark<Tdim>::compute_stress_strain(
    unsigned phase, bool pressure_smoothing) {

  // Iterate over each particle to calculate strain using nodal displacement
  mesh_->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::compute_strain_newmark, std::placeholders::_1));

  // Iterate over each particle to update particle volume
  mesh_->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::update_volume, std::placeholders::_1));

  // Pressure smoothing
  if (pressure_smoothing) this->pressure_smoothing(phase);

  // Iterate over each particle to compute stress
  mesh_->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::compute_stress, std::placeholders::_1));
}

//! Precompute stresses and strains
template <unsigned Tdim>
inline void mpm::MPMSchemeNewmark<Tdim>::precompute_stress_strain(
    unsigned phase, bool pressure_smoothing) {}

//! Postcompute stresses and strains
template <unsigned Tdim>
inline void mpm::MPMSchemeNewmark<Tdim>::postcompute_stress_strain(
    unsigned phase, bool pressure_smoothing) {
  this->compute_stress_strain(phase, pressure_smoothing);
}

// Compute forces
template <unsigned Tdim>
inline void mpm::MPMSchemeNewmark<Tdim>::compute_forces(
    const Eigen::Matrix<double, Tdim, 1>& gravity, unsigned phase,
    unsigned step, bool concentrated_nodal_forces) {
  // Spawn a task for external force
#pragma omp parallel sections
  {
#pragma omp section
    {
      // Iterate over each particle to compute nodal body force
      mesh_->iterate_over_particles(
          std::bind(&mpm::ParticleBase<Tdim>::map_body_force,
                    std::placeholders::_1, gravity));

      // Iterate over each particle to compute nodal inertial force
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::map_inertial_force, std::placeholders::_1));

      // Apply particle traction and map to nodes
      mesh_->apply_traction_on_particles(step * dt_);

      // Iterate over each node to add concentrated node force to external
      // force
      if (concentrated_nodal_forces)
        mesh_->iterate_over_nodes(
            std::bind(&mpm::NodeBase<Tdim>::apply_concentrated_force,
                      std::placeholders::_1, phase, (step * dt_)));
    }

#pragma omp section
    {
      // Spawn a task for internal force
      // Iterate over each particle to compute nodal internal force
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::map_internal_force, std::placeholders::_1));
    }
  }  // Wait for tasks to finish

#ifdef USE_MPI
  // Run if there is more than a single MPI task
  if (mpi_size_ > 1) {
    // MPI all reduce external force
    mesh_->template nodal_halo_exchange<Eigen::Matrix<double, Tdim, 1>, Tdim>(
        std::bind(&mpm::NodeBase<Tdim>::external_force, std::placeholders::_1,
                  phase),
        std::bind(&mpm::NodeBase<Tdim>::update_external_force,
                  std::placeholders::_1, false, phase, std::placeholders::_2));
    // MPI all reduce internal force
    mesh_->template nodal_halo_exchange<Eigen::Matrix<double, Tdim, 1>, Tdim>(
        std::bind(&mpm::NodeBase<Tdim>::internal_force, std::placeholders::_1,
                  phase),
        std::bind(&mpm::NodeBase<Tdim>::update_internal_force,
                  std::placeholders::_1, false, phase, std::placeholders::_2));
  }
#endif
}

// Compute particle kinematics
template <unsigned Tdim>
inline void mpm::MPMSchemeNewmark<Tdim>::compute_particle_kinematics(
    bool velocity_update, unsigned phase, const std::string& damping_type,
    double damping_factor) {

  // Iterate over each particle to compute updated position
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::compute_updated_position_newmark,
                std::placeholders::_1, dt_, velocity_update));
}

//! Postcompute nodal kinematics - map mass and momentum to nodes
template <unsigned Tdim>
inline void mpm::MPMSchemeNewmark<Tdim>::postcompute_nodal_kinematics(
    unsigned phase) {}

//! Stress update scheme
template <unsigned Tdim>
inline std::string mpm::MPMSchemeNewmark<Tdim>::scheme() const {
  return "Newmark";
}

// Compute mass of particle
namespace mpm {
namespace particle {

// Compute particle mass
template <unsigned Tdim>
void compute_mass(std::shared_ptr<mpm::ParticleBase<Tdim>> particle) noexcept {
  // Check if particle volume is set and material ptr is valid
  assert(particle->volume() != std::numeric_limits<double>::max() &&
         particle->material() != nullptr);

  // Mass = volume of particle * mass_density
  particle->update_scalar_property(
      mpm::properties::Scalar::MassDensity, false,
      particle->material()->template property<double>(std::string("density")));

  // Update particle mass
  particle->update_scalar_property(
      mpm::properties::Scalar::Mass, false,
      particle->volume() * particle->mass_density());
}

// Compute volume of the particle
template <unsigned Tdim>
void compute_volume(
    std::shared_ptr<mpm::ParticleBase<Tdim>> particle) noexcept {
  // Check if particle has a valid cell ptr
  assert(particle->cell_ptr());

  // Volume of the cell / # of particles
  particle->assign_volume(particle->cell()->volume() /
                          particle->cell()->nparticles());
}

// Update volume based on the central strain rate
template <unsigned Tdim>
void update_volume(std::shared_ptr<mpm::ParticleBase<Tdim>> particle) noexcept {
  // Check if particle has a valid cell ptr and a valid volume
  assert(particle->cell_ptr() &&
         particle->volume() != std::numeric_limits<double>::max());

  // Compute at centroid
  // Strain rate for reduced integration
  particle->update_scalar_property(
      mpm::properties::Scalar::Volume, false,
      (particle->volume() * (1. + particle->dvolumetric_strain())));
  particle->update_scalar_property(
      mpm::properties::Scalar::MassDensity, false,
      (particle->mass_density() / (1. + particle->dvolumetric_strain())));
}

//! Map particle mass and momentum to nodes
template <unsigned Tdim>
void map_mass_momentum_to_nodes(
    std::shared_ptr<mpm::ParticleBase<Tdim>> particle) noexcept {
  // Check if particle mass is set
  assert(particle->mass() != std::numeric_limits<double>::max());

  // Map mass and momentum to nodes
  particle->map_scalar_property_to_nodes(mpm::properties::Scalar::Mass, true,
                                         mpm::ParticlePhase::Solid);
  particle->map_vector_property_to_nodes(
      mpm::properties::Vector::Momentum, true, mpm::ParticlePhase::Solid,
      particle->mass() * particle->velocity());
}

//! Map particle pressure to nodes
template <unsigned Tdim>
void map_mass_pressure_to_nodes(
    std::shared_ptr<mpm::ParticleBase<Tdim>> particle) noexcept {
  // Mass is initialized
  assert(particle->mass() != std::numeric_limits<double>::max());

  // Check if state variable pressure is found
  if (particle->pressure() != std::numeric_limits<double>::quiet_NaN()) {
    // Map particle pressure to nodes
    particle->map_scalar_property_to_nodes(
        mpm::properties::Scalar::MassPressure, true, mpm::ParticlePhase::Solid,
        particle->mass() * particle->pressure());
  }
}

//! Map body force to nodes
template <unsigned Tdim>
void map_body_force(std::shared_ptr<mpm::ParticleBase<Tdim>> particle,
                    const Eigen::Matrix<double, Tdim, 1>& pgravity) noexcept {
  // Compute nodal body forces
  particle->map_vector_property_to_nodes(mpm::properties::Vector::ExternalForce,
                                         true, mpm::ParticlePhase::Solid,
                                         pgravity * particle->mass());
}

//! Map traction force
template <unsigned Tdim>
void map_traction_force(
    std::shared_ptr<mpm::ParticleBase<Tdim>> particle) noexcept {
  if (particle->set_traction()) {
    // Map particle traction forces to nodes
    particle->map_vector_property_to_nodes(
        mpm::properties::Vector::ExternalForce, true, mpm::ParticlePhase::Solid,
        particle->traction());
  }
}

// Compute pressure smoothing of the particle based on nodal pressure
template <unsigned Tdim>
void compute_pressure_smoothing(
    std::shared_ptr<mpm::ParticleBase<Tdim>> particle) noexcept {
  // Assert
  assert(particle->cell_ptr());

  // Check if particle has pressure
  if (particle->pressure() != std::numeric_limits<double>::quiet_NaN()) {
    double pressure = particle->interpolate_scalar_property_from_nodes(
        mpm::properties::Scalar::Pressure, mpm::ParticlePhase::Solid);
    particle->assign_state_variable("pressure", pressure);
  }
}

}  // namespace particle
}  // namespace mpm

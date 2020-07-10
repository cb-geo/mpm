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

//! Map particle mass and momentum to nodes
template <unsigned Tdim>
void map_mass_momentum_to_nodes(
    std::shared_ptr<mpm::ParticleBase<Tdim>> particle) noexcept {
  // Check if particle mass is set
  assert(particle->mass() != std::numeric_limits<double>::max());

  // Map mass and momentum to nodes
  particle->map_scalar_property_nodes(mpm::properties::Scalar::Mass, true,
                                      mpm::ParticlePhase::Solid);
  particle->map_vector_property_nodes(mpm::properties::Vector::Momentum, true,
                                      mpm::ParticlePhase::Solid,
                                      particle->mass() * particle->velocity());
}

}  // namespace particle
}  // namespace mpm

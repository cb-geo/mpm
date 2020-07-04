// Compute mass of particle
namespace mpm {
namespace particle {
template <unsigned Tdim>
void compute_mass(std::shared_ptr<mpm::ParticleBase<Tdim>> particle) noexcept {
  // Check if particle volume is set and material ptr is valid
  assert(particle->volume() != std::numeric_limits<double>::max() &&
         particle->material() != nullptr);
  // Mass = volume of particle * mass_density
  auto density =
      (particle->material())->template property<double>(std::string("density"));
  // Update particle mass
  particle->update_scalar_property(mpm::properties::Scalar::Mass, false,
                                   particle->volume() * density);
}
}  // namespace particle
}  // namespace mpm

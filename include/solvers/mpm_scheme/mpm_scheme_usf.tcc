//! Constructor
template <unsigned Tdim>
mpm::MPMSchemeUSF<Tdim>::MPMSchemeUSF(
    const std::shared_ptr<mpm::Mesh<Tdim>>& mesh, double dt)
    : mpm::MPMScheme<Tdim>(mesh, dt) {}

//! Precompute stresses and strains
template <unsigned Tdim>
inline void mpm::MPMSchemeUSF<Tdim>::precompute_stress_strain(
    unsigned phase, bool pressure_smoothing) {
  this->compute_stress_strain(phase, pressure_smoothing);
}

//! Postcompute stresses and strains
template <unsigned Tdim>
inline void mpm::MPMSchemeUSF<Tdim>::postcompute_stress_strain(
    unsigned phase, bool pressure_smoothing) {}

//! Recompute nodal kinematics and compute particle position (only for musl)
template <unsigned Tdim>
inline void mpm::MPMSchemeUSF<Tdim>::compute_particle_updated_position(
  bool velocity_update, unsigned phase) {}

//! Stress update scheme
template <unsigned Tdim>
inline std::string mpm::MPMSchemeUSF<Tdim>::scheme() const {
  return "USF";
}

//! Constructor
template <unsigned Tdim>
mpm::StressUpdateUSF<Tdim>::StressUpdateUSF(
    const std::shared_ptr<mpm::Mesh<Tdim>>& mesh, double dt)
    : mpm::StressUpdate<Tdim>(mesh, dt) {
  mesh_ = mesh;
  dt_ = dt;
}

//! Precompute stresses and strains
template <unsigned Tdim>
inline void mpm::StressUpdateUSF<Tdim>::precompute_stress_strain(
    unsigned phase, bool pressure_smoothing) {
  this->compute_stress_strain(phase, pressure_smoothing);
}

//! Postcompute stresses and strains
template <unsigned Tdim>
inline void mpm::StressUpdateUSF<Tdim>::postcompute_stress_strain(
    unsigned phase, bool pressure_smoothing) {}

//! Stress update scheme
template <unsigned Tdim>
inline std::string mpm::StressUpdateUSF<Tdim>::scheme() const {
  return "USF";
}

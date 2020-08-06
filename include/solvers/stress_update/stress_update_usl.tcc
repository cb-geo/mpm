//! Constructor
template <unsigned Tdim>
mpm::StressUpdateUSL<Tdim>::StressUpdateUSL(
    const std::shared_ptr<mpm::Mesh<Tdim>>& mesh, double dt)
    : mpm::StressUpdate<Tdim>(mesh, dt) {
  mesh_ = mesh;
  dt_ = dt;
}

//! Precompute stresses and strains
template <unsigned Tdim>
void mpm::StressUpdateUSL<Tdim>::precompute_stress_strain(
    unsigned phase, bool pressure_smoothing) {}

//! Postcompute stresses and strains
template <unsigned Tdim>
void mpm::StressUpdateUSL<Tdim>::postcompute_stress_strain(
    unsigned phase, bool pressure_smoothing) {
  mpm::StressUpdate<Tdim>::compute_stress_strain(phase, pressure_smoothing);
}

//! Stress update scheme
template <unsigned Tdim>
std::string mpm::StressUpdateUSL<Tdim>::scheme() const {
  return "USL";
}

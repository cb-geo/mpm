//! Constructor with id and coordinates
//! \param[in] id Particle id
//! \param[in] coord coordinates of the particle
//! \tparam Tdim Dimension
//! \tparam Tnphases Number of phases
template <unsigned Tdim, unsigned Tnphases>
mpm::Particle<Tdim, Tnphases>::Particle(Index id, const VectorDim& coord)
    : mpm::ParticleBase<Tdim>(id, coord) {}

//! Constructor with id, coordinates and status
//! \param[in] id Particle id
//! \param[in] coord coordinates of the particle
//! \param[in] status Particle status (active / inactive)
//! \tparam Tdim Dimension
//! \tparam Tnphases Number of phases
template <unsigned Tdim, unsigned Tnphases>
mpm::Particle<Tdim, Tnphases>::Particle(Index id, const VectorDim& coord, bool status)
    : mpm::ParticleBase<Tdim>(id, coord, status) {}

// Assign a cell to particle
//! \param[in] cellptr Pointer to a cell
//! \tparam Tdim Dimension
//! \tparam Tnphases Number of phases
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::assign_cell(
    const std::shared_ptr<Cell<Tdim>>& cellptr) {
  cell_ = cellptr;
  return cell_->add_particle_id(this->id());
}

//! Constructor with id and coordinates
//! \param[in] id Particle id
//! \param[in] coord coordinates of the particle
//! \tparam Tdim Dimension
template <unsigned Tdim>
mpm::Particle<Tdim>::Particle(Index id, const VectorDim& coord) : id_{id} {
  // Check if the dimension is between 1 & 3
  static_assert((Tdim >= 1 && Tdim <= 3), "Invalid global dimension");
  coordinates_ = coord;
  status_ = true;
}

//! Constructor with id, coordinates and status
//! \param[in] id Particle id
//! \param[in] coord coordinates of the particle
//! \param[in] status Particle status (active / inactive)
//! \tparam Tdim Dimension
template <unsigned Tdim>
mpm::Particle<Tdim>::Particle(Index id, const VectorDim& coord, bool status)
    : mpm::Particle<Tdim>::Particle(id, coord) {
  status_ = status;
}

// Assign a cell to particle
//! \param[in] cellptr Pointer to a cell
//! \tparam Tdim Dimension
template <unsigned Tdim>
bool mpm::Particle<Tdim>::assign_cell(
    const std::shared_ptr<Cell<Tdim>>& cellptr) {
  cell_ = cellptr;
  return cell_->add_particle_id(this->id());
}

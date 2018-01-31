//! Constructor with id and coordinates
//! \param[in] id Particle id
//! \param[in] coord coordinates of the particle
//! \tparam Tdim Dimension
template <unsigned Tdim>
mpm::ParticleBase<Tdim>::ParticleBase(Index id, const VectorDim& coord) : id_{id} {
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
mpm::ParticleBase<Tdim>::ParticleBase(Index id, const VectorDim& coord, bool status)
    : mpm::ParticleBase<Tdim>::ParticleBase(id, coord) {
  status_ = status;
}

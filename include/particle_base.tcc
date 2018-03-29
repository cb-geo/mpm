//! Constructor with id and coordinates
template <unsigned Tdim>
mpm::ParticleBase<Tdim>::ParticleBase(Index id, const VectorDim& coord)
    : id_{id} {
  // Check if the dimension is 2 or 3
  static_assert((Tdim == 2 || Tdim == 3), "Invalid global dimension");
  coordinates_ = coord;
  status_ = true;
}

//! Constructor with id, coordinates and status
template <unsigned Tdim>
mpm::ParticleBase<Tdim>::ParticleBase(Index id, const VectorDim& coord,
                                      bool status)
    : mpm::ParticleBase<Tdim>::ParticleBase(id, coord) {
  status_ = status;
}

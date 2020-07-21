//! Constructor with id and coordinates
template <unsigned Tdim>
mpm::ParticleBase<Tdim>::ParticleBase(Index id, const VectorDim& coord)
    : id_{id} {
  // Check if the dimension is between 1 & 3
  static_assert((Tdim >= 1 && Tdim <= 3), "Invalid global dimension");
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

//! Update scalar property at particle
template <unsigned Tdim>
void mpm::ParticleBase<Tdim>::update_scalar_property(
    mpm::properties::Scalar property, bool update, double value) noexcept {
  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;
  scalar_properties_.at(property) =
      scalar_properties_.at(property) * factor + value;
}

//! Update scalar property at particle
template <unsigned Tdim>
double mpm::ParticleBase<Tdim>::scalar_property(
    mpm::properties::Scalar property) const {
  return scalar_properties_.at(property);
}

//! Map scalar property to nodes
template <unsigned Tdim>
void mpm::ParticleBase<Tdim>::map_scalar_property_to_nodes(
    mpm::properties::Scalar property, bool update, unsigned phase) noexcept {
  // Check if particle property is set
  assert(scalar_properties_.at(property) != std::numeric_limits<double>::max());

  // Map scalar property to nodes
  for (unsigned i = 0; i < nodes_.size(); ++i)
    nodes_[i]->update_scalar_property(
        property, update, phase, scalar_properties_.at(property) * shapefn_[i]);
}

//! Map an arbitrary scalar value to nodal scalar property
template <unsigned Tdim>
void mpm::ParticleBase<Tdim>::map_scalar_property_to_nodes(
    mpm::properties::Scalar property, bool update, unsigned phase,
    double value) noexcept {
  // Map scalar value to nodes
  for (unsigned i = 0; i < nodes_.size(); ++i)
    nodes_[i]->update_scalar_property(property, update, phase,
                                      value * shapefn_[i]);
}

//! Interpolate scalar property from nodes
template <unsigned Tdim>
double mpm::ParticleBase<Tdim>::interpolate_scalar_property_from_nodes(
    mpm::properties::Scalar property, unsigned phase) const {
  double value = 0.;
  // Interpolate scalar property from nodes
  for (unsigned i = 0; i < nodes_.size(); ++i)
    value += nodes_[i]->scalar_property(property, phase) * shapefn_[i];
  return value;
}

//! Update vector property at particle
template <unsigned Tdim>
void mpm::ParticleBase<Tdim>::update_vector_property(
    mpm::properties::Vector property, bool update,
    const Eigen::Matrix<double, Tdim, 1>& value) noexcept {
  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;
  vector_properties_.at(property) =
      vector_properties_.at(property) * factor + value;
}

//! Update vector property at particle
template <unsigned Tdim>
Eigen::Matrix<double, Tdim, 1> mpm::ParticleBase<Tdim>::vector_property(
    mpm::properties::Vector property) const {
  return vector_properties_.at(property);
}

//! Map vector property to nodes
template <unsigned Tdim>
void mpm::ParticleBase<Tdim>::map_vector_property_to_nodes(
    mpm::properties::Vector property, bool update, unsigned phase) noexcept {
  // Map vector property to nodes
  for (unsigned i = 0; i < nodes_.size(); ++i)
    nodes_[i]->update_vector_property(
        property, update, phase, vector_properties_.at(property) * shapefn_[i]);
}

//! Map an arbitrary vector value to nodal vector property
template <unsigned Tdim>
void mpm::ParticleBase<Tdim>::map_vector_property_to_nodes(
    mpm::properties::Vector property, bool update, unsigned phase,
    const Eigen::Matrix<double, Tdim, 1>& value) noexcept {
  // Map vector property to nodes
  for (unsigned i = 0; i < nodes_.size(); ++i)
    nodes_[i]->update_vector_property(property, update, phase,
                                      value * shapefn_[i]);
}

//! Interpolate vector property from nodes
template <unsigned Tdim>
Eigen::Matrix<double, Tdim, 1>
    mpm::ParticleBase<Tdim>::interpolate_vector_property_from_nodes(
        mpm::properties::Vector property, unsigned phase) const {
  Eigen::Matrix<double, Tdim, 1> value = Eigen::Matrix<double, Tdim, 1>::Zero();
  // Interpolate vector property from nodes
  for (unsigned i = 0; i < nodes_.size(); ++i)
    value += nodes_[i]->vector_property(property, phase) * shapefn_[i];
  return value;
}

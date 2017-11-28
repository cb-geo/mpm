// Assign a cell to particle
//! \param[in] cellptr Pointer to a cell
//! \tparam Tdim Dimension
template <unsigned Tdim>
bool mpm::Particle<Tdim>::assign_cell(
    const std::shared_ptr<Cell<Tdim>>& cellptr) {
  cell_ = cellptr;
  return cell_->add_particle_id(this->id());
}

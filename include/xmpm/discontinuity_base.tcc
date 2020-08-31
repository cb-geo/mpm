//! Constructor
template <unsigned Tdim>
mpm::DiscontinuityBase<Tdim>::DiscontinuityBase(
    unsigned id, const Json& discontinuity_props) {

  friction_coef_ = 0;

  std::string logger = "discontinuity::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

//! create points from file
template <unsigned Tdim>
bool mpm::DiscontinuityBase<Tdim>::create_points(
    const std::vector<VectorDim>& coordinates) {
  bool status = true;
  try {
    // Check if point coordinates is empty
    if (coordinates.empty())
      throw std::runtime_error("List of coordinates is empty");
    // Iterate over all coordinates
    for (const auto& point_coordinates : coordinates) {

      // Add point
      mpm::discontinuity_point<Tdim> point(point_coordinates);

      points_.emplace_back(point);  //
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}
//! Locate points in a cell
template <unsigned Tdim>
void mpm::DiscontinuityBase<Tdim>::locate_discontinuity_mesh(
    const Vector<Cell<Tdim>>& cells,
    const Map<Cell<Tdim>>& map_cells) noexcept {
  for (auto& point : this->points_)
    point.locate_discontinuity_mesh(cells, map_cells);
}

// Compute updated position of the particle
template <unsigned Tdim>
void mpm::DiscontinuityBase<Tdim>::compute_updated_position(
    double dt) noexcept {
  for (auto& point : this->points_) point.compute_updated_position(dt);
}

// Compute updated position of the particle
template <unsigned Tdim>
void mpm::DiscontinuityBase<Tdim>::compute_shapefn() noexcept {
  for (auto& point : this->points_) point.compute_shapefn();
}
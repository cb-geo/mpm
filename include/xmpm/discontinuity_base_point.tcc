// Assign a cell to particle
template <unsigned Tdim>
bool mpm::discontinuity_point<Tdim>::assign_cell_xi(
    const std::shared_ptr<Cell<Tdim>>& cellptr,
    const Eigen::Matrix<double, Tdim, 1>& xi) {
  bool status = true;
  try {
    // Assign cell to the new cell ptr, if point can be found in new cell
    if (cellptr != nullptr) {

      cell_ = cellptr;
      cell_id_ = cellptr->id();
      nodes_ = cell_->nodes();
      // assign discontinuity_enrich
      for (unsigned i = 0; i < nodes_.size(); ++i)
        nodes_[i]->assign_discontinuity_enrich(true);
      // Assign the reference location of particle
      bool xi_nan = false;

      // Check if point is within the cell
      for (unsigned i = 0; i < xi.size(); ++i)
        if (xi(i) < -1. || xi(i) > 1. || std::isnan(xi(i))) xi_nan = true;

      if (xi_nan == false)
        this->xi_ = xi;
      else
        return false;
    } else {
      console_->warn("Points of discontinuity cannot be found in cell!");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Assign a cell to point
template <unsigned Tdim>
bool mpm::discontinuity_point<Tdim>::assign_cell(
    const std::shared_ptr<Cell<Tdim>>& cellptr) {
  bool status = true;
  try {
    Eigen::Matrix<double, Tdim, 1> xi;
    // Assign cell to the new cell ptr, if point can be found in new cell
    if (cellptr->is_point_in_cell(this->coordinates_, &xi)) {

      cell_ = cellptr;
      cell_id_ = cellptr->id();
      nodes_ = cell_->nodes();
      // assign discontinuity_enrich
      for (unsigned i = 0; i < nodes_.size(); ++i)
        nodes_[i]->assign_discontinuity_enrich(true);
    } else {
      console_->warn("Points of discontinuity cannot be found in cell!");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Compute reference location cell to particle
template <unsigned Tdim>
bool mpm::discontinuity_point<Tdim>::compute_reference_location() noexcept {
  // Set status of compute reference location
  bool status = false;
  // Compute local coordinates
  Eigen::Matrix<double, Tdim, 1> xi;
  // Check if the point is in cell
  if (cell_ != nullptr && cell_->is_point_in_cell(this->coordinates_, &xi)) {
    this->xi_ = xi;
    status = true;
  }

  return status;
}

//! Locate points in a cell
template <unsigned Tdim>
void mpm::discontinuity_point<Tdim>::locate_discontinuity_mesh(const
    Vector<Cell<Tdim>>& cells, const Map<Cell<Tdim>>& map_cells) noexcept {
  // Check the current cell if it is not invalid
  if (cell_id() != std::numeric_limits<mpm::Index>::max()) {
    // If a cell id is present, but not a cell locate the cell from map
    if (!cell_ptr()) assign_cell(map_cells[cell_id()]);
    if (compute_reference_location()) return;

    // Check if discontinuity point is in any of its nearest neighbours
    const auto neighbours = map_cells[cell_id()]->neighbours();
    Eigen::Matrix<double, Tdim, 1> xi;
    for (auto neighbour : neighbours) {
      if (map_cells[neighbour]->is_point_in_cell(coordinates_, &xi)) {
        assign_cell_xi(map_cells[neighbour], xi);
        return;
      }
    }
  }
#pragma omp parallel for schedule(runtime)
  for (auto citr = cells.cbegin(); citr != cells.cend(); ++citr) {
    // Check if particle is already found, if so don't run for other cells
    // Check if co-ordinates is within the cell, if true
    // add particle to cell
    Eigen::Matrix<double, Tdim, 1> xi;
    if ((*citr)->is_point_in_cell(coordinates(), &xi)) {
      assign_cell_xi(*citr, xi);
    }
  }
}

// Compute updated position of the particle
template <unsigned Tdim>
void mpm::discontinuity_point<Tdim>::compute_updated_position(
 const double dt) noexcept {
  // Check if point has a valid cell ptr
  if (cell_ == nullptr) return;
  // Get interpolated nodal velocity
  Eigen::Matrix<double, Tdim, 1> nodal_velocity =
      Eigen::Matrix<double, Tdim, 1>::Zero();
  const double tolerance = 1.E-16;
  unsigned int phase = 0;
  // need to do, points move with which side
  int move_direction = -1;
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    if (nodes_[i]->discontinuity_enrich()) {
      double nodal_mass =
          nodes_[i]->mass(phase) -
          nodes_[i]->discontinuity_property("mass_enrich", 1)(0, 0);
      if (nodal_mass < tolerance) continue;

      nodal_velocity +=
          shapefn_[i] *
          (nodes_[i]->momentum(phase) -
           nodes_[i]->discontinuity_property("momenta_enrich", 3)) /
          nodal_mass;
    } else {
      double nodal_mass = nodes_[i]->mass(phase);
      if (nodal_mass < tolerance) continue;
      nodal_velocity += shapefn_[i] * nodes_[i]->momentum(phase) / nodal_mass;
    }
  }
  // New position  current position + velocity * dt
  this->coordinates_ += nodal_velocity * dt;
}

// Compute updated position of the particle
template <unsigned Tdim>
void mpm::discontinuity_point<Tdim>::compute_shapefn() noexcept {
  // Check if point has a valid cell ptr
  if (cell_ == nullptr) return;
  // Get element ptr of a cell
  const auto element = cell_->element_ptr();

  // Zero matrix
  Eigen::Matrix<double, Tdim, 1> zero = Eigen::Matrix<double, Tdim, 1>::Zero();

  // Compute shape function of the point
  //! Size of particle in natural coordinates
  Eigen::Matrix<double, 1, Tdim> natural_size_;
  natural_size_.setZero();
  shapefn_ = element->shapefn(this->xi_, natural_size_, zero);
}
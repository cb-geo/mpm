// Constructor with id and number of nodes
//! \param[in] id Global cell id
//! \param[in] nnodes Number of nodes per cell
//! \tparam Tdim Dimension
template <unsigned Tdim>
mpm::Cell<Tdim>::Cell(Index id, unsigned nnodes) : id_{id}, nnodes_{nnodes} {
  // Check if the dimension is between 1 & 3
  static_assert((Tdim >= 1 && Tdim <= 3), "Invalid global dimension");
  try {
    // Number of nodes should be greater than dimensions
    if (!(nnodes > Tdim)) {
      throw std::runtime_error(
          "Specified number of nodes for a cell is too low");
    }
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
    std::abort();
  }
}

// Constructor with id, coordinates and shapefn
//! \param[in] id Global cell id
//! \param[in] nnodes Number of nodes per cell
//! \param[in] shapefnptr Pointer to a shape function
//! \tparam Tdim Dimension
template <unsigned Tdim>
mpm::Cell<Tdim>::Cell(Index id, unsigned nnodes,
                      const std::shared_ptr<ShapeFn<Tdim>>& shapefnptr)
    : mpm::Cell<Tdim>::Cell(id, nnodes) {

  try {
    if (shapefnptr->nfunctions() >= this->nnodes_) {
      shapefn_ = shapefnptr;
    } else {
      throw std::runtime_error(
          "Specified number of shape functions is not defined");
    }
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
    std::abort();
  }
}

// Assign a shape function to cell
//! \param[in] shapefnptr Pointer to a shape function
//! \tparam Tdim Dimension
template <unsigned Tdim>
bool mpm::Cell<Tdim>::shapefn(
    const std::shared_ptr<ShapeFn<Tdim>>& shapefnptr) {
  bool status = false;
  try {
    if (shapefnptr->nfunctions() >= this->nnodes_) {
      shapefn_ = shapefnptr;
      status = true;
    } else {
      throw std::runtime_error(
          "Specified number of shape functions is not defined");
    }
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
  }
  return status;
}

//! Add a node pointer
//! \param[in] local_id local id of the node
//! \param[in] ptr A shared pointer
//! \retval insertion_status Return the successful addition of a node
//! \tparam Tdim Dimension
template <unsigned Tdim>
bool mpm::Cell<Tdim>::add_node(
    unsigned local_id, const std::shared_ptr<mpm::NodeBase<Tdim>>& node_ptr) {
  bool insertion_status = false;
  try {
    // If number of node ptrs in a cell is less than the maximum number of nodes
    // per cell
    // The local id should be between 0 and maximum number of nodes
    if (nodes_.size() < this->nnodes_ &&
        (local_id >= 0 && local_id < this->nnodes_)) {
      insertion_status = nodes_.insert(local_id, node_ptr);
    } else {
      throw std::runtime_error(
          "Number nodes in a cell exceeds the maximum allowed per cell");
    }
  } catch (std::exception& exception) {
    std::cerr << exception.what() << "\n";
  }
  return insertion_status;
}

//! Add a neighbour cell
//! \param[in] local_id local id of the cell
//! \param[in] ptr A shared pointer to cell
//! \retval insertion_status Return the successful addition of a node
//! \tparam Tdim Dimension
template <unsigned Tdim>
bool mpm::Cell<Tdim>::add_neighbour(
    unsigned local_id, const std::shared_ptr<mpm::Cell<Tdim>>& cell_ptr) {
  bool insertion_status = false;
  try {
    // If number of cell ptrs id is not the current cell id
    if (cell_ptr->id() != this->id()) {
      insertion_status = neighbour_cells_.insert(local_id, cell_ptr);
    } else {
      throw std::runtime_error("Invalid local id of a cell neighbour");
    }
  } catch (std::exception& exception) {
    std::cerr << exception.what() << "\n";
  }
  return insertion_status;
}

//! Add a particle id
//! \param[in] id Global id of a particle
//! \retval status Return the successful addition of a particle id
//! \tparam Tdim Dimension
template <unsigned Tdim>
bool mpm::Cell<Tdim>::add_particle_id(Index id) {
  bool status = false;
  // Check if it is found in the container
  auto itr = std::find(particles_.begin(), particles_.end(), id);

  if (itr == particles_.end()) {
    particles_.emplace_back(id);
    status = true;
  }

  return status;
}

//! Remove a particle id
//! \param[in] id Global id of a particle
//! \retval status Return the successful removal of a particle id
//! \tparam Tdim Dimension
template <unsigned Tdim>
void mpm::Cell<Tdim>::remove_particle_id(Index id) {
  particles_.erase(std::remove(particles_.begin(), particles_.end(), id),
                   particles_.end());
}

//! Assign mass to nodes
//! param[in] xi local coordinates of particle
//! param[in] pmass mass of particle
//! \tparam Tdim Dimension
template <unsigned Tdim>
void mpm::Cell<Tdim>::assign_mass_to_nodes(const VectorDim& xi, const Eigen::VectorXd& pmass) {
Eigen::MatrixXd shapefns = shapefn_->shapefn(xi);
}

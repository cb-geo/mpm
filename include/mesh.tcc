// Constructor with id
//! \param[in] id Global mesh id
//! \tparam Tdim Dimension
template <unsigned Tdim>
mpm::Mesh<Tdim>::Mesh(unsigned id) : id_{id} {
  // Check if the dimension is between 1 & 3
  static_assert((Tdim >= 1 && Tdim <= 3), "Invalid global dimension");
  particles_.clear();
}

//! Add a neighbour mesh
//! \param[in] local_id local id of the mesh
//! \param[in] ptr A shared pointer to mesh
//! \retval insertion_status Return the successful addition of a node
//! \tparam Tdim Dimension
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::add_neighbour(
    unsigned local_id, const std::shared_ptr<mpm::Mesh<Tdim>>& mesh) {
  bool insertion_status = false;
  try {
    // If number of mesh ptrs id is not the current mesh id
    if (mesh->id() != this->id()) {
      insertion_status = neighbour_meshes_.insert(local_id, mesh);
    } else {
      throw std::runtime_error("Invalid local id of a mesh neighbour");
    }
  } catch (std::exception& exception) {
    std::cerr << exception.what() << "\n";
  }
  return insertion_status;
}

//! Add a particle
//! \param[in] particle A shared pointer to particle
//! \retval insertion_status Return the successful addition of a particle
//! \tparam Tdim Dimension
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::add_particle(
    const std::shared_ptr<mpm::ParticleBase<Tdim>>& particle) {
  bool insertion_status = particles_.add(particle);
  return insertion_status;
}

//! Remove a particle
//! \param[in] particle A shared pointer to particle
//! \retval insertion_status Return the successful addition of a particle
//! \tparam Tdim Dimension
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::remove_particle(
    const std::shared_ptr<mpm::ParticleBase<Tdim>>& particle) {
  // Remove a particle if found in the container
  bool status = particles_.remove(particle);
  return status;
}

//! Add a node
//! \param[in] node A shared pointer to node
//! \retval insertion_status Return the successful addition of a node
//! \tparam Tdim Dimension
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::add_node(
    const std::shared_ptr<mpm::NodeBase<Tdim>>& node) {
  bool insertion_status = nodes_.add(node);
  return insertion_status;
}

//! Remove a node
//! \param[in] node A shared pointer to node
//! \retval insertion_status Return the successful addition of a node
//! \tparam Tdim Dimension
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::remove_node(
    const std::shared_ptr<mpm::NodeBase<Tdim>>& node) {
  // Remove a node if found in the container
  bool status = nodes_.remove(node);
  return status;
}

//! Add a cell
//! \param[in] cell A shared pointer to cell
//! \retval insertion_status Return the successful addition of a cell
//! \tparam Tdim Dimension
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::add_cell(const std::shared_ptr<mpm::Cell<Tdim>>& cell) {
  bool insertion_status = cells_.add(cell);
  return insertion_status;
}

//! Remove a cell
//! \param[in] cell A shared pointer to cell
//! \retval insertion_status Return the successful addition of a cell
//! \tparam Tdim Dimension
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::remove_cell(
    const std::shared_ptr<mpm::Cell<Tdim>>& cell) {
  // Remove a cell if found in the container
  bool status = cells_.remove(cell);
  return status;
}

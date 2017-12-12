// Constructor with id
//! \param[in] id Global mesh id
//! \tparam Tdim Dimension
template <unsigned Tdim>
mpm::Mesh<Tdim>::Mesh(unsigned id) : id_{id} {
  // Check if the dimension is between 1 & 3
  static_assert((Tdim >= 1 && Tdim <= 3), "Invalid global dimension");
}

//! Add a neighbour mesh
//! \param[in] local_id local id of the mesh
//! \param[in] ptr A shared pointer to mesh
//! \retval insertion_status Return the successful addition of a node
//! \tparam Tdim Dimension
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::add_neighbour(
    unsigned local_id, const std::shared_ptr<mpm::Mesh<Tdim>>& mesh_ptr) {
  bool insertion_status = false;
  try {
    // If number of mesh ptrs id is not the current mesh id
    if (mesh_ptr->id() != this->id()) {
      insertion_status = neighbour_meshes_.insert(local_id, mesh_ptr);
    } else {
      throw std::runtime_error("Invalid local id of a mesh neighbour");
    }
  } catch (std::exception& exception) {
    std::cerr << exception.what() << "\n";
  }
  return insertion_status;
}

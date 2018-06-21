// Constructor with id
template <unsigned Tdim>
mpm::Mesh<Tdim>::Mesh(unsigned id) : id_{id} {
  // Check if the dimension is between 1 & 3
  static_assert((Tdim >= 1 && Tdim <= 3), "Invalid global dimension");
  particles_.clear();
}

//! Create nodes from coordinates
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::create_nodes(mpm::Index gnid, const std::string& ntype,
                                   const std::vector<VectorDim>& coordinates) {
  bool status = false;
  // Global node id
  try {
    // Check if nodal coordinates is not empty
    if (!coordinates.empty()) {
      for (const auto& node_coordinates : coordinates) {
        // Add node to mesh and check
        bool insert_status = this->nodes_.add(
            Factory<mpm::NodeBase<Tdim>, mpm::Index,
                    const Eigen::Matrix<double, Tdim, 1>&>::instance()
                ->create(ntype, std::move(gnid), node_coordinates));

        if (insert_status) {
          // Increament node id
          ++gnid;
        } else {
          // When addition of node fails
          throw std::runtime_error("Addition of node to mesh failed!");
        }
      }
      // When successful
      status = true;
    } else {
      // If the coordinates vector is empty
      throw std::runtime_error("List of coordinates is empty");
    }
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
  }
  return status;
}

//! Add a node to the mesh
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::add_node(
    const std::shared_ptr<mpm::NodeBase<Tdim>>& node) {
  bool insertion_status = nodes_.add(node);
  return insertion_status;
}

//! Remove a node from the mesh
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::remove_node(
    const std::shared_ptr<mpm::NodeBase<Tdim>>& node) {
  // Remove a node if found in the container
  bool status = nodes_.remove(node);
  return status;
}

//! Iterate over nodes
template <unsigned Tdim>
template <typename Toper>
void mpm::Mesh<Tdim>::iterate_over_nodes(Toper oper) {
  tbb::parallel_for_each(nodes_.cbegin(), nodes_.cend(), oper);
}

//! Add a cell to the mesh
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::add_cell(const std::shared_ptr<mpm::Cell<Tdim>>& cell) {
  bool insertion_status = cells_.add(cell);
  return insertion_status;
}

//! Remove a cell from the mesh
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::remove_cell(
    const std::shared_ptr<mpm::Cell<Tdim>>& cell) {
  // Remove a cell if found in the container
  bool status = cells_.remove(cell);
  return status;
}

//! Iterate over cells
template <unsigned Tdim>
template <typename Toper>
void mpm::Mesh<Tdim>::iterate_over_cells(Toper oper) {
  tbb::parallel_for_each(cells_.cbegin(), cells_.cend(), oper);
}

//! Add a particle pointer to the mesh
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::add_particle(
    const std::shared_ptr<mpm::ParticleBase<Tdim>>& particle) {
  bool insertion_status = particles_.add(particle);
  return insertion_status;
}

//! Remove a particle pointer from the mesh
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::remove_particle(
    const std::shared_ptr<mpm::ParticleBase<Tdim>>& particle) {
  // Remove a particle if found in the container
  bool status = particles_.remove(particle);
  return status;
}

//! Locate particles in a cell
template <unsigned Tdim>
void mpm::Mesh<Tdim>::locate_particles_mesh() {
  // Iterate through each particle and
  for (auto pitr = particles_.cbegin(); pitr != particles_.cend(); ++pitr) {
    for (auto citr = cells_.cbegin(); citr != cells_.cend(); ++citr) {
      // Check if co-ordinates lie within the cell, if true add particle to cell
      if ((*citr)->point_in_cell((*pitr)->coordinates()))
        (*pitr)->assign_cell(*citr);
    }
  }
}

//! Iterate over particles
template <unsigned Tdim>
template <typename Toper>
void mpm::Mesh<Tdim>::iterate_over_particles(Toper oper) {
  tbb::parallel_for_each(particles_.cbegin(), particles_.cend(), oper);
}

//! Add a neighbour mesh, using the local id of the mesh and a mesh pointer
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::add_neighbour(
    unsigned local_id, const std::shared_ptr<mpm::Mesh<Tdim>>& mesh) {
  bool insertion_status = false;
  try {
    // If the mesh id is not the current mesh id
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

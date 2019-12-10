// Constructor with id
template <unsigned Tdim>
mpm::Mesh<Tdim>::Mesh(unsigned id, bool isoparametric)
    : id_{id}, isoparametric_{isoparametric} {
  // Check if the dimension is between 1 & 3
  static_assert((Tdim >= 1 && Tdim <= 3), "Invalid global dimension");
  //! Logger
  std::string logger = "mesh::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);

  particles_.clear();
}

//! Create nodes from coordinates
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::create_nodes(mpm::Index gnid,
                                   const std::string& node_type,
                                   const std::vector<VectorDim>& coordinates,
                                   bool check_duplicates) {
  bool status = true;
  try {
    // Check if nodal coordinates is empty
    if (coordinates.empty())
      throw std::runtime_error("List of coordinates is empty");
    // Iterate over all coordinates
    for (const auto& node_coordinates : coordinates) {
      // Add node to mesh and check
      bool insert_status = this->add_node(
          // Create a node of particular
          Factory<mpm::NodeBase<Tdim>, mpm::Index,
                  const Eigen::Matrix<double, Tdim, 1>&>::instance()
              ->create(node_type, static_cast<mpm::Index>(gnid),
                       node_coordinates),
          check_duplicates);

      // Increment node id
      if (insert_status) ++gnid;
      // When addition of node fails
      else
        throw std::runtime_error("Addition of node to mesh failed!");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Add a node to the mesh
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::add_node(const std::shared_ptr<mpm::NodeBase<Tdim>>& node,
                               bool check_duplicates) {
  bool insertion_status = nodes_.add(node, check_duplicates);
  // Add node to map
  if (insertion_status) map_nodes_.insert(node->id(), node);
  return insertion_status;
}

//! Remove a node from the mesh
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::remove_node(
    const std::shared_ptr<mpm::NodeBase<Tdim>>& node) {
  const mpm::Index id = node->id();
  // Remove a node if found in the container
  return (nodes_.remove(node) && map_nodes_.remove(id));
}

//! Iterate over nodes
template <unsigned Tdim>
template <typename Toper>
void mpm::Mesh<Tdim>::iterate_over_nodes(Toper oper) {
  tbb::parallel_for_each(nodes_.cbegin(), nodes_.cend(), oper);
}

//! Iterate over nodes
template <unsigned Tdim>
template <typename Toper, typename Tpred>
void mpm::Mesh<Tdim>::iterate_over_nodes_predicate(Toper oper, Tpred pred) {
  tbb::parallel_for_each(nodes_.cbegin(), nodes_.cend(),
                         [=](std::shared_ptr<mpm::NodeBase<Tdim>> node) {
                           if (pred(node)) oper(node);
                         });
}

//! Create a list of active nodes in mesh
template <unsigned Tdim>
void mpm::Mesh<Tdim>::find_active_nodes() {
  // Clear existing list of active nodes
  this->active_nodes_.clear();

  for (auto nitr = nodes_.cbegin(); nitr != nodes_.cend(); ++nitr)
    if ((*nitr)->status()) this->active_nodes_.add(*nitr);
}

//! Iterate over active nodes
template <unsigned Tdim>
template <typename Toper>
void mpm::Mesh<Tdim>::iterate_over_active_nodes(Toper oper) {
  tbb::parallel_for_each(active_nodes_.cbegin(), active_nodes_.cend(), oper);
}

#ifdef USE_MPI
//! All reduce over nodal scalar property
template <unsigned Tdim>
template <typename Tgetfunctor, typename Tsetfunctor>
void mpm::Mesh<Tdim>::allreduce_nodal_scalar_property(Tgetfunctor getter,
                                                      Tsetfunctor setter) {
  // Create vector of nodal scalars
  mpm::Index nnodes = this->nodes_.size();
  std::vector<double> prop_get(nnodes), prop_set(nnodes);

  tbb::parallel_for_each(
      nodes_.cbegin(), nodes_.cend(),
      [=, &prop_get](std::shared_ptr<mpm::NodeBase<Tdim>> node) {
        prop_get.at(node->id()) = getter(node);
      });

  MPI_Allreduce(prop_get.data(), prop_set.data(), nnodes, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);

  tbb::parallel_for_each(
      nodes_.cbegin(), nodes_.cend(),
      [=, &prop_set](std::shared_ptr<mpm::NodeBase<Tdim>> node) {
        setter(node, prop_set.at(node->id()));
      });
}
#endif

#ifdef USE_MPI
//! All reduce over nodal vector property
template <unsigned Tdim>
template <typename Tgetfunctor, typename Tsetfunctor>
void mpm::Mesh<Tdim>::allreduce_nodal_vector_property(Tgetfunctor getter,
                                                      Tsetfunctor setter) {
  // Create vector of nodal vectors
  mpm::Index nnodes = this->nodes_.size();
  std::vector<Eigen::Matrix<double, Tdim, 1>> prop_get(nnodes),
      prop_set(nnodes);

  tbb::parallel_for_each(
      nodes_.cbegin(), nodes_.cend(),
      [=, &prop_get](std::shared_ptr<mpm::NodeBase<Tdim>> node) {
        prop_get.at(node->id()) = getter(node);
      });

  MPI_Allreduce(prop_get.data(), prop_set.data(), nnodes * Tdim, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);

  tbb::parallel_for_each(
      nodes_.cbegin(), nodes_.cend(),
      [=, &prop_set](std::shared_ptr<mpm::NodeBase<Tdim>> node) {
        setter(node, prop_set.at(node->id()));
      });
}
#endif

//! Create cells from node lists
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::create_cells(
    mpm::Index gcid, const std::shared_ptr<mpm::Element<Tdim>>& element,
    const std::vector<std::vector<mpm::Index>>& cells, bool check_duplicates) {
  bool status = true;
  try {
    // Check if nodes in cell list is not empty
    if (cells.empty())
      throw std::runtime_error("List of nodes of cells is empty");

    for (const auto& nodes : cells) {
      // Create cell with element
      auto cell = std::make_shared<mpm::Cell<Tdim>>(gcid, nodes.size(), element,
                                                    this->isoparametric_);

      // Cell local node id
      unsigned local_nid = 0;
      // For nodeids in a given cell
      for (auto nid : nodes) {
        cell->add_node(local_nid, map_nodes_[nid]);
        ++local_nid;
      }

      // Add cell to mesh
      bool insert_cell = false;
      // Check if cell has all nodes before inserting to mesh
      if (cell->nnodes() == nodes.size()) {
        // Initialise cell before insertion
        cell->initialise();
        // If cell is initialised insert to mesh
        if (cell->is_initialised())
          insert_cell = this->add_cell(cell, check_duplicates);
      } else
        throw std::runtime_error("Invalid node ids for cell!");

      // Increment global cell id
      if (insert_cell) ++gcid;
      // When addition of cell fails
      else
        throw std::runtime_error("Addition of cell to mesh failed!");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Add a cell to the mesh
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::add_cell(const std::shared_ptr<mpm::Cell<Tdim>>& cell,
                               bool check_duplicates) {
  bool insertion_status = cells_.add(cell, check_duplicates);
  // Add cell to map
  if (insertion_status) map_cells_.insert(cell->id(), cell);
  return insertion_status;
}

//! Remove a cell from the mesh
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::remove_cell(
    const std::shared_ptr<mpm::Cell<Tdim>>& cell) {
  const mpm::Index id = cell->id();
  // Remove a cell if found in the container
  return (cells_.remove(cell) && map_cells_.remove(id));
}

//! Iterate over cells
template <unsigned Tdim>
template <typename Toper>
void mpm::Mesh<Tdim>::iterate_over_cells(Toper oper) {
  tbb::parallel_for_each(cells_.cbegin(), cells_.cend(), oper);
}

//! Create cells from node lists
template <unsigned Tdim>
void mpm::Mesh<Tdim>::compute_cell_neighbours() {
  // Initialize and compute node cell map
  tsl::robin_map<mpm::Index, std::set<mpm::Index>> node_cell_map;
  for (auto citr = cells_.cbegin(); citr != cells_.cend(); ++citr) {
    // Get cell id and nodes id
    auto cell_id = (*citr)->id();
    const auto nodes_id_list = (*citr)->nodes_id();
    // Populate node_cell_map with the node_id and multiple cell_id
    for (const auto& id : nodes_id_list) node_cell_map[id].insert(cell_id);
  }

  // Assign neighbour to cells
  for (auto citr = cells_.cbegin(); citr != cells_.cend(); ++citr) {
    // Initiate set of neighbouring cells
    std::set<mpm::Index> neighbouring_cell_sets;

    // Loop over the current cell nodes and add ids of the initiated set
    const auto nodes_id_list = (*citr)->nodes_id();
    for (const auto& id : nodes_id_list)
      neighbouring_cell_sets.insert(node_cell_map[id].begin(),
                                    node_cell_map[id].end());

    for (const auto& neighbour_id : neighbouring_cell_sets)
      if (neighbour_id != (*citr)->id())
        map_cells_[(*citr)->id()]->add_neighbour(neighbour_id);
  }
}

//! Create cells from node lists
template <unsigned Tdim>
std::vector<Eigen::Matrix<double, Tdim, 1>>
    mpm::Mesh<Tdim>::generate_material_points(unsigned nquadratures) {
  std::vector<VectorDim> points;
  try {
    if (cells_.size() > 0) {
      points.reserve(cells_.size() * std::pow(nquadratures, Tdim));

      // Generate points
      for (auto citr = cells_.cbegin(); citr != cells_.cend(); ++citr) {
        (*citr)->assign_quadrature(nquadratures);
        const auto cpoints = (*citr)->generate_points();
        points.insert(std::end(points), std::begin(cpoints), std::end(cpoints));
      }
      console_->info(
          "Generate points:\n# of cells: {}\nExpected # of points: {}\n"
          "# of points generated: {}",
          cells_.size(), cells_.size() * std::pow(nquadratures, Tdim),
          points.size());
    } else
      throw std::runtime_error("No cells are found in the mesh!");
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    points.clear();
  }
  return points;
}

//! Create particles from coordinates
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::create_particles(
    const std::vector<mpm::Index>& gp_ids, const std::string& particle_type,
    const std::vector<VectorDim>& coordinates, bool check_duplicates) {
  bool status = true;
  try {
    unsigned gpid = 0;
    // Check if particle coordinates is empty
    if (coordinates.empty())
      throw std::runtime_error("List of coordinates is empty");
    // Iterate over particle coordinates
    for (const auto& particle_coordinates : coordinates) {
      // Add particle to mesh and check
      bool insert_status = this->add_particle(
          Factory<mpm::ParticleBase<Tdim>, mpm::Index,
                  const Eigen::Matrix<double, Tdim, 1>&>::instance()
              ->create(particle_type, static_cast<mpm::Index>(gp_ids.at(gpid)),
                       particle_coordinates),
          check_duplicates);

      // Increment particle id
      if (insert_status) ++gpid;
      // When addition of particle fails
      else
        throw std::runtime_error("Addition of particle to mesh failed!");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Add a particle pointer to the mesh
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::add_particle(
    const std::shared_ptr<mpm::ParticleBase<Tdim>>& particle, bool checks) {
  bool status = false;
  try {
    if (checks) {
      // Add only if particle can be located in any cell of the mesh
      if (this->locate_particle_cells(particle)) {
        status = particles_.add(particle, checks);
        particles_cell_ids_.insert(std::pair<mpm::Index, mpm::Index>(
            particle->id(), particle->cell_id()));
        map_particles_.insert(particle->id(), particle);
      } else {
        throw std::runtime_error("Particle not found in mesh");
      }
    } else {
      status = particles_.add(particle, checks);
      particles_cell_ids_.insert(std::pair<mpm::Index, mpm::Index>(
          particle->id(), particle->cell_id()));
      map_particles_.insert(particle->id(), particle);
    }
    if (!status) throw std::runtime_error("Particle addition failed");
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Remove a particle pointer from the mesh
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::remove_particle(
    const std::shared_ptr<mpm::ParticleBase<Tdim>>& particle) {
  const mpm::Index id = particle->id();
  // Remove associated cell for the particle
  map_particles_[id]->remove_cell();
  // Remove a particle if found in the container and map
  return (particles_.remove(particle) && map_particles_.remove(id));
}

//! Remove a particle by id
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::remove_particle_by_id(mpm::Index id) {
  // Remove associated cell for the particle
  map_particles_[id]->remove_cell();
  bool result = particles_.remove(map_particles_[id]);
  return (result && map_particles_.remove(id));
}

//! Remove all particles in a cell given cell id
template <unsigned Tdim>
void mpm::Mesh<Tdim>::remove_all_nonrank_particles(unsigned rank) {
  // Remove associated cell for the particle
  for (auto citr = this->cells_.cbegin(); citr != this->cells_.cend(); ++citr) {
    // If cell is non empty
    if ((*citr)->particles().size() != 0 && (*citr)->rank() != rank) {
      auto particle_ids = (*citr)->particles();
      for (auto& id : particle_ids) {
        map_particles_[id]->remove_cell();
        particles_.remove(map_particles_[id]);
        map_particles_.remove(id);
      }
      (*citr)->clear_particle_ids();
    }
  }
}

//! Locate particles in a cell
template <unsigned Tdim>
std::vector<std::shared_ptr<mpm::ParticleBase<Tdim>>>
    mpm::Mesh<Tdim>::locate_particles_mesh() {

  std::vector<std::shared_ptr<mpm::ParticleBase<Tdim>>> particles;

  std::for_each(particles_.cbegin(), particles_.cend(),
                [=, &particles](
                    const std::shared_ptr<mpm::ParticleBase<Tdim>>& particle) {
                  // If particle is not found in mesh add to a list of particles
                  if (!this->locate_particle_cells(particle))
                    particles.emplace_back(particle);
                });

  return particles;
}

//! Locate particles in a cell
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::locate_particle_cells(
    const std::shared_ptr<mpm::ParticleBase<Tdim>>& particle) {
  // Check the current cell if it is not invalid
  if (particle->cell_id() != std::numeric_limits<mpm::Index>::max()) {
    // If a cell id is present, but not a cell locate the cell from map
    if (!particle->cell_ptr())
      particle->assign_cell(map_cells_[particle->cell_id()]);
    if (particle->compute_reference_location()) return true;

    // Check if material point is in any of its nearest neighbours
    const auto neighbours = map_cells_[particle->cell_id()]->neighbours();
    Eigen::Matrix<double, Tdim, 1> xi;
    Eigen::Matrix<double, Tdim, 1> coordinates = particle->coordinates();
    for (auto neighbour : neighbours) {
      if (map_cells_[neighbour]->is_point_in_cell(coordinates, &xi)) {
        particle->assign_cell_xi(map_cells_[neighbour], xi);
        return true;
      }
    }
  }

  bool status = false;
  tbb::parallel_for_each(
      cells_.cbegin(), cells_.cend(),
      [=, &status](const std::shared_ptr<mpm::Cell<Tdim>>& cell) {
        // Check if particle is already found, if so don't run for other cells
        // Check if co-ordinates is within the cell, if true
        // add particle to cell
        Eigen::Matrix<double, Tdim, 1> xi;
        if (!status && cell->is_point_in_cell(particle->coordinates(), &xi)) {
          particle->assign_cell_xi(cell, xi);
          status = true;
        }
      });

  return status;
}

//! Iterate over particles
template <unsigned Tdim>
template <typename Toper>
void mpm::Mesh<Tdim>::iterate_over_particles(Toper oper) {
  tbb::parallel_for_each(particles_.cbegin(), particles_.cend(), oper);
}

//! Iterate over particle set
template <unsigned Tdim>
template <typename Toper>
void mpm::Mesh<Tdim>::iterate_over_particle_set(unsigned set_id, Toper oper) {
  tbb::parallel_for_each(this->particle_sets_.at(set_id).cbegin(),
                         this->particle_sets_.at(set_id).cend(), oper);
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
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return insertion_status;
}

//! Return particle coordinates
template <unsigned Tdim>
std::vector<Eigen::Matrix<double, 3, 1>>
    mpm::Mesh<Tdim>::particle_coordinates() {
  std::vector<Eigen::Matrix<double, 3, 1>> particle_coordinates;
  for (auto pitr = particles_.cbegin(); pitr != particles_.cend(); ++pitr) {
    Eigen::Vector3d coordinates;
    coordinates.setZero();
    auto pcoords = (*pitr)->coordinates();
    // Fill coordinates to the size of dimensions
    for (unsigned i = 0; i < Tdim; ++i) coordinates(i) = pcoords(i);
    particle_coordinates.emplace_back(coordinates);
  }
  return particle_coordinates;
}

//! Return particle vector data
template <unsigned Tdim>
std::vector<Eigen::Matrix<double, 3, 1>> mpm::Mesh<Tdim>::particles_vector_data(
    const std::string& attribute, unsigned phase) {
  std::vector<Eigen::Matrix<double, 3, 1>> vector_data;
  try {
    // Iterate over particles
    for (auto pitr = particles_.cbegin(); pitr != particles_.cend(); ++pitr) {
      Eigen::Vector3d data;
      data.setZero();
      auto pdata = (*pitr)->vector_data(attribute);
      // Fill stresses to the size of dimensions
      for (unsigned i = 0; i < Tdim; ++i) data(i) = pdata(i);

      // Add to a vector of data
      vector_data.emplace_back(data);
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {} {}\n", __FILE__, __LINE__, exception.what(),
                    attribute);
    vector_data.clear();
  }
  return vector_data;
}

//! Assign velocity constraints to nodes
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::assign_velocity_constraints(
    const std::vector<std::tuple<mpm::Index, unsigned, double>>&
        velocity_constraints) {
  bool status = false;
  try {
    if (!nodes_.size())
      throw std::runtime_error(
          "No nodes have been assigned in mesh, cannot assign velocity "
          "constraints");

    for (const auto& velocity_constraint : velocity_constraints) {
      // Node id
      mpm::Index nid = std::get<0>(velocity_constraint);
      // Direction
      unsigned dir = std::get<1>(velocity_constraint);
      // Velocity
      double velocity = std::get<2>(velocity_constraint);

      // Apply constraint
      status = map_nodes_[nid]->assign_velocity_constraint(dir, velocity);

      if (!status)
        throw std::runtime_error("Node or velocity constraint is invalid");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assign friction constraints to nodes
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::assign_friction_constraints(
    const std::vector<std::tuple<mpm::Index, unsigned, int, double>>&
        friction_constraints) {
  bool status = false;
  try {
    if (!nodes_.size())
      throw std::runtime_error(
          "No nodes have been assigned in mesh, cannot assign friction "
          "constraints");

    for (const auto& friction_constraint : friction_constraints) {
      // Node id
      mpm::Index nid = std::get<0>(friction_constraint);
      // Direction
      unsigned dir = std::get<1>(friction_constraint);
      // Sign
      int sign = std::get<2>(friction_constraint);
      // Friction
      double friction = std::get<3>(friction_constraint);

      // Apply constraint
      status = map_nodes_[nid]->assign_friction_constraint(dir, sign, friction);

      if (!status)
        throw std::runtime_error("Node or friction constraint is invalid");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assign pressure constraints to nodes
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::assign_pressure_constraints(
    const unsigned phase,
    const std::vector<std::tuple<mpm::Index, double>>& pressure_constraints) {
  bool status = false;
  try {
    if (!nodes_.size())
      throw std::runtime_error(
          "No nodes have been assigned in mesh, cannot assign pressure "
          "constraints");

    for (const auto& pressure_constraint : pressure_constraints) {
      // Node id
      mpm::Index nid = std::get<0>(pressure_constraint);
      // Pressure
      double pressure = std::get<1>(pressure_constraint);

      // Apply constraint
      status = map_nodes_[nid]->assign_pressure_constraint(phase, pressure);

      if (!status)
        throw std::runtime_error("Node or pressure constraint is invalid");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assign particles volumes
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::assign_particles_volumes(
    const std::vector<std::tuple<mpm::Index, double>>& particle_volumes) {
  bool status = true;
  const unsigned phase = 0;
  try {
    if (!particles_.size())
      throw std::runtime_error(
          "No particles have been assigned in mesh, cannot assign volume");

    for (const auto& particle_volume : particle_volumes) {
      // Particle id
      mpm::Index pid = std::get<0>(particle_volume);
      // Volume
      double volume = std::get<1>(particle_volume);

      if (map_particles_.find(pid) != map_particles_.end())
        status = map_particles_[pid]->assign_volume(volume);

      if (!status)
        throw std::runtime_error("Cannot assign invalid particle volume");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Compute and assign rotation matrix to nodes
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::compute_nodal_rotation_matrices(
    const std::map<mpm::Index, Eigen::Matrix<double, Tdim, 1>>& euler_angles) {
  bool status = false;
  try {
    if (!nodes_.size())
      throw std::runtime_error(
          "No nodes have been assigned in mesh, cannot assign rotation "
          "matrix");

    // Loop through nodal_euler_angles of different nodes
    for (const auto& nodal_euler_angles : euler_angles) {
      // Node id
      mpm::Index nid = nodal_euler_angles.first;
      // Euler angles
      Eigen::Matrix<double, Tdim, 1> angles = nodal_euler_angles.second;
      // Compute rotation matrix
      const auto rotation_matrix = mpm::geometry::rotation_matrix(angles);

      // Apply rotation matrix to nodes
      map_nodes_[nid]->assign_rotation_matrix(rotation_matrix);
      status = true;
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assign particle tractions
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::assign_particles_tractions(
    const std::vector<std::tuple<mpm::Index, unsigned, double>>&
        particle_tractions) {
  bool status = true;
  // TODO: Remove phase
  const unsigned phase = 0;
  try {
    if (!particles_.size())
      throw std::runtime_error(
          "No particles have been assigned in mesh, cannot assign traction");
    for (const auto& particle_traction : particle_tractions) {
      // Particle id
      mpm::Index pid = std::get<0>(particle_traction);
      // Direction
      unsigned dir = std::get<1>(particle_traction);
      // Traction
      double traction = std::get<2>(particle_traction);

      if (map_particles_.find(pid) != map_particles_.end())
        status = map_particles_[pid]->assign_traction(dir, traction);

      if (!status) throw std::runtime_error("Traction is invalid for particle");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assign particles velocity constraints
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::assign_particles_velocity_constraints(
    const std::vector<std::tuple<mpm::Index, unsigned, double>>&
        particle_velocity_constraints) {
  bool status = true;
  // TODO: Remove phase
  const unsigned phase = 0;
  try {
    if (!particles_.size())
      throw std::runtime_error(
          "No particles have been assigned in mesh, cannot assign velocity");
    for (const auto& particle_velocity_constraint :
         particle_velocity_constraints) {
      // Particle id
      mpm::Index pid = std::get<0>(particle_velocity_constraint);
      // Direction
      unsigned dir = std::get<1>(particle_velocity_constraint);
      // Velocity
      double velocity = std::get<2>(particle_velocity_constraint);

      if (map_particles_.find(pid) != map_particles_.end())
        status = map_particles_[pid]->assign_particle_velocity_constraint(
            dir, velocity);

      if (!status)
        throw std::runtime_error("Velocity constraint is invalid for particle");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assign particles pressure constraints
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::assign_particles_pressure_constraints(
    const unsigned phase, const std::vector<std::tuple<mpm::Index, double>>&
                              particle_pressure_constraints) {
  bool status = true;

  try {
    if (!particles_.size())
      throw std::runtime_error(
          "No particles have been assigned in mesh, cannot assign pressure");
    for (const auto& particle_pressure_constraint :
         particle_pressure_constraints) {
      // Particle id
      mpm::Index pid = std::get<0>(particle_pressure_constraint);
      // Pressure
      double pressure = std::get<1>(particle_pressure_constraint);

      if (map_particles_.find(pid) != map_particles_.end())
        status = map_particles_[pid]->assign_particle_pressure_constraint(
            phase, pressure);

      if (!status)
        throw std::runtime_error("Pressure constraint is invalid for particle");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assign node tractions
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::assign_nodal_tractions(
    const std::vector<std::tuple<mpm::Index, unsigned, double>>&
        node_tractions) {
  bool status = true;
  // TODO: Remove phase
  const unsigned phase = 0;
  try {
    if (!nodes_.size())
      throw std::runtime_error(
          "No nodes have been assigned in mesh, cannot assign traction");
    for (const auto& node_traction : node_tractions) {
      // Node id
      mpm::Index pid = std::get<0>(node_traction);
      // Direction
      unsigned dir = std::get<1>(node_traction);
      // Traction
      double traction = std::get<2>(node_traction);

      if (map_nodes_.find(pid) != map_nodes_.end())
        status = map_nodes_[pid]->assign_traction_force(phase, dir, traction);

      if (!status) throw std::runtime_error("Traction is invalid for node");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assign particle stresses
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::assign_particles_stresses(
    const std::vector<Eigen::Matrix<double, 6, 1>>& particle_stresses) {
  bool status = true;
  // TODO: Remove phase
  const unsigned phase = 0;
  try {
    if (!particles_.size())
      throw std::runtime_error(
          "No particles have been assigned in mesh, cannot assign stresses");

    if (particles_.size() != particle_stresses.size())
      throw std::runtime_error(
          "Number of particles in mesh and initial stresses don't match");

    unsigned i = 0;
    for (auto pitr = particles_.cbegin(); pitr != particles_.cend(); ++pitr) {
      (*pitr)->initial_stress(particle_stresses.at(i));
      ++i;
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assign particle pore pressures
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::assign_particles_pore_pressures(
    const unsigned pore_fluid,
    const std::vector<double>& particle_pore_pressure) {
  bool status = true;

  try {
    if (!particles_.size())
      throw std::runtime_error(
          "No particles have been assigned in mesh, cannot assign pore "
          "pressures");

    if (particles_.size() != particle_pore_pressure.size())
      throw std::runtime_error(
          "Number of particles in mesh and initial pore pressures don't match");

    unsigned i = 0;
    for (auto pitr = particles_.cbegin(); pitr != particles_.cend(); ++pitr) {
      (*pitr)->initial_pore_pressure(pore_fluid, particle_pore_pressure.at(i));
      ++i;
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assign particle cells
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::assign_particles_cells(
    const std::vector<std::array<mpm::Index, 2>>& particles_cells) {
  bool status = true;
  try {
    if (!particles_.size())
      throw std::runtime_error(
          "No particles have been assigned in mesh, cannot assign cells");
    for (const auto& particle_cell : particles_cells) {
      // Particle id
      mpm::Index pid = particle_cell[0];
      // Cell id
      mpm::Index cid = particle_cell[1];

      map_particles_[pid]->assign_cell_id(cid);
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Return particle cells
template <unsigned Tdim>
std::vector<std::array<mpm::Index, 2>> mpm::Mesh<Tdim>::particles_cells()
    const {
  std::vector<std::array<mpm::Index, 2>> particles_cells;
  try {
    if (!particles_.size())
      throw std::runtime_error(
          "No particles have been assigned in mesh, cannot write cells");
    for (auto pitr = particles_.cbegin(); pitr != particles_.cend(); ++pitr) {
      if ((*pitr)->cell_id() != std::numeric_limits<mpm::Index>::max())
        particles_cells.emplace_back(
            std::array<mpm::Index, 2>({(*pitr)->id(), (*pitr)->cell_id()}));
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    particles_cells.clear();
  }
  return particles_cells;
}

//! Assign velocity constraints to cells
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::assign_cell_velocity_constraints(
    const std::vector<std::tuple<mpm::Index, unsigned, unsigned, double>>&
        velocity_constraints) {
  bool status = false;
  try {
    // Check the number of nodes
    if (!nodes_.size()) {
      throw std::runtime_error(
          "No cells have been assigned in mesh, cannot assign velocity "
          "constraints");
    }
    // Loop through all the velocity constraints
    for (const auto& velocity_constraint : velocity_constraints) {
      // Cell id
      const mpm::Index cell_id = std::get<0>(velocity_constraint);
      // Face id
      const unsigned face_id = std::get<1>(velocity_constraint);
      // Direction of the local coordinate system of the face
      // Tdim = 2, Normal is y local axis, dir = 1
      // Tdim = 3, Normal is z local axis, dir = 2
      const unsigned dir = std::get<2>(velocity_constraint);
      // Velocity
      const double velocity = std::get<3>(velocity_constraint);

      // Apply constraint
      for (auto citr = cells_.cbegin(); citr != cells_.cend(); ++citr) {
        if ((*citr)->id() == cell_id) {
          status = (*citr)->assign_velocity_constraint(face_id, dir, velocity);
          break;
        }
      }

      if (!status)
        throw std::runtime_error(
            "Cell or face or velocity constraint is invalid");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Write particles to HDF5
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::write_particles_hdf5(unsigned phase,
                                           const std::string& filename) {
  const unsigned nparticles = this->nparticles();

  std::vector<HDF5Particle> particle_data;  // = new HDF5Particle[nparticles];
  particle_data.reserve(nparticles);

  for (auto pitr = particles_.cbegin(); pitr != particles_.cend(); ++pitr)
    particle_data.emplace_back((*pitr)->hdf5());

  // Calculate the size and the offsets of our struct members in memory
  const hsize_t NRECORDS = nparticles;

  const hsize_t NFIELDS = mpm::hdf5::particle::NFIELDS;

  hid_t string_type;
  hid_t file_id;
  hsize_t chunk_size = 10000;
  int* fill_data = NULL;
  int compress = 0;

  // Create a new file using default properties.
  file_id =
      H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  // make a table
  H5TBmake_table(
      "Table Title", file_id, "table", NFIELDS, NRECORDS,
      mpm::hdf5::particle::dst_size, mpm::hdf5::particle::field_names,
      mpm::hdf5::particle::dst_offset, mpm::hdf5::particle::field_type,
      chunk_size, fill_data, compress, particle_data.data());

  H5Fclose(file_id);
  return true;
}

//! Write particles to HDF5
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::read_particles_hdf5(unsigned phase,
                                          const std::string& filename) {

  // Create a new file using default properties.
  hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  // Throw an error if file can't be found
  if (file_id < 0) throw std::runtime_error("HDF5 particle file is not found");

  // Calculate the size and the offsets of our struct members in memory
  const unsigned nparticles = this->nparticles();
  const hsize_t NRECORDS = nparticles;

  const hsize_t NFIELDS = mpm::hdf5::particle::NFIELDS;

  std::vector<HDF5Particle> dst_buf;
  dst_buf.reserve(nparticles);
  // Read the table
  H5TBread_table(file_id, "table", mpm::hdf5::particle::dst_size,
                 mpm::hdf5::particle::dst_offset,
                 mpm::hdf5::particle::dst_sizes, dst_buf.data());

  unsigned i = 0;
  for (auto pitr = particles_.cbegin(); pitr != particles_.cend(); ++pitr) {
    HDF5Particle particle = dst_buf[i];
    // Get particle's material from list of materials
    auto material = materials_.at(particle.material_id);
    // Initialise particle with HDF5 data
    (*pitr)->initialise_particle(particle, material);
    ++i;
  }
  // close the file
  H5Fclose(file_id);
  return true;
}

//! Nodal coordinates
template <unsigned Tdim>
std::vector<Eigen::Matrix<double, 3, 1>> mpm::Mesh<Tdim>::nodal_coordinates()
    const {

  // Nodal coordinates
  std::vector<Eigen::Matrix<double, 3, 1>> coordinates;
  coordinates.reserve(nodes_.size());

  try {
    if (nodes_.size() == 0)
      throw std::runtime_error("No nodes have been initialised!");

    // Fill nodal coordinates
    for (auto nitr = nodes_.cbegin(); nitr != nodes_.cend(); ++nitr) {
      // initialise coordinates
      Eigen::Matrix<double, 3, 1> node;
      node.setZero();
      auto coords = (*nitr)->coordinates();

      for (unsigned i = 0; i < coords.size(); ++i) node(i) = coords(i);

      coordinates.emplace_back(node);
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    coordinates.clear();
  }
  return coordinates;
}

//! Cell node pairs
template <unsigned Tdim>
std::vector<std::array<mpm::Index, 2>> mpm::Mesh<Tdim>::node_pairs() const {
  // Vector of node_pairs
  std::vector<std::array<mpm::Index, 2>> node_pairs;

  try {
    if (cells_.size() == 0)
      throw std::runtime_error("No cells have been initialised!");

    for (auto citr = cells_.cbegin(); citr != cells_.cend(); ++citr) {
      const auto pairs = (*citr)->side_node_pairs();
      node_pairs.insert(std::end(node_pairs), std::begin(pairs),
                        std::end(pairs));
    }

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    node_pairs.clear();
  }
  return node_pairs;
}

//! Create map of container of particles in sets
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::create_particle_sets(
    const tsl::robin_map<mpm::Index, std::vector<mpm::Index>>& particle_sets,
    bool check_duplicates) {
  bool status = false;
  try {
    // Create container for each particle set
    for (auto sitr = particle_sets.begin(); sitr != particle_sets.end();
         ++sitr) {
      // Create a container for the set
      Container<ParticleBase<Tdim>> particles;
      // Reserve the size of the container
      particles.reserve((sitr->second).size());
      // Add particles to the container
      for (auto pid : sitr->second) {
        bool insertion_status =
            particles.add(map_particles_[pid], check_duplicates);
      }
      // Create the map of the container
      status = this->particle_sets_
                   .insert(std::pair<mpm::Index, Container<ParticleBase<Tdim>>>(
                       sitr->first, particles))
                   .second;
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return status;
}

template <unsigned Tdim>
mpm::Container<mpm::Cell<Tdim>> mpm::Mesh<Tdim>::cells() {
  return this->cells_;
}

//! return particle_ptr
template <unsigned Tdim>
std::map<mpm::Index, mpm::Index>* mpm::Mesh<Tdim>::particles_cell_ids() {
  return &(this->particles_cell_ids_);
}

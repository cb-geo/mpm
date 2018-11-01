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

  // Create a local variable to pass as lambda
  Container<NodeBase<Tdim>> active_nodes;

  tbb::parallel_for_each(
      nodes_.cbegin(), nodes_.cend(),
      [=, &active_nodes](std::shared_ptr<mpm::NodeBase<Tdim>> node) {
        // If node is active add to a list of active nodes
        std::lock_guard<std::mutex> guard(mesh_mutex_);
        if (node->status()) {
          active_nodes.add(node);
        }
      });

  this->active_nodes_ = active_nodes;
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
      if (this->locate_particle_cells(particle))
        status = particles_.add(particle, checks);
      else
        throw std::runtime_error("Particle not found in mesh");
    } else {
      status = particles_.add(particle, checks);
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
  // Remove a particle if found in the container
  bool status = particles_.remove(particle);
  return status;
}

//! Locate particles in a cell
template <unsigned Tdim>
std::vector<std::shared_ptr<mpm::ParticleBase<Tdim>>>
    mpm::Mesh<Tdim>::locate_particles_mesh() {

  std::vector<std::shared_ptr<mpm::ParticleBase<Tdim>>> particles;

  tbb::parallel_for_each(
      particles_.cbegin(), particles_.cend(),
      [=,
       &particles](const std::shared_ptr<mpm::ParticleBase<Tdim>>& particle) {
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
  if (particle->cell_id() != std::numeric_limits<mpm::Index>::max())
    if (particle->compute_reference_location()) return true;

  bool status = false;
  std::for_each(
      cells_.cbegin(), cells_.cend(),
      [=, &status](const std::shared_ptr<mpm::Cell<Tdim>>& cell) {
        // Check if particle is already found, if so don't run for other cells
        // Check if co-ordinates is within the cell, if true
        // add particle to cell
        if (!status && cell->is_point_in_cell(particle->coordinates())) {
          particle->assign_cell(cell);
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
      // Stresses
      if (attribute == "stresses") {
        auto pdata = (*pitr)->stress(phase);
        // Fill stresses to the size of dimensions
        for (unsigned i = 0; i < Tdim; ++i) data(i) = pdata(i);
      }
      // Strains
      else if (attribute == "strains") {
        auto pdata = (*pitr)->strain(phase);
        // Fill stresses to the size of dimensions
        for (unsigned i = 0; i < Tdim; ++i) data(i) = pdata(i);
      }
      // Velocities
      else if (attribute == "velocities") {
        auto pdata = (*pitr)->velocity(phase);
        // Fill stresses to the size of dimensions
        for (unsigned i = 0; i < Tdim; ++i) data(i) = pdata(i);
      }
      // Error
      else
        throw std::runtime_error("Invalid particle vector data attribute: !");
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

//! Assign particle tractions
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::assign_particles_tractions(
    const std::vector<std::tuple<mpm::Index, unsigned, double>>&
        particle_tractions) {
  bool status = false;
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

      // Apply traction
      for (auto pitr = particles_.cbegin(); pitr != particles_.cend(); ++pitr) {
        if ((*pitr)->id() == pid) {
          status = (*pitr)->assign_traction(phase, dir, traction);
          break;
        }
      }
      if (!status)
        throw std::runtime_error("Particle not found / traction is invalid");
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
      (*pitr)->initial_stress(phase, particle_stresses.at(i));
      ++i;
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
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

  mpm::Index i = 0;
  for (auto pitr = particles_.cbegin(); pitr != particles_.cend(); ++pitr) {

    Eigen::Vector3d coordinates;
    coordinates.setZero();
    Eigen::VectorXd coords = (*pitr)->coordinates();
    for (unsigned j = 0; j < Tdim; ++j) coordinates[j] = coords[j];

    Eigen::Vector3d velocity;
    velocity.setZero();
    for (unsigned j = 0; j < Tdim; ++j)
      velocity[j] = (*pitr)->velocity(phase)[j];

    // Particle local size
    Eigen::Vector3d nsize;
    nsize.setZero();
    Eigen::VectorXd size = (*pitr)->natural_size();
    for (unsigned j = 0; j < Tdim; ++j) nsize[j] = size[j];

    Eigen::Matrix<double, 6, 1> stress = (*pitr)->stress(phase);

    Eigen::Matrix<double, 6, 1> strain = (*pitr)->strain(phase);

    particle_data[i].id = (*pitr)->id();
    particle_data[i].mass = (*pitr)->mass(phase);
    particle_data[i].volume = (*pitr)->volume(phase);
    particle_data[i].pressure = (*pitr)->pressure(phase);

    particle_data[i].coord_x = coordinates[0];
    particle_data[i].coord_y = coordinates[1];
    particle_data[i].coord_z = coordinates[2];

    particle_data[i].nsize_x = nsize[0];
    particle_data[i].nsize_y = nsize[1];
    particle_data[i].nsize_z = nsize[2];

    particle_data[i].velocity_x = velocity[0];
    particle_data[i].velocity_y = velocity[1];
    particle_data[i].velocity_z = velocity[2];

    particle_data[i].stress_xx = stress[0];
    particle_data[i].stress_yy = stress[1];
    particle_data[i].stress_zz = stress[2];
    particle_data[i].tau_xy = stress[3];
    particle_data[i].tau_yz = stress[4];
    particle_data[i].tau_xz = stress[5];

    particle_data[i].strain_xx = strain[0];
    particle_data[i].strain_yy = strain[1];
    particle_data[i].strain_zz = strain[2];
    particle_data[i].gamma_xy = strain[3];
    particle_data[i].gamma_yz = strain[4];
    particle_data[i].gamma_xz = strain[5];

    particle_data[i].epsilon_v = (*pitr)->volumetric_strain_centroid(phase);

    particle_data[i].status = (*pitr)->status();

    // Counter
    ++i;
  }
  // Calculate the size and the offsets of our struct members in memory
  const hsize_t NRECORDS = nparticles;

  const hsize_t NFIELDS = 27;

  size_t dst_size = sizeof(HDF5Particle);
  size_t dst_offset[NFIELDS] = {
      HOFFSET(HDF5Particle, id),         HOFFSET(HDF5Particle, mass),
      HOFFSET(HDF5Particle, volume),     HOFFSET(HDF5Particle, pressure),
      HOFFSET(HDF5Particle, coord_x),    HOFFSET(HDF5Particle, coord_y),
      HOFFSET(HDF5Particle, coord_z),    HOFFSET(HDF5Particle, nsize_x),
      HOFFSET(HDF5Particle, nsize_y),    HOFFSET(HDF5Particle, nsize_z),
      HOFFSET(HDF5Particle, velocity_x), HOFFSET(HDF5Particle, velocity_y),
      HOFFSET(HDF5Particle, velocity_z), HOFFSET(HDF5Particle, stress_xx),
      HOFFSET(HDF5Particle, stress_yy),  HOFFSET(HDF5Particle, stress_zz),
      HOFFSET(HDF5Particle, tau_xy),     HOFFSET(HDF5Particle, tau_yz),
      HOFFSET(HDF5Particle, tau_xz),     HOFFSET(HDF5Particle, strain_xx),
      HOFFSET(HDF5Particle, strain_yy),  HOFFSET(HDF5Particle, strain_zz),
      HOFFSET(HDF5Particle, gamma_xy),   HOFFSET(HDF5Particle, gamma_yz),
      HOFFSET(HDF5Particle, gamma_xz),   HOFFSET(HDF5Particle, epsilon_v),
      HOFFSET(HDF5Particle, status),
  };

  size_t dst_sizes[NFIELDS] = {
      sizeof(particle_data[0].id),         sizeof(particle_data[0].mass),
      sizeof(particle_data[0].volume),     sizeof(particle_data[0].pressure),
      sizeof(particle_data[0].coord_x),    sizeof(particle_data[0].coord_y),
      sizeof(particle_data[0].coord_z),    sizeof(particle_data[0].nsize_x),
      sizeof(particle_data[0].nsize_y),    sizeof(particle_data[0].nsize_z),
      sizeof(particle_data[0].velocity_x), sizeof(particle_data[0].velocity_y),
      sizeof(particle_data[0].velocity_z), sizeof(particle_data[0].stress_xx),
      sizeof(particle_data[0].stress_yy),  sizeof(particle_data[0].stress_zz),
      sizeof(particle_data[0].tau_xy),     sizeof(particle_data[0].tau_yz),
      sizeof(particle_data[0].tau_xz),     sizeof(particle_data[0].strain_xx),
      sizeof(particle_data[0].strain_yy),  sizeof(particle_data[0].strain_zz),
      sizeof(particle_data[0].gamma_xy),   sizeof(particle_data[0].gamma_yz),
      sizeof(particle_data[0].gamma_xz),   sizeof(particle_data[0].epsilon_v),
      sizeof(particle_data[0].status),
  };

  // Define particle field information
  const char* field_names[NFIELDS] = {
      "id",         "mass",       "volume",     "pressure",  "coord_x",
      "coord_y",    "coord_z",    "nsize_x",    "nsize_y",   "nsize_z",
      "velocity_x", "velocity_y", "velocity_z", "stress_xx", "stress_yy",
      "stress_zz",  "tau_xy",     "tau_yz",     "tau_xz",    "strain_xx",
      "strain_yy",  "strain_zz",  "gamma_xy",   "gamma_yz",  "gamma_xz",
      "epsilon_v",  "status"};

  hid_t field_type[NFIELDS];
  hid_t string_type;
  hid_t file_id;
  hsize_t chunk_size = 10000;
  int* fill_data = NULL;
  int compress = 0;

  // Initialize the field_type
  field_type[0] = H5T_NATIVE_LLONG;
  field_type[1] = H5T_NATIVE_DOUBLE;
  field_type[2] = H5T_NATIVE_DOUBLE;
  field_type[3] = H5T_NATIVE_DOUBLE;
  field_type[4] = H5T_NATIVE_DOUBLE;
  field_type[5] = H5T_NATIVE_DOUBLE;
  field_type[6] = H5T_NATIVE_DOUBLE;
  field_type[7] = H5T_NATIVE_DOUBLE;
  field_type[8] = H5T_NATIVE_DOUBLE;
  field_type[9] = H5T_NATIVE_DOUBLE;
  field_type[10] = H5T_NATIVE_DOUBLE;
  field_type[11] = H5T_NATIVE_DOUBLE;
  field_type[12] = H5T_NATIVE_DOUBLE;
  field_type[13] = H5T_NATIVE_DOUBLE;
  field_type[14] = H5T_NATIVE_DOUBLE;
  field_type[15] = H5T_NATIVE_DOUBLE;
  field_type[16] = H5T_NATIVE_DOUBLE;
  field_type[17] = H5T_NATIVE_DOUBLE;
  field_type[18] = H5T_NATIVE_DOUBLE;
  field_type[19] = H5T_NATIVE_DOUBLE;
  field_type[20] = H5T_NATIVE_DOUBLE;
  field_type[21] = H5T_NATIVE_DOUBLE;
  field_type[22] = H5T_NATIVE_DOUBLE;
  field_type[23] = H5T_NATIVE_DOUBLE;
  field_type[24] = H5T_NATIVE_DOUBLE;
  field_type[25] = H5T_NATIVE_DOUBLE;
  field_type[26] = H5T_NATIVE_HBOOL;

  // Create a new file using default properties.
  file_id =
      H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  // make a table
  H5TBmake_table("Table Title", file_id, "table", NFIELDS, NRECORDS, dst_size,
                 field_names, dst_offset, field_type, chunk_size, fill_data,
                 compress, particle_data.data());

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

  const hsize_t NFIELDS = 27;

  size_t dst_size = sizeof(HDF5Particle);
  size_t dst_offset[NFIELDS] = {
      HOFFSET(HDF5Particle, id),         HOFFSET(HDF5Particle, mass),
      HOFFSET(HDF5Particle, volume),     HOFFSET(HDF5Particle, pressure),
      HOFFSET(HDF5Particle, coord_x),    HOFFSET(HDF5Particle, coord_y),
      HOFFSET(HDF5Particle, coord_z),    HOFFSET(HDF5Particle, nsize_x),
      HOFFSET(HDF5Particle, nsize_y),    HOFFSET(HDF5Particle, nsize_z),
      HOFFSET(HDF5Particle, velocity_x), HOFFSET(HDF5Particle, velocity_y),
      HOFFSET(HDF5Particle, velocity_z), HOFFSET(HDF5Particle, stress_xx),
      HOFFSET(HDF5Particle, stress_yy),  HOFFSET(HDF5Particle, stress_zz),
      HOFFSET(HDF5Particle, tau_xy),     HOFFSET(HDF5Particle, tau_yz),
      HOFFSET(HDF5Particle, tau_xz),     HOFFSET(HDF5Particle, strain_xx),
      HOFFSET(HDF5Particle, strain_yy),  HOFFSET(HDF5Particle, strain_zz),
      HOFFSET(HDF5Particle, gamma_xy),   HOFFSET(HDF5Particle, gamma_yz),
      HOFFSET(HDF5Particle, gamma_xz),   HOFFSET(HDF5Particle, epsilon_v),
      HOFFSET(HDF5Particle, status),
  };

  // To get size
  HDF5Particle particle;

  size_t dst_sizes[NFIELDS] = {
      sizeof(particle.id),         sizeof(particle.mass),
      sizeof(particle.volume),     sizeof(particle.pressure),
      sizeof(particle.coord_x),    sizeof(particle.coord_y),
      sizeof(particle.coord_z),    sizeof(particle.nsize_x),
      sizeof(particle.nsize_y),    sizeof(particle.nsize_z),
      sizeof(particle.velocity_x), sizeof(particle.velocity_y),
      sizeof(particle.velocity_z), sizeof(particle.stress_xx),
      sizeof(particle.stress_yy),  sizeof(particle.stress_zz),
      sizeof(particle.tau_xy),     sizeof(particle.tau_yz),
      sizeof(particle.tau_xz),     sizeof(particle.strain_xx),
      sizeof(particle.strain_yy),  sizeof(particle.strain_zz),
      sizeof(particle.gamma_xy),   sizeof(particle.gamma_yz),
      sizeof(particle.gamma_xz),   sizeof(particle.epsilon_v),
      sizeof(particle.status),
  };

  std::vector<HDF5Particle> dst_buf;
  dst_buf.reserve(nparticles);
  // Read the table
  H5TBread_table(file_id, "table", dst_size, dst_offset, dst_sizes,
                 dst_buf.data());

  unsigned i = 0;
  for (auto pitr = particles_.cbegin(); pitr != particles_.cend(); ++pitr) {
    particle = dst_buf[i];
    // Initialise particle with HDF5 data
    (*pitr)->initialise_particle(particle);
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

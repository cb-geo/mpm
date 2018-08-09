// Constructor with id
template <unsigned Tdim>
mpm::Mesh<Tdim>::Mesh(unsigned id) : id_{id} {
  // Check if the dimension is between 1 & 3
  static_assert((Tdim >= 1 && Tdim <= 3), "Invalid global dimension");
  particles_.clear();
}

//! Create nodes from coordinates
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::create_nodes(mpm::Index gnid,
                                   const std::string& node_type,
                                   const std::vector<VectorDim>& coordinates) {
  bool status = true;
  try {
    // Check if nodal coordinates is not empty
    if (!coordinates.empty()) {
      for (const auto& node_coordinates : coordinates) {
        // Add node to mesh and check
        bool insert_status = this->add_node(
            // Create a node of particular
            Factory<mpm::NodeBase<Tdim>, mpm::Index,
                    const Eigen::Matrix<double, Tdim, 1>&>::instance()
                ->create(node_type, static_cast<mpm::Index>(gnid),
                         node_coordinates));

        // Increament node id
        if (insert_status) ++gnid;
        // When addition of node fails
        else
          throw std::runtime_error("Addition of node to mesh failed!");
      }
    } else
      // If the coordinates vector is empty
      throw std::runtime_error("List of coordinates is empty");
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
    status = false;
  }
  return status;
}

//! Add a node to the mesh
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::add_node(
    const std::shared_ptr<mpm::NodeBase<Tdim>>& node) {
  bool insertion_status = nodes_.add(node);
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
  //  tbb::parallel_for_each(nodes_.cbegin(), nodes_.cend(), oper);
  for (auto itr = nodes_.cbegin(); itr != nodes_.cend(); ++itr) {
    oper(*itr);
  }

}

//! Iterate over nodes
template <unsigned Tdim>
template <typename Toper, typename Tpred>
void mpm::Mesh<Tdim>::iterate_over_nodes_predicate(Toper oper, Tpred pred) {
  for (auto itr = nodes_.cbegin(); itr != nodes_.cend(); ++itr) {
    if (pred(*itr)) oper(*itr);
  }
}

//! Create cells from node lists
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::create_cells(
    mpm::Index gcid, const std::shared_ptr<mpm::ShapeFn<Tdim>>& shapefn,
    const std::vector<std::vector<mpm::Index>>& cells) {
  bool status = true;
  try {
    // Check if node id list is not empty
    if (!cells.empty()) {
      for (const auto& nodes : cells) {
        // Create cell with shapefn
        auto cell =
            std::make_shared<mpm::Cell<Tdim>>(gcid, nodes.size(), shapefn);

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
          if (cell->is_initialised()) insert_cell = this->add_cell(cell);
        } else
          throw std::runtime_error("Invalid node ids for cell!");

        // Increament global cell id
        if (insert_cell) ++gcid;
        // When addition of cell fails
        else
          throw std::runtime_error("Addition of cell to mesh failed!");
      }
    } else {
      // If the coordinates vector is empty
      throw std::runtime_error("List of nodes of cells is empty");
    }
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
    status = false;
  }
  return status;
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

//! Create particles from coordinates
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::create_particles(
    mpm::Index gpid, const std::string& particle_type,
    const std::vector<VectorDim>& coordinates) {
  bool status = true;
  try {
    // Check if particle coordinates is not empty
    if (!coordinates.empty()) {
      for (const auto& particle_coordinates : coordinates) {
        // Add particle to mesh and check
        bool insert_status = this->add_particle(
            Factory<mpm::ParticleBase<Tdim>, mpm::Index,
                    const Eigen::Matrix<double, Tdim, 1>&>::instance()
                ->create(particle_type, static_cast<mpm::Index>(gpid),
                         particle_coordinates));

        // Increament particle id
        if (insert_status) ++gpid;
        // When addition of particle fails
        else
          throw std::runtime_error("Addition of particle to mesh failed!");
      }
    } else {
      // If the coordinates vector is empty
      throw std::runtime_error("List of coordinates is empty");
    }
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
    status = false;
  }
  return status;
}

//! Add a particle pointer to the mesh
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::add_particle(
    const std::shared_ptr<mpm::ParticleBase<Tdim>>& particle) {
  bool status = false;
  try {
    // Add only if particle can be located in any cell of the mesh
    if (this->locate_particle_cells(particle))
      status = particles_.add(particle);
    else
      throw std::runtime_error("Particle not found in mesh");
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
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
      [=, &particles](std::shared_ptr<mpm::ParticleBase<Tdim>> particle) {
        // If particle is not found in mesh add to a list of particles
        if (!this->locate_particle_cells(particle))
          // Needs a lock guard here
          particles.emplace_back(particle);
      });

  return particles;
}

//! Locate particles in a cell
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::locate_particle_cells(
    std::shared_ptr<mpm::ParticleBase<Tdim>> particle) {
  // Check the current cell if it is not invalid
  if (particle->cell_id() != std::numeric_limits<mpm::Index>::max())
    if (particle->compute_reference_location()) return true;

  bool status = false;
  tbb::parallel_for_each(
      cells_.cbegin(), cells_.cend(),
      [=, &status](std::shared_ptr<mpm::Cell<Tdim>> cell) {
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
  //  tbb::parallel_for_each(particles_.cbegin(), particles_.cend(), oper);
  for (auto itr = particles_.cbegin(); itr != particles_.cend(); ++itr) {
    oper(*itr);
  }

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

//! Return particle stresses
template <unsigned Tdim>
std::vector<Eigen::Matrix<double, 3, 1>> mpm::Mesh<Tdim>::particle_stresses(
    unsigned phase) {
  std::vector<Eigen::Matrix<double, 3, 1>> particle_stresses;
  for (auto pitr = particles_.cbegin(); pitr != particles_.cend(); ++pitr) {
    Eigen::Vector3d stresses;
    stresses.setZero();
    auto pstress = (*pitr)->stress(phase);
    // Fill stresses to the size of dimensions
    for (unsigned i = 0; i < Tdim; ++i) stresses(i) = pstress(i);
    particle_stresses.emplace_back(stresses);
  }
  return particle_stresses;
}

//! Assign velocity constraints
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::assign_velocity_constraints(
    const std::vector<std::tuple<mpm::Index, unsigned, double>>&
        velocity_constraints) {
  bool status = false;
  try {
    if (nodes_.size()) {
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
    } else {
      throw std::runtime_error(
          "No nodes have been assigned in mesh, cannot assign velocity "
          "constraints");
    }
  } catch (std::exception& exception) {
    std::cerr << exception.what() << "\n";
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

    Eigen::Matrix<double, 6, 1> stress = (*pitr)->stress(phase);

    particle_data[i].id = (*pitr)->id();

    particle_data[i].coord_x = coordinates[0];
    particle_data[i].coord_y = coordinates[1];
    particle_data[i].coord_z = coordinates[2];

    particle_data[i].velocity_x = velocity[0];
    particle_data[i].velocity_y = velocity[1];
    particle_data[i].velocity_z = velocity[2];

    particle_data[i].stress_xx = stress[0];
    particle_data[i].stress_yy = stress[1];
    particle_data[i].stress_zz = stress[2];
    particle_data[i].tau_xy = stress[3];
    particle_data[i].tau_yz = stress[4];
    particle_data[i].tau_xz = stress[5];

    particle_data[i].status = (*pitr)->status();

    // Counter
    ++i;
  }
  // Calculate the size and the offsets of our struct members in memory
  const hsize_t NRECORDS = nparticles;

  const hsize_t NFIELDS = 14;

  size_t dst_size = sizeof(HDF5Particle);
  size_t dst_offset[NFIELDS] = {
      HOFFSET(HDF5Particle, id),         HOFFSET(HDF5Particle, coord_x),
      HOFFSET(HDF5Particle, coord_y),    HOFFSET(HDF5Particle, coord_z),
      HOFFSET(HDF5Particle, velocity_x), HOFFSET(HDF5Particle, velocity_y),
      HOFFSET(HDF5Particle, velocity_z), HOFFSET(HDF5Particle, stress_xx),
      HOFFSET(HDF5Particle, stress_yy),  HOFFSET(HDF5Particle, stress_zz),
      HOFFSET(HDF5Particle, tau_xy),     HOFFSET(HDF5Particle, tau_yz),
      HOFFSET(HDF5Particle, tau_xz),     HOFFSET(HDF5Particle, status),
  };

  size_t dst_sizes[NFIELDS] = {
      sizeof(particle_data[0].id),         sizeof(particle_data[0].coord_x),
      sizeof(particle_data[0].coord_y),    sizeof(particle_data[0].coord_z),
      sizeof(particle_data[0].velocity_x), sizeof(particle_data[0].velocity_y),
      sizeof(particle_data[0].velocity_z), sizeof(particle_data[0].stress_xx),
      sizeof(particle_data[0].stress_yy),  sizeof(particle_data[0].stress_zz),
      sizeof(particle_data[0].tau_xy),     sizeof(particle_data[0].tau_yz),
      sizeof(particle_data[0].tau_xz),     sizeof(particle_data[0].status),
  };

  // Define particle field information
  const char* field_names[NFIELDS] = {
      "id",         "coord_x",    "coord_y",   "coord_z",   "velocity_x",
      "velocity_y", "velocity_z", "stress_xx", "stress_yy", "stress_zz",
      "tau_xy",     "tau_yz",     "tau_xz",    "status"};

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
  field_type[13] = H5T_NATIVE_HBOOL;

  // Create a new file using default properties.
  file_id =
      H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  // make a table
  H5TBmake_table("Table Title", file_id, "table", NFIELDS, NRECORDS, dst_size,
                 field_names, dst_offset, field_type, chunk_size, fill_data,
                 compress, particle_data.data());

#ifdef DEBUG
  std::vector<HDF5Particle> dst_buf;
  dst_buf.reserve(nparticles);
  // Read the table
  H5TBread_table(file_id, "table", dst_size, dst_offset, dst_sizes,
                 dst_buf.data());

  // print it by rows
  std::cout << "Printing HDF5 data: \n";
  for (unsigned i = 0; i < this->nparticles(); ++i) {
    std::cout << dst_buf[i].id << '\t' << dst_buf[i].coord_x << '\t'
              << dst_buf[i].coord_y << '\t' << dst_buf[i].coord_z << '\t'
              << dst_buf[i].velocity_x << '\t' << dst_buf[i].velocity_y << '\t'
              << dst_buf[i].velocity_z << '\t' << dst_buf[i].stress_xx << '\t'
              << dst_buf[i].stress_yy << '\t' << dst_buf[i].stress_zz << '\t'
              << dst_buf[i].tau_xy << '\t' << dst_buf[i].tau_yz << '\t'
              << dst_buf[i].tau_xz << '\t' << dst_buf[i].status << '\n';
  }
  std::cout << "End of HDF5 data\n";
#endif

  // close the file
  H5Fclose(file_id);
  return true;
}

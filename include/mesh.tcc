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
  // Check if node_cell_maps is empty
  if (node_cell_maps_.empty()) {
    compute_node_cell_maps();
  }

  for (auto citr = cells_.cbegin(); citr != cells_.cend(); ++citr) {
    // Initiate set of neighbouring cells
    std::set<mpm::Index> neighbouring_cell_sets;

    // Loop over the current cell nodes and add ids of the initiated set
    const auto nodes_id_list = (*citr)->nodes_id();
    for (const auto& id : nodes_id_list)
      neighbouring_cell_sets.insert(node_cell_maps_[id].begin(),
                                    node_cell_maps_[id].end());

    for (const auto& neighbour_id : neighbouring_cell_sets)
      if (neighbour_id != (*citr)->id())
        map_cells_[(*citr)->id()]->add_neighbour(neighbour_id);
  }
}

<<<<<<< HEAD
=======
template <unsigned Tdim>
void mpm::Mesh<Tdim>::compute_node_cell_maps() {
  // Loop over the cells
  for (auto citr = cells_.cbegin(); citr != cells_.cend(); ++citr) {
    // Get cell id and nodes id
    const auto cell_id = (*citr)->id();
    const auto nodes_id_list = (*citr)->nodes_id();
    // Populate node_cell_maps_ with the node_id and multiple cell_id
    for (const auto& id : nodes_id_list) node_cell_maps_[id].insert(cell_id);
  }
}

template <unsigned Tdim>
bool mpm::Mesh<Tdim>::assign_particle_neighbours() {
  bool status = true;

  try {
    // Check if node_cell_maps is empty, this should only be performed at the
    // beginning of the simulation
    if (node_cell_maps_.empty()) {
      compute_node_cell_maps();
    }

    // Loop over cells
    for (auto citr = cells_.cbegin(); citr != cells_.cend(); ++citr) {
      // Initiate set of neighbouring cells
      std::set<mpm::Index> neighbouring_cell_sets;
      neighbouring_cell_sets.insert((*citr)->id());

      // Loop over the current cell nodes and add ids of the initiated set
      const auto nodes_id_list = (*citr)->nodes_id();
      for (const auto& id : nodes_id_list)
        neighbouring_cell_sets.insert(node_cell_maps_[id].begin(),
                                      node_cell_maps_[id].end());

      // Loop over the neighbouring cell particles
      std::vector<mpm::Index> neighbouring_particle_sets;
      for (const auto& neighbour_cell_id : neighbouring_cell_sets)
        neighbouring_particle_sets.push_back(
            map_cells_[neighbour_cell_id]->particles());

      // Loop over the particles in the current cells and assign particle
      // neighbours
      const auto particles_id_list = (*citr)->particles();
      for (const auto& id : particles_id_list) {
        status = map_particles_[id]->assign_particle_neighbours(
            neighbouring_particle_sets);
        if (!status)
          throw std::runtime_error("Cannot assign valid particle neighbours");
      }
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

template <unsigned Tdim>
bool mpm::Mesh<Tdim>::assign_particle_neighbours(const mpm::Index particle_id) {
  bool status = true;

  try {
    std::shared_ptr<mpm::ParticleBase<Tdim>>& particle =
        map_particles_[particle_id];

    // Check if node_cell_maps is empty, this should only be performed at the
    // beginning of the simulation
    if (node_cell_maps_.empty()) {
      compute_node_cell_maps();
    }

    // Get the current particle's cell
    const auto current_cell_id = (*particle)->cell_id();
    const auto& current_cell = map_cells_[current_cell_id];

    // Initiate set of neighbouring cells
    std::set<mpm::Index> neighbouring_cell_sets;
    neighbouring_cell_sets.insert(current_cell_id);

    // Loop over the current cell nodes and add ids of the initiated set
    const auto nodes_id_list = (*current_cell)->nodes_id();
    for (const auto& id : nodes_id_list)
      neighbouring_cell_sets.insert(node_cell_maps_[id].begin(),
                                    node_cell_maps_[id].end());

    // Loop over the neighbouring cell particles
    std::vector<mpm::Index> neighbouring_particle_sets;
    for (const auto& neighbour_cell_id : neighbouring_cell_sets)
      neighbouring_particle_sets.push_back(
          map_cells_[neighbour_cell_id]->particles());

    // Loop over the particles in the current cells and assign particle
    // neighbours
    status =
        (*particle)->assign_particle_neighbours(neighbouring_particle_sets);

    if (!status)
      throw std::runtime_error("Cannot assign valid particle neighbours");

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

template <unsigned Tdim>
bool mpm::Mesh<Tdim>::compute_particle_normal(const mpm::Index particle_id) {
  bool status = true;

  try {
    std::shared_ptr<mpm::ParticleBase<Tdim>>& particle =
        map_particles_[particle_id];

    // Check whether neighbour is empty
    if ((*particle)->nneighbour_particles() == 0)
      throw std::runtime_error(
          "No neighbour particles have been assigned to particle, cannot "
          "compute particle normal");

    Eigen::Matrix<double, Tdim, 1> normal_vector;
    normal_vector.setZero();
    const auto& current_coordinates = (*particle)->coordinates();
    const auto& neighbour_particles = (*particle)->neighbour_particles();
    for (const auto& neighbour_id : neighbour_particles) {
      normal_vector -=
          (map_particles_[neighbour_id]->coordinates() - current_coordinates);
    }
    normal_vector = normal_vector.normalized();

    // Assign normal_vector to particles
    status = (*particle)->assign_particle_normal(normal_vector);

    if (!status)
      throw std::runtime_error("Cannot assign valid particle normal vector");

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Clone particle quantities
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::interpolate_particle_quantities(
    const mpm::Index particle_id) {
  status = true;
  try {
    std::shared_ptr<mpm::ParticleBase<Tdim>>& particle =
        map_particles_[particle_id];

    // Check whether neighbour is empty
    if ((*particle)->nneighbour_particles() == 0)
      throw std::runtime_error(
          "No neighbour particles have been assigned to particle, cannot "
          "interpolate particle quantities");

    const auto& current_coordinates = (*particle)->coordinates();
    const auto& mean_cell_length =
        map_cells_[(*particle)->cell_id()]->mean_length();
    const auto& neighbour_particles = (*particle)->neighbour_particles();

    // Current particle data
    mpm::HDF5Particle particle_data = (*particle)->hdf5();

    // Reinitialise variables which necessary for interpolation
    particle_data.pressure = 0.0;

    particle_data.displacement_x = 0.0;
    particle_data.displacement_y = 0.0;
    particle_data.displacement_z = 0.0;

    particle_data.velocity_x = 0.0;
    particle_data.velocity_y = 0.0;
    particle_data.velocity_z = 0.0;

    particle_data.stress_xx = 0.0;
    particle_data.stress_yy = 0.0;
    particle_data.stress_zz = 0.0;
    particle_data.tau_xy = 0.0;
    particle_data.tau_yz = 0.0;
    particle_data.tau_xz = 0.0;

    particle_data.strain_xx = 0.0;
    particle_data.strain_yy = 0.0;
    particle_data.strain_zz = 0.0;
    particle_data.gamma_xy = 0.0;
    particle_data.gamma_yz = 0.0;
    particle_data.gamma_xz = 0.0;

    particle_data.epsilon_v = 0.0;

    for (const auto& neighbour_id : neighbour_particles) {
      // Get Basis function
      const auto relative_coordinates =
          (map_particles_[neighbour_id]->coordinates() - current_coordinates);
      const double norm_distance = relative_coordinates.squaredNorm();
      const double RBF3 =
          cubic_radial_basis_function(mean_cell_length, norm_distance);

      // Get Neighbour volume
      const double neighbour_volume = map_particles_[neighbour_id]->volume();

      // Get neighbours quantities
      const auto& neighbour_particle_data =
          map_particles_[neighbour_id]->hdf5();

      // Interpolate necessary quantities
      particle_data.pressure =
          neighbour_particle_data.pressure * neighbour_volume * RBF3;

      particle_data.displacement_x +=
          neighbour_particle_data.displacement_x * neighbour_volume * RBF3;
      particle_data.displacement_y +=
          neighbour_particle_data.displacement_y * neighbour_volume * RBF3;
      particle_data.displacement_z +=
          neighbour_particle_data.displacement_z * neighbour_volume * RBF3;

      particle_data.velocity_x +=
          neighbour_particle_data.velocity_x * neighbour_volume * RBF3;
      particle_data.velocity_y +=
          neighbour_particle_data.velocity_y * neighbour_volume * RBF3;
      particle_data.velocity_z +=
          neighbour_particle_data.velocity_z * neighbour_volume * RBF3;

      particle_data.stress_xx +=
          neighbour_particle_data.stress_xx * neighbour_volume * RBF3;
      particle_data.stress_yy +=
          neighbour_particle_data.stress_yy * neighbour_volume * RBF3;
      particle_data.stress_zz +=
          neighbour_particle_data.stress_zz * neighbour_volume * RBF3;
      particle_data.tau_xy +=
          neighbour_particle_data.tau_xy * neighbour_volume * RBF3;
      particle_data.tau_yz +=
          neighbour_particle_data.tau_yz * neighbour_volume * RBF3;
      particle_data.tau_xz +=
          neighbour_particle_data.tau_xz * neighbour_volume * RBF3;

      particle_data.strain_xx +=
          neighbour_particle_data.strain_xx * neighbour_volume * RBF3;
      particle_data.strain_yy +=
          neighbour_particle_data.strain_yy * neighbour_volume * RBF3;
      particle_data.strain_zz +=
          neighbour_particle_data.strain_zz * neighbour_volume * RBF3;
      particle_data.gamma_xy +=
          neighbour_particle_data.gamma_xy * neighbour_volume * RBF3;
      particle_data.gamma_yz +=
          neighbour_particle_data.gamma_yz * neighbour_volume * RBF3;
      particle_data.gamma_xz +=
          neighbour_particle_data.gamma_xz * neighbour_volume * RBF3;

      particle_data.epsilon_v +=
          neighbour_particle_data.epsilon_v * neighbour_volume * RBF3;
    }

    // Assign particle quantities to particle
    status = (*particle)->initialise_particle(particle_data);

    if (!status)
      throw std::runtime_error("Cannot assign valid particle state variables");

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Clone material and state variables
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::interpolate_particle_material(
    const mpm::Index particle_id) {
  status = true;
  try {
    std::shared_ptr<mpm::ParticleBase<Tdim>>& particle =
        map_particles_[particle_id];

    // Check whether neighbour is empty
    if ((*particle)->nneighbour_particles() == 0)
      throw std::runtime_error(
          "No neighbour particles have been assigned to particle, cannot "
          "interpolate particle material");

    const auto& current_coordinates = (*particle)->coordinates();
    const auto& mean_cell_length =
        map_cells_[(*particle)->cell_id()]->mean_length();
    const auto& neighbour_particles = (*particle)->neighbour_particles();

    // Check whether all material id are the same otherwise throw error
    // if all the same, then clone from begin()
    const unsigned material_id =
        map_particles_[(*neighbour_particles.cbegin())]->material()->id();
    for (const auto& neighbour_id : neighbour_particles) {
      const auto neighbour_mat_id =
          map_particles_[neighbour_id]->material()->id();
      if (material_id != neighbour_mat_id) {
        std::string error_mssg =
            "Neighbour particles have different material id. ";
        error_mssg += "Currently interpolate_particle_state_variable ";
        error_mssg += "only support single material id.";
        throw std::runtime_error(error_mssg);
      }
    }

    // Assign particle material (clone it from the first neighbour)
    status = (*particle)->assign_material(
        map_particles_[(*neighbour_particles.cbegin())]->material());
    auto particle_state_vars =
        (*particle)->material()->initialise_state_variables();

    if (!status) throw std::runtime_error("Cannot assign valid material");

    for (const auto& neighbour_id : neighbour_particles) {
      // Get Basis function
      const auto relative_coordinates =
          (map_particles_[neighbour_id]->coordinates() - current_coordinates);
      const double norm_distance = relative_coordinates.squaredNorm();
      const double RBF3 =
          cubic_radial_basis_function(mean_cell_length, norm_distance);

      // Get Neighbour volume
      const double neighbour_volume = map_particles_[neighbour_id]->volume();

      // Get neighbours state variables
      const auto& neighbour_state_vars =
          map_particles_[neighbour_id]->state_variables();

      // Loop over neighbours state variables
      for (const auto& state_vars : neighbour_state_vars) {
        // Interpolate value for given key
        const auto var_name = state_vars.key();
        const auto var_value = state_vars.value();
        particle_state_vars.at(var_name) += var_value * neighbour_volume * RBF3;
      }
    }

    // Assign state variable to particle
    status = (*particle)->assign_state_variables(particle_state_vars);

    if (!status)
      throw std::runtime_error("Cannot assign valid particle state variables");

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// FIXME: Move to a separate file
//! Temporary Radial Basis Function
template <unsigned Tdim>
double mpm::Mesh<Tdim>::cubic_radial_basis_function(
    const double& smoothing_length, const double& norm_distance,
    const double& scaling_factor) {

  // Assign multiplier depends on dimension
  double multiplier;
  if (Tdim == 2)
    multiplier = 15.0 / (7.0 * M_PI * std::pow(smoothing_length, 2));
  else if (Tdim == 3)
    multiplier = 3.0 / (2.0 * M_PI * std::pow(smoothing_length, 3));
  else
    throw std::runtime_error("Tdim is invalid");

  // Compute basis function
  double basis_function = multiplier;
  const double radius = norm_distance / smoothing_length;
  if (radius >= 0.0 && radius < 1.0)
    basis_function *=
        2.0 / 3.0 - std::pow(radius, 2) + 0.5 * std::pow(radius, 3);
  else if (radius >= 1.0 && radius < 2.0)
    basis_function *= 1.0 / 6.0 * std::pow((2.0 - radius), 3);
  else
    basis_function = 0.0;

  return basis_function;
}

// FIXME: Remove this temporary function
template <unsigned Tdim>
unsigned mpm::Mesh<Tdim>::u_power(unsigned x, unsigned p) {
  if (p == 0) return 1;
  if (p == 1) return x;
  return x * this->u_power(x, p - 1);
}

// FIXME: Currently only working for quad4N and quad8N elements
template <unsigned Tdim>
unsigned mpm::Mesh<Tdim>::compute_cell_zone(
    const Eigen::Matrix<double, Tdim, 1>& xi) {
  unsigned zone_id = 0;
  if (Tdim == 2) {
    Eigen::Matrix<double, Tdim, 1> max, min, median;
    max << 1.0, 1.0;
    min << -1.0, -1.0;
    median = 0.5 * (max + min);

    std::vector<bool> exponent(4, false);
    if (xi[0] > median[0]) exponent[0] = true;
    if (xi[1] > median[1]) exponent[1] = true;

    max = median;
    if (exponent[0]) {
      max[0] += 1.0;
      min[0] += 1.0;
    }
    if (exponent[1]) {
      max[1] += 1.0;
      min[1] += 1.0;
    }

    median = 0.5 * (max + min);
    if (xi[0] > median[0]) exponent[2] = true;
    if (xi[1] > median[1]) exponent[3] = true;

    // Loop over exponents and check
    for (unsigned i = 0; i < 4; i++)
      if (exponent[i]) zone_id += u_power(2, (3 - i));

  } else {
    throw std::runtime_error("Tdim is invalid");
  }
  return zone_id;
}

// FIXME: Currently only working for quad4N and quad8N elements (nzones==16)
template <unsigned Tdim>
Eigen::Matrix<double, Tdim, 1> mpm::Mesh<Tdim>::cell_zone_center(
    const unsigned zone_id) {
  Eigen::Matrix<double, Tdim, 1> center_xi;
  switch (zone_id) {
    // Quadrant 0
    case (0):
      center_xi << -0.75, -0.75;
      break;
    case (1):
      center_xi << -0.25, -0.75;
      break;
    case (2):
      center_xi << -0.75, -0.25;
      break;
    case (3):
      center_xi << -0.25, -0.25;
      break;
    // Quadrant 1
    case (4):
      center_xi << 0.25, -0.75;
      break;
    case (5):
      center_xi << 0.75, -0.75;
      break;
    case (6):
      center_xi << 0.25, -0.25;
      break;
    case (7):
      center_xi << 0.75, -0.25;
      break;
    // Quadrant 2
    case (8):
      center_xi << -0.75, 0.25;
      break;
    case (9):
      center_xi << -0.25, 0.25;
      break;
    case (10):
      center_xi << -0.75, 0.75;
      break;
    case (11):
      center_xi << -0.25, 0.75;
      break;
    // Quadrant 3
    case (12):
      center_xi << 0.25, 0.25;
      break;
    case (13):
      center_xi << 0.75, 0.25;
      break;
    case (14):
      center_xi << 0.25, 0.75;
      break;
    case (15):
      center_xi << 0.75, 0.75;
      break;
  }
  return center_xi;
}

template <unsigned Tdim>
bool mpm::Mesh<Tdim>::populate_cells_with_material_points(
    const mpm::Index cell_id) {
  status = true;
  try {
    std::shared_ptr<mpm::Cell<Tdim>>& cell = map_cells_[cell_id];

    // Check whether we need refinement or not
    bool refine = (*cell)->particle_refinement();
    if (refine) {
      // Partition elements (natural coordinates) to smaller zones
      const auto& element = (*cell)->element_ptr();
      const auto nzone = (*element)->nzones();
      std::vector<bool> zones(nzone, true);

      // Loop over particles inside cell
      for (const auto particle_id : (*cell)->particles()) {
        // Locate particle location and color zones
        const auto& particle = map_particles_[particle_id];
        Eigen::Matrix<double, Tdim, 1> xi;
        xi.setZero();
        if ((*particle)->compute_reference_location())
          xi = (*particle)->reference_location();

        // Check where particle lies in
        const unsigned zone_id = this->compute_cell_zone(xi);
        zones[zone_id] = false;
      }

      // Prepare particle neighbour set for Jacobi-style interpolation
      auto neighbouring_particle_sets =
          map_particles_[(*(*cell)->particles().cbegin())]
              ->neighbour_particles();
      neighbouring_particle_sets.push_back((*(*cell)->particles().cbegin()));

      // Loop over zones
      for (unsigned i = 0; i < nzone; i++) {
        // Check zone colors
        if (zones[i]) {
          // Create particles at the center of the zone
          const auto center_xi = this->cell_zone_center(i);
          // Add particle to mesh and check
          mpm::Index particle_id = this->nparticles();
          bool insert_status = this->add_particle(
              Factory<mpm::ParticleBase<Tdim>, mpm::Index,
                      const Eigen::Matrix<double, Tdim, 1>&>::instance()
                  ->create("P2D", particle_id, center_xi),
              true);

          // Assign particle neighbour
          auto& particle = map_particles_[particle_id];
          status = (*particle)->assign_particle_neighbours(
              neighbouring_particle_sets);

          // Interpolate particle quantities
          this->interpolate_particle_quantities(particle_id);

          // Interpolate particle material
          this->interpolate_particle_material(particle_id);
        }
      }
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

>>>>>>> 5142e31... rename cell maps
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

      // if (map_particles_.find(pid) != map_particles_.end())
      status = map_particles_[pid]->assign_particle_velocity_constraint(
          dir, velocity);

      if (!status) throw std::runtime_error("Velocity is invalid for particle");
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

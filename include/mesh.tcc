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
#pragma omp parallel for schedule(runtime)
  for (auto nitr = nodes_.cbegin(); nitr != nodes_.cend(); ++nitr) oper(*nitr);
}

//! Iterate over nodes
template <unsigned Tdim>
template <typename Toper, typename Tpred>
void mpm::Mesh<Tdim>::iterate_over_nodes_predicate(Toper oper, Tpred pred) {
#pragma omp parallel for schedule(runtime)
  for (auto nitr = nodes_.cbegin(); nitr != nodes_.cend(); ++nitr) {
    if (pred(*nitr)) oper(*nitr);
  }
}

//! Create a list of active nodes in mesh
template <unsigned Tdim>
void mpm::Mesh<Tdim>::find_active_nodes() {
  // Clear existing list of active nodes
  this->active_nodes_.clear();

  for (auto nitr = nodes_.cbegin(); nitr != nodes_.cend(); ++nitr)
    if ((*nitr)->status()) this->active_nodes_.add(*nitr);
}

//! Create a list of active nodes in mesh and assign active node id (rank-wise)
template <unsigned Tdim>
unsigned mpm::Mesh<Tdim>::assign_active_nodes_id() {
  // Clear existing list of active nodes
  this->active_nodes_.clear();
  Index active_id = 0;

  for (auto nitr = nodes_.cbegin(); nitr != nodes_.cend(); ++nitr) {
    if ((*nitr)->status()) {
      this->active_nodes_.add(*nitr);
      (*nitr)->assign_active_id(active_id);
      active_id++;
    } else {
      (*nitr)->assign_active_id(std::numeric_limits<Index>::max());
    }
  }

  return active_id;
}

//! Return global node indices
template <unsigned Tdim>
std::vector<Eigen::VectorXi> mpm::Mesh<Tdim>::global_node_indices() const {
  // Vector of node_pairs
  std::vector<Eigen::VectorXi> node_indices;
  try {
    // Iterate over cells
    for (auto citr = cells_.cbegin(); citr != cells_.cend(); ++citr) {
      if ((*citr)->status()) {
        node_indices.emplace_back((*citr)->local_node_indices());
      }
    }

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return node_indices;
}

//! Iterate over active nodes
template <unsigned Tdim>
template <typename Toper>
void mpm::Mesh<Tdim>::iterate_over_active_nodes(Toper oper) {
#pragma omp parallel for schedule(runtime)
  for (auto nitr = active_nodes_.cbegin(); nitr != active_nodes_.cend(); ++nitr)
    oper(*nitr);
}

#ifdef USE_MPI
#ifdef USE_HALO_EXCHANGE
//! Nodal halo exchange
template <unsigned Tdim>
template <typename Ttype, unsigned Tnparam, typename Tgetfunctor,
          typename Tsetfunctor>
void mpm::Mesh<Tdim>::nodal_halo_exchange(Tgetfunctor getter,
                                          Tsetfunctor setter) {
  // Create vector of nodal vectors
  unsigned nnodes = this->domain_shared_nodes_.size();

  // Get number of MPI ranks
  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  if (mpi_size > 1) {
    // Vector of send requests
    std::vector<MPI_Request> send_requests;
    send_requests.reserve(ncomms_);

    unsigned j = 0;
    // Non-blocking send
    for (unsigned i = 0; i < nnodes; ++i) {
      Ttype property = getter(domain_shared_nodes_[i]);
      std::set<unsigned> node_mpi_ranks = domain_shared_nodes_[i]->mpi_ranks();
      for (auto& node_rank : node_mpi_ranks) {
        if (node_rank != mpi_rank) {
          MPI_Isend(&property, Tnparam, MPI_DOUBLE, node_rank, 0,
                    MPI_COMM_WORLD, &send_requests[j]);
          ++j;
        }
      }
    }

    // send complete
    for (unsigned i = 0; i < ncomms_; ++i)
      MPI_Wait(&send_requests[i], MPI_STATUS_IGNORE);

    for (unsigned i = 0; i < nnodes; ++i) {
      // Get value at current node
      Ttype property = getter(domain_shared_nodes_[i]);

      std::set<unsigned> node_mpi_ranks = domain_shared_nodes_[i]->mpi_ranks();
      // Receive from all shared ranks
      for (auto& node_rank : node_mpi_ranks) {
        if (node_rank != mpi_rank) {
          Ttype value;
          MPI_Recv(&value, Tnparam, MPI_DOUBLE, node_rank, 0, MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);
          property += value;
        }
      }
      setter(domain_shared_nodes_[i], property);
    }
  }
}

#else
//! All reduce over nodal scalar property
template <unsigned Tdim>
template <typename Ttype, unsigned Tnparam, typename Tgetfunctor,
          typename Tsetfunctor>
void mpm::Mesh<Tdim>::nodal_halo_exchange(Tgetfunctor getter,
                                          Tsetfunctor setter) {
  // Create vector of nodal scalars
  std::vector<Ttype> prop_get(nhalo_nodes_, mpm::zero<Ttype>());
  std::vector<Ttype> prop_set(nhalo_nodes_, mpm::zero<Ttype>());

#pragma omp parallel for schedule(runtime) shared(prop_get)
  for (auto nitr = domain_shared_nodes_.cbegin();
       nitr != domain_shared_nodes_.cend(); ++nitr)
    prop_get.at((*nitr)->ghost_id()) = getter((*nitr));

  MPI_Allreduce(prop_get.data(), prop_set.data(), nhalo_nodes_ * Tnparam,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#pragma omp parallel for schedule(runtime)
  for (auto nitr = domain_shared_nodes_.cbegin();
       nitr != domain_shared_nodes_.cend(); ++nitr)
    setter((*nitr), prop_set.at((*nitr)->ghost_id()));
}
#endif
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
#pragma omp parallel for schedule(runtime)
  for (auto citr = cells_.cbegin(); citr != cells_.cend(); ++citr) oper(*citr);
}

//! Create cells from node lists
template <unsigned Tdim>
void mpm::Mesh<Tdim>::find_cell_neighbours() {
  // Initialize and compute node cell map
  tsl::robin_map<mpm::Index, std::set<mpm::Index>> node_cell_map;
  for (auto citr = cells_.cbegin(); citr != cells_.cend(); ++citr) {
    // Populate node_cell_map with the node_id and multiple cell_id
    auto cell_id = (*citr)->id();
    for (auto id : (*citr)->nodes_id()) node_cell_map[id].insert(cell_id);
  }

#pragma omp parallel for schedule(runtime)
  for (auto citr = cells_.cbegin(); citr != cells_.cend(); ++citr) {
    // Iterate over each node in current cell
    for (auto id : (*citr)->nodes_id()) {
      auto cell_id = (*citr)->id();
      // Get the cells associated with each node
      for (auto neighbour_id : node_cell_map[id])
        if (neighbour_id != cell_id) (*citr)->add_neighbour(neighbour_id);
    }
  }
}

//! Find global number of particles across MPI ranks / cell
template <unsigned Tdim>
void mpm::Mesh<Tdim>::find_nglobal_particles_cells() {
  int mpi_rank = 0;
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  for (auto citr = cells_.cbegin(); citr != cells_.cend(); ++citr) {
    int nparticles;
    // Determine the rank of the broadcast emitter process
    if ((*citr)->rank() == mpi_rank) nparticles = (*citr)->nparticles();
    MPI_Bcast(&nparticles, 1, MPI_INT, (*citr)->rank(), MPI_COMM_WORLD);
    // Receive broadcast and update on all ranks
    (*citr)->nglobal_particles(nparticles);
  }
#endif
}

//! Find particle neighbours for all particle
template <unsigned Tdim>
void mpm::Mesh<Tdim>::find_particle_neighbours() {
  for (auto citr = cells_.cbegin(); citr != cells_.cend(); ++citr)
    this->find_particle_neighbours(*citr);
}

//! Find particle neighbours for specific cell particle
template <unsigned Tdim>
void mpm::Mesh<Tdim>::find_particle_neighbours(
    const std::shared_ptr<mpm::Cell<Tdim>>& cell) {
  int mpi_rank = 0;
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif

  // Particles in current cell
  std::vector<mpm::Index> neighbour_particles = cell->particles();
  // Loop over all neighboring cells, and append particle ids from each cell
  for (const auto& neighbour_cell_id : cell->neighbours()) {
    // Get the MPI rank of the neighbour cell
    int neighbour_cell_rank = map_cells_[neighbour_cell_id]->rank();
    if (neighbour_cell_rank != cell->rank()) {
#ifdef USE_MPI
      // Send particle ids
      if (neighbour_cell_rank == mpi_rank) {
        // Get particle ids from each cell
        auto send_particle_ids = map_cells_[neighbour_cell_id]->particles();
        // Get size of the particle ids
        int pid_size = send_particle_ids.size();
        // Send the size of the particles in cell
        MPI_Send(&pid_size, 1, MPI_INT, cell->rank(), neighbour_cell_id,
                 MPI_COMM_WORLD);

        // Send particle ids if it is not empty
        if (pid_size > 0)
          MPI_Send(send_particle_ids.data(), pid_size, MPI_UNSIGNED_LONG_LONG,
                   cell->rank(), neighbour_cell_id, MPI_COMM_WORLD);
      }
      // Receive particle ids in the current MPI rank
      if (cell->rank() == mpi_rank) {
        // Particle ids at local cell MPI rank
        std::vector<mpm::Index> received_particle_ids;
        int nparticles = 0;
        MPI_Recv(&nparticles, 1, MPI_INT, neighbour_cell_rank,
                 neighbour_cell_id, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (nparticles > 0) {
          received_particle_ids.resize(nparticles);
          MPI_Recv(received_particle_ids.data(), nparticles,
                   MPI_UNSIGNED_LONG_LONG, neighbour_cell_rank,
                   neighbour_cell_id, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        neighbour_particles.insert(neighbour_particles.end(),
                                   received_particle_ids.begin(),
                                   received_particle_ids.end());
      }
#endif
    } else {
      const auto& particle_ids = map_cells_[neighbour_cell_id]->particles();
      neighbour_particles.insert(neighbour_particles.end(),
                                 particle_ids.begin(), particle_ids.end());
    }
  }

  // Assign neighbouring particle ids to particles in the current cell
  for (auto particle_id : cell->particles())
    map_particles_[particle_id]->assign_neighbours(neighbour_particles);
}

//! Find ghost cell neighbours
template <unsigned Tdim>
void mpm::Mesh<Tdim>::find_ghost_boundary_cells() {
#ifdef USE_MPI
  // Get number of MPI ranks
  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  bool check_duplicates = true;
  if (mpi_size > 1) {
    ghost_cells_.clear();
    local_ghost_cells_.clear();
    ghost_cells_neighbour_ranks_.clear();
    // Iterate through cells
    for (auto citr = cells_.cbegin(); citr != cells_.cend(); ++citr) {
      std::set<unsigned> neighbour_ranks;
      // If cell rank is the current MPI rank
      if ((*citr)->rank() == mpi_rank) {
        // Iterate through the neighbours of a cell
        auto neighbours = (*citr)->neighbours();
        for (auto neighbour : neighbours) {
          // If the neighbour is in a different MPI rank
          if (map_cells_[neighbour]->rank() != mpi_rank) {
            ghost_cells_.add(map_cells_[neighbour], check_duplicates);
            // Add mpi rank to set
            neighbour_ranks.insert(map_cells_[neighbour]->rank());
          }
        }
      }
      // Set the number of different MPI rank neighbours to a ghost cell
      if (neighbour_ranks.size() > 0) {
        // Also add the current cell, as this would be a receiver
        local_ghost_cells_.add(*citr, check_duplicates);

        // Update the neighbouring ranks of the local ghost cell
        std::vector<unsigned> mpi_neighbours;
        for (auto rank : neighbour_ranks) mpi_neighbours.emplace_back(rank);
        ghost_cells_neighbour_ranks_[(*citr)->id()] = mpi_neighbours;
      }
    }
  }
#endif
}

//! Find ncells in rank
template <unsigned Tdim>
mpm::Index mpm::Mesh<Tdim>::ncells_rank(bool active_cells) {
  unsigned ncells_rank = 0;

  int mpi_rank = 0;
  int mpi_size = 1;
#ifdef USE_MPI
  // Get number of MPI ranks
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif

  if (active_cells) {
    for (auto citr = cells_.cbegin(); citr != cells_.cend(); ++citr)
      if ((*citr)->rank() == mpi_rank && (*citr)->nparticles() > 0)
        ncells_rank += 1;
  } else {
    for (auto citr = cells_.cbegin(); citr != cells_.cend(); ++citr)
      if ((*citr)->rank() == mpi_rank) ncells_rank += 1;
  }
  return ncells_rank;
}

//! Find nnodes in rank
template <unsigned Tdim>
mpm::Index mpm::Mesh<Tdim>::nnodes_rank() {
  unsigned nnodes_rank = 0;

  int mpi_rank = 0;
#ifdef USE_MPI
  // Get number of MPI ranks
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif
  for (auto nitr = nodes_.cbegin(); nitr != nodes_.cend(); ++nitr) {
    // Get MPI ranks for the node
    auto mpi_ranks = (*nitr)->mpi_ranks();
    // Check if the local rank is in the list of ranks for the node
    const bool local_node = mpi_ranks.find(mpi_rank) != mpi_ranks.end();
    if (local_node) nnodes_rank += 1;
  }
  return nnodes_rank;
}

//! Create cells from node lists
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::generate_material_points(
    unsigned nquadratures, const std::string& particle_type,
    const std::vector<unsigned>& material_ids, int cset_id, unsigned pset_id) {
  bool status = true;
  try {
    if (cells_.size() > 0) {
      // Particle ids
      std::vector<mpm::Index> pids;
      unsigned before_generation = this->nparticles();
      bool checks = false;
      // Get material
      std::vector<std::shared_ptr<mpm::Material<Tdim>>> materials;
      for (auto m_id : material_ids)
        materials.emplace_back(materials_.at(m_id));

      // If set id is -1, use all cells
      auto cset = (cset_id == -1) ? this->cells_ : cell_sets_.at(cset_id);
      // Iterate over each cell to generate points
      for (auto citr = cset.cbegin(); citr != cset.cend(); ++citr) {
        (*citr)->assign_quadrature(nquadratures);
        // Genereate particles at the Gauss points
        const auto cpoints = (*citr)->generate_points();
        // Iterate over each coordinate to generate material points
        for (const auto& coordinates : cpoints) {
          // Particle id
          mpm::Index pid = particles_.size();
          // Create particle
          auto particle =
              Factory<mpm::ParticleBase<Tdim>, mpm::Index,
                      const Eigen::Matrix<double, Tdim, 1>&>::instance()
                  ->create(particle_type, static_cast<mpm::Index>(pid),
                           coordinates);

          // Add particle to mesh
          status = this->add_particle(particle, checks);
          if (status) {
            map_particles_[pid]->assign_cell(*citr);
            for (unsigned phase = 0; phase < materials.size(); phase++)
              map_particles_[pid]->assign_material(materials[phase], phase);
            pids.emplace_back(pid);
          } else
            throw std::runtime_error("Generate particles in mesh failed");
        }
      }
      if (before_generation == this->nparticles())
        throw std::runtime_error("No particles were generated!");

      // Add particles to set
      status = this->particle_sets_
                   .insert(std::pair<mpm::Index, std::vector<mpm::Index>>(
                       pset_id, pids))
                   .second;
      if (!status) throw std::runtime_error("Particle set creation failed");

      console_->info(
          "Generate points:\n# of cells: {}\nExpected # of points: {}\n"
          "# of points generated: {}",
          cells_.size(), cells_.size() * std::pow(nquadratures, Tdim),
          particles_.size());
    } else
      throw std::runtime_error("No cells are found in the mesh!");
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Create particles from coordinates
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::create_particles(
    const std::string& particle_type, const std::vector<VectorDim>& coordinates,
    const std::vector<unsigned>& material_ids, unsigned pset_id,
    bool check_duplicates) {
  bool status = true;
  try {
    // Particle ids
    std::vector<mpm::Index> pids;
    // Get material
    std::vector<std::shared_ptr<mpm::Material<Tdim>>> materials;
    for (auto m_id : material_ids) materials.emplace_back(materials_.at(m_id));
    // Check if particle coordinates is empty
    if (coordinates.empty())
      throw std::runtime_error("List of coordinates is empty");
    // Iterate over particle coordinates
    for (const auto& particle_coordinates : coordinates) {
      // Particle id
      mpm::Index pid = particles_.size();
      // Create particle
      auto particle = Factory<mpm::ParticleBase<Tdim>, mpm::Index,
                              const Eigen::Matrix<double, Tdim, 1>&>::instance()
                          ->create(particle_type, static_cast<mpm::Index>(pid),
                                   particle_coordinates);

      // Add particle to mesh and check
      bool insert_status = this->add_particle(particle, check_duplicates);

      // If insertion is successful
      if (insert_status) {
        for (unsigned phase = 0; phase < materials.size(); phase++)
          map_particles_[pid]->assign_material(materials[phase], phase);
        pids.emplace_back(pid);
      } else
        throw std::runtime_error("Addition of particle to mesh failed!");
    }
    // Add particles to set
    status = this->particle_sets_
                 .insert(std::pair<mpm::Index, std::vector<mpm::Index>>(pset_id,
                                                                        pids))
                 .second;
    if (!status) throw std::runtime_error("Particle set creation failed");
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

//! Remove a particle by id
template <unsigned Tdim>
void mpm::Mesh<Tdim>::remove_particles(const std::vector<mpm::Index>& pids) {
  if (!pids.empty()) {
    // Get MPI rank
    int mpi_size = 1;
#ifdef USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif
    for (auto& id : pids) {
      map_particles_[id]->remove_cell();
      map_particles_.remove(id);
    }

    // Get number of particles to reserve size
    unsigned nparticles = this->nparticles();
    // Clear particles and start a new element of particles
    particles_.clear();
    particles_.reserve(static_cast<int>(nparticles / mpi_size));
    // Iterate over the map of particles and add them to container
    for (auto& particle : map_particles_)
      particles_.add(particle.second, false);
  }
}

//! Remove all particles in a cell given cell id
template <unsigned Tdim>
void mpm::Mesh<Tdim>::remove_all_nonrank_particles() {
  // Get MPI rank
  int mpi_rank = 0;
  int mpi_size = 1;
#ifdef USE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif

  // Remove associated cell for the particle
  for (auto citr = this->cells_.cbegin(); citr != this->cells_.cend(); ++citr) {
    // If cell is non empty
    if ((*citr)->particles().size() != 0 && (*citr)->rank() != mpi_rank) {
      auto pids = (*citr)->particles();
      // Remove particles from map
      for (auto& id : pids) {
        map_particles_[id]->remove_cell();
        map_particles_.remove(id);
      }
      (*citr)->clear_particle_ids();
    }
  }

  // Get number of particles to reserve size
  unsigned nparticles = this->nparticles();
  // Clear particles and start a new element of particles
  particles_.clear();
  particles_.reserve(static_cast<int>(nparticles / mpi_size));
  // Iterate over the map of particles and add them to container
  for (auto& particle : map_particles_) particles_.add(particle.second, false);
}

//! Transfer all particles in cells that are not in local rank
template <unsigned Tdim>
void mpm::Mesh<Tdim>::transfer_halo_particles() {
#ifdef USE_MPI
  // Get number of MPI ranks
  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  if (mpi_size > 1) {
    std::vector<MPI_Request> send_requests;
    send_requests.reserve(ghost_cells_.size());
    unsigned i = 0;

    std::vector<mpm::Index> remove_pids;
    // Iterate through the ghost cells and send particles
    for (auto citr = this->ghost_cells_.cbegin();
         citr != this->ghost_cells_.cend(); ++citr, ++i) {
      // Create a vector of h5_particles
      std::vector<mpm::HDF5Particle> h5_particles;
      auto particle_ids = (*citr)->particles();
      // Create a vector of HDF5 data of particles to send
      // delete particle
      for (auto& id : particle_ids) {
        // Append to vector of particles
        h5_particles.emplace_back(map_particles_[id]->hdf5());
        // Particles to be removed from the current rank
        remove_pids.emplace_back(id);
      }
      (*citr)->clear_particle_ids();

      // Send number of particles to receiver rank
      unsigned nparticles = h5_particles.size();
      MPI_Isend(&nparticles, 1, MPI_UNSIGNED, (*citr)->rank(), 0,
                MPI_COMM_WORLD, &send_requests[i]);
      if (nparticles != 0) {
        mpm::HDF5Particle h5_particle;
        // Initialize MPI datatypes and send vector of particles
        MPI_Datatype particle_type =
            mpm::register_mpi_particle_type(h5_particle);
        MPI_Send(h5_particles.data(), nparticles, particle_type,
                 (*citr)->rank(), 0, MPI_COMM_WORLD);
        mpm::deregister_mpi_particle_type(particle_type);
      }
      h5_particles.clear();
    }
    // Remove all sent particles
    this->remove_particles(remove_pids);
    // Send complete
    for (unsigned i = 0; i < this->ghost_cells_.size(); ++i)
      MPI_Wait(&send_requests[i], MPI_STATUS_IGNORE);

    // Iterate through the local ghost cells and receive particles
    for (auto citr = this->local_ghost_cells_.cbegin();
         citr != this->local_ghost_cells_.cend(); ++citr) {
      std::vector<unsigned> neighbour_ranks =
          ghost_cells_neighbour_ranks_[(*citr)->id()];

      for (unsigned i = 0; i < neighbour_ranks.size(); ++i) {
        // Receive number of particles
        unsigned nrecv_particles;
        MPI_Recv(&nrecv_particles, 1, MPI_UNSIGNED, neighbour_ranks[i], 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (nrecv_particles != 0) {
          std::vector<mpm::HDF5Particle> recv_particles;
          recv_particles.resize(nrecv_particles);
          // Receive the vector of particles
          mpm::HDF5Particle received;
          MPI_Datatype particle_type =
              mpm::register_mpi_particle_type(received);
          MPI_Recv(recv_particles.data(), nrecv_particles, particle_type,
                   neighbour_ranks[i], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          mpm::deregister_mpi_particle_type(particle_type);

          // Iterate through n number of received particles
          for (const auto& rparticle : recv_particles) {
            mpm::Index id = 0;
            // Initial particle coordinates
            Eigen::Matrix<double, Tdim, 1> pcoordinates;
            pcoordinates.setZero();

            // Received particle
            auto received_particle =
                std::make_shared<mpm::Particle<Tdim>>(id, pcoordinates);
            // Get material
            auto material = materials_.at(rparticle.material_id);
            // Reinitialise particle from HDF5 data
            received_particle->initialise_particle(rparticle, material);

            // Add particle to mesh
            this->add_particle(received_particle, true);
          }
        }
      }
    }
  }
#endif
}

//! Transfer all particles in cells that are not in local rank
template <unsigned Tdim>
void mpm::Mesh<Tdim>::transfer_nonrank_particles(
    const std::vector<mpm::Index>& exchange_cells) {
#ifdef USE_MPI
  // Get number of MPI ranks
  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  if (mpi_size > 1) {
    std::vector<MPI_Request> send_requests;
    send_requests.reserve(exchange_cells.size());
    unsigned nsend_requests = 0;

    std::vector<mpm::Index> remove_pids;
    // Iterate through the ghost cells and send particles
    for (auto cid : exchange_cells) {
      // Get cell pointer
      auto cell = map_cells_[cid];
      // If the previous rank of cell is the current MPI rank,
      // then send all particles
      if ((cell->rank() != cell->previous_mpirank()) &&
          (cell->previous_mpirank() == mpi_rank)) {
        // Create a vector of h5_particles
        std::vector<mpm::HDF5Particle> h5_particles;
        auto particle_ids = cell->particles();
        // Create a vector of HDF5 data of particles to send
        // delete particle
        for (auto& id : particle_ids) {
          // Append to vector of particles
          h5_particles.emplace_back(map_particles_[id]->hdf5());
          // Particles to be removed from the current rank
          remove_pids.emplace_back(id);
        }
        cell->clear_particle_ids();

        // Send number of particles to receiver rank
        unsigned nparticles = particle_ids.size();
        MPI_Isend(&nparticles, 1, MPI_UNSIGNED, cell->rank(), 0, MPI_COMM_WORLD,
                  &send_requests[nsend_requests]);
        if (nparticles > 0) {
          mpm::HDF5Particle h5_particle;
          // Initialize MPI datatypes and send vector of particles
          MPI_Datatype particle_type =
              mpm::register_mpi_particle_type(h5_particle);
          MPI_Send(h5_particles.data(), nparticles, particle_type, cell->rank(),
                   0, MPI_COMM_WORLD);
          mpm::deregister_mpi_particle_type(particle_type);
        }
        h5_particles.clear();
        ++nsend_requests;
      }
    }
    // Remove all sent particles
    this->remove_particles(remove_pids);
    // Send complete iterate only upto valid send requests
    for (unsigned i = 0; i < nsend_requests; ++i)
      MPI_Wait(&send_requests[i], MPI_STATUS_IGNORE);

    // Iterate through the ghost cells and send particles
    for (auto cid : exchange_cells) {
      // Get cell pointer
      auto cell = map_cells_[cid];
      // If the current rank is the MPI rank receive particles
      if ((cell->rank() != cell->previous_mpirank()) &&
          (cell->rank() == mpi_rank)) {
        // Receive number of particles
        unsigned nrecv_particles;
        MPI_Recv(&nrecv_particles, 1, MPI_UNSIGNED, cell->previous_mpirank(), 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (nrecv_particles != 0) {
          std::vector<mpm::HDF5Particle> recv_particles;
          recv_particles.resize(nrecv_particles);
          // Receive the vector of particles
          mpm::HDF5Particle received;
          MPI_Datatype particle_type =
              mpm::register_mpi_particle_type(received);
          MPI_Recv(recv_particles.data(), nrecv_particles, particle_type,
                   cell->previous_mpirank(), 0, MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);
          mpm::deregister_mpi_particle_type(particle_type);

          // Iterate through n number of received particles
          for (const auto& rparticle : recv_particles) {
            mpm::Index id = 0;
            // Initial particle coordinates
            Eigen::Matrix<double, Tdim, 1> pcoordinates;
            pcoordinates.setZero();

            // Received particle
            auto received_particle =
                std::make_shared<mpm::Particle<Tdim>>(id, pcoordinates);
            // Get material
            auto material = materials_.at(rparticle.material_id);
            // Reinitialise particle from HDF5 data
            received_particle->initialise_particle(rparticle, material);

            // Add particle to mesh
            this->add_particle(received_particle, true);
          }
        }
      }
    }
  }
#endif
}

//! Find shared nodes across MPI domains
template <unsigned Tdim>
void mpm::Mesh<Tdim>::find_domain_shared_nodes() {
  // Clear MPI rank at the nodes
#pragma omp parallel for schedule(runtime)
  for (auto nitr = nodes_.cbegin(); nitr != nodes_.cend(); ++nitr)
    (*nitr)->clear_mpi_ranks();

  // Get MPI rank
  int mpi_rank = 0;
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif

#pragma omp parallel for schedule(runtime)
  // Assign MPI rank to nodes of cell
  for (auto citr = cells_.cbegin(); citr != cells_.cend(); ++citr)
    (*citr)->assign_mpi_rank_to_nodes();

  this->domain_shared_nodes_.clear();

#ifdef USE_HALO_EXCHANGE
  ncomms_ = 0;
  for (auto nitr = nodes_.cbegin(); nitr != nodes_.cend(); ++nitr) {
    // If node has more than 1 MPI rank
    std::set<unsigned> nodal_mpi_ranks = (*nitr)->mpi_ranks();
    const unsigned nodal_mpi_ranks_size = nodal_mpi_ranks.size();
    if (nodal_mpi_ranks_size > 1) {
      if (nodal_mpi_ranks.find(mpi_rank) != nodal_mpi_ranks.end()) {
        // Create Ghost ID
        (*nitr)->ghost_id(ncomms_);
        // Add to list of shared nodes on local rank
        domain_shared_nodes_.add(*nitr);
        ncomms_ += nodal_mpi_ranks_size - 1;
      }
    }
  }
#else
  nhalo_nodes_ = 0;
  for (auto nitr = nodes_.cbegin(); nitr != nodes_.cend(); ++nitr) {
    std::set<unsigned> nodal_mpi_ranks = (*nitr)->mpi_ranks();
    // If node has more than 1 MPI rank
    if (nodal_mpi_ranks.size() > 1) {
      (*nitr)->ghost_id(nhalo_nodes_);
      nhalo_nodes_ += 1;
      // Add to domain shared nodes only if active on current MPI rank
      if (nodal_mpi_ranks.find(mpi_rank) != nodal_mpi_ranks.end())
        domain_shared_nodes_.add(*nitr);
    }
  }
#endif
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
#pragma omp parallel for schedule(runtime)
  for (auto citr = cells_.cbegin(); citr != cells_.cend(); ++citr) {
    // Check if particle is already found, if so don't run for other cells
    // Check if co-ordinates is within the cell, if true
    // add particle to cell
    Eigen::Matrix<double, Tdim, 1> xi;
    if (!status && (*citr)->is_point_in_cell(particle->coordinates(), &xi)) {
      particle->assign_cell_xi(*citr, xi);
      status = true;
    }
  }

  return status;
}

//! Iterate over particles
template <unsigned Tdim>
template <typename Toper>
void mpm::Mesh<Tdim>::iterate_over_particles(Toper oper) {
#pragma omp parallel for schedule(runtime)
  for (auto pitr = particles_.cbegin(); pitr != particles_.cend(); ++pitr)
    oper(*pitr);
}

//! Iterate over particles
template <unsigned Tdim>
template <typename Toper, typename Tpred>
void mpm::Mesh<Tdim>::iterate_over_particles_predicate(Toper oper, Tpred pred) {
#pragma omp parallel for schedule(runtime)
  for (auto pitr = particles_.cbegin(); pitr != particles_.cend(); ++pitr) {
    if (pred(*pitr)) oper(*pitr);
  }
}

//! Iterate over particle set
template <unsigned Tdim>
template <typename Toper>
void mpm::Mesh<Tdim>::iterate_over_particle_set(int set_id, Toper oper) {
  // If set id is -1, use all particles
  if (set_id == -1) {
    this->iterate_over_particles(oper);
  } else {
    // Iterate over the particle set
    auto set = particle_sets_.at(set_id);
#pragma omp parallel for schedule(runtime)
    for (auto sitr = set.begin(); sitr != set.cend(); ++sitr) {
      unsigned pid = (*sitr);
      if (map_particles_.find(pid) != map_particles_.end())
        oper(map_particles_[pid]);
    }
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

//! Return particle tensor data
template <unsigned Tdim>
template <unsigned Tsize>
std::vector<Eigen::Matrix<double, Tsize, 1>>
    mpm::Mesh<Tdim>::particles_tensor_data(const std::string& attribute) {
  std::vector<Eigen::Matrix<double, Tsize, 1>> tensor_data;
  try {
    // Iterate over particles
    for (auto pitr = particles_.cbegin(); pitr != particles_.cend(); ++pitr) {
      Eigen::Matrix<double, Tsize, 1> data;
      data.setZero();
      auto pdata = (*pitr)->tensor_data(attribute);
      // Fill stresses to the size of dimensions
      for (unsigned i = 0; i < pdata.size(); ++i) data(i) = pdata(i);

      // Add to a tensor of data
      tensor_data.emplace_back(data);
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {} {}\n", __FILE__, __LINE__, exception.what(),
                    attribute);
    tensor_data.clear();
  }
  return tensor_data;
}

//! Return particle scalar data
template <unsigned Tdim>
std::vector<double> mpm::Mesh<Tdim>::particles_statevars_data(
    const std::string& attribute) {
  std::vector<double> scalar_data;
  scalar_data.reserve(particles_.size());
  // Iterate over particles and add scalar value to data
  for (auto pitr = particles_.cbegin(); pitr != particles_.cend(); ++pitr)
    scalar_data.emplace_back((*pitr)->state_variable(attribute));
  return scalar_data;
}

//! Assign particles volumes
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::assign_particles_volumes(
    const std::vector<std::tuple<mpm::Index, double>>& particle_volumes) {
  bool status = true;
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

//! Create particle tractions
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::create_particles_tractions(
    const std::shared_ptr<FunctionBase>& mfunction, int set_id, unsigned dir,
    double traction) {
  bool status = true;
  try {
    if (set_id == -1 || particle_sets_.find(set_id) != particle_sets_.end())
      // Create a particle traction load
      particle_tractions_.emplace_back(
          std::make_shared<mpm::Traction>(set_id, mfunction, dir, traction));
    else
      throw std::runtime_error("No particle set found to assign traction");

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Apply particle tractions
template <unsigned Tdim>
void mpm::Mesh<Tdim>::apply_traction_on_particles(double current_time) {
  // Iterate over all particle tractions
  for (const auto& ptraction : particle_tractions_) {
    int set_id = ptraction->setid();
    unsigned dir = ptraction->dir();
    double traction = ptraction->traction(current_time);
    this->iterate_over_particle_set(
        set_id, std::bind(&mpm::ParticleBase<Tdim>::assign_traction,
                          std::placeholders::_1, dir, traction));
  }
  if (!particle_tractions_.empty()) {
    this->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::map_traction_force, std::placeholders::_1));
  }
}

//! Create particle velocity constraints
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::create_particle_velocity_constraint(
    int set_id, const std::shared_ptr<mpm::VelocityConstraint>& constraint) {
  bool status = true;
  try {
    if (set_id == -1 || particle_sets_.find(set_id) != particle_sets_.end()) {
      // Create a particle velocity constraint
      if (constraint->dir() < Tdim)
        particle_velocity_constraints_.emplace_back(constraint);
      else
        throw std::runtime_error("Invalid direction of velocity constraint");
    } else
      throw std::runtime_error(
          "No particle set found to assign velocity constraint");

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Apply particle tractions
template <unsigned Tdim>
void mpm::Mesh<Tdim>::apply_particle_velocity_constraints() {
  // Iterate over all particle velocity constraints
  for (const auto& pvelocity : particle_velocity_constraints_) {
    // If set id is -1, use all particles
    int set_id = pvelocity->setid();
    unsigned dir = pvelocity->dir();
    double velocity = pvelocity->velocity();

    this->iterate_over_particle_set(
        set_id,
        std::bind(&mpm::ParticleBase<Tdim>::apply_particle_velocity_constraints,
                  std::placeholders::_1, dir, velocity));
  }
}

//! Assign node tractions
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::assign_nodal_concentrated_forces(
    const std::shared_ptr<FunctionBase>& mfunction, int set_id, unsigned dir,
    double concentrated_force) {
  bool status = true;
  // TODO: Remove phase
  const unsigned phase = 0;
  try {
    if (!nodes_.size())
      throw std::runtime_error(
          "No nodes have been assigned in mesh, cannot assign concentrated "
          "force");

    // Set id of -1, is all nodes
    Vector<NodeBase<Tdim>> nodes =
        (set_id == -1) ? this->nodes_ : node_sets_.at(set_id);

#pragma omp parallel for schedule(runtime)
    for (auto nitr = nodes.cbegin(); nitr != nodes.cend(); ++nitr) {
      if (!(*nitr)->assign_concentrated_force(phase, dir, concentrated_force,
                                              mfunction))
        throw std::runtime_error("Setting concentrated force failed");
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
  hsize_t nrecords = 0;
  hsize_t nfields = 0;
  H5TBget_table_info(file_id, "table", &nfields, &nrecords);

  if (nfields != mpm::hdf5::particle::NFIELDS)
    throw std::runtime_error("HDF5 table has incorrect number of fields");

  std::vector<HDF5Particle> dst_buf;
  dst_buf.reserve(nrecords);
  // Read the table
  H5TBread_table(file_id, "table", mpm::hdf5::particle::dst_size,
                 mpm::hdf5::particle::dst_offset,
                 mpm::hdf5::particle::dst_sizes, dst_buf.data());

  // Vector of particles
  Vector<ParticleBase<Tdim>> particles;

  // Clear map of particles
  map_particles_.clear();

  unsigned i = 0;
  for (auto pitr = particles_.cbegin(); pitr != particles_.cend(); ++pitr) {
    if (i < nrecords) {
      HDF5Particle particle = dst_buf[i];
      // Get particle's material from list of materials
      auto material = materials_.at(particle.material_id);
      // Initialise particle with HDF5 data
      (*pitr)->initialise_particle(particle, material);
      // Add particle to map
      map_particles_.insert(particle.id, *pitr);
      particles.add(*pitr);
      ++i;
    }
  }
  // close the file
  H5Fclose(file_id);

  // Overwrite particles container
  this->particles_ = particles;

  // Remove associated cell for the particle
  for (auto citr = this->cells_.cbegin(); citr != this->cells_.cend(); ++citr)
    (*citr)->clear_particle_ids();

  return true;
}

//! Write particles to HDF5
template <unsigned Tdim>
std::vector<mpm::HDF5Particle> mpm::Mesh<Tdim>::particles_hdf5() const {
  const unsigned nparticles = this->nparticles();

  std::vector<mpm::HDF5Particle> particles_hdf5;
  particles_hdf5.reserve(nparticles);

  for (auto pitr = particles_.cbegin(); pitr != particles_.cend(); ++pitr)
    particles_hdf5.emplace_back((*pitr)->hdf5());

  return particles_hdf5;
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
std::vector<std::array<mpm::Index, 2>> mpm::Mesh<Tdim>::node_pairs(
    bool active) const {
  // Vector of node_pairs
  std::vector<std::array<mpm::Index, 2>> node_pairs;

  try {
    int mpi_rank = 0;
#ifdef USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif
    if (cells_.size() == 0)
      throw std::runtime_error("No cells have been initialised!");

    for (auto citr = cells_.cbegin(); citr != cells_.cend(); ++citr) {
      // If node pairs are only requested for active nodes
      bool get_pairs = (active == true) ? ((*citr)->rank() == mpi_rank) : true;
      if (get_pairs) {
        const auto pairs = (*citr)->side_node_pairs();
        node_pairs.insert(std::end(node_pairs), std::begin(pairs),
                          std::end(pairs));
      }
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
      std::vector<mpm::Index> particles((sitr->second).begin(),
                                        (sitr->second).end());

      // Create the map of the container
      status = this->particle_sets_
                   .insert(std::pair<mpm::Index, std::vector<mpm::Index>>(
                       sitr->first, particles))
                   .second;
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return status;
}

//! Create map of container of nodes in sets
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::create_node_sets(
    const tsl::robin_map<mpm::Index, std::vector<mpm::Index>>& node_sets,
    bool check_duplicates) {
  bool status = false;
  try {
    // Create container for each node set
    for (auto sitr = node_sets.begin(); sitr != node_sets.end(); ++sitr) {
      // Create a vector for the set
      Vector<NodeBase<Tdim>> nodes;
      // Reserve the size of the container
      nodes.reserve((sitr->second).size());
      // Add nodes to the container
      for (auto pid : sitr->second) {
        nodes.add(map_nodes_[pid], check_duplicates);
      }

      // Create the map of the vector
      status = this->node_sets_
                   .insert(std::pair<mpm::Index, Vector<NodeBase<Tdim>>>(
                       sitr->first, nodes))
                   .second;
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return status;
}

// Return cells
template <unsigned Tdim>
mpm::Vector<mpm::Cell<Tdim>> mpm::Mesh<Tdim>::cells() {
  return this->cells_;
}

//! Create map of container of cells in sets
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::create_cell_sets(
    const tsl::robin_map<mpm::Index, std::vector<mpm::Index>>& cell_sets,
    bool check_duplicates) {
  bool status = false;
  try {
    // Create container for each cell set
    for (auto sitr = cell_sets.begin(); sitr != cell_sets.end(); ++sitr) {
      // Create a container for the set
      Vector<Cell<Tdim>> cells;
      // Reserve the size of the container
      cells.reserve((sitr->second).size());
      // Add cells to the container
      for (auto pid : sitr->second) {
        cells.add(map_cells_[pid], check_duplicates);
      }

      // Create the map of the container
      status = this->cell_sets_
                   .insert(std::pair<mpm::Index, Vector<Cell<Tdim>>>(
                       sitr->first, cells))
                   .second;
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return status;
}

//! return particle_ptr
template <unsigned Tdim>
std::map<mpm::Index, mpm::Index>* mpm::Mesh<Tdim>::particles_cell_ids() {
  return &(this->particles_cell_ids_);
}

//! Generate particles
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::generate_particles(const std::shared_ptr<mpm::IO>& io,
                                         const Json& generator) {
  bool status = true;
  try {
    // Particle generator
    const auto generator_type = generator["type"].template get<std::string>();

    // Generate particles from file
    if (generator_type == "file") {
      // Particle set id
      unsigned pset_id = generator["pset_id"].template get<unsigned>();
      status = this->read_particles_file(io, generator, pset_id);
    }

    // Generate material points at the Gauss location in all cells
    else if (generator_type == "gauss") {
      // Number of particles per dir
      unsigned nparticles_dir =
          generator["nparticles_per_dir"].template get<unsigned>();
      // Particle type
      auto particle_type =
          generator["particle_type"].template get<std::string>();
      // Material id
      std::vector<unsigned> material_ids;
      if (generator.at("material_id").is_array())
        material_ids =
            generator["material_id"].template get<std::vector<unsigned>>();
      else
        material_ids.emplace_back(
            generator["material_id"].template get<unsigned>());
      // Cell set id
      int cset_id = generator["cset_id"].template get<int>();
      // Particle set id
      unsigned pset_id = generator["pset_id"].template get<unsigned>();
      status = this->generate_material_points(nparticles_dir, particle_type,
                                              material_ids, cset_id, pset_id);
    }

    // Generate material points at the Gauss location in all cells
    else if (generator_type == "inject") {
      mpm::Injection inject;
      // Number of particles per dir
      inject.nparticles_dir =
          generator["nparticles_per_dir"].template get<unsigned>();
      // Particle type
      inject.particle_type =
          generator["particle_type"].template get<std::string>();
      // Material id
      if (generator.at("material_id").is_array())
        inject.material_ids =
            generator["material_id"].template get<std::vector<unsigned>>();
      else
        inject.material_ids.emplace_back(
            generator["material_id"].template get<unsigned>());
      // Cell set id
      inject.cell_set_id = generator["cset_id"].template get<int>();
      // Duration of injection
      if (generator.contains("duration") && generator["duration"].is_array() &&
          generator["duration"].size() == 2) {
        inject.start_time = generator["duration"].at(0);
        inject.end_time = generator["duration"].at(1);
      }

      // Velocity
      inject.velocity.resize(Tdim, 0.);
      if (generator["velocity"].is_array() &&
          generator["velocity"].size() == Tdim) {
        for (unsigned i = 0; i < Tdim; ++i)
          inject.velocity[i] = generator["velocity"].at(i);
      }
      // Add to particle injections
      particle_injections_.emplace_back(inject);
    }

    else
      throw std::runtime_error(
          "Particle generator type is not properly specified");

  } catch (std::exception& exception) {
    console_->error("{}: #{} Generating particle failed", __FILE__, __LINE__);
    status = false;
  }
  return status;
}

//! Generate particles
template <unsigned Tdim>
void mpm::Mesh<Tdim>::inject_particles(double current_time) {
  int mpi_rank = 0;
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif
  // Container of new injected particles
  std::vector<std::shared_ptr<ParticleBase<Tdim>>> injected_particles;
  // Iterate over all injection cells
  for (auto injection : particle_injections_) {
    unsigned pid = this->nparticles();
    bool checks = false;
    // Get material
    std::vector<std::shared_ptr<mpm::Material<Tdim>>> materials;
    for (auto m_id : injection.material_ids)
      materials.emplace_back(materials_.at(m_id));

    // Check if duration is within the current time
    if (injection.start_time <= current_time &&
        injection.end_time > current_time) {
      // If set id is -1, use all cells
      auto cset = (injection.cell_set_id == -1)
                      ? this->cells_
                      : cell_sets_.at(injection.cell_set_id);
      // Iterate over each cell to generate points
      for (auto citr = cset.cbegin(); citr != cset.cend(); ++citr) {
        if ((*citr)->rank() == mpi_rank && (*citr)->nparticles() == 0) {
          // Assign quadratures based on number of particles
          (*citr)->assign_quadrature(injection.nparticles_dir);

          // Genereate particles at the Gauss points
          const auto cpoints = (*citr)->generate_points();
          // Iterate over each coordinate to generate material points
          for (const auto& coordinates : cpoints) {
            // Create particle
            auto particle =
                Factory<mpm::ParticleBase<Tdim>, mpm::Index,
                        const Eigen::Matrix<double, Tdim, 1>&>::instance()
                    ->create(injection.particle_type,
                             static_cast<mpm::Index>(pid), coordinates);

            // particle velocity
            Eigen::Matrix<double, Tdim, 1> pvelocity(injection.velocity.data());
            particle->assign_velocity(pvelocity);

            // Add particle to mesh
            unsigned status = this->add_particle(particle, checks);
            if (status) {
              map_particles_[pid]->assign_cell(*citr);
              for (unsigned phase = 0; phase < materials.size(); phase++)
                map_particles_[pid]->assign_material(materials[phase], phase);
              ++pid;
              injected_particles.emplace_back(particle);
            }
          }
        }
      }
    }
    for (auto particle : injected_particles) {
      particle->compute_volume();
      particle->compute_mass();
    }
  }
}

// Read particles file
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::read_particles_file(const std::shared_ptr<mpm::IO>& io,
                                          const Json& generator,
                                          unsigned pset_id) {
  // Particle type
  auto particle_type = generator["particle_type"].template get<std::string>();

  // File location
  auto file_loc =
      io->file_name(generator["location"].template get<std::string>());

  // Check duplicates
  bool check_duplicates = generator["check_duplicates"].template get<bool>();

  // Material id
  std::vector<unsigned> material_ids;
  if (generator.at("material_id").is_array())
    material_ids =
        generator["material_id"].template get<std::vector<unsigned>>();
  else
    material_ids.emplace_back(
        generator["material_id"].template get<unsigned>());

  const std::string reader = generator["io_type"].template get<std::string>();

  // Create a particle reader
  auto particle_io = Factory<mpm::IOMesh<Tdim>>::instance()->create(reader);

  // Get coordinates
  auto coords = particle_io->read_particles(file_loc);

  // Create particles from coordinates
  bool status = this->create_particles(particle_type, coords, material_ids,
                                       pset_id, check_duplicates);

  if (!status) throw std::runtime_error("Addition of particles to mesh failed");

  return status;
}

//! Assign nodal concentrated force
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::assign_nodal_concentrated_forces(
    const std::vector<std::tuple<mpm::Index, unsigned, double>>& nodal_forces) {
  bool status = true;
  // TODO: Remove phase
  const unsigned phase = 0;
  try {
    if (!nodes_.size())
      throw std::runtime_error(
          "No nodes have been assigned in mesh, cannot assign traction");
    for (const auto& nodal_force : nodal_forces) {
      // Node id
      mpm::Index pid = std::get<0>(nodal_force);
      // Direction
      unsigned dir = std::get<1>(nodal_force);
      // Force
      double force = std::get<2>(nodal_force);

      if (map_nodes_.find(pid) != map_nodes_.end())
        status = map_nodes_[pid]->assign_concentrated_force(phase, dir, force,
                                                            nullptr);

      if (!status) throw std::runtime_error("Force is invalid for node");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Create the nodal properties' map
template <unsigned Tdim>
void mpm::Mesh<Tdim>::create_nodal_properties() {
  // Initialise the shared pointer to nodal properties
  nodal_properties_ = std::make_shared<mpm::NodalProperties>();

  // Check if nodes_ and materials_is empty and throw runtime error if they are
  if (nodes_.size() != 0 && materials_.size() != 0) {
    // Compute number of rows in nodal properties for vector entities
    const unsigned nrows = nodes_.size() * Tdim;
    // Create pool data for each property in the nodal properties struct
    // object. Properties must be named in the plural form
    nodal_properties_->create_property("masses", nodes_.size(),
                                       materials_.size());
    nodal_properties_->create_property("momenta", nrows, materials_.size());
    nodal_properties_->create_property("change_in_momenta", nrows,
                                       materials_.size());
    nodal_properties_->create_property("displacements", nrows,
                                       materials_.size());
    nodal_properties_->create_property("separation_vectors", nrows,
                                       materials_.size());
    nodal_properties_->create_property("domain_gradients", nrows,
                                       materials_.size());
    nodal_properties_->create_property("normal_unit_vectors", nrows,
                                       materials_.size());

    // Iterate over all nodes to initialise the property handle in each node
    // and assign its node id as the prop id in the nodal property data pool
    for (auto nitr = nodes_.cbegin(); nitr != nodes_.cend(); ++nitr)
      (*nitr)->initialise_property_handle((*nitr)->id(), nodal_properties_);
  } else {
    throw std::runtime_error("Number of nodes or number of materials is zero");
  }
}

//! Compute nodal correction force
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::compute_nodal_correction_force(
    const Eigen::SparseMatrix<double>& correction_matrix,
    const Eigen::VectorXd& pressure_increment, double dt) {
  bool status = true;
  try {
    //! Active node size
    const auto nactive_node = active_nodes_.size();

    // Part of nodal correction force of one direction
    Eigen::MatrixXd correction_force;
    correction_force.resize(nactive_node, Tdim);

    // Iterate over each direction
    for (unsigned i = 0; i < Tdim; ++i) {
      correction_force.block(0, i, nactive_node, 1) =
          -correction_matrix.block(0, nactive_node * i, nactive_node,
                                   nactive_node) *
          pressure_increment;
    }

    // Iterate over each active node
    VectorDim nodal_correction_force;
    for (auto nitr = active_nodes_.cbegin(); nitr != active_nodes_.cend();
         ++nitr) {
      unsigned active_id = (*nitr)->active_id();
      nodal_correction_force = (correction_force.row(active_id)).transpose();

      // Compute correction force for each node
      map_nodes_[(*nitr)->id()]->compute_nodal_correction_force(
          nodal_correction_force);
    }

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return status;
};

//! Compute free surface cells, nodes, and particles
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::compute_free_surface(double tolerance) {
  bool status = true;
  try {
    // Reset free surface cell and particles
    this->iterate_over_cells(std::bind(&mpm::Cell<Tdim>::assign_free_surface,
                                       std::placeholders::_1, false));

    VectorDim temp_normal;
    temp_normal.setZero();
    this->iterate_over_particles_predicate(
        std::bind(&mpm::ParticleBase<Tdim>::assign_normal,
                  std::placeholders::_1, temp_normal),
        std::bind(&mpm::ParticleBase<Tdim>::free_surface,
                  std::placeholders::_1));

    this->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::assign_free_surface,
                  std::placeholders::_1, false));

    // Reset volume fraction
    this->iterate_over_cells(std::bind(&mpm::Cell<Tdim>::assign_volume_fraction,
                                       std::placeholders::_1, 0.0));

    // Compute and assign volume fraction to each cell
    for (auto citr = this->cells_.cbegin(); citr != this->cells_.cend();
         ++citr) {
      if ((*citr)->status()) {
        // Compute volume fraction
        double cell_volume_fraction = 0.0;
        for (const auto p_id : (*citr)->particles())
          cell_volume_fraction += map_particles_[p_id]->volume();

        cell_volume_fraction = cell_volume_fraction / (*citr)->volume();
        (*citr)->assign_volume_fraction(cell_volume_fraction);
      }
    }

    // First, we detect the cell with possible free surfaces
    // Compute boundary cells and nodes based on geometry
    std::set<mpm::Index> free_surface_candidate_cells;
    std::set<mpm::Index> free_surface_candidate_nodes;
    for (auto citr = this->cells_.cbegin(); citr != this->cells_.cend();
         ++citr) {
      // Cell contains particles
      if ((*citr)->status()) {
        bool candidate_cell = false;
        const auto& node_id = (*citr)->nodes_id();
        if ((*citr)->volume_fraction() < tolerance) {
          candidate_cell = true;
          for (const auto id : node_id) {
            map_nodes_[id]->assign_free_surface(true);
          }
        } else {
          // Loop over neighbouring cells
          for (const auto n_id : (*citr)->neighbours()) {
            if (!map_cells_[n_id]->status()) {
              candidate_cell = true;
              const auto& n_node_id = map_cells_[n_id]->nodes_id();

              // Detect common node id
              std::set<mpm::Index> common_node_id;
              std::set_intersection(
                  node_id.begin(), node_id.end(), n_node_id.begin(),
                  n_node_id.end(),
                  std::inserter(common_node_id, common_node_id.begin()));

              // Assign free surface nodes
              if (!common_node_id.empty()) {
                for (const auto common_id : common_node_id) {
                  map_nodes_[common_id]->assign_free_surface(true);
                }
              }
            }
          }
        }

        // Assign free surface cell
        if (candidate_cell) {
          (*citr)->assign_free_surface(true);
          free_surface_candidate_cells.insert((*citr)->id());
          free_surface_candidate_nodes.insert(node_id.begin(), node_id.end());
        }
      }
    }

    // Compute particle neighbours for particles at candidate cells
    std::vector<mpm::Index> free_surface_candidate_particles;
    for (const auto cell_id : free_surface_candidate_cells) {
      this->find_particle_neighbours(map_cells_[cell_id]);
      const auto& particle_ids = map_cells_[cell_id]->particles();
      free_surface_candidate_particles.insert(
          free_surface_candidate_particles.end(), particle_ids.begin(),
          particle_ids.end());
    }

    // Compute boundary particles based on density function
    // Lump cell volume to nodes
    this->iterate_over_cells(std::bind(
        &mpm::Cell<Tdim>::map_cell_volume_to_nodes, std::placeholders::_1, 0));

    // Compute nodal value of mass density
    this->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_density, std::placeholders::_1),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    std::set<mpm::Index> free_surface_candidate_particles_second;
    for (const auto p_id : free_surface_candidate_particles) {
      const auto& particle = map_particles_[p_id];
      bool status = particle->compute_free_surface();
      if (status) free_surface_candidate_particles_second.insert(p_id);
    }

    // Find free surface particles through geometry
    std::set<mpm::Index> free_surface_particles;
    for (const auto p_id : free_surface_candidate_particles_second) {
      // Initialize renormalization matrix
      Eigen::Matrix<double, Tdim, Tdim> renormalization_matrix_inv;
      renormalization_matrix_inv.setZero();

      // Loop over neighbours
      const auto& particle = map_particles_[p_id];
      const auto& p_coord = particle->coordinates();
      const auto& neighbour_particles = particle->neighbours();
      const double smoothing_length = 1.33 * particle->diameter();
      for (const auto n_id : neighbour_particles) {
        const auto& n_coord = map_particles_[n_id]->coordinates();
        const VectorDim rel_coord = n_coord - p_coord;

        // Compute kernel gradient
        const VectorDim kernel_gradient =
            mpm::RadialBasisFunction::gradient<Tdim>(smoothing_length,
                                                     -rel_coord, "gaussian");

        // Inverse of renormalization matrix B
        renormalization_matrix_inv +=
            (particle->mass() / particle->mass_density()) * kernel_gradient *
            rel_coord.transpose();
      }

      // Compute lambda: minimum eigenvalue of B_inverse
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(
          renormalization_matrix_inv);
      double lambda = es.eigenvalues().minCoeff();

      // Categorize particle based on lambda
      bool free_surface = false;
      bool secondary_check = false;
      bool interior = false;
      if (lambda <= 0.2)
        free_surface = true;
      else if (lambda > 0.2 && lambda <= 0.75)
        secondary_check = true;
      else
        interior = true;

      // Compute numerical normal vector
      VectorDim normal;
      normal.setZero();
      if (!interior) {
        VectorDim temporary_vec;
        temporary_vec.setZero();
        for (const auto n_id : neighbour_particles) {
          const auto& n_coord = map_particles_[n_id]->coordinates();
          const VectorDim rel_coord = n_coord - p_coord;

          // Compute kernel gradient
          const VectorDim kernel_gradient =
              mpm::RadialBasisFunction::gradient<Tdim>(smoothing_length,
                                                       -rel_coord, "gaussian");

          // Sum of kernel by volume
          temporary_vec +=
              (particle->mass() / particle->mass_density()) * kernel_gradient;
        }
        normal = -renormalization_matrix_inv.inverse() * temporary_vec;
        if (normal.norm() > std::numeric_limits<double>::epsilon())
          normal.normalize();
        else
          normal.setZero();
      }

      // If secondary check is needed
      if (secondary_check) {
        // Construct scanning region
        // TODO: spacing distance should be a function of porosity
        const double spacing_distance = smoothing_length;
        VectorDim t_coord = p_coord + spacing_distance * normal;

        // Check all neighbours
        for (const auto n_id : neighbour_particles) {
          const auto& n_coord = map_particles_[n_id]->coordinates();
          const VectorDim rel_coord_np = n_coord - p_coord;
          const double distance_np = rel_coord_np.norm();
          const VectorDim rel_coord_nt = n_coord - t_coord;
          const double distance_nt = rel_coord_nt.norm();

          free_surface = true;
          if (distance_np < std::sqrt(2) * spacing_distance) {
            if (std::acos(normal.dot(rel_coord_np) / distance_np) < M_PI / 4) {
              free_surface = false;
              break;
            }
          } else {
            if (distance_nt < spacing_distance) {
              free_surface = false;
              break;
            }
          }
        }
      }

      // Assign normal only to validated free surface
      if (free_surface) {
        particle->assign_free_surface(true);
        particle->assign_normal(normal);
        free_surface_particles.insert(p_id);
      }
    }

    // Compute node sign distance function
    for (const auto node_id : free_surface_candidate_nodes) {
      const auto& node_coord = map_nodes_[node_id]->coordinates();
      double closest_distance = std::numeric_limits<double>::max();
      double signed_distance = std::numeric_limits<double>::max();
      for (const auto fs_id : free_surface_particles) {
        const auto& fs_particle = map_particles_[fs_id];
        const VectorDim fs_coord =
            fs_particle->coordinates() +
            0.5 * fs_particle->diameter() * fs_particle->normal();
        const VectorDim rel_coord = fs_coord - node_coord;
        const double distance = rel_coord.norm();
        if (distance < closest_distance) {
          closest_distance = distance;
          signed_distance = (rel_coord).dot(fs_particle->normal());
        }
      }

      // Assign signed distance to node
      map_nodes_[node_id]->assign_signed_distance(signed_distance);
    }

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return status;
}

// //! Compute free surface cells, nodes, and particles
// template <unsigned Tdim>
// bool mpm::Mesh<Tdim>::compute_free_surface(double tolerance) {
//   bool status = true;
//   try {
//     // Reset free surface cell
//     this->iterate_over_cells(std::bind(&mpm::Cell<Tdim>::assign_free_surface,
//                                        std::placeholders::_1, false));

//     // Reset volume fraction
//     this->iterate_over_cells(std::bind(&mpm::Cell<Tdim>::assign_volume_fraction,
//                                        std::placeholders::_1, 0.0));

//     // Reset free surface particle
//     this->iterate_over_particles(
//         std::bind(&mpm::ParticleBase<Tdim>::assign_free_surface,
//                   std::placeholders::_1, false));

//     // Compute and assign volume fraction to each cell
//     for (auto citr = this->cells_.cbegin(); citr != this->cells_.cend();
//          ++citr) {
//       if ((*citr)->status()) {
//         // Compute volume fraction
//         double cell_volume_fraction = 0.0;
//         for (const auto p_id : (*citr)->particles())
//           cell_volume_fraction += map_particles_[p_id]->volume();

//         cell_volume_fraction = cell_volume_fraction / (*citr)->volume();
//         (*citr)->assign_volume_fraction(cell_volume_fraction);
//       }
//     }

//     // Compute boundary cells and nodes based on geometry
//     std::set<mpm::Index> boundary_cells;
//     std::set<mpm::Index> boundary_nodes;
//     for (auto citr = this->cells_.cbegin(); citr != this->cells_.cend();
//          ++citr) {

//       if ((*citr)->status()) {
//         bool cell_at_interface = false;
//         const auto& node_id = (*citr)->nodes_id();
//         bool internal = true;

//         //! Check internal cell
//         for (const auto c_id : (*citr)->neighbours()) {
//           if (!map_cells_[c_id]->status()) {
//             internal = false;
//             break;
//           }
//         }

//         //! Check volume fraction only for boundary cell
//         if (!internal) {
//           if ((*citr)->volume_fraction() < tolerance) {
//             cell_at_interface = true;
//             for (const auto id : node_id) {
//               map_nodes_[id]->assign_free_surface(cell_at_interface);
//               boundary_nodes.insert(id);
//             }
//           } else {
//             for (const auto n_id : (*citr)->neighbours()) {
//               if (map_cells_[n_id]->volume_fraction() < tolerance) {
//                 cell_at_interface = true;
//                 const auto& n_node_id = map_cells_[n_id]->nodes_id();

//                 // Detect common node id
//                 std::set<mpm::Index> common_node_id;
//                 std::set_intersection(
//                     node_id.begin(), node_id.end(), n_node_id.begin(),
//                     n_node_id.end(),
//                     std::inserter(common_node_id, common_node_id.begin()));

//                 // Assign free surface nodes
//                 if (!common_node_id.empty()) {
//                   for (const auto common_id : common_node_id) {
//                     map_nodes_[common_id]->assign_free_surface(
//                         cell_at_interface);
//                     boundary_nodes.insert(common_id);
//                   }
//                 }
//               }
//             }
//           }

//           // Assign free surface cell
//           if (cell_at_interface) {
//             (*citr)->assign_free_surface(cell_at_interface);
//             boundary_cells.insert((*citr)->id());
//           }
//         }
//       }
//     }

//     // Compute boundary particles based on density function
//     // Lump cell volume to nodes
//     this->iterate_over_cells(std::bind(
//         &mpm::Cell<Tdim>::map_cell_volume_to_nodes, std::placeholders::_1,
//         0));

//     // Compute nodal value of mass density
//     this->iterate_over_nodes_predicate(
//         std::bind(&mpm::NodeBase<Tdim>::compute_density,
//         std::placeholders::_1), std::bind(&mpm::NodeBase<Tdim>::status,
//         std::placeholders::_1));

//     // Evaluate free surface particles
//     std::set<mpm::Index> boundary_particles;
//     for (auto pitr = this->particles_.cbegin(); pitr !=
//     this->particles_.cend();
//          ++pitr) {
//       bool status = (*pitr)->compute_free_surface();
//       if (status) {
//         (*pitr)->assign_free_surface(status);
//         boundary_particles.insert((*pitr)->id());
//       }
//     }

//   } catch (std::exception& exception) {
//     console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
//   }
//   return status;
// }

//! Get free surface node set
template <unsigned Tdim>
std::set<mpm::Index> mpm::Mesh<Tdim>::free_surface_nodes() {
  std::set<mpm::Index> id_set;
  for (auto nitr = this->nodes_.cbegin(); nitr != this->nodes_.cend(); ++nitr)
    if ((*nitr)->free_surface()) id_set.insert((*nitr)->id());
  return id_set;
}

//! Get free surface cell set
template <unsigned Tdim>
std::set<mpm::Index> mpm::Mesh<Tdim>::free_surface_cells() {
  std::set<mpm::Index> id_set;
  for (auto citr = this->cells_.cbegin(); citr != this->cells_.cend(); ++citr)
    if ((*citr)->free_surface()) id_set.insert((*citr)->id());
  return id_set;
}

//! Get free surface particle set
template <unsigned Tdim>
std::set<mpm::Index> mpm::Mesh<Tdim>::free_surface_particles() {
  std::set<mpm::Index> id_set;
  for (auto pitr = this->particles_.cbegin(); pitr != this->particles_.cend();
       ++pitr)
    if ((*pitr)->free_surface()) id_set.insert((*pitr)->id());
  return id_set;
}

// Initialise the nodal properties' map
template <unsigned Tdim>
void mpm::Mesh<Tdim>::initialise_nodal_properties() {
  // Call initialise_properties function from the nodal properties
  nodal_properties_->initialise_nodal_properties();
}

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
    unsigned np = 0;
    std::vector<mpm::Index> remove_pids;
    // Iterate through the ghost cells and send particles
    for (auto citr = this->ghost_cells_.cbegin();
         citr != this->ghost_cells_.cend(); ++citr, ++i) {

      // Send number of particles to receiver rank
      auto particle_ids = (*citr)->particles();
      unsigned nparticles = particle_ids.size();
      MPI_Isend(&nparticles, 1, MPI_UNSIGNED, (*citr)->rank(), 1,
                MPI_COMM_WORLD, &send_requests[i]);
    }

    // Iterate through the ghost cells and send particles
    for (auto citr = this->ghost_cells_.cbegin();
         citr != this->ghost_cells_.cend(); ++citr, ++i) {
      // Send number of particles to receiver rank
      auto particle_ids = (*citr)->particles();
      for (auto& id : particle_ids) {
        // Create a vector of serialized particle
        std::vector<uint8_t> buffer = map_particles_[id]->serialize();
        MPI_Send(buffer.data(), buffer.size(), MPI_UINT8_T, (*citr)->rank(), 0,
                 MPI_COMM_WORLD);
        ++np;
        // Particles to be removed from the current rank
        remove_pids.emplace_back(id);
      }
      (*citr)->clear_particle_ids();
    }
    // Remove all sent particles
    this->remove_particles(remove_pids);
    // Send complete
    for (unsigned i = 0; i < this->ghost_cells_.size(); ++i)
      MPI_Wait(&send_requests[i], MPI_STATUS_IGNORE);

    // Particle id
    mpm::Index pid = 0;
    // Initial particle coordinates
    Eigen::Matrix<double, Tdim, 1> pcoordinates =
        Eigen::Matrix<double, Tdim, 1>::Zero();

    // Iterate through the local ghost cells and receive particles
    for (auto citr = this->local_ghost_cells_.cbegin();
         citr != this->local_ghost_cells_.cend(); ++citr) {
      std::vector<unsigned> neighbour_ranks =
          ghost_cells_neighbour_ranks_[(*citr)->id()];
      // Total number of particles
      std::vector<unsigned> nrank_particles(neighbour_ranks.size(), 0);
      for (unsigned i = 0; i < neighbour_ranks.size(); ++i)
        MPI_Recv(&nrank_particles[i], 1, MPI_UNSIGNED, neighbour_ranks[i], 1,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      // Receive number of particles
      unsigned nrecv_particles =
          std::accumulate(nrank_particles.begin(), nrank_particles.end(), 0);

      for (unsigned j = 0; j < nrecv_particles; ++j) {
        // Retrieve information about the incoming message
        MPI_Status status;
        MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

        // Get buffer size
        int size;
        MPI_Get_count(&status, MPI_UINT8_T, &size);

        // Allocate the buffer now that we know how many elements there are
        std::vector<uint8_t> buffer;
        buffer.resize(size);

        // Finally receive the message
        MPI_Recv(buffer.data(), size, MPI_UINT8_T, MPI_ANY_SOURCE, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        uint8_t* bufptr = const_cast<uint8_t*>(&buffer[0]);
        int position = 0;

        // Get particle type
        int ptype;
        MPI_Unpack(bufptr, buffer.size(), &position, &ptype, 1, MPI_INT,
                   MPI_COMM_WORLD);
        std::string particle_type = mpm::ParticleTypeName.at(ptype);

        // Get materials material id
        int nmaterials = 0;
        MPI_Unpack(bufptr, buffer.size(), &position, &nmaterials, 1,
                   MPI_UNSIGNED, MPI_COMM_WORLD);
        // Vector of materials
        std::vector<std::shared_ptr<mpm::Material<Tdim>>> materials;
        materials.reserve(nmaterials);
        for (unsigned k = 0; k < nmaterials; ++k) {
          int mat_id;
          MPI_Unpack(bufptr, buffer.size(), &position, &mat_id, 1, MPI_UNSIGNED,
                     MPI_COMM_WORLD);
          materials.emplace_back(materials_.at(mat_id));
        }

        // Create particle
        auto particle =
            Factory<mpm::ParticleBase<Tdim>, mpm::Index,
                    const Eigen::Matrix<double, Tdim, 1>&>::instance()
                ->create(particle_type, static_cast<mpm::Index>(pid),
                         pcoordinates);
        particle->deserialize(buffer, materials);
        // Add particle to mesh
        this->add_particle(particle, true);
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
    std::vector<MPI_Request> send_particle_requests;
    send_particle_requests.reserve(exchange_cells.size() * 100);
    unsigned nsend_requests = 0;
    unsigned np = 0;
    std::vector<mpm::Index> remove_pids;
    // Iterate through the ghost cells and send particles
    for (auto cid : exchange_cells) {
      // Get cell pointer
      auto cell = map_cells_[cid];
      // If the previous rank of cell is the current MPI rank,
      // then send all particles
      if ((cell->rank() != cell->previous_mpirank()) &&
          (cell->previous_mpirank() == mpi_rank)) {

        // Send number of particles to receiver rank
        unsigned nparticles = cell->nparticles();
        MPI_Isend(&nparticles, 1, MPI_UNSIGNED, cell->rank(), 0, MPI_COMM_WORLD,
                  &send_requests[nsend_requests]);

        auto particle_ids = cell->particles();
        for (auto& id : particle_ids) {
          // Create a vector of serialized particle
          std::vector<uint8_t> buffer = map_particles_[id]->serialize();
          MPI_Isend(buffer.data(), buffer.size(), MPI_UINT8_T, cell->rank(), 0,
                    MPI_COMM_WORLD, &send_particle_requests[np]);
          ++np;

          // Particles to be removed from the current rank
          remove_pids.emplace_back(id);
        }
        cell->clear_particle_ids();
        ++nsend_requests;
      }
    }
    // Remove all sent particles
    this->remove_particles(remove_pids);
    // Send complete iterate only upto valid send requests
    for (unsigned i = 0; i < nsend_requests; ++i)
      MPI_Wait(&send_requests[i], MPI_STATUS_IGNORE);
    // Send particles complete
    for (unsigned i = 0; i < np; ++i)
      MPI_Wait(&send_particle_requests[i], MPI_STATUS_IGNORE);

    // Particle id
    mpm::Index pid = 0;
    // Initial particle coordinates
    Eigen::Matrix<double, Tdim, 1> pcoordinates =
        Eigen::Matrix<double, Tdim, 1>::Zero();

    // Iterate through the ghost cells and receive particles
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

        for (unsigned j = 0; j < nrecv_particles; ++j) {
          // Retrieve information about the incoming message
          MPI_Status status;
          MPI_Probe(cell->previous_mpirank(), 0, MPI_COMM_WORLD, &status);

          // Get buffer size
          int size;
          MPI_Get_count(&status, MPI_UINT8_T, &size);

          // Allocate the buffer now that we know how many elements there are
          std::vector<uint8_t> buffer;
          buffer.resize(size);

          // Finally receive the message
          MPI_Recv(buffer.data(), size, MPI_UINT8_T, cell->previous_mpirank(),
                   0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

          uint8_t* bufptr = const_cast<uint8_t*>(&buffer[0]);
          int position = 0;

          // Get particle type
          int ptype;
          MPI_Unpack(bufptr, buffer.size(), &position, &ptype, 1, MPI_INT,
                     MPI_COMM_WORLD);
          std::string particle_type = mpm::ParticleTypeName.at(ptype);

          // Get materials material id
          int nmaterials = 0;
          MPI_Unpack(bufptr, buffer.size(), &position, &nmaterials, 1,
                     MPI_UNSIGNED, MPI_COMM_WORLD);
          std::vector<std::shared_ptr<mpm::Material<Tdim>>> materials;
          materials.reserve(nmaterials);
          for (unsigned k = 0; k < nmaterials; ++k) {
            int mat_id;
            MPI_Unpack(bufptr, buffer.size(), &position, &mat_id, 1,
                       MPI_UNSIGNED, MPI_COMM_WORLD);
            // Get material
            materials.emplace_back(materials_.at(mat_id));
          }

          // Create particle
          auto particle =
              Factory<mpm::ParticleBase<Tdim>, mpm::Index,
                      const Eigen::Matrix<double, Tdim, 1>&>::instance()
                  ->create(particle_type, static_cast<mpm::Index>(pid),
                           pcoordinates);
          particle->deserialize(buffer, materials);
          // Add particle to mesh
          this->add_particle(particle, true);
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

//! Return particle scalar data
template <unsigned Tdim>
std::vector<double> mpm::Mesh<Tdim>::particles_scalar_data(
    const std::string& attribute) const {
  std::vector<double> scalar_data;
  scalar_data.reserve(particles_.size());
  // Iterate over particles and add scalar value to data
  for (auto pitr = particles_.cbegin(); pitr != particles_.cend(); ++pitr)
    scalar_data.emplace_back((*pitr)->scalar_data(attribute));
  return scalar_data;
}

//! Return particle vector data
template <unsigned Tdim>
std::vector<Eigen::Matrix<double, 3, 1>> mpm::Mesh<Tdim>::particles_vector_data(
    const std::string& attribute) const {
  std::vector<Eigen::Matrix<double, 3, 1>> vector_data;
  // Iterate over particles
  for (auto pitr = particles_.cbegin(); pitr != particles_.cend(); ++pitr) {
    Eigen::Matrix<double, 3, 1> data;
    data.setZero();
    auto pdata = (*pitr)->vector_data(attribute);
    // Fill vector_data to the size of dimensions
    for (unsigned i = 0; i < pdata.size(); ++i) data(i) = pdata(i);

    // Add to a tensor of data
    vector_data.emplace_back(data);
  }
  return vector_data;
}

//! Return particle tensor data
template <unsigned Tdim>
template <unsigned Tsize>
std::vector<Eigen::Matrix<double, Tsize, 1>>
    mpm::Mesh<Tdim>::particles_tensor_data(const std::string& attribute) const {
  std::vector<Eigen::Matrix<double, Tsize, 1>> tensor_data;
  // Iterate over particles
  for (auto pitr = particles_.cbegin(); pitr != particles_.cend(); ++pitr) {
    Eigen::Matrix<double, Tsize, 1> data;
    data.setZero();
    auto pdata = (*pitr)->tensor_data(attribute);
    // Fill tensor_data to the size of dimensions
    for (unsigned i = 0; i < pdata.size(); ++i) data(i) = pdata(i);

    // Add to a tensor of data
    tensor_data.emplace_back(data);
  }
  return tensor_data;
}

//! Return particle state variable data
template <unsigned Tdim>
std::vector<double> mpm::Mesh<Tdim>::particles_statevars_data(
    const std::string& attribute, unsigned phase) {
  std::vector<double> statevars_data;
  statevars_data.reserve(particles_.size());
  // Iterate over particles and add scalar value to data
  for (auto pitr = particles_.cbegin(); pitr != particles_.cend(); ++pitr)
    statevars_data.emplace_back((*pitr)->state_variable(attribute, phase));
  return statevars_data;
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

// Initialise the nodal properties' map
template <unsigned Tdim>
void mpm::Mesh<Tdim>::initialise_nodal_properties() {
  // Call initialise_properties function from the nodal properties
  nodal_properties_->initialise_nodal_properties();
}

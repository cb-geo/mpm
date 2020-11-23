//! Compute free surface cells, nodes, and particles
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::compute_free_surface(const std::string& method,
                                           double volume_tolerance) {
  if (method == "density") {
    return this->compute_free_surface_by_density(volume_tolerance);
  } else if (method == "geometry") {
    return this->compute_free_surface_by_geometry(volume_tolerance);
  } else {
    console_->info(
        "The selected free-surface computation method: {}\n is not "
        "available. "
        "Using density approach as default method.",
        method);
    return this->compute_free_surface_by_density(volume_tolerance);
  }
}

//! Compute free surface cells, nodes, and particles by density and geometry
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::compute_free_surface_by_geometry(
    double volume_tolerance) {
  bool status = true;

  // Reset free surface cell and particles
  this->iterate_over_cells(std::bind(&mpm::Cell<Tdim>::assign_free_surface,
                                     std::placeholders::_1, false));

  VectorDim temp_normal;
  temp_normal.setZero();
  this->iterate_over_particles_predicate(
      std::bind(&mpm::ParticleBase<Tdim>::assign_normal, std::placeholders::_1,
                temp_normal),
      std::bind(&mpm::ParticleBase<Tdim>::free_surface, std::placeholders::_1));

  this->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::assign_free_surface,
                std::placeholders::_1, false));

  // Reset volume fraction
  this->iterate_over_cells(std::bind(&mpm::Cell<Tdim>::assign_volume_fraction,
                                     std::placeholders::_1, 0.0));

  // Compute and assign volume fraction to each cell
  this->compute_cell_vol_fraction();

  // First, we detect the cell with possible free surfaces
  // Compute boundary cells and nodes based on geometry
  std::set<mpm::Index> free_surface_candidate_cells;
  for (auto citr = this->cells_.cbegin(); citr != this->cells_.cend(); ++citr) {
    // Cell contains particles
    if ((*citr)->status()) {
      bool candidate_cell = false;
      const auto& node_id = (*citr)->nodes_id();
      if ((*citr)->volume_fraction() < volume_tolerance) {
        candidate_cell = true;
        for (const auto id : node_id) {
          map_nodes_[id]->assign_free_surface(true);
        }
      } else {
        // Loop over neighbouring cells
        for (const auto neighbour_cell_id : (*citr)->neighbours()) {
          if (!map_cells_[neighbour_cell_id]->status()) {
            candidate_cell = true;
            const auto& n_node_id = map_cells_[neighbour_cell_id]->nodes_id();

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
      }
    }
  }

  // Compute particle neighbours for particles at candidate cells
  std::vector<mpm::Index> free_surface_candidate_particles_first;
  for (const auto cell_id : free_surface_candidate_cells) {
    this->find_particle_neighbours(map_cells_[cell_id]);
    const auto& particle_ids = map_cells_[cell_id]->particles();
    free_surface_candidate_particles_first.insert(
        free_surface_candidate_particles_first.end(), particle_ids.begin(),
        particle_ids.end());
  }

  // Compute boundary particles based on density function
  // Lump cell volume to nodes
  this->iterate_over_cells(std::bind(&mpm::Cell<Tdim>::map_cell_volume_to_nodes,
                                     std::placeholders::_1, 0));

  // Compute nodal value of mass density
  this->iterate_over_nodes_predicate(
      std::bind(&mpm::NodeBase<Tdim>::compute_density, std::placeholders::_1),
      std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

  std::set<mpm::Index> free_surface_candidate_particles_second;
  for (const auto p_id : free_surface_candidate_particles_first) {
    const auto& particle = map_particles_[p_id];
    status = particle->compute_free_surface_by_density();
    if (status) free_surface_candidate_particles_second.insert(p_id);
  }

  // Find free surface particles through geometry
  for (const auto p_id : free_surface_candidate_particles_second) {
    // Initialize renormalization matrix
    Eigen::Matrix<double, Tdim, Tdim> renormalization_matrix_inv;
    renormalization_matrix_inv.setZero();

    // Loop over neighbours
    const auto& particle = map_particles_[p_id];
    const auto& p_coord = particle->coordinates();
    const auto& neighbour_particles = particle->neighbours();
    const double smoothing_length = 1.33 * particle->diameter();
    for (const auto neighbour_particle_id : neighbour_particles) {
      const auto& n_coord =
          map_particles_[neighbour_particle_id]->coordinates();
      const VectorDim rel_coord = n_coord - p_coord;

      // Compute kernel gradient
      const VectorDim kernel_gradient =
          mpm::RadialBasisFunction::gradient<Tdim>(smoothing_length, -rel_coord,
                                                   "gaussian");

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
      for (const auto neighbour_particle_id : neighbour_particles) {
        const auto& n_coord =
            map_particles_[neighbour_particle_id]->coordinates();
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
      for (const auto neighbour_particle_id : neighbour_particles) {
        const auto& n_coord =
            map_particles_[neighbour_particle_id]->coordinates();
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
    }
  }

  return status;
}

//! Compute free surface cells, nodes, and particles by density
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::compute_free_surface_by_density(double volume_tolerance) {
  bool status = true;

  // Get number of MPI ranks
  int mpi_rank = 0;
  int mpi_size = 1;

#ifdef USE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  // Initialize pointer of booleans for send and receive
  bool* send_cell_solving_status = new bool[ncells()];
  memset(send_cell_solving_status, 0, ncells() * sizeof(bool));
  bool* receive_cell_solving_status = new bool[ncells()];

  for (auto citr = cells_.cbegin(); citr != cells_.cend(); ++citr)
    if ((*citr)->status())
      // Assign solving status for MPI solver
      send_cell_solving_status[(*citr)->id()] = true;
    else
      send_cell_solving_status[(*citr)->id()] = false;

  MPI_Allreduce(send_cell_solving_status, receive_cell_solving_status, ncells(),
                MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);

  for (auto citr = cells_.cbegin(); citr != cells_.cend(); ++citr) {
    // Assign solving status for MPI solver
    (*citr)->assign_solving_status(receive_cell_solving_status[(*citr)->id()]);
  }

  delete[] send_cell_solving_status;
  delete[] receive_cell_solving_status;
#endif

  // Reset free surface cell
  this->iterate_over_cells(std::bind(&mpm::Cell<Tdim>::assign_free_surface,
                                     std::placeholders::_1, false));

  // Reset volume fraction
  this->iterate_over_cells(std::bind(&mpm::Cell<Tdim>::assign_volume_fraction,
                                     std::placeholders::_1, 0.0));

  // Reset free surface particle
  this->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::assign_free_surface,
                std::placeholders::_1, false));

  // Compute and assign volume fraction to each cell
  this->compute_cell_vol_fraction();

#ifdef USE_MPI
  // Initialize vector of double for send and receive
  std::vector<double> send_cell_vol_fraction;
  send_cell_vol_fraction.resize(ncells());
  std::vector<double> receive_cell_vol_fraction;
  receive_cell_vol_fraction.resize(ncells());

  for (auto citr = cells_.cbegin(); citr != cells_.cend(); ++citr)
    if ((*citr)->status())
      // Assign volume_fraction for MPI solver
      send_cell_vol_fraction[(*citr)->id()] = (*citr)->volume_fraction();
    else
      send_cell_vol_fraction[(*citr)->id()] = 0.0;

  MPI_Allreduce(send_cell_vol_fraction.data(), receive_cell_vol_fraction.data(),
                ncells(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for (auto citr = cells_.cbegin(); citr != cells_.cend(); ++citr) {
    // Assign volume_fraction for MPI solver
    (*citr)->assign_volume_fraction(receive_cell_vol_fraction[(*citr)->id()]);
  }
#endif

  // Compute boundary cells and nodes based on geometry
  for (auto citr = this->cells_.cbegin(); citr != this->cells_.cend(); ++citr) {

    if ((*citr)->status()) {
      bool cell_at_interface = false;
      const auto& node_id = (*citr)->nodes_id();
      bool internal = true;

      //! Check internal cell
      for (const auto neighbour_cell_id : (*citr)->neighbours()) {
#if USE_MPI
        if (!map_cells_[neighbour_cell_id]->solving_status()) {
          internal = false;
          break;
        }
#else
        if (!map_cells_[neighbour_cell_id]->status()) {
          internal = false;
          break;
        }
#endif
      }

      //! Check volume fraction only for boundary cell
      if (!internal) {
        if ((*citr)->volume_fraction() < volume_tolerance) {
          cell_at_interface = true;
          for (const auto id : node_id) {
            map_nodes_[id]->assign_free_surface(cell_at_interface);
          }
        } else {
          for (const auto neighbour_cell_id : (*citr)->neighbours()) {
            if (map_cells_[neighbour_cell_id]->volume_fraction() <
                volume_tolerance) {
              cell_at_interface = true;
              const auto& n_node_id = map_cells_[neighbour_cell_id]->nodes_id();

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
        if (cell_at_interface) {
          (*citr)->assign_free_surface(cell_at_interface);
        }
      }
    }
  }

  // Compute boundary particles based on density function
  // Lump cell volume to nodes
  this->iterate_over_cells(std::bind(&mpm::Cell<Tdim>::map_cell_volume_to_nodes,
                                     std::placeholders::_1,
                                     mpm::ParticlePhase::SinglePhase));

#ifdef USE_MPI
  // Run if there is more than a single MPI task
  if (mpi_size > 1) {
    // MPI all reduce nodal volume
    this->template nodal_halo_exchange<double, 1>(
        std::bind(&mpm::NodeBase<Tdim>::volume, std::placeholders::_1,
                  mpm::ParticlePhase::SinglePhase),
        std::bind(&mpm::NodeBase<Tdim>::update_volume, std::placeholders::_1,
                  false, mpm::ParticlePhase::SinglePhase,
                  std::placeholders::_2));
  }
#endif

  // Compute nodal value of mass density
  this->iterate_over_nodes_predicate(
      std::bind(&mpm::NodeBase<Tdim>::compute_density, std::placeholders::_1),
      std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

  // Evaluate free surface particles
  this->iterate_over_particles(
      [](std::shared_ptr<mpm::ParticleBase<Tdim>> ptr) {
        bool status = ptr->compute_free_surface_by_density();
        if (status) {
          return ptr->assign_free_surface(status);
        }
      });

  return status;
}

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

//! Compute cell volume fraction
template <unsigned Tdim>
void mpm::Mesh<Tdim>::compute_cell_vol_fraction() {
  this->iterate_over_cells([&map_particles = map_particles_](
                               std::shared_ptr<mpm::Cell<Tdim>> c_ptr) {
    if (c_ptr->status()) {
      // Compute volume fraction
      double cell_volume_fraction = 0.0;
      for (const auto p_id : c_ptr->particles())
        cell_volume_fraction += map_particles[p_id]->volume();

      cell_volume_fraction = cell_volume_fraction / c_ptr->volume();
      return c_ptr->assign_volume_fraction(cell_volume_fraction);
    }
  });
}

//! Assign particle pore pressures
template <unsigned Tdim>
bool mpm::Mesh<Tdim>::assign_particles_pore_pressures(
    const std::vector<std::tuple<mpm::Index, double>>&
        particle_pore_pressures) {
  bool status = true;

  try {
    if (!particles_.size())
      throw std::runtime_error(
          "No particles have been assigned in mesh, cannot assign pore "
          "pressures");

    for (const auto& particle_pore_pressure : particle_pore_pressures) {
      // Particle id
      mpm::Index pid = std::get<0>(particle_pore_pressure);
      // Pore pressure
      double pore_pressure = std::get<1>(particle_pore_pressure);

      if (map_particles_.find(pid) != map_particles_.end())
        map_particles_[pid]->assign_pressure(pore_pressure,
                                             mpm::ParticlePhase::Liquid);
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Create a list of active nodes in mesh and assign active node id (rank-wise)
template <unsigned Tdim>
unsigned mpm::Mesh<Tdim>::assign_active_nodes_id() {
  // Clear existing list of active nodes
  this->active_nodes_.clear();
  Index active_id = 0;

#ifdef USE_MPI
  // Initialize pointer of booleans for send and receive
  bool* send_nodal_solving_status = new bool[nnodes()];
  memset(send_nodal_solving_status, 0, nnodes() * sizeof(bool));
  bool* receive_nodal_solving_status = new bool[nnodes()];

  for (auto nitr = nodes_.cbegin(); nitr != nodes_.cend(); ++nitr) {
    if ((*nitr)->status()) {
      this->active_nodes_.add(*nitr);
      (*nitr)->assign_active_id(active_id);
      active_id++;

      // Assign solving status for MPI solver
      send_nodal_solving_status[(*nitr)->id()] = true;

    } else {
      (*nitr)->assign_active_id(std::numeric_limits<Index>::max());
    }
  }

  MPI_Allreduce(send_nodal_solving_status, receive_nodal_solving_status,
                nnodes(), MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);

  for (auto nitr = nodes_.cbegin(); nitr != nodes_.cend(); ++nitr) {
    if (receive_nodal_solving_status[(*nitr)->id()]) {
      // Assign solving status for MPI solver
      (*nitr)->assign_solving_status(true);
    }
  }

  delete[] send_nodal_solving_status;
  delete[] receive_nodal_solving_status;

#else
  for (auto nitr = nodes_.cbegin(); nitr != nodes_.cend(); ++nitr) {
    if ((*nitr)->status()) {
      this->active_nodes_.add(*nitr);
      (*nitr)->assign_active_id(active_id);
      active_id++;

    } else {
      (*nitr)->assign_active_id(std::numeric_limits<Index>::max());
    }
  }
#endif

  return active_id;
}

//! Assign active node id (globally in All MPI ranks)
template <unsigned Tdim>
unsigned mpm::Mesh<Tdim>::assign_global_active_nodes_id() {
  // Clear existing list of active nodes
  Index global_active_id = 0;

#ifdef USE_MPI
  for (auto nitr = nodes_.cbegin(); nitr != nodes_.cend(); ++nitr) {
    if ((*nitr)->solving_status()) {
      (*nitr)->assign_global_active_id(global_active_id);
      global_active_id++;
    } else {
      (*nitr)->assign_global_active_id(std::numeric_limits<Index>::max());
    }
  }
#endif

  return global_active_id;
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
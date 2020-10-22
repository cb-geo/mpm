//! Construct a semi-implicit eigen matrix assembler
template <unsigned Tdim>
mpm::AssemblerEigenSemiImplicitTwoPhase<
    Tdim>::AssemblerEigenSemiImplicitTwoPhase()
    : mpm::AssemblerBase<Tdim>() {
  //! Logger
  console_ = spdlog::stdout_color_mt("AssemblerEigenSemiImplicitTwoPhase");
}

//! Assign global node indices
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhase<Tdim>::assign_global_node_indices(
    unsigned nactive_node, unsigned nglobal_active_node) {
  bool status = true;
  try {
    // Total number of active node (in a rank) and (rank) node indices
    active_dof_ = nactive_node;
    global_node_indices_ = mesh_->global_node_indices();

#ifdef USE_MPI
    // Total number of active node (in all rank)
    global_active_dof_ = nglobal_active_node;

    // Initialise mapping vector
    rank_global_mapper_.resize(active_dof_);

    // Nodes container
    const auto& nodes = mesh_->active_nodes();
    for (int counter = 0; counter < nodes.size(); counter++) {
      // Assign get nodal global index
      rank_global_mapper_[counter] = nodes[counter]->global_active_id();
    }
#endif

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assemble stiffness matrix (semi-implicit)
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhase<Tdim>::assemble_stiffness_matrix(
    unsigned dir, double dt) {
  bool status = true;
  try {
    // Initialise stiffness_matrix
    Eigen::SparseMatrix<double> stiffness_matrix;
    // Resize stiffness matrix
    stiffness_matrix.resize(2 * active_dof_, 2 * active_dof_);
    stiffness_matrix.setZero();

    // Reserve storage for sparse matrix
    switch (Tdim) {
      // For 2d: 10 entries /column
      case (2): {
        stiffness_matrix.reserve(
            Eigen::VectorXi::Constant(2 * active_dof_, 10));
        break;
      }
      // For 3d: 30 entries /column
      case (3): {
        stiffness_matrix.reserve(
            Eigen::VectorXi::Constant(2 * active_dof_, 30));
        break;
      }
    }

    // Cell pointer
    const auto& cells = mesh_->cells();
    // Iterate over cells for drag force coefficient
    mpm::Index cid = 0;
    for (auto cell_itr = cells.cbegin(); cell_itr != cells.cend(); ++cell_itr) {
      if ((*cell_itr)->status()) {
        // Node ids in each cell
        const auto nids = global_node_indices_.at(cid);
        // Local stiffness matrix
        auto k_inter_element = (*cell_itr)->K_inter_element();
        // Assemble global stiffness matrix
        for (unsigned i = 0; i < nids.size(); ++i) {
          for (unsigned j = 0; j < nids.size(); ++j) {
            stiffness_matrix.coeffRef(nids(i) + active_dof_, nids(j)) +=
                -k_inter_element(i, j) * dt;
            stiffness_matrix.coeffRef(nids(i) + active_dof_,
                                      nids(j) + active_dof_) +=
                k_inter_element(i, j) * dt;
          }
        }
        ++cid;
      }
    }

    // Active nodes pointer
    const auto& nodes = mesh_->active_nodes();
    // Iterate over cells for mass coefficient
    for (auto node_itr = nodes.cbegin(); node_itr != nodes.cend(); ++node_itr) {
      // Id for active node
      auto active_id = (*node_itr)->active_id();
      // Assemble global stiffness matrix for solid mass
      stiffness_matrix.coeffRef(active_id, active_id) +=
          (*node_itr)->mass(mpm::NodePhase::nSolid);
      // Assemble global stiffness matrix for liquid mass
      stiffness_matrix.coeffRef(active_id + active_dof_,
                                active_id + active_dof_) +=
          (*node_itr)->mass(mpm::NodePhase::nLiquid);
      stiffness_matrix.coeffRef(active_id, active_id + active_dof_) +=
          (*node_itr)->mass(mpm::NodePhase::nLiquid);
    }

    // Add stiffness matrix to map
    stiffness_matrix_.insert(
        std::make_pair<unsigned, Eigen::SparseMatrix<double>>(
            static_cast<unsigned>(dir),
            static_cast<Eigen::SparseMatrix<double>>(stiffness_matrix)));

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assemble force vector (semi-implicit)
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhase<Tdim>::assemble_force_vector(
    double dt) {
  bool status = true;
  try {
    // Resize force vector
    force_inter_.resize(active_dof_ * 2, Tdim);
    force_inter_.setZero();
    // Resize intermediate velocity vector
    acceleration_inter_.resize(active_dof_ * 2, Tdim);
    acceleration_inter_.setZero();

    // Initialise relative velocity matrix
    Eigen::MatrixXd relative_velocity;
    // Initialise drag force matrix
    Eigen::MatrixXd drag_force;

    // Cell pointer
    const auto& cells = mesh_->cells();
    for (auto cell_itr = cells.cbegin(); cell_itr != cells.cend(); ++cell_itr) {
      if ((*cell_itr)->status()) {
        // Nodes pointer of the cell
        auto nodes = (*cell_itr)->nodes();
        // Reset relative velocity matrix
        relative_velocity.resize(nodes.size(), Tdim);
        relative_velocity.setZero();
        // Compute relative velocity matrix
        for (unsigned i = 0; i < nodes.size(); ++i)
          relative_velocity.row(i) =
              (nodes.at(i)->velocity(mpm::NodePhase::nLiquid) -
               nodes.at(i)->velocity(mpm::NodePhase::nSolid))
                  .transpose();

        // Reset drag force matrix
        drag_force.resize(nodes.size(), Tdim);
        drag_force.setZero();
        // Local stiffness matrix
        drag_force = (*cell_itr)->K_inter_element() * relative_velocity;
        // Update nodal drag force
        for (unsigned i = 0; i < nodes.size(); ++i) {
          nodes.at(i)->update_drag_force(drag_force.row(i).transpose());
        }
      }
    }

    // Active nodes pointer
    const auto& nodes = mesh_->active_nodes();
    // Iterate over nodes
    mpm::Index nid = 0;
    for (auto node_itr = nodes.cbegin(); node_itr != nodes.cend(); ++node_itr) {
      // Compute nodal intermediate force
      (*node_itr)->compute_intermediate_force(dt);
      // Assemble intermediate force vector
      force_inter_.row(nid) = (*node_itr)->force_total_inter().transpose();
      force_inter_.row(nid + active_dof_) =
          (*node_itr)->force_fluid_inter().transpose();
      ++nid;
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Apply velocity constraints to force vector
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhase<
    Tdim>::apply_velocity_constraints() {
  bool status = false;
  try {
    // Modify the force vector(b = b - A * bc)
    for (unsigned i = 0; i < Tdim; i++) {
      force_inter_.col(i) -=
          stiffness_matrix_.at(i) * velocity_constraints_.col(i);

      // Iterate over velocity constraints (non-zero elements)
      for (unsigned j = 0; j < velocity_constraints_.outerSize(); ++j) {
        for (Eigen::SparseMatrix<double>::InnerIterator itr(
                 velocity_constraints_, j);
             itr; ++itr) {
          // Check direction
          if (itr.col() == i) {
            // Assign 0 to specified column
            stiffness_matrix_.at(i).col(itr.row()) *= 0;
            // Assign 0 to specified row
            stiffness_matrix_.at(i).row(itr.row()) *= 0;
            // Assign 1 to diagnal element
            stiffness_matrix_.at(i).coeffRef(itr.row(), itr.row()) = 1;
            // Assign 0 to specified component
            force_inter_(itr.row(), itr.col()) = 0;
          }
        }
      }
    }

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return status;
}

//! Assemble Laplacian matrix
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhase<Tdim>::assemble_laplacian_matrix(
    double dt) {
  bool status = true;
  try {
    // Initialise Laplacian matrix
    laplacian_matrix_.resize(active_dof_, active_dof_);
    laplacian_matrix_.setZero();

    // Reserve storage for sparse matrix
    switch (Tdim) {
      // For 2d: 10 entries /column
      case (2): {
        laplacian_matrix_.reserve(Eigen::VectorXi::Constant(active_dof_, 10));
        break;
      }
      // For 3d: 30 entries /column
      case (3): {
        laplacian_matrix_.reserve(Eigen::VectorXi::Constant(active_dof_, 30));
        break;
      }
    }

    // Cell pointer
    const auto& cells = mesh_->cells();

    // Iterate over cells
    mpm::Index cid = 0;
    for (auto cell_itr = cells.cbegin(); cell_itr != cells.cend(); ++cell_itr) {
      if ((*cell_itr)->status()) {
        // Node ids in each cell
        const auto nids = global_node_indices_.at(cid);

        // Laplacian element of cell
        const auto cell_laplacian = (*cell_itr)->laplacian_matrix();

        // Assemble global laplacian matrix
        for (unsigned i = 0; i < nids.size(); ++i) {
          for (unsigned j = 0; j < nids.size(); ++j) {
            laplacian_matrix_.coeffRef(global_node_indices_.at(cid)(i),
                                       global_node_indices_.at(cid)(j)) +=
                cell_laplacian(i, j);
          }
        }

        ++cid;
      }
    }

    laplacian_matrix_ *= dt;

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Assemble poisson right vector
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhase<Tdim>::assemble_poisson_right(
    double dt) {
  bool status = true;
  try {
    // Initialise Poisson RHS matrix
    Eigen::SparseMatrix<double> poisson_right_matrix, poisson_right_matrix_m;
    // Resize poisson right matrix for solid
    poisson_right_matrix.resize(active_dof_, active_dof_ * Tdim);
    poisson_right_matrix.setZero();
    // Resize poisson right matrix for mixture
    poisson_right_matrix_m.resize(active_dof_, active_dof_ * Tdim);
    poisson_right_matrix_m.setZero();

    // Reserve storage for sparse matrix
    switch (Tdim) {
      // For 2d: 10 entries /column
      case (2): {
        poisson_right_matrix.reserve(
            Eigen::VectorXi::Constant(active_dof_ * Tdim, 10));
        poisson_right_matrix_m.reserve(
            Eigen::VectorXi::Constant(active_dof_ * Tdim, 10));
        break;
      }
      // For 3d: 30 entries /column
      case (3): {
        poisson_right_matrix.reserve(
            Eigen::VectorXi::Constant(active_dof_ * Tdim, 30));
        poisson_right_matrix_m.reserve(
            Eigen::VectorXi::Constant(active_dof_ * Tdim, 30));
        break;
      }
    }

    // Cell pointer
    const auto& cells = mesh_->cells();

    // Iterate over cells
    mpm::Index cid = 0;
    for (auto cell_itr = cells.cbegin(); cell_itr != cells.cend(); ++cell_itr) {
      if ((*cell_itr)->status()) {
        // Node ids in each cell
        const auto nids = global_node_indices_.at(cid);

        // Local Poisson RHS matrix for solid
        auto cell_poisson_right = (*cell_itr)->poisson_right_matrix();
        // Local Poisson RHS matrix for mixture
        auto cell_poisson_right_m = (*cell_itr)->poisson_right_matrix_m();

        // Assemble global poisson RHS matrix
        for (unsigned i = 0; i < nids.size(); ++i) {
          for (unsigned j = 0; j < nids.size(); ++j) {
            for (unsigned k = 0; k < Tdim; ++k) {
              poisson_right_matrix.coeffRef(
                  global_node_indices_.at(cid)(i),
                  global_node_indices_.at(cid)(j) + k * active_dof_) +=
                  cell_poisson_right(i, j + k * nids.size());
              poisson_right_matrix_m.coeffRef(
                  global_node_indices_.at(cid)(i),
                  global_node_indices_.at(cid)(j) + k * active_dof_) +=
                  cell_poisson_right_m(i, j + k * nids.size());
            }
          }
        }
        cid++;
      }
    }

    // Resize poisson right vector
    poisson_rhs_vector_.resize(active_dof_);
    poisson_rhs_vector_.setZero();

    // Compute velocity
    Eigen::MatrixXd solid_velocity, relative_velocity;
    solid_velocity.resize(active_dof_, Tdim);
    solid_velocity.setZero();
    relative_velocity.resize(active_dof_, Tdim);
    relative_velocity.setZero();

    // Active nodes
    const auto& active_nodes = mesh_->active_nodes();
    unsigned node_index = 0;

    for (auto node_itr = active_nodes.cbegin(); node_itr != active_nodes.cend();
         ++node_itr) {
      // Compute nodal intermediate force
      solid_velocity.row(node_index) =
          (*node_itr)
              ->intermediate_velocity(mpm::ParticlePhase::Solid)
              .transpose();
      relative_velocity.row(node_index) =
          ((*node_itr)->intermediate_velocity(mpm::ParticlePhase::Liquid) -
           (*node_itr)->intermediate_velocity(mpm::ParticlePhase::Solid))
              .transpose();
      node_index++;
    }

    // Resize velocity vectors
    solid_velocity.resize(active_dof_ * Tdim, 1);
    relative_velocity.resize(active_dof_ * Tdim, 1);

    // Compute poisson RHS vector
    poisson_rhs_vector_ = -poisson_right_matrix * solid_velocity +
                          poisson_right_matrix_m * relative_velocity;

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Apply pressure constraints vector
template <unsigned Tdim>
void mpm::AssemblerEigenSemiImplicitTwoPhase<
    Tdim>::apply_pressure_constraints() {
  try {
    // Modify poisson_rhs_vector_
    poisson_rhs_vector_ -= laplacian_matrix_ * pressure_constraints_;

    // Apply free surface
    if (!free_surface_.empty()) {
      const auto& nodes = mesh_->nodes();
      for (const auto& free_node : free_surface_) {
        const auto column_index = nodes[free_node]->active_id();
        // Modify poisson_rhs_vector
        poisson_rhs_vector_(column_index) = 0;
        // Modify laplacian_matrix
        laplacian_matrix_.row(column_index) *= 0;
        laplacian_matrix_.col(column_index) *= 0;
        laplacian_matrix_.coeffRef(column_index, column_index) = 1;
      }
      // Clear the vector
      free_surface_.clear();
    }

    // Apply pressure constraints
    for (Eigen::SparseVector<double>::InnerIterator it(pressure_constraints_);
         it; ++it) {
      // Modify poisson_rhs_vector
      poisson_rhs_vector_(it.index()) = it.value();
      // Modify laplacian_matrix
      laplacian_matrix_.row(it.index()) *= 0;
      laplacian_matrix_.col(it.index()) *= 0;
      laplacian_matrix_.coeffRef(it.index(), it.index()) = 1;
    }

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
}
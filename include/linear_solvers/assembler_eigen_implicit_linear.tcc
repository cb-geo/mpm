//! Construct a implicit eigen matrix assembler
template <unsigned Tdim>
mpm::AssemblerEigenImplicitLinear<Tdim>::AssemblerEigenImplicitLinear(
    unsigned node_neighbourhood)
    : mpm::AssemblerBase<Tdim>(node_neighbourhood) {
  //! Logger
  console_ = spdlog::stdout_color_mt("AssemblerEigenImplicitLinear");
}

//! Assign global node indices
template <unsigned Tdim>
bool mpm::AssemblerEigenImplicitLinear<Tdim>::assign_global_node_indices(
    unsigned nactive_node, unsigned nglobal_active_node) {
  bool status = true;
  try {
    // Total number of active node (in a rank) and (rank) node indices
    active_dof_ = nactive_node;
    global_node_indices_ = mesh_->global_node_indices();

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assemble stiffness matrix
template <unsigned Tdim>
bool mpm::AssemblerEigenImplicitLinear<Tdim>::assemble_stiffness_matrix() {
  bool status = true;
  try {
    // Initialise stiffness matrix
    stiffness_matrix_.resize(active_dof_ * Tdim, Tdim * active_dof_ * Tdim);
    stiffness_matrix_.setZero();

    // Reserve storage for sparse matrix
    stiffness_matrix_.reserve(
        Eigen::VectorXi::Constant(active_dof_ * Tdim, sparse_row_size_));

    // Cell pointer
    const auto& cells = mesh_->cells();

    // Active nodes
    const auto& active_nodes = mesh_->active_nodes();
    const unsigned nactive_node = active_nodes.size();

    // Iterate over cells
    mpm::Index cid = 0;
    for (auto cell_itr = cells.cbegin(); cell_itr != cells.cend(); ++cell_itr) {
      if ((*cell_itr)->status()) {
        // Node ids in each cell
        const auto nids = global_node_indices_.at(cid);

        // Element stiffness of cell
        const auto cell_stiffness = (*cell_itr)->stiffness_matrix();

        // Assemble global laplacian matrix
        for (unsigned i = 0; i < nids.size(); ++i) {
          for (unsigned j = 0; j < nids.size(); ++j) {
            for (unsigned k = 0; k < Tdim; ++k) {
              for (unsigned l = 0; l < Tdim; ++l) {
                stiffness_matrix_.coeffRef(
                    nactive_node * k + global_node_indices_.at(cid)(i),
                    global_node_indices_.at(cid)(j)) +=
                    cell_stiffness(nids.size() * i + k, nids.size() * j + l);
              }
            }
          }
        }
        ++cid;
      }
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Assemble residual force right vector
template <unsigned Tdim>
bool mpm::AssemblerEigenImplicitLinear<Tdim>::assemble_residual_force_right() {
  bool status = true;
  try {
    // Initialise residual force RHS vector
    residual_force_rhs_vector_.resize(active_dof_ * Tdim);
    residual_force_rhs_vector_.setZero();

    // Active nodes
    const auto& active_nodes = mesh_->active_nodes();
    const unsigned nactive_node = active_nodes.size();
    const unsigned solid = mpm::ParticlePhase::SinglePhase;
    unsigned node_index = 0;

    // Iterate over nodes
    for (auto node_itr = active_nodes.cbegin(); node_itr != active_nodes.cend();
         ++node_itr) {

      for (unsigned i = 0; i < Tdim; ++i) {
        // Nodal residual force
        residual_force_rhs_vector_(nactive_node * i + node_index) =
            (*node_itr)->external_force(solid)(i) -
            (*node_itr)->internal_force(solid)(i);
      }
      node_index++;
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! TODO: Assign displacement constraints
template <unsigned Tdim>
bool mpm::AssemblerEigenImplicitLinear<Tdim>::assign_displacement_constraints(
    double current_time) {
  bool status = false;
  try {
    // // Resize pressure constraints vector
    // pressure_constraints_.resize(active_dof_);
    // pressure_constraints_.reserve(int(0.5 * active_dof_));

    // // Nodes container
    // const auto& nodes = mesh_->active_nodes();
    // // Iterate over nodes to get pressure constraints
    // for (auto node = nodes.cbegin(); node != nodes.cend(); ++node) {
    //   // Assign total pressure constraint
    //   const double pressure_constraint =
    //       (*node)->pressure_constraint(0, current_time);

    //   // Check if there is a pressure constraint
    //   if (pressure_constraint != std::numeric_limits<double>::max()) {
    //     // Insert the pressure constraints
    //     pressure_constraints_.insert((*node)->active_id()) =
    //         pressure_constraint;
    //   }
    // }
    status = true;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return status;
}

//! TODO: Apply displacement constraints vector
template <unsigned Tdim>
void mpm::AssemblerEigenImplicitLinear<Tdim>::apply_displacement_constraints() {
  //   try {
  //     // Modify poisson_rhs_vector_
  //     poisson_rhs_vector_ -= laplacian_matrix_ * pressure_constraints_;

  //     // Apply free surface
  //     if (!free_surface_.empty()) {
  //       const auto& nodes = mesh_->nodes();
  //       for (const auto& free_node : free_surface_) {
  //         const auto column_index = nodes[free_node]->active_id();
  //         // Modify poisson_rhs_vector
  //         poisson_rhs_vector_(column_index) = 0;
  //         // Modify laplacian_matrix
  //         laplacian_matrix_.row(column_index) *= 0;
  //         laplacian_matrix_.col(column_index) *= 0;
  //         laplacian_matrix_.coeffRef(column_index, column_index) = 1;
  //       }
  //       // Clear the vector
  //       free_surface_.clear();
  //     }

  //     // Apply pressure constraints
  //     for (Eigen::SparseVector<double>::InnerIterator
  //     it(pressure_constraints_);
  //          it; ++it) {
  //       // Modify poisson_rhs_vector
  //       poisson_rhs_vector_(it.index()) = it.value();
  //       // Modify laplacian_matrix
  //       laplacian_matrix_.row(it.index()) *= 0;
  //       laplacian_matrix_.col(it.index()) *= 0;
  //       laplacian_matrix_.coeffRef(it.index(), it.index()) = 1;
  //     }

  //   } catch (std::exception& exception) {
  //     console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  //   }
}
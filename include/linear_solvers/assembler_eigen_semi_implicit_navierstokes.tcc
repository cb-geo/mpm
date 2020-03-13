//! Construct a semi-implicit eigen matrix assembler
template <unsigned Tdim>
mpm::AssemblerEigenSemiImplicitNavierStokes<
    Tdim>::AssemblerEigenSemiImplicitNavierStokes()
    : mpm::AssemblerBase<Tdim>() {
  //! Logger
  console_ = spdlog::stdout_color_mt("AssemblerEigenSemiImplicitNavierStokes");
}

//! Assign global node indices
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitNavierStokes<
    Tdim>::assign_global_node_indices(unsigned nactive_node) {
  bool status = true;
  try {
    // Total number of active node
    active_dof_ = nactive_node;

    global_node_indices_ = mesh_->global_node_indices();

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assemble Laplacian matrix
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitNavierStokes<
    Tdim>::assemble_laplacian_matrix(double dt) {
  bool status = true;
  // try {
  //   // Initialise Laplacian matrix
  //   laplacian_matrix_.resize(active_dof_, active_dof_);
  //   laplacian_matrix_.setZero();
  //   // TODO: Make the storage being able to adaptively resize
  //   // Reserve storage for sparse matrix
  //   switch (Tdim) {
  //     // For 2d: 10 entries /column
  //     case (2): {
  //       laplacian_matrix_.reserve(Eigen::VectorXi::Constant(active_dof_,
  //       10)); break;
  //     }
  //     // For 3d: 30 entries /column
  //     case (3): {
  //       laplacian_matrix_.reserve(Eigen::VectorXi::Constant(active_dof_,
  //       30)); break;
  //     }
  //   }
  //   // Cell pointer
  //   const auto& cells = mesh_->cells();
  //   // Iterate over cells
  //   // active cell id
  //   mpm::Index cid = 0;
  //   for (auto cell_itr = cells.cbegin(); cell_itr != cells.cend();
  //   ++cell_itr) {
  //     if ((*cell_itr)->status()) {
  //       // Node ids in each cell
  //       const auto nids = global_node_indices_.at(cid);
  //       // Laplacian element of cell
  //       const auto L_element = (*cell_itr)->L_element();
  //       // Compute Laplacian elements
  //       for (unsigned i = 0; i < nids.size(); ++i) {
  //         for (unsigned j = 0; j < nids.size(); ++j) {
  //           laplacian_matrix_.coeffRef(global_node_indices_.at(cid)(i),
  //                                      global_node_indices_.at(cid)(j)) +=
  //               L_element(i, j);
  //         }
  //       }
  //       ++cid;
  //     }
  //   }
  //   laplacian_matrix_ *= dt;
  // } catch (std::exception& exception) {
  //   console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  //   status = false;
  // }
  return status;
}

template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitNavierStokes<Tdim>::assemble_poisson_right(
    std::shared_ptr<mpm::Mesh<Tdim>>& mesh_, double dt) {
  bool status = true;
  // try {
  //   // Initialise Fs & Fm matrix
  //   Eigen::SparseMatrix<double> F_s_matrix;
  //   // Resize Fs matrix
  //   F_s_matrix.resize(active_dof_, active_dof_ * Tdim);
  //   F_s_matrix.setZero();
  //   // Reserve storage for sparse matrix
  //   switch (Tdim) {
  //     // For 2d: 10 entries /column
  //     case (2): {
  //       F_s_matrix.reserve(Eigen::VectorXi::Constant(active_dof_ * Tdim,
  //       10)); break;
  //     }
  //     // For 3d: 30 entries /column
  //     case (3): {
  //       F_s_matrix.reserve(Eigen::VectorXi::Constant(active_dof_ * Tdim,
  //       30)); break;
  //     }
  //   }
  //   // Cell pointer
  //   const auto& cells = mesh_->cells();
  //   // Iterate over cells
  //   // active cell id
  //   mpm::Index cid = 0;
  //   for (auto cell_itr = cells.cbegin(); cell_itr != cells.cend();
  //   ++cell_itr) {
  //     // Cell id
  //     if ((*cell_itr)->status()) {
  //       // Node ids in each cell
  //       const auto nids = global_node_indices_.at(cid);
  //       // Fs element of cell
  //       auto F_s_element = (*cell_itr)->F_s_element();
  //       // Compute Fs & Fm
  //       for (unsigned i = 0; i < nids.size(); ++i) {
  //         for (unsigned j = 0; j < nids.size(); ++j) {
  //           for (unsigned k = 0; k < Tdim; ++k) {
  //             F_s_matrix.coeffRef(
  //                 global_node_indices_.at(cid)(i),
  //                 global_node_indices_.at(cid)(j) + k * active_dof_) +=
  //                 F_s_element(i, j + k * nids.size());
  //           }
  //         }
  //       }
  //       cid++;
  //     }
  //   }

  //   // Resize poisson right matrix
  //   force_laplacian_matrix_.resize(active_dof_);
  //   force_laplacian_matrix_.setZero();
  //   // Compute velocity
  //   Eigen::MatrixXd fluid_velocity;
  //   fluid_velocity.resize(active_dof_, Tdim);
  //   fluid_velocity.setZero();

  //   // Active nodes
  //   const auto& active_nodes = mesh_->active_nodes();
  //   const unsigned fluid = mpm::ParticlePhase::SinglePhase;
  //   unsigned node_index = 0;

  //   for (auto node_itr = active_nodes.cbegin(); node_itr !=
  //   active_nodes.cend();
  //        ++node_itr) {
  //     // Compute nodal intermediate force
  //     fluid_velocity.row(node_index) =
  //     (*node_itr)->velocity(fluid).transpose(); node_index++;
  //   }

  //   fluid_velocity.resize(active_dof_ * Tdim, 1);

  //   force_laplacian_matrix_ = -F_s_matrix * fluid_velocity;
  // } catch (std::exception& exception) {
  //   console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  //   status = false;
  // }
  return status;
}

//! Assemble K_cor_matrix
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitNavierStokes<Tdim>::assemble_K_cor_matrix(
    std::shared_ptr<mpm::Mesh<Tdim>>& mesh_, double dt) {
  bool status = true;
  // try {
  //   K_cor_matrix_.resize(active_dof_, active_dof_ * Tdim);
  //   K_cor_matrix_.setZero();

  //   K_cor_matrix_.reserve(Eigen::VectorXi::Constant(active_dof_ * Tdim, 20));

  //   unsigned nnodes_per_cell = global_node_indices_.at(0).size();

  //   const auto& cell = mesh_->cells();

  //   unsigned cid = 0;
  //   for (auto cell_itr = cell.cbegin(); cell_itr != cell.cend(); ++cell_itr)
  //   {
  //     if ((*cell_itr)->status()) {
  //       auto k_cor_element_water = (*cell_itr)->K_cor_w_element();
  //       for (unsigned k = 0; k < Tdim; k++) {
  //         for (unsigned i = 0; i < nnodes_per_cell; i++) {
  //           for (unsigned j = 0; j < nnodes_per_cell; j++) {
  //             // Fluid
  //             K_cor_matrix_.coeffRef(
  //                 global_node_indices_.at(cid)(i),
  //                 k * active_dof_ + global_node_indices_.at(cid)(j)) +=
  //                 k_cor_element_water(i, j + k * nnodes_per_cell);
  //           }
  //         }
  //       }
  //       cid++;
  //     }
  //   }
  // } catch (std::exception& exception) {
  //   console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  //   status = false;
  // }
  return status;
}

//! Apply pressure constraints vector
template <unsigned Tdim>
void mpm::AssemblerEigenSemiImplicitNavierStokes<
    Tdim>::apply_pressure_constraints() {
  // try {
  //   // Modify force_laplacian_matrix_
  //   force_laplacian_matrix_ -= laplacian_matrix_ * pressure_constraints_;

  //   // Apply free surface
  //   if (!free_surface_.empty()) {
  //     const auto& nodes = mesh_->nodes();
  //     for (const auto& free_node : free_surface_) {
  //       const auto column_index = nodes[free_node]->active_id();
  //       // Modify force_laplacian_matrix
  //       force_laplacian_matrix_(column_index) = 0;
  //       // Modify laplacian_matrix
  //       laplacian_matrix_.row(column_index) *= 0;
  //       laplacian_matrix_.col(column_index) *= 0;
  //       laplacian_matrix_.coeffRef(column_index, column_index) = 1;
  //     }
  //     // Clear the vector
  //     free_surface_.clear();
  //   }

  //   // Iterate over pressure constraints
  //   for (Eigen::SparseVector<double>::InnerIterator
  //   it(pressure_constraints_);
  //        it; ++it) {
  //     // Modify force_laplacian_matrix
  //     force_laplacian_matrix_(it.index()) = it.value();
  //     // Modify laplacian_matrix
  //     laplacian_matrix_.row(it.index()) *= 0;
  //     laplacian_matrix_.col(it.index()) *= 0;
  //     laplacian_matrix_.coeffRef(it.index(), it.index()) = 1;
  //   }
  // } catch (std::exception& exception) {
  //   console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  // }
}

//! Assign pressure constraints
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitNavierStokes<
    Tdim>::assign_pressure_constraints(double beta, const double current_time) {
  bool status = false;
  try {
    // Resize pressure constraints vector
    pressure_constraints_.resize(active_dof_);
    pressure_constraints_.reserve(int(0.5 * active_dof_));

    // Nodes container
    const auto& nodes = mesh_->active_nodes();
    // Iterate over nodes to get pressure constraints
    for (auto node = nodes.cbegin(); node != nodes.cend(); ++node) {
      // Assign total pressure constraint
      const double pressure_constraint =
          (*node)->pressure_constraint(1, current_time);
      // Check if there is a pressure constraint
      if (pressure_constraint != std::numeric_limits<double>::max()) {
        // Insert the pressure constraints
        pressure_constraints_.insert((*node)->active_id()) =
            pressure_constraint;
      }
    }
    status = true;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return status;
}

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

template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitNavierStokes<Tdim>::assemble_poisson_right(
    double dt) {
  bool status = true;
  try {
    // Initialise Poisson RHS matrix
    Eigen::SparseMatrix<double> poisson_right_matrix;
    poisson_right_matrix.resize(active_dof_, active_dof_ * Tdim);
    poisson_right_matrix.setZero();

    // Reserve storage for sparse matrix
    switch (Tdim) {
      // For 2d: 10 entries /column
      case (2): {
        poisson_right_matrix.reserve(
            Eigen::VectorXi::Constant(active_dof_ * Tdim, 10));
        break;
      }
      // For 3d: 30 entries /column
      case (3): {
        poisson_right_matrix.reserve(
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

        // Local Poisson RHS matrix
        auto cell_poisson_right = (*cell_itr)->poisson_right_matrix();

        // Assemble global poisson RHS matrix
        for (unsigned i = 0; i < nids.size(); ++i) {
          for (unsigned j = 0; j < nids.size(); ++j) {
            for (unsigned k = 0; k < Tdim; ++k) {
              poisson_right_matrix.coeffRef(
                  global_node_indices_.at(cid)(i),
                  global_node_indices_.at(cid)(j) + k * active_dof_) +=
                  cell_poisson_right(i, j + k * nids.size());
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
    Eigen::MatrixXd fluid_velocity;
    fluid_velocity.resize(active_dof_, Tdim);
    fluid_velocity.setZero();

    // Active nodes
    const auto& active_nodes = mesh_->active_nodes();
    const unsigned fluid = mpm::ParticlePhase::SinglePhase;
    unsigned node_index = 0;

    for (auto node_itr = active_nodes.cbegin(); node_itr != active_nodes.cend();
         ++node_itr) {
      // Compute nodal intermediate force
      fluid_velocity.row(node_index) = (*node_itr)->velocity(fluid).transpose();
      node_index++;
    }

    // Resize fluid velocity
    fluid_velocity.resize(active_dof_ * Tdim, 1);

    // Compute poisson RHS vector
    poisson_rhs_vector_ = -poisson_right_matrix * fluid_velocity;

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assemble corrector right matrix
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitNavierStokes<
    Tdim>::assemble_corrector_right(double dt) {
  bool status = true;
  try {
    // Resize correction matrix
    correction_matrix_.resize(active_dof_, active_dof_ * Tdim);
    correction_matrix_.setZero();

    // Reserve storage for sparse matrix
    switch (Tdim) {
      // For 2d: 10 entries /column
      case (2): {
        correction_matrix_.reserve(
            Eigen::VectorXi::Constant(active_dof_ * Tdim, 10));
        break;
      }
      // For 3d: 30 entries /column
      case (3): {
        correction_matrix_.reserve(
            Eigen::VectorXi::Constant(active_dof_ * Tdim, 30));
        break;
      }
    }

    // Cell pointer
    unsigned nnodes_per_cell = global_node_indices_.at(0).size();
    const auto& cells = mesh_->cells();

    // Iterate over cells
    unsigned cid = 0;
    for (auto cell_itr = cells.cbegin(); cell_itr != cells.cend(); ++cell_itr) {
      if ((*cell_itr)->status()) {
        auto cell_correction_matrix = (*cell_itr)->correction_matrix();
        for (unsigned k = 0; k < Tdim; k++) {
          for (unsigned i = 0; i < nnodes_per_cell; i++) {
            for (unsigned j = 0; j < nnodes_per_cell; j++) {
              // Fluid
              correction_matrix_.coeffRef(
                  global_node_indices_.at(cid)(i),
                  k * active_dof_ + global_node_indices_.at(cid)(j)) +=
                  cell_correction_matrix(i, j + k * nnodes_per_cell);
            }
          }
        }
        cid++;
      }
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Apply pressure constraints vector
template <unsigned Tdim>
void mpm::AssemblerEigenSemiImplicitNavierStokes<
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

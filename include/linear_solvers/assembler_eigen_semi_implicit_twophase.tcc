//! Construct a semi-implicit eigen matrix assembler
template <unsigned Tdim>
mpm::AssemblerEigenSemiImplicitTwoPhase<
    Tdim>::AssemblerEigenSemiImplicitTwoPhase()
    : mpm::AssemblerEigenSemiImplicitNavierStokes<Tdim>() {
  //! Logger
  console_ = spdlog::stdout_color_mt("AssemblerEigenSemiImplicitTwoPhase");
}

//! Assemble coefficient matrix for two-phase predictor
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhase<Tdim>::assemble_predictor_left(
    double dt) {
  bool status = true;
  try {
    // Loop over three direction
    for (unsigned dir = 0; dir < Tdim; dir++) {
      // Initialise coefficient_matrix
      Eigen::SparseMatrix<double> coefficient_matrix;
      coefficient_matrix.setZero();

      // Resize coefficient matrix
      coefficient_matrix.resize(2 * active_dof_, 2 * active_dof_);

      // Reserve storage for sparse matrix
      switch (Tdim) {
        // For 2d: 10 entries /column
        case (2): {
          coefficient_matrix.reserve(
              Eigen::VectorXi::Constant(2 * active_dof_, 2 * 10));
          break;
        }
        // For 3d: 30 entries /column
        case (3): {
          coefficient_matrix.reserve(
              Eigen::VectorXi::Constant(2 * active_dof_, 2 * 30));
          break;
        }
      }

      // Cell pointer
      const auto& cells = mesh_->cells();

      // Iterate over cells for drag force coefficient
      mpm::Index cid = 0;
      for (auto cell_itr = cells.cbegin(); cell_itr != cells.cend();
           ++cell_itr) {
        if ((*cell_itr)->status()) {
          // Node ids in each cell
          const auto nids = global_node_indices_.at(cid);
          // Local drag matrix
          auto cell_drag_matrix = (*cell_itr)->drag_matrix(dir);
          // Assemble global coefficient matrix
          for (unsigned i = 0; i < nids.size(); ++i) {
            for (unsigned j = 0; j < nids.size(); ++j) {
              coefficient_matrix.coeffRef(nids(i) + active_dof_, nids(j)) +=
                  -cell_drag_matrix(i, j) * dt;
              coefficient_matrix.coeffRef(nids(i) + active_dof_,
                                          nids(j) + active_dof_) +=
                  cell_drag_matrix(i, j) * dt;
            }
          }
          ++cid;
        }
      }

      // Active nodes pointer
      const auto& nodes = mesh_->active_nodes();
      // Iterate over cells for mass coefficient
      for (auto node_itr = nodes.cbegin(); node_itr != nodes.cend();
           ++node_itr) {
        // Id for active node
        auto active_id = (*node_itr)->active_id();
        // Assemble global coefficient matrix for solid mass
        coefficient_matrix.coeffRef(active_id, active_id) +=
            (*node_itr)->mass(mpm::NodePhase::NSolid);
        // Assemble global coefficient matrix for liquid mass
        coefficient_matrix.coeffRef(active_id + active_dof_,
                                    active_id + active_dof_) +=
            (*node_itr)->mass(mpm::NodePhase::NLiquid);
        coefficient_matrix.coeffRef(active_id, active_id + active_dof_) +=
            (*node_itr)->mass(mpm::NodePhase::NLiquid);
      }

      // Add coefficient matrix to map
      if (predictor_lhs_matrix_.find(dir) != predictor_lhs_matrix_.end())
        predictor_lhs_matrix_.erase(dir);

      predictor_lhs_matrix_.insert(
          std::make_pair<unsigned, Eigen::SparseMatrix<double>>(
              static_cast<unsigned>(dir),
              static_cast<Eigen::SparseMatrix<double>>(coefficient_matrix)));
    }

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assemble predictor RHS force vector
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhase<Tdim>::assemble_predictor_right(
    double dt) {
  bool status = true;
  try {
    // Resize force vector
    predictor_rhs_vector_.resize(active_dof_ * 2, Tdim);
    predictor_rhs_vector_.setZero();

    // Resize intermediate velocity vector
    intermediate_acceleration_.resize(active_dof_ * 2, Tdim);
    intermediate_acceleration_.setZero();

    // Initialise relative velocity matrix
    Eigen::MatrixXd relative_velocity;
    // Initialise cell drag force matrix
    Eigen::MatrixXd cell_drag_force;

    // Cell pointer
    const auto& cells = mesh_->cells();

    for (auto cell_itr = cells.cbegin(); cell_itr != cells.cend(); ++cell_itr) {
      if ((*cell_itr)->status()) {
        // Nodes pointer of the cell
        const auto& nodes = (*cell_itr)->nodes();

        // Reset relative velocity matrix
        relative_velocity.resize(nodes.size(), Tdim);
        relative_velocity.setZero();
        // Compute relative velocity matrix
        for (unsigned i = 0; i < nodes.size(); ++i)
          relative_velocity.row(i) =
              (nodes.at(i)->velocity(mpm::NodePhase::NLiquid) -
               nodes.at(i)->velocity(mpm::NodePhase::NSolid))
                  .transpose();

        // Reset drag force matrix
        cell_drag_force.resize(nodes.size(), Tdim);
        cell_drag_force.setZero();

        // Local drag matrix
        for (unsigned dir = 0; dir < Tdim; dir++)
          cell_drag_force.col(dir) =
              (*cell_itr)->drag_matrix(dir) * relative_velocity.col(dir);

        // Update nodal drag force
        for (unsigned i = 0; i < nodes.size(); ++i) {
          nodes.at(i)->update_drag_force(cell_drag_force.row(i).transpose());
        }
      }
    }

    // Active nodes pointer
    const auto& nodes = mesh_->active_nodes();
    // Iterate over nodes
    mpm::Index nid = 0;
    for (auto node_itr = nodes.cbegin(); node_itr != nodes.cend(); ++node_itr) {
      // Compute nodal intermediate force
      (*node_itr)->compute_intermediate_force();

      // Assemble intermediate force vector
      predictor_rhs_vector_.row(nid) =
          (*node_itr)->force_total_inter().transpose();
      predictor_rhs_vector_.row(nid + active_dof_) =
          (*node_itr)->force_fluid_inter().transpose();
      ++nid;
    }
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
    Eigen::SparseMatrix<double> solid_poisson_right_matrix,
        liquid_poisson_right_matrix;

    // Resize poisson right matrix for solid
    solid_poisson_right_matrix.resize(active_dof_, active_dof_ * Tdim);
    solid_poisson_right_matrix.setZero();
    // Resize poisson right matrix for liquid
    liquid_poisson_right_matrix.resize(active_dof_, active_dof_ * Tdim);
    liquid_poisson_right_matrix.setZero();

    // Reserve storage for sparse matrix
    switch (Tdim) {
      // For 2d: 10 entries /column
      case (2): {
        solid_poisson_right_matrix.reserve(
            Eigen::VectorXi::Constant(active_dof_ * Tdim, 10));
        liquid_poisson_right_matrix.reserve(
            Eigen::VectorXi::Constant(active_dof_ * Tdim, 10));
        break;
      }
      // For 3d: 30 entries /column
      case (3): {
        solid_poisson_right_matrix.reserve(
            Eigen::VectorXi::Constant(active_dof_ * Tdim, 30));
        liquid_poisson_right_matrix.reserve(
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
        auto cell_poisson_right_solid =
            (*cell_itr)->poisson_right_matrix(mpm::NodePhase::NSolid);
        // Local Poisson RHS matrix for liquid
        auto cell_poisson_right_liquid =
            (*cell_itr)->poisson_right_matrix(mpm::NodePhase::NLiquid);

        // Assemble global poisson RHS matrix
        for (unsigned i = 0; i < nids.size(); ++i) {
          for (unsigned j = 0; j < nids.size(); ++j) {
            for (unsigned k = 0; k < Tdim; ++k) {
              solid_poisson_right_matrix.coeffRef(
                  global_node_indices_.at(cid)(i),
                  global_node_indices_.at(cid)(j) + k * active_dof_) +=
                  cell_poisson_right_solid(i, j + k * nids.size());
              liquid_poisson_right_matrix.coeffRef(
                  global_node_indices_.at(cid)(i),
                  global_node_indices_.at(cid)(j) + k * active_dof_) +=
                  cell_poisson_right_liquid(i, j + k * nids.size());
            }
          }
        }
        cid++;
      }
    }

    // Resize poisson right vector
    poisson_rhs_vector_.resize(active_dof_);
    poisson_rhs_vector_.setZero();

    // Compute intermediate solid and liquid velocity
    Eigen::MatrixXd solid_velocity, liquid_velocity;
    solid_velocity.resize(active_dof_, Tdim);
    solid_velocity.setZero();
    liquid_velocity.resize(active_dof_, Tdim);
    liquid_velocity.setZero();

    // Active nodes
    const auto& active_nodes = mesh_->active_nodes();
    unsigned node_index = 0;

    for (auto node_itr = active_nodes.cbegin(); node_itr != active_nodes.cend();
         ++node_itr) {
      // Compute nodal intermediate force
      solid_velocity.row(node_index) =
          (*node_itr)
              ->intermediate_velocity(mpm::NodePhase::NSolid)
              .transpose();
      liquid_velocity.row(node_index) =
          (*node_itr)
              ->intermediate_velocity(mpm::NodePhase::NLiquid)
              .transpose();
      node_index++;
    }

    // Resize velocity vectors
    solid_velocity.resize(active_dof_ * Tdim, 1);
    liquid_velocity.resize(active_dof_ * Tdim, 1);

    // Compute poisson RHS vector
    poisson_rhs_vector_ = -(solid_poisson_right_matrix * solid_velocity) -
                          (liquid_poisson_right_matrix * liquid_velocity);

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assemble corrector right matrix
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhase<Tdim>::assemble_corrector_right(
    double dt) {
  bool status = true;
  try {
    // Resize correction matrix
    correction_matrix_.resize(2 * active_dof_, active_dof_ * Tdim);
    correction_matrix_.setZero();

    // Reserve storage for sparse matrix
    switch (Tdim) {
      // For 2d: 10 entries /column x 2 phase
      case (2): {
        correction_matrix_.reserve(
            Eigen::VectorXi::Constant(active_dof_ * Tdim, 2 * 10));
        break;
      }
      // For 3d: 30 entries /column x 2 phase
      case (3): {
        correction_matrix_.reserve(
            Eigen::VectorXi::Constant(active_dof_ * Tdim, 2 * 30));
        break;
      }
    }

    // Cell pointer
    const auto& cells = mesh_->cells();

    // Iterate over cells
    unsigned cid = 0;
    for (auto cell_itr = cells.cbegin(); cell_itr != cells.cend(); ++cell_itr) {
      if ((*cell_itr)->status()) {
        // Number of nodes in cell
        unsigned nnodes_per_cell = global_node_indices_.at(cid).size();
        // Local correction matrix for solid
        auto correction_matrix_solid =
            (*cell_itr)->correction_matrix(mpm::NodePhase::NSolid);
        // Local correction matrix for liquid
        auto coefficient_matrix_liquid =
            (*cell_itr)->correction_matrix(mpm::NodePhase::NLiquid);
        for (unsigned k = 0; k < Tdim; k++) {
          for (unsigned i = 0; i < nnodes_per_cell; i++) {
            for (unsigned j = 0; j < nnodes_per_cell; j++) {
              // Solid phase
              correction_matrix_.coeffRef(
                  global_node_indices_.at(cid)(i),
                  k * active_dof_ + global_node_indices_.at(cid)(j)) +=
                  correction_matrix_solid(i, j + k * nnodes_per_cell);
              // Liquid phase
              correction_matrix_.coeffRef(
                  global_node_indices_.at(cid)(i) + active_dof_,
                  k * active_dof_ + global_node_indices_.at(cid)(j)) +=
                  coefficient_matrix_liquid(i, j + k * nnodes_per_cell);
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

//! Assign pressure constraints
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhase<Tdim>::assign_pressure_constraints(
    double beta, double current_time) {
  bool status = false;
  try {
    // Resize pressure constraints vector
    pressure_constraints_.setZero();
    pressure_constraints_.data().squeeze();
    pressure_constraints_.resize(active_dof_);
    pressure_constraints_.reserve(int(0.5 * active_dof_));

    // Nodes container
    const auto& nodes = mesh_->active_nodes();
    // Iterate over nodes to get pressure constraints
    for (auto node = nodes.cbegin(); node != nodes.cend(); ++node) {
      // Assign total pressure constraint
      const double pressure_constraint =
          (*node)->pressure_constraint(mpm::NodePhase::NLiquid, current_time);

      // Check if there is a pressure constraint
      if (pressure_constraint != std::numeric_limits<double>::max()) {
        // Insert the pressure constraints
        pressure_constraints_.insert((*node)->active_id()) =
            (1 - beta) * pressure_constraint;
      }
    }
    status = true;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return status;
}

//! Assign velocity constraints
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhase<
    Tdim>::assign_velocity_constraints() {
  bool status = false;
  try {
    // Initialise constraints matrix from triplet
    std::vector<Eigen::Triplet<double>> triplet_list;
    // Nodes container
    const auto& nodes = mesh_->active_nodes();
    // Iterate over nodes
    for (auto node = nodes.cbegin(); node != nodes.cend(); ++node) {
      // Get velocity constraints
      const auto& velocity_constraints = (*node)->velocity_constraints();
      // Assign constraints matrix
      for (const auto constraint : velocity_constraints) {
        // Insert constraint to the matrix
        triplet_list.push_back(Eigen::Triplet<double>(
            (constraint).first / Tdim * active_dof_ + (*node)->active_id(),
            (constraint).first % Tdim, (constraint).second));
      }
    }
    // Reserve the storage for the velocity constraints matrix
    velocity_constraints_.setZero();
    velocity_constraints_.data().squeeze();
    velocity_constraints_.resize(active_dof_ * 2, Tdim);
    velocity_constraints_.reserve(
        Eigen::VectorXi::Constant(Tdim, triplet_list.size() + 10));
    // Assemble the velocity constraints matrix
    velocity_constraints_.setFromTriplets(triplet_list.begin(),
                                          triplet_list.end());

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return status;
}

//! Apply velocity constraints to preditor RHS vector
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhase<
    Tdim>::apply_velocity_constraints() {
  bool status = false;
  try {
    // Modify the force vector(b = b - A * bc)
    for (unsigned dir = 0; dir < Tdim; dir++) {
      predictor_rhs_vector_.col(dir) -=
          predictor_lhs_matrix_.at(dir) * velocity_constraints_.col(dir);

      // Iterate over velocity constraints (non-zero elements)
      for (unsigned j = 0; j < velocity_constraints_.outerSize(); ++j) {
        for (Eigen::SparseMatrix<double>::InnerIterator itr(
                 velocity_constraints_, j);
             itr; ++itr) {
          // Check direction
          if (itr.col() == dir) {
            // Assign 0 to specified column
            predictor_lhs_matrix_.at(dir).col(itr.row()) *= 0;
            // Assign 0 to specified row
            predictor_lhs_matrix_.at(dir).row(itr.row()) *= 0;
            // Assign 1  to diagnal element
            predictor_lhs_matrix_.at(dir).coeffRef(itr.row(), itr.row()) = 1.0;

            predictor_rhs_vector_(itr.row(), itr.col()) = 0.;
          }
        }
      }
    }

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return status;
}
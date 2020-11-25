//! Map cell volume to nodes
template <unsigned Tdim>
void mpm::Cell<Tdim>::map_cell_volume_to_nodes(unsigned phase) {
  if (this->status()) {
    // Check if cell volume is set
    if (volume_ <= std::numeric_limits<double>::lowest())
      this->compute_volume();

    for (unsigned i = 0; i < nodes_.size(); ++i) {
      nodes_[i]->update_volume(true, phase, volume_ / nnodes_);
    }
  }
}

//! Return local node indices
template <unsigned Tdim>
Eigen::VectorXi mpm::Cell<Tdim>::local_node_indices() {
  Eigen::VectorXi indices;
  try {
    indices.resize(nodes_.size());
    indices.setZero();
    unsigned node_idx = 0;
    for (auto node_itr = nodes_.cbegin(); node_itr != nodes_.cend();
         ++node_itr) {
      indices(node_idx) = (*node_itr)->active_id();
      node_idx++;
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return indices;
}

//! Initialise element matrix
template <unsigned Tdim>
bool mpm::Cell<Tdim>::initialise_element_matrix() {
  bool status = true;
  if (this->status()) {
    try {
      // Initialse Laplacian matrix (NxN)
      laplacian_matrix_.resize(nnodes_, nnodes_);
      laplacian_matrix_.setZero();

      // Initialse poisson RHS matrix (Nx(N*Tdim))
      poisson_right_matrix_.resize(nnodes_, nnodes_ * Tdim);
      poisson_right_matrix_.setZero();

      // Initialse correction RHS matrix (NxTdim)
      correction_matrix_.resize(nnodes_, nnodes_ * Tdim);
      correction_matrix_.setZero();

    } catch (std::exception& exception) {
      console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
      status = false;
    }
  }
  return status;
}

//! Initialize local elemental matrices for two-phase one-point solver
template <unsigned Tdim>
bool mpm::Cell<Tdim>::initialise_element_matrix_twophase() {
  bool status = true;
  if (this->status()) {
    try {
      // Initialise drag matrix
      drag_matrix_.resize(Tdim);
      for (unsigned dir = 0; dir < Tdim; dir++) {
        drag_matrix_[dir].resize(nnodes_, nnodes_);
        drag_matrix_[dir].setZero();
      }

      // Initialse Laplacian matrix (NxN)
      laplacian_matrix_.resize(nnodes_, nnodes_);
      laplacian_matrix_.setZero();

      // Initialse poisson RHS matrix (Nx(N*Tdim))
      poisson_right_matrix_twophase_.resize(2);
      for (unsigned phase = 0; phase < poisson_right_matrix_twophase_.size();
           phase++) {
        poisson_right_matrix_twophase_[phase].resize(nnodes_, nnodes_ * Tdim);
        poisson_right_matrix_twophase_[phase].setZero();
      }

      // Initialse correction RHS matrix (NxTdim)
      correction_matrix_twophase_.resize(2);
      for (unsigned phase = 0; phase < correction_matrix_twophase_.size();
           phase++) {
        correction_matrix_twophase_[phase].resize(nnodes_, nnodes_ * Tdim);
        correction_matrix_twophase_[phase].setZero();
      }

    } catch (std::exception& exception) {
      console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
      status = false;
    }
  }
  return status;
}

//! Compute local matrix for drag force coupling
template <unsigned Tdim>
void mpm::Cell<Tdim>::compute_local_drag_matrix(
    const Eigen::VectorXd& shapefn, double pvolume,
    const VectorDim& multiplier) noexcept {

  assert(drag_matrix_.size() == Tdim);

  // Lock the storage
  std::lock_guard<std::mutex> guard(cell_mutex_);
  // Compute local drag matrix
  for (unsigned dir = 0; dir < Tdim; dir++)
    drag_matrix_[dir] +=
        multiplier(dir) * pvolume * shapefn * (shapefn.transpose());
}

//! Compute local matrix of laplacian
template <unsigned Tdim>
void mpm::Cell<Tdim>::compute_local_laplacian(
    const Eigen::MatrixXd& grad_shapefn, double pvolume,
    double multiplier) noexcept {

  std::lock_guard<std::mutex> guard(cell_mutex_);
  laplacian_matrix_ +=
      grad_shapefn * grad_shapefn.transpose() * multiplier * pvolume;
}

//! Compute local poisson RHS matrix
//! Used in poisson equation RHS for Navier Stokes solver
template <unsigned Tdim>
void mpm::Cell<Tdim>::compute_local_poisson_right(
    const Eigen::VectorXd& shapefn, const Eigen::MatrixXd& grad_shapefn,
    double pvolume, double multiplier) noexcept {

  std::lock_guard<std::mutex> guard(cell_mutex_);
  for (unsigned i = 0; i < Tdim; i++) {
    poisson_right_matrix_.block(0, i * nnodes_, nnodes_, nnodes_) +=
        shapefn * grad_shapefn.col(i).transpose() * multiplier * pvolume;
  }
}

//! Compute local poisson RHS matrix
//! Used in poisson equation RHS for TwoPhase solver
template <unsigned Tdim>
void mpm::Cell<Tdim>::compute_local_poisson_right_twophase(
    unsigned phase, const Eigen::VectorXd& shapefn,
    const Eigen::MatrixXd& grad_shapefn, double pvolume,
    double multiplier) noexcept {

  assert(phase < poisson_right_matrix_twophase_.size() &&
         poisson_right_matrix_twophase_.size() == 2);

  std::lock_guard<std::mutex> guard(cell_mutex_);
  for (unsigned i = 0; i < Tdim; i++) {
    poisson_right_matrix_twophase_[phase].block(0, i * nnodes_, nnodes_,
                                                nnodes_) +=
        shapefn * grad_shapefn.col(i).transpose() * multiplier * pvolume;
  }
}

//! Compute local correction matrix
//! Used to compute corrector of nodal velocity for Navier Stokes solver
template <unsigned Tdim>
void mpm::Cell<Tdim>::compute_local_correction_matrix(
    const Eigen::VectorXd& shapefn, const Eigen::MatrixXd& grad_shapefn,
    double pvolume) noexcept {

  std::lock_guard<std::mutex> guard(cell_mutex_);
  for (unsigned i = 0; i < Tdim; i++) {
    correction_matrix_.block(0, i * nnodes_, nnodes_, nnodes_) +=
        shapefn * grad_shapefn.col(i).transpose() * pvolume;
  }
}

//! Compute local correction matrix
//! Used to compute corrector of nodal velocity for Two Phase solver
template <unsigned Tdim>
void mpm::Cell<Tdim>::compute_local_correction_matrix_twophase(
    unsigned phase, const Eigen::VectorXd& shapefn,
    const Eigen::MatrixXd& grad_shapefn, double pvolume,
    double multiplier) noexcept {

  assert(phase < correction_matrix_twophase_.size() &&
         correction_matrix_twophase_.size() == 2);

  std::lock_guard<std::mutex> guard(cell_mutex_);
  for (unsigned i = 0; i < Tdim; i++) {
    correction_matrix_twophase_[phase].block(0, i * nnodes_, nnodes_,
                                             nnodes_) +=
        shapefn * grad_shapefn.col(i).transpose() * multiplier * pvolume;
  }
}
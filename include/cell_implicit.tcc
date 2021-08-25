//! Initialise element stiffness matrix
template <unsigned Tdim>
bool mpm::Cell<Tdim>::initialise_element_stiffness_matrix() {
  bool status = true;
  if (this->status()) {
    try {
      // Initialse stiffness matrix ((N*Tdim)x(N*Tdim))
      stiffness_matrix_.resize(nnodes_ * Tdim, nnodes_ * Tdim);
      stiffness_matrix_.setZero();

    } catch (std::exception& exception) {
      console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
      status = false;
    }
  }
  return status;
}

//! Compute local material stiffness matrix
template <unsigned Tdim>
void mpm::Cell<Tdim>::compute_local_material_stiffness_matrix(
    const Eigen::MatrixXd& bmatrix, const Eigen::MatrixXd& dmatrix,
    double pvolume, double multiplier) noexcept {

  std::lock_guard<std::mutex> guard(cell_mutex_);
  stiffness_matrix_ +=
      bmatrix.transpose() * dmatrix * bmatrix * multiplier * pvolume;
}

//! Compute local mass matrix
template <unsigned Tdim>
inline void mpm::Cell<Tdim>::compute_local_mass_matrix(
    const Eigen::VectorXd& shapefn, double pvolume,
    double multiplier) noexcept {

  std::lock_guard<std::mutex> guard(cell_mutex_);
  for (unsigned i = 0; i < this->nnodes_; ++i) {
    for (unsigned j = 0; j < this->nnodes_; ++j) {
      for (unsigned k = 0; k < Tdim; ++k) {
        stiffness_matrix_(Tdim * i + k, Tdim * j + k) +=
            shapefn(i) * shapefn(j) * multiplier * pvolume;
      }
    }
  }
}
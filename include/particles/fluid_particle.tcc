//! Construct a two phase particle with id and coordinates
template <unsigned Tdim>
mpm::FluidParticle<Tdim>::FluidParticle(Index id, const VectorDim& coord)
    : mpm::Particle<Tdim>(id, coord) {

  // Initialize vector data properties
  this->properties_["pressure"] = [&]() {
    Eigen::VectorXd vec_pressure = Eigen::VectorXd::Zero(3);
    vec_pressure[0] = this->pressure();
    // FIXME: This is to check free surface particles
    // TODO: To be removed somewhere
    vec_pressure[1] = this->free_surface();
    return vec_pressure;
  };

  // Logger
  std::string logger =
      "FluidParticle" + std::to_string(Tdim) + "d::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

// Compute stress
template <unsigned Tdim>
void mpm::FluidParticle<Tdim>::compute_stress() noexcept {
  // Run particle compute stress
  mpm::Particle<Tdim>::compute_stress();

  // Calculate fluid turbulent stress
  this->stress_ += this->compute_turbulent_stress();
}

// Compute turbulent stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1>
    mpm::FluidParticle<Tdim>::compute_turbulent_stress() {
  // Compute turbulent stress depends on the model
  Eigen::Matrix<double, 6, 1> tstress;
  tstress.setZero();

  // TODO: To be implemented

  return tstress;
}

//! Map internal force
template <>
inline void mpm::FluidParticle<1>::map_internal_force() noexcept {
  // initialise a vector of total stress (deviatoric + turbulent - pressure)
  Eigen::Matrix<double, 6, 1> total_stress = this->stress_;
  total_stress(0) -= this->projection_param_ * state_variables_.at("pressure");
  // Compute nodal internal forces
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    // Compute force: -pstress * volume
    Eigen::Matrix<double, 1, 1> force;
    force[0] = dn_dx_(i, 0) * total_stress[0];
    force *= -1 * this->volume_;

    nodes_[i]->update_internal_force(true, mpm::ParticlePhase::SinglePhase,
                                     force);
  }
}

//! Map internal force
template <>
inline void mpm::FluidParticle<2>::map_internal_force() noexcept {
  // initialise a vector of total stress (deviatoric + turbulent - pressure)
  Eigen::Matrix<double, 6, 1> total_stress = this->stress_;
  total_stress(0) -= this->projection_param_ * state_variables_.at("pressure");
  total_stress(1) -= this->projection_param_ * state_variables_.at("pressure");

  // Compute nodal internal forces
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    // Compute force: -pstress * volume
    Eigen::Matrix<double, 2, 1> force;
    force[0] = dn_dx_(i, 0) * total_stress[0] + dn_dx_(i, 1) * total_stress[3];
    force[1] = dn_dx_(i, 1) * total_stress[1] + dn_dx_(i, 0) * total_stress[3];

    force *= -1. * this->volume_;

    nodes_[i]->update_internal_force(true, mpm::ParticlePhase::SinglePhase,
                                     force);
  }
}

//! Map internal force
template <>
inline void mpm::FluidParticle<3>::map_internal_force() noexcept {
  // initialise a vector of total stress (deviatoric + turbulent - pressure)
  Eigen::Matrix<double, 6, 1> total_stress = this->stress_;
  total_stress(0) -= this->projection_param_ * state_variables_.at("pressure");
  total_stress(1) -= this->projection_param_ * state_variables_.at("pressure");
  total_stress(2) -= this->projection_param_ * state_variables_.at("pressure");

  // Compute nodal internal forces
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    // Compute force: -pstress * volume
    Eigen::Matrix<double, 3, 1> force;
    force[0] = dn_dx_(i, 0) * total_stress[0] + dn_dx_(i, 1) * total_stress[3] +
               dn_dx_(i, 2) * total_stress[5];

    force[1] = dn_dx_(i, 1) * total_stress[1] + dn_dx_(i, 0) * total_stress[3] +
               dn_dx_(i, 2) * total_stress[4];

    force[2] = dn_dx_(i, 2) * total_stress[2] + dn_dx_(i, 1) * total_stress[4] +
               dn_dx_(i, 0) * total_stress[5];

    force *= -1. * this->volume_;

    nodes_[i]->update_internal_force(true, mpm::ParticlePhase::SinglePhase,
                                     force);
  }
}

//! Map laplacian element matrix to cell (used in poisson equation LHS)
template <unsigned Tdim>
bool mpm::FluidParticle<Tdim>::map_laplacian_to_cell() {
  bool status = true;
  try {
    // Compute local matrix of Laplacian
    cell_->compute_local_laplacian(dn_dx_, volume_);
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Map poisson rhs element matrix to cell (used in poisson equation RHS)
template <unsigned Tdim>
bool mpm::FluidParticle<Tdim>::map_poisson_right_to_cell() {
  bool status = true;
  try {
    // Compute local poisson rhs matrix
    cell_->compute_local_poisson_right(
        shapefn_, dn_dx_,
        material_->template property<double>(std::string("density")) * volume_);
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Map correction matrix element matrix to cell (used to correct velocity)
template <unsigned Tdim>
bool mpm::FluidParticle<Tdim>::map_correction_matrix_to_cell() {
  bool status = true;
  try {
    // TODO: Tobe uncomment once cell is implemented
    // cell_->compute_K_cor_element(shapefn_, dn_dx_, volume_);

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return status;
}

// Compute updated pore pressure of the particle based on nodal pressure
template <unsigned Tdim>
bool mpm::FluidParticle<Tdim>::compute_updated_pressure() {
  bool status = true;

  try {
    double pressure_increment = 0;
    for (unsigned i = 0; i < nodes_.size(); ++i) {
      pressure_increment += shapefn_(i) * nodes_[i]->pressure_increment();
    }
    // Get interpolated nodal pore pressure
    state_variables_.at("pressure") =
        state_variables_.at("pressure") * projection_param_ +
        pressure_increment;
    // Apply free surface
    if (this->free_surface()) state_variables_.at("pressure") = 0.0;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}
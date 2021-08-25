//! Map particle mass, momentum and inertia to nodes
template <unsigned Tdim>
void mpm::Particle<Tdim>::map_mass_momentum_inertia_to_nodes() noexcept {
  // Map mass and momentum to nodes
  this->map_mass_momentum_to_nodes();

  // Map inertia to nodes
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    nodes_[i]->update_inertia(true, mpm::ParticlePhase::Solid,
                              mass_ * shapefn_[i] * acceleration_);
  }
}

// Initialise particle displacement and strain rate at the beginning of time
// step
template <unsigned Tdim>
void mpm::Particle<Tdim>::initialise_displacement_strain_rate() {
  displacement_.setZero();
  strain_rate_.setZero();
}

//! Map inertial force
template <unsigned Tdim>
void mpm::Particle<Tdim>::map_inertial_force() noexcept {
  // Check if particle has a valid cell ptr
  assert(cell_ != nullptr);
  // Get interpolated nodal acceleration
  Eigen::Matrix<double, Tdim, 1> nodal_acceleration =
      Eigen::Matrix<double, Tdim, 1>::Zero();

  for (unsigned i = 0; i < nodes_.size(); ++i)
    nodal_acceleration +=
        shapefn_[i] * nodes_[i]->acceleration(mpm::ParticlePhase::Solid);

  // Compute nodal inertial forces
  for (unsigned i = 0; i < nodes_.size(); ++i)
    nodes_[i]->update_external_force(
        true, mpm::ParticlePhase::Solid,
        (-1. * nodal_acceleration * mass_ * shapefn_(i)));
}

//! Map material stiffness matrix to cell (used in poisson equation LHS)
template <unsigned Tdim>
inline bool mpm::Particle<Tdim>::map_material_stiffness_matrix_to_cell() {
  bool status = true;
  try {
    // Check if material ptr is valid
    assert(this->material() != nullptr);
    // Calculate constitutive relations matrix
    Eigen::MatrixXd dmatrix;
    dmatrix =
        (this->material())
            ->compute_dmatrix(stress_, dstrain_, this,
                              &state_variables_[mpm::ParticlePhase::Solid]);

    // Reduce constitutive relations matrix depending on the dimension
    Eigen::MatrixXd reduced_dmatrix;
    reduced_dmatrix = this->reduce_dmatrix(dmatrix);

    // Calculate B matrix
    Eigen::MatrixXd bmatrix;
    bmatrix = this->compute_bmatrix();

    // Compute local material stiffness matrix
    cell_->compute_local_material_stiffness_matrix(bmatrix, reduced_dmatrix,
                                                   volume_);
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Reduce constitutive relations matrix depending on the dimension
template <>
inline Eigen::MatrixXd mpm::Particle<1>::reduce_dmatrix(
    const Eigen::MatrixXd& dmatrix) noexcept {

  // Convert to 1x1 matrix in 1D
  Eigen::MatrixXd dmatrix1x1;
  dmatrix1x1.resize(1, 1);
  dmatrix1x1(0, 0) = dmatrix(0, 0);

  return dmatrix1x1;
}

//! Reduce constitutive relations matrix depending on the dimension
template <>
inline Eigen::MatrixXd mpm::Particle<2>::reduce_dmatrix(
    const Eigen::MatrixXd& dmatrix) noexcept {

  // Convert to 3x3 matrix in 2D
  Eigen::MatrixXd dmatrix3x3;
  dmatrix3x3.resize(3, 3);
  dmatrix3x3(0, 0) = dmatrix(0, 0);
  dmatrix3x3(0, 1) = dmatrix(0, 1);
  dmatrix3x3(0, 2) = 0.;
  dmatrix3x3(1, 0) = dmatrix(1, 0);
  dmatrix3x3(1, 1) = dmatrix(1, 1);
  dmatrix3x3(1, 2) = 0.;
  dmatrix3x3(2, 0) = 0.;
  dmatrix3x3(2, 1) = 0.;
  dmatrix3x3(2, 2) = dmatrix(4, 4);

  return dmatrix3x3;
}

//! Reduce constitutive relations matrix depending on the dimension
template <>
inline Eigen::MatrixXd mpm::Particle<3>::reduce_dmatrix(
    const Eigen::MatrixXd& dmatrix) noexcept {
  return dmatrix;
}

//! Map mass matrix to cell (used in poisson equation LHS)
template <unsigned Tdim>
inline bool mpm::Particle<Tdim>::map_mass_matrix_to_cell(double multiplier) {
  bool status = true;
  try {
    // Check if material ptr is valid
    assert(this->material() != nullptr);

    // Compute local mass matrix
    cell_->compute_local_mass_matrix(shapefn_, volume_,
                                     mass_density_ * multiplier);
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Compute updated position of the particle by Newmark scheme
template <unsigned Tdim>
void mpm::Particle<Tdim>::compute_updated_position_newmark(
    double dt, bool velocity_update) noexcept {
  // Check if particle has a valid cell ptr
  assert(cell_ != nullptr);
  // Get interpolated nodal displacement and acceleration
  Eigen::Matrix<double, Tdim, 1> nodal_displacement =
      Eigen::Matrix<double, Tdim, 1>::Zero();
  Eigen::Matrix<double, Tdim, 1> nodal_acceleration =
      Eigen::Matrix<double, Tdim, 1>::Zero();
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    nodal_displacement +=
        shapefn_[i] * nodes_[i]->displacement(mpm::ParticlePhase::Solid);
    nodal_acceleration +=
        shapefn_[i] * nodes_[i]->acceleration(mpm::ParticlePhase::Solid);
  }

  // Acceleration update of velocity
  if (!velocity_update) {
    // Update particle velocity from interpolated nodal acceleration
    this->velocity_ += 0.5 * (this->acceleration_ + nodal_acceleration) * dt;
  }
  // Update particle velocity using interpolated nodal velocity
  else {
    // Get interpolated nodal velocity
    Eigen::Matrix<double, Tdim, 1> nodal_velocity =
        Eigen::Matrix<double, Tdim, 1>::Zero();
    for (unsigned i = 0; i < nodes_.size(); ++i)
      nodal_velocity +=
          shapefn_[i] * nodes_[i]->velocity(mpm::ParticlePhase::Solid);

    this->velocity_ = nodal_velocity;
  }

  // Update acceleration
  this->acceleration_ = nodal_acceleration;

  // New position  current position + (displacement - previous displacement)
  this->coordinates_ += nodal_displacement - this->displacement_;
  // Update displacement (displacement is initialized from zero)
  this->displacement_ = nodal_displacement;
}

// Compute strain increment of the particle
template <>
inline Eigen::Matrix<double, 6, 1> mpm::Particle<1>::compute_strain_increment(
    const Eigen::MatrixXd& dn_dx, unsigned phase) noexcept {
  // Define strain rincrement
  Eigen::Matrix<double, 6, 1> strain_increment =
      Eigen::Matrix<double, 6, 1>::Zero();

  for (unsigned i = 0; i < this->nodes_.size(); ++i) {
    Eigen::Matrix<double, 1, 1> dsp = nodes_[i]->displacement(phase);
    strain_increment[0] += dn_dx(i, 0) * dsp[0];
  }

  if (std::fabs(strain_increment(0)) < 1.E-15) strain_increment[0] = 0.;
  return strain_increment;
}

// Compute strain increment of the particle
template <>
inline Eigen::Matrix<double, 6, 1> mpm::Particle<2>::compute_strain_increment(
    const Eigen::MatrixXd& dn_dx, unsigned phase) noexcept {
  // Define strain increment
  Eigen::Matrix<double, 6, 1> strain_increment =
      Eigen::Matrix<double, 6, 1>::Zero();

  for (unsigned i = 0; i < this->nodes_.size(); ++i) {
    Eigen::Matrix<double, 2, 1> dsp = nodes_[i]->displacement(phase);
    strain_increment[0] += dn_dx(i, 0) * dsp[0];
    strain_increment[1] += dn_dx(i, 1) * dsp[1];
    strain_increment[3] += dn_dx(i, 1) * dsp[0] + dn_dx(i, 0) * dsp[1];
  }

  if (std::fabs(strain_increment[0]) < 1.E-15) strain_increment[0] = 0.;
  if (std::fabs(strain_increment[1]) < 1.E-15) strain_increment[1] = 0.;
  if (std::fabs(strain_increment[3]) < 1.E-15) strain_increment[3] = 0.;
  return strain_increment;
}

// Compute strain increment of the particle
template <>
inline Eigen::Matrix<double, 6, 1> mpm::Particle<3>::compute_strain_increment(
    const Eigen::MatrixXd& dn_dx, unsigned phase) noexcept {
  // Define strain increment
  Eigen::Matrix<double, 6, 1> strain_increment =
      Eigen::Matrix<double, 6, 1>::Zero();

  for (unsigned i = 0; i < this->nodes_.size(); ++i) {
    Eigen::Matrix<double, 3, 1> dsp = nodes_[i]->displacement(phase);
    strain_increment[0] += dn_dx(i, 0) * dsp[0];
    strain_increment[1] += dn_dx(i, 1) * dsp[1];
    strain_increment[2] += dn_dx(i, 2) * dsp[2];
    strain_increment[3] += dn_dx(i, 1) * dsp[0] + dn_dx(i, 0) * dsp[1];
    strain_increment[4] += dn_dx(i, 2) * dsp[1] + dn_dx(i, 1) * dsp[2];
    strain_increment[5] += dn_dx(i, 2) * dsp[0] + dn_dx(i, 0) * dsp[2];
  }

  for (unsigned i = 0; i < strain_increment.size(); ++i)
    if (std::fabs(strain_increment[i]) < 1.E-15) strain_increment[i] = 0.;
  return strain_increment;
}

// Compute strain of the particle using nodal displacement
template <unsigned Tdim>
void mpm::Particle<Tdim>::compute_strain_newmark() noexcept {
  // Compute total strain increment within time step
  const Eigen::Matrix<double, 6, 1> strain_increment =
      this->compute_strain_increment(dn_dx_, mpm::ParticlePhase::Solid);
  // Update dstrain = latest total strain increment - previous total strain
  // increment
  dstrain_ = strain_increment - strain_rate_;
  // Update total strain increment within time step
  strain_rate_ = strain_increment;
  // Update strain += dstrain
  strain_ += dstrain_;

  // Compute at centroid
  // Strain rate for reduced integration
  const Eigen::Matrix<double, 6, 1> strain_increment_centroid =
      this->compute_strain_increment(dn_dx_centroid_,
                                     mpm::ParticlePhase::Solid);

  // Assign volumetric strain at centroid
  dvolumetric_strain_ = strain_increment_centroid.head(Tdim).sum();
  volumetric_strain_centroid_ += dvolumetric_strain_;
}
//! Update drag force coefficient
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_drag_force_coefficient(
    bool update, const Eigen::Matrix<double, Tdim, 1>& drag_force_coefficient) {

  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign drag force coefficient
  node_mutex_.lock();
  drag_force_coefficient_ =
      drag_force_coefficient_ * factor + drag_force_coefficient;
  node_mutex_.unlock();
}

//! Compute acceleration and velocity for two phase
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::
    compute_acceleration_velocity_twophase_explicit(double dt) noexcept {
  bool status = false;
  const double tolerance = 1.0E-15;
  if (this->mass(mpm::NodePhase::nSolid) > tolerance &&
      this->mass(mpm::NodePhase::nLiquid) > tolerance) {
    // Compute drag force
    VectorDim drag_force = drag_force_coefficient_.cwiseProduct(
        velocity_.col(mpm::NodePhase::nLiquid) -
        velocity_.col(mpm::NodePhase::nSolid));

    // Acceleration of pore fluid (momentume balance of fluid phase)
    this->acceleration_.col(mpm::NodePhase::nLiquid) =
        (this->external_force_.col(mpm::NodePhase::nLiquid) +
         this->internal_force_.col(mpm::NodePhase::nLiquid) - drag_force) /
        this->mass_(mpm::NodePhase::nLiquid);

    // Acceleration of solid skeleton (momentume balance of mixture)
    this->acceleration_.col(mpm::NodePhase::nSolid) =
        (this->external_force_.col(mpm::NodePhase::nMixture) +
         this->internal_force_.col(mpm::NodePhase::nMixture) -
         this->mass_(mpm::NodePhase::nLiquid) *
             this->acceleration_.col(mpm::NodePhase::nLiquid)) /
        this->mass_(mpm::NodePhase::nSolid);

    // Apply friction constraints
    this->apply_friction_constraints(dt);

    // Velocity += acceleration * dt
    this->velocity_ += this->acceleration_ * dt;

    // Apply velocity constraints, which also sets acceleration to 0,
    // when velocity is set.
    this->apply_velocity_constraints();

    // Set a threshold
    for (unsigned i = 0; i < Tdim; ++i) {
      if (std::abs(velocity_.col(mpm::NodePhase::nSolid)(i)) < tolerance)
        velocity_.col(mpm::NodePhase::nSolid)(i) = 0.;
      if (std::abs(acceleration_.col(mpm::NodePhase::nSolid)(i)) < tolerance)
        acceleration_.col(mpm::NodePhase::nSolid)(i) = 0.;
      if (std::abs(velocity_.col(mpm::NodePhase::nLiquid)(i)) < tolerance)
        velocity_.col(mpm::NodePhase::nLiquid)(i) = 0.;
      if (std::abs(acceleration_.col(mpm::NodePhase::nLiquid)(i)) < tolerance)
        acceleration_.col(mpm::NodePhase::nLiquid)(i) = 0.;
    }
    status = true;
  }
  return status;
}

//! Compute acceleration and velocity for two phase with damping
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::
    compute_acceleration_velocity_twophase_explicit_cundall(
        double dt, double damping_factor) noexcept {
  bool status = false;
  const double tolerance = 1.0E-15;

  if (this->mass(mpm::NodePhase::nSolid) > tolerance &&
      this->mass(mpm::NodePhase::nLiquid) > tolerance) {
    // Compute drag force
    VectorDim drag_force = drag_force_coefficient_.cwiseProduct(
        velocity_.col(mpm::NodePhase::nLiquid) -
        velocity_.col(mpm::NodePhase::nSolid));

    // Unbalanced force of liquid phase
    auto unbalanced_force_liquid =
        this->external_force_.col(mpm::NodePhase::nLiquid) +
        this->internal_force_.col(mpm::NodePhase::nLiquid) - drag_force;
    // Acceleration of liquid phase (momentume balance of fluid phase)
    this->acceleration_.col(mpm::NodePhase::nLiquid) =
        (unbalanced_force_liquid -
         damping_factor * unbalanced_force_liquid.norm() *
             this->velocity_.col(mpm::NodePhase::nLiquid).cwiseSign()) /
        this->mass(mpm::NodePhase::nLiquid);

    // Unbalanced force of solid phase
    auto unbalanced_force_solid =
        this->external_force_.col(mpm::NodePhase::nMixture) +
        this->internal_force_.col(mpm::NodePhase::nMixture) -
        this->mass_(mpm::NodePhase::nLiquid) *
            this->acceleration_.col(mpm::NodePhase::nLiquid);
    // Acceleration of solid phase (momentume balance of mixture)
    this->acceleration_.col(mpm::NodePhase::nSolid) =
        (unbalanced_force_solid -
         damping_factor * unbalanced_force_solid.norm() *
             this->velocity_.col(mpm::NodePhase::nSolid).cwiseSign()) /
        this->mass(mpm::NodePhase::nSolid);

    // Apply friction constraints
    this->apply_friction_constraints(dt);

    // Velocity += acceleration * dt
    this->velocity_ += this->acceleration_ * dt;

    // Apply velocity constraints, which also sets acceleration to 0,
    // when velocity is set.
    this->apply_velocity_constraints();

    // Set a threshold
    for (unsigned i = 0; i < Tdim; ++i) {
      if (std::abs(velocity_.col(mpm::NodePhase::nSolid)(i)) < tolerance)
        velocity_.col(mpm::NodePhase::nSolid)(i) = 0.;
      if (std::abs(acceleration_.col(mpm::NodePhase::nSolid)(i)) < tolerance)
        acceleration_.col(mpm::NodePhase::nSolid)(i) = 0.;
      if (std::abs(velocity_.col(mpm::NodePhase::nLiquid)(i)) < tolerance)
        velocity_.col(mpm::NodePhase::nLiquid)(i) = 0.;
      if (std::abs(acceleration_.col(mpm::NodePhase::nLiquid)(i)) < tolerance)
        acceleration_.col(mpm::NodePhase::nLiquid)(i) = 0.;
    }
    status = true;
  }
  return status;
}

//! Compute intermediate force
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::compute_intermediate_force(
    const double dt) {
  bool status = true;

  try {
    node_mutex_.lock();

    // Total force matrix
    auto force_total = internal_force_ + external_force_;
    // Force vector for mixture
    force_total_inter_ = force_total.col(mpm::NodePhase::nMixture);
    // Force vector for liquid
    force_fluid_inter_ = force_total.col(mpm::NodePhase::nLiquid) - drag_force_;

    node_mutex_.unlock();
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Compute semi-implicit acceleration and velocity
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::
    compute_acceleration_velocity_twophase_semi_implicit(unsigned phase,
                                                         double dt) {
  bool status = true;
  const double tolerance = std::numeric_limits<double>::min();
  try {

    // Semi-implicit solver
    Eigen::Matrix<double, Tdim, 1> acceleration_corrected =
        correction_force_.col(phase) / mass_(phase);

    // Acceleration
    this->acceleration_.col(phase) = acceleration_corrected;

    // Update velocity
    velocity_.col(phase) += acceleration_corrected * dt;

    // Apply friction constraints
    this->apply_friction_constraints(dt);

    // Apply velocity constraints, which also sets acceleration to 0,
    // when velocity is set.
    this->apply_velocity_constraints();

    // Set a threshold
    for (unsigned i = 0; i < Tdim; ++i) {
      if (std::abs(velocity_.col(phase)(i)) < tolerance)
        velocity_.col(phase)(i) = 0.;
      if (std::abs(acceleration_.col(phase)(i)) < tolerance)
        acceleration_.col(phase)(i) = 0.;
    }

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Compute semi-implicit acceleration and velocity  with cundall damping factor
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::
    compute_acceleration_velocity_twophase_semi_implicit_cundall(
        unsigned phase, double dt, double damping_factor) {
  bool status = true;
  const double tolerance = std::numeric_limits<double>::min();
  try {
    if (mass_(phase) > tolerance) {
      // Unbalance force
      auto unbalanced_force = correction_force_.col(phase);

      // Semi-implicit solver
      auto acceleration_corrected =
          (unbalanced_force - damping_factor * unbalanced_force.norm() *
                                  this->velocity_.col(phase).cwiseSign()) /
          this->mass_(phase);

      // Acceleration
      this->acceleration_.col(phase) = acceleration_corrected;

      // Update velocity
      velocity_.col(phase) += acceleration_corrected * dt;

      // Apply friction constraints
      this->apply_friction_constraints(dt);

      // Apply velocity constraints, which also sets acceleration to 0,
      // when velocity is set.
      this->apply_velocity_constraints();

      // Set a threshold
      for (unsigned i = 0; i < Tdim; ++i) {
        if (std::abs(velocity_.col(phase)(i)) < tolerance)
          velocity_.col(phase)(i) = 0.;
        if (std::abs(acceleration_.col(phase)(i)) < tolerance)
          acceleration_.col(phase)(i) = 0.;
      }
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

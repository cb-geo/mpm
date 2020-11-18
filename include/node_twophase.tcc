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
  if (this->mass(mpm::NodePhase::NSolid) > tolerance &&
      this->mass(mpm::NodePhase::NLiquid) > tolerance) {
    // Compute drag force
    VectorDim drag_force = drag_force_coefficient_.cwiseProduct(
        velocity_.col(mpm::NodePhase::NLiquid) -
        velocity_.col(mpm::NodePhase::NSolid));

    // Acceleration of pore fluid (momentume balance of fluid phase)
    this->acceleration_.col(mpm::NodePhase::NLiquid) =
        (this->external_force_.col(mpm::NodePhase::NLiquid) +
         this->internal_force_.col(mpm::NodePhase::NLiquid) - drag_force) /
        this->mass_(mpm::NodePhase::NLiquid);

    // Acceleration of solid skeleton (momentume balance of mixture)
    this->acceleration_.col(mpm::NodePhase::NSolid) =
        (this->external_force_.col(mpm::NodePhase::NMixture) +
         this->internal_force_.col(mpm::NodePhase::NMixture) -
         this->mass_(mpm::NodePhase::NLiquid) *
             this->acceleration_.col(mpm::NodePhase::NLiquid)) /
        this->mass_(mpm::NodePhase::NSolid);

    // Apply friction constraints
    this->apply_friction_constraints(dt);

    // Velocity += acceleration * dt
    this->velocity_ += this->acceleration_ * dt;

    // Apply velocity constraints, which also sets acceleration to 0,
    // when velocity is set.
    this->apply_velocity_constraints();

    // Set a threshold
    for (unsigned i = 0; i < Tdim; ++i) {
      if (std::abs(velocity_.col(mpm::NodePhase::NSolid)(i)) < tolerance)
        velocity_.col(mpm::NodePhase::NSolid)(i) = 0.;
      if (std::abs(acceleration_.col(mpm::NodePhase::NSolid)(i)) < tolerance)
        acceleration_.col(mpm::NodePhase::NSolid)(i) = 0.;
      if (std::abs(velocity_.col(mpm::NodePhase::NLiquid)(i)) < tolerance)
        velocity_.col(mpm::NodePhase::NLiquid)(i) = 0.;
      if (std::abs(acceleration_.col(mpm::NodePhase::NLiquid)(i)) < tolerance)
        acceleration_.col(mpm::NodePhase::NLiquid)(i) = 0.;
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

  if (this->mass(mpm::NodePhase::NSolid) > tolerance &&
      this->mass(mpm::NodePhase::NLiquid) > tolerance) {
    // Compute drag force
    VectorDim drag_force = drag_force_coefficient_.cwiseProduct(
        velocity_.col(mpm::NodePhase::NLiquid) -
        velocity_.col(mpm::NodePhase::NSolid));

    // Unbalanced force of liquid phase
    auto unbalanced_force_liquid =
        this->external_force_.col(mpm::NodePhase::NLiquid) +
        this->internal_force_.col(mpm::NodePhase::NLiquid) - drag_force;
    // Acceleration of liquid phase (momentume balance of fluid phase)
    this->acceleration_.col(mpm::NodePhase::NLiquid) =
        (unbalanced_force_liquid -
         damping_factor * unbalanced_force_liquid.norm() *
             this->velocity_.col(mpm::NodePhase::NLiquid).cwiseSign()) /
        this->mass(mpm::NodePhase::NLiquid);

    // Unbalanced force of solid phase
    auto unbalanced_force_solid =
        this->external_force_.col(mpm::NodePhase::NMixture) +
        this->internal_force_.col(mpm::NodePhase::NMixture) -
        this->mass_(mpm::NodePhase::NLiquid) *
            this->acceleration_.col(mpm::NodePhase::NLiquid);
    // Acceleration of solid phase (momentume balance of mixture)
    this->acceleration_.col(mpm::NodePhase::NSolid) =
        (unbalanced_force_solid -
         damping_factor * unbalanced_force_solid.norm() *
             this->velocity_.col(mpm::NodePhase::NSolid).cwiseSign()) /
        this->mass(mpm::NodePhase::NSolid);

    // Apply friction constraints
    this->apply_friction_constraints(dt);

    // Velocity += acceleration * dt
    this->velocity_ += this->acceleration_ * dt;

    // Apply velocity constraints, which also sets acceleration to 0,
    // when velocity is set.
    this->apply_velocity_constraints();

    // Set a threshold
    for (unsigned i = 0; i < Tdim; ++i) {
      if (std::abs(velocity_.col(mpm::NodePhase::NSolid)(i)) < tolerance)
        velocity_.col(mpm::NodePhase::NSolid)(i) = 0.;
      if (std::abs(acceleration_.col(mpm::NodePhase::NSolid)(i)) < tolerance)
        acceleration_.col(mpm::NodePhase::NSolid)(i) = 0.;
      if (std::abs(velocity_.col(mpm::NodePhase::NLiquid)(i)) < tolerance)
        velocity_.col(mpm::NodePhase::NLiquid)(i) = 0.;
      if (std::abs(acceleration_.col(mpm::NodePhase::NLiquid)(i)) < tolerance)
        acceleration_.col(mpm::NodePhase::NLiquid)(i) = 0.;
    }
    status = true;
  }
  return status;
}
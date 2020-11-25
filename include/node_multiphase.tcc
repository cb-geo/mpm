//! Compute mass density (Z. Wiezckowski, 2004)
//! density = mass / lumped volume
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::compute_density() {
  const double tolerance = 1.E-16;  // std::numeric_limits<double>::lowest();

  for (unsigned phase = 0; phase < Tnphases; ++phase) {
    if (mass_(phase) > tolerance) {
      if (volume_(phase) > tolerance)
        density_(phase) = mass_(phase) / volume_(phase);

      // Check to see if value is below threshold
      if (std::abs(density_(phase)) < tolerance) density_(phase) = 0.;
    }
  }
}

//! Initialise two-phase nodal properties
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::initialise_twophase() noexcept {
  this->initialise();
  // Specific variables for two phase
  drag_force_coefficient_.setZero();
  drag_force_.setZero();
  force_total_inter_.setZero();
  force_fluid_inter_.setZero();
  velocity_inter_.setZero();
  acceleration_inter_.setZero();
}

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

//! Compute semi-implicit acceleration and velocity
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::
    compute_acceleration_velocity_semi_implicit_corrector(unsigned phase,
                                                          double dt) {
  bool status = false;
  const double tolerance = std::numeric_limits<double>::min();
  if (mass_(phase) > tolerance) {
    Eigen::Matrix<double, Tdim, 1> acceleration_corrected =
        correction_force_.col(phase) / mass_(phase);

    // Acceleration
    this->acceleration_.col(phase) += acceleration_corrected;

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

    status = true;
  }

  return status;
}

//! Compute semi-implicit acceleration and velocity  with cundall damping
//! factor
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::
    compute_acceleration_velocity_semi_implicit_corrector_cundall(
        unsigned phase, double dt, double damping_factor) {
  bool status = false;
  const double tolerance = std::numeric_limits<double>::min();

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

    status = true;
  }

  return status;
}

//! Update pressure increment at the node
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_pressure_increment(
    const Eigen::VectorXd& pressure_increment, unsigned phase,
    double current_time) {
  this->pressure_increment_ = pressure_increment(active_id_);

  // If pressure boundary, increment is zero
  if (pressure_constraints_.find(phase) != pressure_constraints_.end() ||
      this->free_surface())
    this->pressure_increment_ = 0;
}

//! Compute intermediate force
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::compute_intermediate_force() {
  // Total force matrix
  const auto force_total = internal_force_ + external_force_;
  // Force vector for mixture
  force_total_inter_ = force_total.col(mpm::NodePhase::NMixture);
  // Force vector for liquid
  force_fluid_inter_ = force_total.col(mpm::NodePhase::NLiquid) - drag_force_;
}

//! Update intermediate acceleration and velocity at the node
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_intermediate_acceleration_velocity(
    const unsigned phase, const Eigen::MatrixXd& acceleration_inter,
    double dt) {
  // Update nodal intermediate acceleration
  acceleration_inter_.col(phase) =
      acceleration_inter.row(active_id_).transpose();

  // Update nodal intermediate velocity
  velocity_inter_.col(phase) =
      velocity_.col(phase) +
      dt * acceleration_inter.row(active_id_).transpose();
}

//! Update correction force
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_correction_force(
    bool update, unsigned phase,
    const Eigen::Matrix<double, Tdim, 1>& force) noexcept {
  // Assert
  assert(phase < Tnphases);

  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign correction force
  node_mutex_.lock();
  correction_force_.col(phase) = correction_force_.col(phase) * factor + force;
  node_mutex_.unlock();
}

//! Compute nodal correction force
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::compute_nodal_correction_force(
    const VectorDim& correction_force) {
  // Compute correction force for fluid phase
  correction_force_.col(0) = correction_force;
}

//! Compute nodal corrected force for two-phase
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::compute_nodal_correction_force(
    const VectorDim& solid_correction_force,
    const VectorDim& liquid_correction_force) {
  // Compute corrected force for solid phase
  correction_force_.col(0) =
      mass_(0) * acceleration_inter_.col(0) + solid_correction_force;
  // Compute corrected force for liquid phase
  correction_force_.col(1) =
      mass_(1) * acceleration_inter_.col(1) + liquid_correction_force;
}
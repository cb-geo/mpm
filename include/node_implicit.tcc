//! Assign nodal inertia
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_inertia(
    bool update, unsigned phase,
    const Eigen::Matrix<double, Tdim, 1>& inertia) noexcept {
  // Assert
  assert(phase < Tnphases);

  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign inertia
  node_mutex_.lock();
  inertia_.col(phase) = inertia_.col(phase) * factor + inertia;
  node_mutex_.unlock();
}

//! Compute velocity and acceleration from momentum and inertia
//! velocity = momentum / mass
//! acceleration = inertia / mass
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::compute_velocity_acceleration() {
  const double tolerance = 1.E-16;
  for (unsigned phase = 0; phase < Tnphases; ++phase) {
    if (mass_(phase) > tolerance) {
      velocity_.col(phase) = momentum_.col(phase) / mass_(phase);
      acceleration_.col(phase) = inertia_.col(phase) / mass_(phase);

      // Check to see if value is below threshold
      for (unsigned i = 0; i < velocity_.rows(); ++i)
        if (std::abs(velocity_.col(phase)(i)) < 1.E-15)
          velocity_.col(phase)(i) = 0.;

      for (unsigned i = 0; i < acceleration_.rows(); ++i)
        if (std::abs(acceleration_.col(phase)(i)) < 1.E-15)
          acceleration_.col(phase)(i) = 0.;
    }
  }
}

//! Update velocity and acceleration by Newmark scheme
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_velocity_acceleration_newmark(
    unsigned phase, double newmark_beta, double newmark_gamma, double dt) {
  const double tolerance = 1.E-16;
  VectorDim velocity_pre;
  VectorDim acceleration_pre;
  //! Compute velocity and acceleration at the previous time step
  if (mass_(phase) > tolerance) {
    velocity_pre = momentum_.col(phase) / mass_(phase);
    acceleration_pre = inertia_.col(phase) / mass_(phase);
  }

  //! Update of velocity and acceleration
  velocity_.col(phase) =
      newmark_gamma / newmark_beta / dt * displacement_.col(phase) -
      (newmark_gamma / newmark_beta - 1.) * velocity_pre -
      0.5 * dt * (newmark_gamma / newmark_beta - 2.) * acceleration_pre;

  acceleration_.col(phase) =
      1. / newmark_beta / dt / dt * displacement_.col(phase) -
      1. / newmark_beta / dt * velocity_pre -
      (0.5 / newmark_beta - 1.) * acceleration_pre;

  // Check to see if value is below threshold
  for (unsigned i = 0; i < velocity_.rows(); ++i)
    if (std::abs(velocity_.col(phase)(i)) < 1.E-15)
      velocity_.col(phase)(i) = 0.;

  for (unsigned i = 0; i < acceleration_.rows(); ++i)
    if (std::abs(acceleration_.col(phase)(i)) < 1.E-15)
      acceleration_.col(phase)(i) = 0.;
}

//! Update displacement increment at the node
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_displacement_increment(
    const Eigen::VectorXd& displacement_increment, unsigned phase,
    const unsigned nactive_node) {

  for (unsigned i = 0; i < Tdim; ++i) {
    displacement_.col(phase)(i) +=
        displacement_increment(nactive_node * i + active_id_);
  }
}
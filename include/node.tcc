//! Constructor with id, coordinates and dof
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
mpm::Node<Tdim, Tdof, Tnphases>::Node(
    Index id, const Eigen::Matrix<double, Tdim, 1>& coord)
    : NodeBase<Tdim>(id, coord) {
  // Check if the dimension is between 1 & 3
  static_assert((Tdim >= 1 && Tdim <= 3), "Invalid global dimension");
  id_ = id;
  coordinates_ = coord;
  dof_ = Tdof;

  //! Logger
  std::string logger =
      "node" + std::to_string(Tdim) + "d::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);

  // Clear any velocity constraints
  velocity_constraints_.clear();
  this->initialise();
}

//! Initialise nodal properties
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::initialise() {
  mass_.setZero();
  volume_.setZero();
  external_force_.setZero();
  internal_force_.setZero();
  concentrated_force_.setZero();
  pressure_.setZero();
  velocity_.setZero();
  momentum_.setZero();
  acceleration_.setZero();
  status_ = false;
  material_ids_.clear();
}

//! Update mass at the nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_mass(bool update, unsigned phase,
                                                  double mass) {
  // Decide to update or assign
  double factor = 1.0;
  if (!update) factor = 0.;

  // Update/assign mass
  std::lock_guard<std::mutex> guard(node_mutex_);
  mass_(phase) = (mass_(phase) * factor) + mass;
}

//! Update volume at the nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_volume(bool update, unsigned phase,
                                                    double volume) {
  // Decide to update or assign
  double factor = 1.0;
  if (!update) factor = 0.;

  // Update/assign volume
  std::lock_guard<std::mutex> guard(node_mutex_);
  volume_(phase) = volume_(phase) * factor + volume;
}

// Assign concentrated force to the node
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::assign_concentrated_force(
    unsigned phase, unsigned direction, double concentrated_force,
    const std::shared_ptr<FunctionBase>& function) {
  bool status = false;
  try {
    if (phase >= Tnphases || direction >= Tdim) {
      throw std::runtime_error(
          "Cannot assign nodal concentrated forcey: Direction / phase is "
          "invalid");
    }
    // Assign concentrated force
    concentrated_force_(direction, phase) = concentrated_force;
    status = true;
    this->force_function_ = function;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Apply concentrated force to the node
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::apply_concentrated_force(
    unsigned phase, double current_time) {
  double scalar = 1.0;
  if (force_function_ != nullptr) scalar = force_function_->value(current_time);
  this->update_external_force(true, phase,
                              scalar * concentrated_force_.col(phase));
}

//! Update external force (body force / traction force)
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::update_external_force(
    bool update, unsigned phase, const Eigen::Matrix<double, Tdim, 1>& force) {
  bool status = false;
  try {
    if (phase >= Tnphases)
      throw std::runtime_error("Nodal external force: Invalid phase");

    // Decide to update or assign
    double factor = 1.0;
    if (!update) factor = 0.;

    // Update/assign external force
    std::lock_guard<std::mutex> guard(node_mutex_);
    external_force_.col(phase) = external_force_.col(phase) * factor + force;
    status = true;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Update internal force (body force / traction force)
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::update_internal_force(
    bool update, unsigned phase, const Eigen::Matrix<double, Tdim, 1>& force) {
  bool status = false;
  try {
    if (phase >= Tnphases)
      throw std::runtime_error("Nodal internal force: Invalid phase");

    // Decide to update or assign
    double factor = 1.0;
    if (!update) factor = 0.;

    // Update/assign internal force
    std::lock_guard<std::mutex> guard(node_mutex_);
    internal_force_.col(phase) = internal_force_.col(phase) * factor + force;
    status = true;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assign nodal momentum
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::update_momentum(
    bool update, unsigned phase,
    const Eigen::Matrix<double, Tdim, 1>& momentum) {
  bool status = false;
  try {
    if (phase >= Tnphases)
      throw std::runtime_error("Nodal momentum: Invalid phase");

    // Decide to update or assign
    double factor = 1.0;
    if (!update) factor = 0.;

    // Update/assign momentum
    std::lock_guard<std::mutex> guard(node_mutex_);
    momentum_.col(phase) = momentum_.col(phase) * factor + momentum;
    status = true;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Update pressure at the nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_mass_pressure(
    unsigned phase, double mass_pressure) {
  try {
    const double tolerance = 1.E-16;

    // Compute pressure from mass*pressure
    if (mass_(phase) > tolerance) {
      std::lock_guard<std::mutex> guard(node_mutex_);
      pressure_(phase) += mass_pressure / mass_(phase);
    } else
      throw std::runtime_error("Nodal mass is zero or below threshold");
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
}

//! Assign pressure at the nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::assign_pressure(unsigned phase,
                                                      double pressure) {
  try {
    const double tolerance = 1.E-16;

    // Compute pressure from mass*pressure
    std::lock_guard<std::mutex> guard(node_mutex_);
    pressure_(phase) = pressure;

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
}

//! Compute velocity from momentum
//! velocity = momentum / mass
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::compute_velocity() {
  try {
    const double tolerance = 1.E-16;  // std::numeric_limits<double>::lowest();

    for (unsigned phase = 0; phase < Tnphases; ++phase) {
      if (mass_(phase) > tolerance) {
        velocity_.col(phase) = momentum_.col(phase) / mass_(phase);

        // Check to see if value is below threshold
        for (unsigned i = 0; i < velocity_.rows(); ++i)
          if (std::abs(velocity_.col(phase)(i)) < 1.E-15)
            velocity_.col(phase)(i) = 0.;
      } else
        throw std::runtime_error("Nodal mass is zero or below threshold");
    }

    // Apply velocity constraints, which also sets acceleration to 0,
    // when velocity is set.
    this->apply_velocity_constraints();
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
}

//! Update nodal acceleration
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::update_acceleration(
    bool update, unsigned phase,
    const Eigen::Matrix<double, Tdim, 1>& acceleration) {
  bool status = false;
  try {
    if (phase >= Tnphases)
      throw std::runtime_error("Nodal acceleration: Invalid phase");

    // Decide to update or assign
    double factor = 1.0;
    if (!update) factor = 0.;

    //! Update/assign acceleration
    std::lock_guard<std::mutex> guard(node_mutex_);
    acceleration_.col(phase) = acceleration_.col(phase) * factor + acceleration;
    status = true;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Compute acceleration and velocity
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::compute_acceleration_velocity(
    unsigned phase, double dt) {
  bool status = true;
  const double tolerance = 1.0E-15;
  try {
    if (mass_(phase) > tolerance) {
      // acceleration (unbalaced force / mass)
      this->acceleration_.col(phase) = (this->external_force_.col(phase) +
                                        this->internal_force_.col(phase)) /
                                       this->mass_(phase);

      // Apply friction constraints
      this->apply_friction_constraints(dt);

      // Velocity += acceleration * dt
      this->velocity_.col(phase) += this->acceleration_.col(phase) * dt;
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
    } else
      throw std::runtime_error("Nodal mass is zero or below threshold");

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assign velocity constraint
//! Constrain directions can take values between 0 and Dim * Nphases
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::assign_velocity_constraint(
    unsigned dir, double velocity) {
  bool status = true;
  try {
    //! Constrain directions can take values between 0 and Dim * Nphases
    if (dir < (Tdim * Tnphases))
      this->velocity_constraints_.insert(std::make_pair<unsigned, double>(
          static_cast<unsigned>(dir), static_cast<double>(velocity)));
    else
      throw std::runtime_error("Constraint direction is out of bounds");

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Apply velocity constraints
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::apply_velocity_constraints() {
  // Set velocity constraint
  for (const auto& constraint : this->velocity_constraints_) {
    // Direction value in the constraint (0, Dim * Nphases)
    const unsigned dir = constraint.first;
    // Direction: dir % Tdim (modulus)
    const auto direction = static_cast<unsigned>(dir % Tdim);
    // Phase: Integer value of division (dir / Tdim)
    const auto phase = static_cast<unsigned>(dir / Tdim);

    if (!generic_boundary_constraints_) {
      // Velocity constraints are applied on Cartesian boundaries
      this->velocity_(direction, phase) = constraint.second;
      // Set acceleration to 0 in direction of velocity constraint
      this->acceleration_(direction, phase) = 0.;
    } else {
      // Velocity constraints on general boundaries
      // Compute inverse rotation matrix
      const Eigen::Matrix<double, Tdim, Tdim> inverse_rotation_matrix =
          rotation_matrix_.inverse();
      // Transform to local coordinate
      Eigen::Matrix<double, Tdim, Tnphases> local_velocity =
          inverse_rotation_matrix * this->velocity_;
      Eigen::Matrix<double, Tdim, Tnphases> local_acceleration =
          inverse_rotation_matrix * this->acceleration_;
      // Apply boundary condition in local coordinate
      local_velocity(direction, phase) = constraint.second;
      local_acceleration(direction, phase) = 0.;
      // Transform back to global coordinate
      this->velocity_ = rotation_matrix_ * local_velocity;
      this->acceleration_ = rotation_matrix_ * local_acceleration;
    }
  }
}

//! Assign friction constraint
//! Constrain directions can take values between 0 and Dim * Nphases
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::assign_friction_constraint(
    unsigned dir, int sign_n, double friction) {
  bool status = true;
  try {
    //! Constrain directions can take values between 0 and Dim * Nphases
    if (dir < Tdim) {
      this->friction_constraint_ =
          std::make_tuple(static_cast<unsigned>(dir), static_cast<int>(sign_n),
                          static_cast<double>(friction));
      this->friction_ = true;
    } else
      throw std::runtime_error("Constraint direction is out of bounds");

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Apply friction constraints
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::apply_friction_constraints(double dt) {
  if (friction_) {
    auto sign = [](double value) { return (value > 0.) ? 1. : -1.; };

    // Set friction constraint
    // Direction value in the constraint (0, Dim)
    const unsigned dir_n = std::get<0>(this->friction_constraint_);

    // Normal direction of friction
    const double sign_dir_n = sign(std::get<1>(this->friction_constraint_));

    // Friction co-efficient
    const double mu = std::get<2>(this->friction_constraint_);

    const unsigned phase = 0;

    // Acceleration and velocity
    double acc_n, acc_t, vel_t;

    if (Tdim == 2) {
      // tangential direction to boundary
      const unsigned dir_t = (Tdim - 1) - dir_n;

      if (!generic_boundary_constraints_) {
        // Cartesian case
        // Normal and tangential acceleration
        acc_n = this->acceleration_(dir_n, phase);
        acc_t = this->acceleration_(dir_t, phase);
        // Velocity tangential
        vel_t = this->velocity_(dir_t, phase);
      } else {
        // General case, transform to local coordinate
        // Compute inverse rotation matrix
        const Eigen::Matrix<double, Tdim, Tdim> inverse_rotation_matrix =
            rotation_matrix_.inverse();
        // Transform to local coordinate
        Eigen::Matrix<double, Tdim, Tnphases> local_acceleration =
            inverse_rotation_matrix * this->acceleration_;
        Eigen::Matrix<double, Tdim, Tnphases> local_velocity =
            inverse_rotation_matrix * this->velocity_;
        // Normal and tangential acceleration
        acc_n = local_acceleration(dir_n, phase);
        acc_t = local_acceleration(dir_t, phase);
        // Velocity tangential
        vel_t = local_velocity(dir_t, phase);
      }

      if ((acc_n * sign_dir_n) > 0.0) {
        if (vel_t != 0.0) {  // kinetic friction
          const double vel_net = dt * acc_t + vel_t;
          const double vel_frictional = dt * mu * std::abs(acc_n);
          if (std::abs(vel_net) <= vel_frictional)
            acc_t = -vel_t / dt;
          else
            acc_t -= sign(vel_net) * mu * std::abs(acc_n);
        } else {  // static friction
          if (std::abs(acc_t) <= mu * std::abs(acc_n))
            acc_t = 0.0;
          else
            acc_t -= sign(acc_t) * mu * std::abs(acc_n);
        }

        if (!generic_boundary_constraints_) {
          // Cartesian case
          this->acceleration_(dir_t, phase) = acc_t;
        } else {
          // Local acceleration in terms of tangential and normal
          Eigen::Matrix<double, Tdim, Tnphases> acc;
          acc(dir_t, phase) = acc_t;
          acc(dir_n, phase) = acc_n;

          // General case, transform to global coordinate
          this->acceleration_.col(phase) = rotation_matrix_ * acc.col(phase);
        }
      }
    } else if (Tdim == 3) {
      Eigen::Matrix<int, 3, 2> dir;
      dir(0, 0) = 1;
      dir(0, 1) = 2;  // tangential directions for dir_n = 0
      dir(1, 0) = 0;
      dir(1, 1) = 2;  // tangential directions for dir_n = 1
      dir(2, 0) = 0;
      dir(2, 1) = 1;  // tangential directions for dir_n = 2

      const unsigned dir_t0 = dir(dir_n, 0);
      const unsigned dir_t1 = dir(dir_n, 1);

      Eigen::Matrix<double, Tdim, 1> acc, vel;
      if (!generic_boundary_constraints_) {
        // Cartesian case
        acc = this->acceleration_.col(phase);
        vel = this->velocity_.col(phase);
      } else {
        // General case, transform to local coordinate
        // Compute inverse rotation matrix
        const Eigen::Matrix<double, Tdim, Tdim> inverse_rotation_matrix =
            rotation_matrix_.inverse();
        // Transform to local coordinate
        acc = inverse_rotation_matrix * this->acceleration_.col(phase);
        vel = inverse_rotation_matrix * this->velocity_.col(phase);
      }

      const auto acc_n = acc(dir_n);
      auto acc_t =
          std::sqrt(acc(dir_t0) * acc(dir_t0) + acc(dir_t1) * acc(dir_t1));
      const auto vel_t =
          std::sqrt(vel(dir_t0) * vel(dir_t0) + vel(dir_t1) * vel(dir_t1));

      if (acc_n * sign_dir_n > 0.0) {
        // kinetic friction
        if (vel_t != 0.0) {
          Eigen::Matrix<double, 2, 1> vel_net;
          // friction is applied opposite to the vel_net
          vel_net(0) = vel(dir_t0) + acc(dir_t0) * dt;
          vel_net(1) = vel(dir_t1) + acc(dir_t1) * dt;
          const double vel_net_t =
              sqrt(vel_net(0) * vel_net(0) + vel_net(1) * vel_net(1));
          const double vel_fricion = mu * std::abs(acc_n) * dt;

          if (vel_net_t <= vel_fricion) {
            acc(dir_t0) = -vel(dir_t0);  // To set particle velocity to zero
            acc(dir_t1) = -vel(dir_t1);
          } else {

            acc(dir_t0) -= mu * std::abs(acc_n) * (vel_net(0) / vel_net_t);
            acc(dir_t1) -= mu * std::abs(acc_n) * (vel_net(1) / vel_net_t);
          }
        } else {                                // static friction
          if (acc_t <= mu * std::abs(acc_n)) {  // since acc_t is positive
            acc(dir_t0) = 0;
            acc(dir_t1) = 0;
          } else {
            acc_t -= mu * std::abs(acc_n);
            acc(dir_t0) -= mu * std::abs(acc_n) * (acc(dir_t0) / acc_t);
            acc(dir_t1) -= mu * std::abs(acc_n) * (acc(dir_t1) / acc_t);
          }
        }

        if (!generic_boundary_constraints_) {
          // Cartesian case
          this->acceleration_.col(phase) = acc;
        } else {
          // General case, transform to global coordinate
          this->acceleration_.col(phase) = rotation_matrix_ * acc;
        }
      }
    }
  }
}

//! Add material id from material points to material_ids_
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::append_material_id(unsigned id) {
  std::lock_guard<std::mutex> guard(node_mutex_);
  material_ids_.emplace(id);
}

// Assign MPI rank to node
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::mpi_rank(unsigned rank) {
  std::lock_guard<std::mutex> guard(node_mutex_);
  auto status = this->mpi_ranks_.insert(rank);
  return status.second;
}

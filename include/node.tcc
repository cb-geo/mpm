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
  pressure_.setZero();
  velocity_.setZero();
  momentum_.setZero();
  acceleration_.setZero();
  coordinates_from_particles_.setZero();
  status_ = false;
}

//! Initialise nodal properties of a subdomain
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::initialise_subdomain(unsigned mid) {
  normal_vector_[mid].setZero();
  tangent_vector_[mid].setZero();
  mass_subdomain_[mid].setZero();
  coordinates_from_particles_subdomain_[mid].setZero();
  external_force_subdomain_[mid].setZero();
  internal_force_subdomain_[mid].setZero();
  velocity_subdomain_[mid].setZero();
  momentum_subdomain_[mid].setZero();
  acceleration_subdomain_[mid].setZero();
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

//! Update mass of a subdomain at the nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_mass_subdomain(bool update,
                                                            unsigned phase,
                                                            double mass,
                                                            unsigned mid) {
  // Decide to update or assign
  double factor = 1.0;
  if (!update) factor = 0.;

  // Update/assign subdomain mass
  std::lock_guard<std::mutex> guard(node_mutex_);
  Eigen::Matrix<double, 1, Tnphases> mass_subdomain =
      this->mass_subdomain_.at(mid);
  mass_subdomain(phase) = (mass_subdomain(phase) * factor) + mass;
  this->mass_subdomain_.at(mid) = mass_subdomain;
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

// Assign traction force to the node
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::assign_traction_force(unsigned phase,
                                                            unsigned direction,
                                                            double traction) {
  bool status = false;
  try {
    if (phase < 0 || phase >= Tnphases || direction < 0 || direction >= Tdim) {
      throw std::runtime_error(
          "Nodal traction property: Direction / phase is invalid");
    }
    // Assign traction
    external_force_(direction, phase) = traction;
    status = true;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Assign traction force of a subdomain to the node
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::assign_traction_force_subdomain(
    unsigned phase, unsigned direction, double traction, unsigned mid) {
  bool status = false;
  try {
    if (phase < 0 || phase >= Tnphases || direction < 0 || direction >= Tdim) {
      throw std::runtime_error(
          "Nodal traction property: Direction / phase is invalid");
    }
    // Assign subdomain traction
    auto exforce = this->external_force_subdomain_.at(mid);
    exforce(direction, phase) = traction;
    this->external_force_subdomain_.at(mid) = exforce;

    status = true;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
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

//! Update external force of a subdomain (body force / traction force)
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::update_external_force_subdomain(
    bool update, unsigned phase, const Eigen::Matrix<double, Tdim, 1>& force,
    unsigned mid) {
  bool status = false;
  try {
    if (phase >= Tnphases)
      throw std::runtime_error("Nodal subdomain external force: Invalid phase");

    // Decide to update or assign
    double factor = 1.0;
    if (!update) factor = 0.;

    // Update/assign external force
    std::lock_guard<std::mutex> guard(node_mutex_);
    Eigen::Matrix<double, Tdim, Tnphases> external_force_subdomain =
        this->external_force_subdomain_.at(mid);
    external_force_subdomain.col(phase) =
        external_force_subdomain.col(phase) * factor + force;
    this->external_force_subdomain_.at(mid) = external_force_subdomain;
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

//! Update internal force of a subdomain (body force / traction force)
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::update_internal_force_subdomain(
    bool update, unsigned phase, const Eigen::Matrix<double, Tdim, 1>& force,
    unsigned mid) {
  bool status = false;
  try {
    if (phase >= Tnphases)
      throw std::runtime_error("Nodal internal force: Invalid phase");

    // Decide to update or assign
    double factor = 1.0;
    if (!update) factor = 0.;

    // Update/assign internal force
    std::lock_guard<std::mutex> guard(node_mutex_);
    Eigen::Matrix<double, Tdim, Tnphases> internal_force_subdomain =
        this->internal_force_subdomain_.at(mid);
    internal_force_subdomain.col(phase) =
        internal_force_subdomain.col(phase) * factor + force;
    this->internal_force_subdomain_.at(mid) = internal_force_subdomain;
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

//! Assign nodal momentum of a subdomain
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::update_momentum_subdomain(
    bool update, unsigned phase, const Eigen::Matrix<double, Tdim, 1>& momentum,
    unsigned mid) {
  bool status = false;
  try {
    if (phase >= Tnphases)
      throw std::runtime_error("Nodal subdomain momentum: Invalid phase");

    // Decide to update or assign
    double factor = 1.0;
    if (!update) factor = 0.;

    // Update/assign subdomain momentum
    std::lock_guard<std::mutex> guard(node_mutex_);
    Eigen::Matrix<double, Tdim, Tnphases> momentum_subdomain =
        this->momentum_subdomain_.at(mid);
    momentum_subdomain.col(phase) =
        momentum_subdomain.col(phase) * factor + momentum;
    this->momentum_subdomain_.at(mid) = momentum_subdomain;
    status = true;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assign nodal coordinates (mapped from particles)
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::update_coordinates(
    bool update, double pmass,
    const Eigen::Matrix<double, Tdim, 1>& shapefn_coordinates) {
  bool status = false;
  try {
    const double tolerance = 1.E-16;

    // Decide to update or assign
    double factor = 1.0;
    if (!update) factor = 0.;

    // Update/assign coordinates
    std::lock_guard<std::mutex> guard(node_mutex_);
    if (mass_(0) > tolerance) {
      coordinates_from_particles_ = coordinates_from_particles_ * factor +
                                    shapefn_coordinates * pmass / mass_(0);
      status = true;
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assign nodal coordinates of a subdomain (mapped from particles)
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::update_coordinates_subdomain(
    bool update, double pmass,
    const Eigen::Matrix<double, Tdim, 1>& shapefn_coordinates, unsigned mid) {
  bool status = false;
  try {
    const double tolerance = 1.E-16;

    // Decide to update or assign
    double factor = 1.0;
    if (!update) factor = 0.;

    // Update/assign subdomain momentum
    std::lock_guard<std::mutex> guard(node_mutex_);
    if (mass_subdomain_.at(mid)(0) > tolerance) {
      coordinates_from_particles_subdomain_.at(mid) =
          coordinates_from_particles_subdomain_.at(mid) * factor +
          shapefn_coordinates * pmass / mass_subdomain_.at(mid)(0);
    }
    status = true;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Update pressure at the nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_pressure(bool update,
                                                      unsigned phase,
                                                      double mass_pressure) {
  try {
    const double tolerance = 1.E-16;

    // Decide to update or assign
    double factor = 1.0;
    if (!update) factor = 0.;

    // Update/assign pressure
    if (mass_(phase) > tolerance) {
      std::lock_guard<std::mutex> guard(node_mutex_);
      pressure_(phase) =
          (pressure_(phase) * factor) + mass_pressure / mass_(phase);
    } else
      throw std::runtime_error("Nodal mass is zero or below threshold");
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

//! Compute velocity of a subdomain from momentum
//! velocity = momentum / mass
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::compute_velocity_subdomain(unsigned mid) {
  try {
    const double tolerance = 1.E-16;  // std::numeric_limits<double>::lowest();

    for (unsigned phase = 0; phase < Tnphases; ++phase) {
      if (mass_subdomain_.at(mid)(phase) > tolerance) {
        Eigen::Matrix<double, Tdim, Tnphases> velocity_subdomain =
            this->velocity_subdomain_.at(mid);
        velocity_subdomain.col(phase) = momentum_subdomain_.at(mid).col(phase) /
                                        mass_subdomain_.at(mid)(phase);
        this->velocity_subdomain_.at(mid) = velocity_subdomain;

        // Check to see if value is below threshold
        for (unsigned i = 0; i < velocity_subdomain.rows(); ++i)
          if (std::abs(velocity_subdomain.col(phase)(i)) < 1.E-15) {
            velocity_subdomain.col(phase)(i) = 0.;
            this->velocity_subdomain_.at(mid) = velocity_subdomain;
          }
      }  // else
         // throw std::runtime_error(
         //"Nodal subdomain mass is zero or below threshold");
    }
    // Apply velocity constraints, which also sets acceleration to 0,
    // when velocity is set.
    this->apply_velocity_constraints_subdomain(mid);
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

//! Compute acceleration and velocity of a subdomain
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::compute_acceleration_velocity_subdomain(
    unsigned phase, double dt, unsigned mid) {
  bool status = true;
  const double tolerance = 1.0E-15;
  try {
    // Get the last acceleration value
    Eigen::Matrix<double, Tdim, Tnphases> acceleration_subdomain =
        this->acceleration_subdomain_.at(mid);
    if (mass_subdomain_.at(mid)(phase) > tolerance) {
      // Compute acceleration (unbalaced force / mass)
      acceleration_subdomain.col(phase) =
          (this->external_force_subdomain_.at(mid).col(phase) +
           this->internal_force_subdomain_.at(mid).col(phase)) /
          this->mass_subdomain_.at(mid)(phase);
      // Update acceleration
      this->acceleration_subdomain_.at(mid) = acceleration_subdomain;

      // Apply friction constraints
      this->apply_friction_constraints_subdomain(dt, mid);

      // Get the last velocity value
      Eigen::Matrix<double, Tdim, Tnphases> velocity_subdomain =
          this->velocity_subdomain_.at(mid);
      // Compute velocity: Velocity += acceleration * dt
      velocity_subdomain.col(phase) =
          velocity_subdomain.col(phase) +
          this->acceleration_subdomain_.at(mid).col(phase) * dt;
      // Update velocity
      this->velocity_subdomain_.at(mid) = velocity_subdomain;
      // Apply velocity constraints, which also sets acceleration to 0,
      // when velocity is set.
      this->apply_velocity_constraints_subdomain(mid);

      // Set a threshold
      for (unsigned i = 0; i < Tdim; ++i) {
        if (std::abs(velocity_subdomain.col(phase)(i)) < tolerance) {
          velocity_subdomain.col(phase)(i) = 0.;
          this->velocity_subdomain_.at(mid) = velocity_subdomain;
        }
        if (std::abs(acceleration_subdomain.col(phase)(i)) < tolerance) {
          acceleration_subdomain.col(phase)(i) = 0.;
          this->acceleration_subdomain_.at(mid) = acceleration_subdomain;
        }
      }
    }  // else
       // throw std::runtime_error(
       //"Nodal subdomain mass is zero or below threshold");
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
    if (dir >= 0 && dir < (Tdim * Tnphases))
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

//! Apply velocity constraints of a subdomain
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::apply_velocity_constraints_subdomain(
    unsigned mid) {
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
      Eigen::Matrix<double, Tdim, Tnphases> velocity_subdomain =
          this->velocity_subdomain_.at(mid);
      velocity_subdomain(direction, phase) = constraint.second;
      this->velocity_subdomain_.at(mid) = velocity_subdomain;
      // Set acceleration to 0 in direction of velocity constraint
      Eigen::Matrix<double, Tdim, Tnphases> acceleration_subdomain =
          this->acceleration_subdomain_.at(mid);
      acceleration_subdomain(direction, phase) = 0.;
      this->acceleration_subdomain_.at(mid) = acceleration_subdomain;
    } else {
      // Velocity constraints on general boundaries
      // Compute inverse rotation matrix
      const Eigen::Matrix<double, Tdim, Tdim> inverse_rotation_matrix =
          rotation_matrix_.inverse();
      // Transform to local coordinate
      Eigen::Matrix<double, Tdim, Tnphases> local_velocity =
          inverse_rotation_matrix * this->velocity_subdomain_.at(mid);
      Eigen::Matrix<double, Tdim, Tnphases> local_acceleration =
          inverse_rotation_matrix * this->acceleration_subdomain_.at(mid);
      // Apply boundary condition in local coordinate
      local_velocity(direction, phase) = constraint.second;
      local_acceleration(direction, phase) = 0.;
      // Transform back to global coordinate
      this->velocity_subdomain_.at(mid) = rotation_matrix_ * local_velocity;
      this->acceleration_subdomain_.at(mid) =
          rotation_matrix_ * local_acceleration;
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
    if (dir >= 0 && dir < Tdim) {
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

//! Apply subdomain friction constraints
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::apply_friction_constraints_subdomain(
    double dt, unsigned mid) {
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
        acc_n = this->acceleration_subdomain_.at(mid)(dir_n, phase);
        acc_t = this->acceleration_subdomain_.at(mid)(dir_t, phase);
        // Velocity tangential
        vel_t = this->velocity_subdomain_.at(mid)(dir_t, phase);
      } else {
        // General case, transform to local coordinate
        // Compute inverse rotation matrix
        const Eigen::Matrix<double, Tdim, Tdim> inverse_rotation_matrix =
            rotation_matrix_.inverse();
        // Transform to local coordinate
        Eigen::Matrix<double, Tdim, Tnphases> local_acceleration =
            inverse_rotation_matrix * this->acceleration_subdomain_.at(mid);
        Eigen::Matrix<double, Tdim, Tnphases> local_velocity =
            inverse_rotation_matrix * this->velocity_subdomain_.at(mid);
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
          this->acceleration_subdomain_.at(mid)(dir_t, phase) = acc_t;
        } else {
          // Local acceleration in terms of tangential and normal
          Eigen::Matrix<double, Tdim, Tnphases> acc;
          acc(dir_t, phase) = acc_t;
          acc(dir_n, phase) = acc_n;

          // General case, transform to global coordinate
          this->acceleration_subdomain_.at(mid).col(phase) =
              rotation_matrix_ * acc.col(phase);
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
        acc = this->acceleration_subdomain_.at(mid).col(phase);
        vel = this->velocity_subdomain_.at(mid).col(phase);
      } else {
        // General case, transform to local coordinate
        // Compute inverse rotation matrix
        const Eigen::Matrix<double, Tdim, Tdim> inverse_rotation_matrix =
            rotation_matrix_.inverse();
        // Transform to local coordinate
        acc = inverse_rotation_matrix *
              this->acceleration_subdomain_.at(mid).col(phase);
        vel = inverse_rotation_matrix *
              this->velocity_subdomain_.at(mid).col(phase);
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
          this->acceleration_subdomain_.at(mid).col(phase) = acc;
        } else {
          // General case, transform to global coordinate
          this->acceleration_subdomain_.at(mid).col(phase) =
              rotation_matrix_ * acc;
        }
      }
    }
  }
}

//! Compute normal vector of contact interface of a subdomain
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::compute_normal_vector(
    const std::string normal_vector_type,
    const Eigen::Matrix<double, Tdim, 1> normal_vector, unsigned mid) {
  const double tolerance = 1.0E-15;
  // Compute the normal vector by "Special Normal (SN)" method
  if (normal_vector_type == "SN") {
    // Compute the length of the normal vector
    double normal_length = 0;
    if (Tdim == 2)
      normal_length = sqrt(normal_vector(0) * normal_vector(0) +
                           normal_vector(1) * normal_vector(1));
    if (Tdim == 3)
      normal_length = sqrt(normal_vector(0) * normal_vector(0) +
                           normal_vector(1) * normal_vector(1) +
                           normal_vector(2) * normal_vector(2));
    // Normalise the normal vector
    if (normal_length > tolerance)
      this->normal_vector_.at(mid) = normal_vector / normal_length;
  }
  // Compute the normal vector by "Maximum Volume Gradient (MVG)" method
  // Compute the normal vector by "Average Volume Gradient (AVG)" method
}

//! Compute contact components
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::compute_contact_components(
    const Eigen::Matrix<double, Tdim, 1> separation_vector,
    double* separation_normal_length, double* separation_tangent_length,
    unsigned mid) {
  const double tolerance = 1.0E-15;
  // Compute the length of the separation vector in normal directions
  (*separation_normal_length) =
      separation_vector.dot(this->normal_vector_.at(mid));
  // Compute the separation vector in normal direction
  Eigen::Matrix<double, Tdim, 1> separation_normal =
      (*separation_normal_length) * this->normal_vector_.at(mid);
  // Compute the separation vector in tangent direction
  Eigen::Matrix<double, Tdim, 1> separation_tangent =
      separation_vector - separation_normal;
  // Compute the length of the separation vector in tangent directions
  if (Tdim == 2)
    (*separation_tangent_length) =
        sqrt(separation_tangent(0) * separation_tangent(0) +
             separation_tangent(1) * separation_tangent(1));
  if (Tdim == 3)
    (*separation_tangent_length) =
        sqrt(separation_tangent(0) * separation_tangent(0) +
             separation_tangent(1) * separation_tangent(1) +
             separation_tangent(2) * separation_tangent(2));
  // Compute the normalized tangent vector
  if ((*separation_tangent_length) > tolerance)
    this->tangent_vector_.at(mid) =
        separation_tangent / (*separation_tangent_length);
}

//! Compute correction momentum
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::compute_correction_momentum(
    const std::string friction_type, const double friction_coefficient,
    const Eigen::Matrix<double, Tdim, 1> momentum_difference,
    const double separation_normal_length, unsigned mid) {
  const double tolerance = 1.0E-15;
  // Initialise correction momentum
  Eigen::Matrix<double, Tdim, 1> momentum_correction = momentum_difference;
  // Check the friction limitation
  if (friction_type == "friction") {
    // Compute normal componenet of momentum_difference
    double momentum_correction_normal =
        momentum_difference.dot(this->normal_vector_.at(mid));
    // Compute tangent componenet of momentum_difference
    double momentum_correction_tangent =
        momentum_difference.dot(this->tangent_vector_.at(mid));
    // Compute (tangent componenet of momentum)/(normal componenet of momentum)
    double mt_mn =
        fabs((momentum_correction_tangent) / (momentum_correction_normal));
    // Check sliding
    if (friction_coefficient < mt_mn) {
      // Modify tangent componenet of momentum_correction
      momentum_correction_tangent =
          friction_coefficient * momentum_correction_normal;
      // Update momentum_correction by friction
      momentum_correction =
          momentum_correction_normal * this->normal_vector_.at(mid) +
          momentum_correction_tangent * this->tangent_vector_.at(mid);
    }
  }
  // Correct the subdomain momentume
  Eigen::Matrix<double, Tdim, Tnphases> momentum_subdomain =
      this->momentum_subdomain_.at(mid);
  momentum_subdomain.col(0) = momentum_subdomain.col(0) + momentum_correction;
  this->momentum_subdomain_.at(mid) = momentum_subdomain;
}

//! Compute contact force
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::compute_contact_force(
    const double separation_normal_length,
    const double separation_tangent_length, const double separation_cut_off,
    const double dc_n, const double dc_t, unsigned mid) {
  Eigen::Matrix<double, Tdim, Tdim> dn_, dt_;
  dn_.setZero();
  dt_.setZero();
  dn_(0, 0) = dc_n;
  dn_(1, 1) = dc_n;
  dt_(0, 0) = dc_t;
  dt_(1, 1) = dc_t;
  if (Tdim == 3) {
    dn_(2, 2) = dc_n;
    dt_(2, 2) = dc_t;
  }
  double ai = 100;
  Eigen::Matrix<double, Tdim, Tnphases> internal_force_subdomain =
      this->internal_force_subdomain_.at(mid);
  internal_force_subdomain.col(0) -=
      (dn_ * (fabs(separation_cut_off - separation_normal_length) *
              this->normal_vector_.at(mid)) +
       dt_ * (separation_tangent_length * this->tangent_vector_.at(mid))) *
      ai;
  this->internal_force_subdomain_.at(mid) = internal_force_subdomain;
}

//! Compute contact interface on node
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::compute_contact_interface(
    bool contact_force, const std::string friction_type,
    const double friction_coefficient, const double separation_cut_off,
    const double dc_n, const double dc_t, unsigned mid) {
  const double tolerance = 1.0E-12;
  // Check if there is any subdomain particle
  if (mass_subdomain_.at(mid)(0) > tolerance &&
      (mass_(0) - mass_subdomain_.at(mid)(0)) > tolerance) {
    // Compute the difference between total momentum and subdomain momentum
    Eigen::Matrix<double, Tdim, 1> momentum_difference =
        velocity_.col(0) * mass_subdomain_.at(mid)(0) -
        momentum_subdomain_.at(mid).col(0);
    // Compute the separation vector
    Eigen::Matrix<double, Tdim, 1> separation_vector =
        mass_(0) / (mass_(0) - mass_subdomain_.at(mid)(0)) *
        (coordinates_from_particles_ -
         coordinates_from_particles_subdomain_.at(mid));
    // Detect contact
    if (momentum_difference.dot(this->normal_vector_.at(mid)) < -1.0E-15 &&
        separation_vector.dot(this->normal_vector_.at(mid)) <
            separation_cut_off) {
      // Normal length of separation vector
      double separation_normal_length = 0;
      // Tangent length of separation vector
      double separation_tangent_length = 0;
      // Compute contact components
      this->compute_contact_components(separation_vector,
                                       &separation_normal_length,
                                       &separation_tangent_length, mid);
      //! Contact_force FALSE: Modify subdomain momentume
      if (!contact_force)
        this->compute_correction_momentum(friction_type, friction_coefficient,
                                          momentum_difference,
                                          separation_normal_length, mid);
      //! Contact_force TRUE: Implement contact force
      if (contact_force)
        this->compute_contact_force(separation_normal_length,
                                    separation_tangent_length,
                                    separation_cut_off, dc_n, dc_t, mid);
    }
  }
}

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
  concentrated_force_.setZero();

  // Initialize scalar properties
  scalar_properties_.emplace(
      std::make_pair(mpm::properties::Scalar::Mass,
                     Eigen::Matrix<double, 1, Tnphases>::Zero()));
  scalar_properties_.emplace(
      std::make_pair(mpm::properties::Scalar::Volume,
                     Eigen::Matrix<double, 1, Tnphases>::Zero()));
  scalar_properties_.emplace(
      std::make_pair(mpm::properties::Scalar::Pressure,
                     Eigen::Matrix<double, 1, Tnphases>::Zero()));

  // Initialize vector properties
  vector_properties_.emplace(
      std::make_pair(mpm::properties::Vector::Velocity,
                     Eigen::Matrix<double, Tdim, Tnphases>::Zero()));
  vector_properties_.emplace(
      std::make_pair(mpm::properties::Vector::Acceleration,
                     Eigen::Matrix<double, Tdim, Tnphases>::Zero()));
  vector_properties_.emplace(
      std::make_pair(mpm::properties::Vector::Momentum,
                     Eigen::Matrix<double, Tdim, Tnphases>::Zero()));
  vector_properties_.emplace(
      std::make_pair(mpm::properties::Vector::ExternalForce,
                     Eigen::Matrix<double, Tdim, Tnphases>::Zero()));
  vector_properties_.emplace(
      std::make_pair(mpm::properties::Vector::InternalForce,
                     Eigen::Matrix<double, Tdim, Tnphases>::Zero()));

  this->initialise();
}

//! Initialise nodal properties
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::initialise() noexcept {
  status_ = false;

  // Initialise nodal scalar properties
  scalar_properties_.at(mpm::properties::Scalar::Mass).setZero();
  scalar_properties_.at(mpm::properties::Scalar::Volume).setZero();
  scalar_properties_.at(mpm::properties::Scalar::Pressure).setZero();

  // Initialise nodal vector properties
  vector_properties_.at(mpm::properties::Vector::Velocity).setZero();
  vector_properties_.at(mpm::properties::Vector::Acceleration).setZero();
  vector_properties_.at(mpm::properties::Vector::Momentum).setZero();
  vector_properties_.at(mpm::properties::Vector::ExternalForce).setZero();
  vector_properties_.at(mpm::properties::Vector::InternalForce).setZero();

  // Initialise variables for contact
  material_ids_.clear();
  contact_displacement_.setZero();
}

//! Initialise shared pointer to nodal properties pool
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::initialise_property_handle(
    unsigned prop_id,
    std::shared_ptr<mpm::NodalProperties> property_handle) noexcept {
  // the property handle and the property id is set in the node
  this->property_handle_ = property_handle;
  this->prop_id_ = prop_id;
}

//! Update scalar property at the nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_scalar_property(
    mpm::properties::Scalar property, bool update, unsigned phase,
    double value) noexcept {
  // Assert phase
  assert(phase < Tnphases);

  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign value
  std::lock_guard<std::mutex> guard(node_mutex_);
  scalar_properties_.at(property)[phase] =
      (scalar_properties_.at(property)[phase] * factor) + value;
}

//! Update scalar property at the nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
double mpm::Node<Tdim, Tdof, Tnphases>::scalar_property(
    mpm::properties::Scalar property, unsigned phase) const {
  // Assert phase
  assert(phase < Tnphases);
  return scalar_properties_.at(property)[phase];
}

//! Update vector property at the nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_vector_property(
    mpm::properties::Vector property, bool update, unsigned phase,
    const Eigen::Matrix<double, Tdim, 1>& value) noexcept {
  // Assert phase
  assert(phase < Tnphases);

  // Decide to update or assign
  const double factor = (update == true) ? 1. : 0.;

  // Update/assign value
  std::lock_guard<std::mutex> guard(node_mutex_);
  Eigen::Matrix<double, Tdim, 1> vecvalue =
      vector_properties_.at(property).col(phase);
  vector_properties_.at(property).col(phase) = (vecvalue * factor) + value;
}

//! Update vector property at the nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
Eigen::Matrix<double, Tdim, 1> mpm::Node<Tdim, Tdof, Tnphases>::vector_property(
    mpm::properties::Vector property, unsigned phase) const {
  // Assert phase
  assert(phase < Tnphases);
  return vector_properties_.at(property).col(phase);
}

//! Update mass at the nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_mass(bool update, unsigned phase,
                                                  double mass) noexcept {
  this->update_scalar_property(mpm::properties::Scalar::Mass, update, phase,
                               mass);
}

//! Update volume at the nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_volume(bool update, unsigned phase,
                                                    double volume) noexcept {
  this->update_scalar_property(mpm::properties::Scalar::Volume, update, phase,
                               volume);
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
  const double scalar =
      (force_function_ != nullptr) ? force_function_->value(current_time) : 1.0;
  this->update_external_force(true, phase,
                              scalar * concentrated_force_.col(phase));
}

//! Update external force (body force / traction force)
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_external_force(
    bool update, unsigned phase,
    const Eigen::Matrix<double, Tdim, 1>& force) noexcept {
  this->update_vector_property(mpm::properties::Vector::ExternalForce, update,
                               phase, force);
}

//! Update internal force (body force / traction force)
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_internal_force(
    bool update, unsigned phase,
    const Eigen::Matrix<double, Tdim, 1>& force) noexcept {
  this->update_vector_property(mpm::properties::Vector::InternalForce, update,
                               phase, force);
}

//! Assign nodal momentum
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_momentum(
    bool update, unsigned phase,
    const Eigen::Matrix<double, Tdim, 1>& momentum) noexcept {
  this->update_vector_property(mpm::properties::Vector::Momentum, update, phase,
                               momentum);
}

//! Update pressure at the nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_mass_pressure(
    unsigned phase, double mass_pressure) noexcept {
  // Assert
  assert(phase < Tnphases);

  const double tolerance = 1.E-16;
  // Compute pressure from mass*pressure
  if (this->mass(phase) > tolerance) {
    std::lock_guard<std::mutex> guard(node_mutex_);
    scalar_properties_.at(mpm::properties::Scalar::Pressure)(phase) +=
        mass_pressure / this->mass(phase);
  }
}

//! Assign pressure at the nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::assign_pressure(unsigned phase,
                                                      double pressure) {
  // Compute pressure from mass*pressure
  std::lock_guard<std::mutex> guard(node_mutex_);
  scalar_properties_.at(mpm::properties::Scalar::Pressure)(phase) = pressure;
}

//! Compute velocity from momentum
//! velocity = momentum / mass
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::compute_velocity() {
  const double tolerance = 1.E-16;
  for (unsigned phase = 0; phase < Tnphases; ++phase) {
    if (this->mass(phase) > tolerance) {
      vector_properties_.at(mpm::properties::Vector::Velocity).col(phase) =
          this->momentum(phase) / this->mass(phase);

      // Check to see if value is below threshold
      for (unsigned i = 0; i < Tdim; ++i)
        if (std::abs(this->velocity(phase)(i)) < 1.E-15)
          vector_properties_.at(mpm::properties::Vector::Velocity)
              .col(phase)(i) = 0.;
    }
  }

  // Apply velocity constraints, which also sets acceleration to 0,
  // when velocity is set.
  this->apply_velocity_constraints();
}

//! Update nodal acceleration
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_acceleration(
    bool update, unsigned phase,
    const Eigen::Matrix<double, Tdim, 1>& acceleration) noexcept {
  this->update_vector_property(mpm::properties::Vector::Acceleration, update,
                               phase, acceleration);
}

//! Compute acceleration and velocity
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::compute_acceleration_velocity(
    unsigned phase, double dt) noexcept {
  bool status = false;
  const double tolerance = 1.0E-15;
  if (this->mass(phase) > tolerance) {
    // acceleration = (unbalaced force / mass)
    vector_properties_.at(mpm::properties::Vector::Acceleration).col(phase) =
        (this->external_force(phase) + this->internal_force(phase)) /
        this->mass(phase);

    // Apply friction constraints
    this->apply_friction_constraints(dt);

    // Velocity += acceleration * dt
    vector_properties_.at(mpm::properties::Vector::Velocity).col(phase) +=
        this->acceleration(phase) * dt;
    // Apply velocity constraints, which also sets acceleration to 0,
    // when velocity is set.
    this->apply_velocity_constraints();

    // Set a threshold
    for (unsigned i = 0; i < Tdim; ++i)
      if (std::abs(this->velocity(phase)(i)) < tolerance)
        vector_properties_.at(mpm::properties::Vector::Velocity).col(phase)(i) =
            0.;
    for (unsigned i = 0; i < Tdim; ++i)
      if (std::abs(this->acceleration(phase)(i)) < tolerance)
        vector_properties_.at(mpm::properties::Vector::Acceleration)
            .col(phase)(i) = 0.;
    status = true;
  }
  return status;
}

//! Compute acceleration and velocity with cundall damping factor
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::compute_acceleration_velocity_cundall(
    unsigned phase, double dt, double damping_factor) noexcept {
  bool status = false;
  const double tolerance = 1.0E-15;
  if (this->mass(phase) > tolerance) {
    // acceleration = (unbalaced force / mass)
    auto unbalanced_force =
        this->external_force(phase) + this->internal_force(phase);
    vector_properties_.at(mpm::properties::Vector::Acceleration).col(phase) =
        (unbalanced_force - damping_factor * unbalanced_force.norm() *
                                this->velocity(phase).cwiseSign()) /
        this->mass(phase);

    // Apply friction constraints
    this->apply_friction_constraints(dt);

    // Velocity += acceleration * dt
    vector_properties_.at(mpm::properties::Vector::Velocity).col(phase) +=
        this->acceleration(phase) * dt;
    // Apply velocity constraints, which also sets acceleration to 0,
    // when velocity is set.
    this->apply_velocity_constraints();

    // Set a threshold
    for (unsigned i = 0; i < Tdim; ++i)
      if (std::abs(this->velocity(phase)(i)) < tolerance)
        vector_properties_.at(mpm::properties::Vector::Velocity).col(phase)(i) =
            0.;
    for (unsigned i = 0; i < Tdim; ++i)
      if (std::abs(this->acceleration(phase)(i)) < tolerance)
        vector_properties_.at(mpm::properties::Vector::Acceleration)
            .col(phase)(i) = 0.;
    status = true;
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
      vector_properties_.at(mpm::properties::Vector::Velocity)(
          direction, phase) = constraint.second;
      // Set acceleration to 0 in direction of velocity constraint
      vector_properties_.at(mpm::properties::Vector::Acceleration)(direction,
                                                                   phase) = 0.;
    } else {
      // Velocity constraints on general boundaries
      // Compute inverse rotation matrix
      const Eigen::Matrix<double, Tdim, Tdim> inverse_rotation_matrix =
          rotation_matrix_.inverse();
      // Transform to local coordinate
      Eigen::Matrix<double, Tdim, Tnphases> local_velocity =
          inverse_rotation_matrix *
          vector_properties_.at(mpm::properties::Vector::Velocity);
      Eigen::Matrix<double, Tdim, Tnphases> local_acceleration =
          inverse_rotation_matrix *
          vector_properties_.at(mpm::properties::Vector::Acceleration);
      // Apply boundary condition in local coordinate
      local_velocity(direction, phase) = constraint.second;
      local_acceleration(direction, phase) = 0.;
      // Transform back to global coordinate
      vector_properties_.at(mpm::properties::Vector::Velocity) =
          rotation_matrix_ * local_velocity;
      vector_properties_.at(mpm::properties::Vector::Acceleration) =
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
        acc_n = vector_properties_.at(mpm::properties::Vector::Acceleration)(
            dir_n, phase);
        acc_t = vector_properties_.at(mpm::properties::Vector::Acceleration)(
            dir_t, phase);
        // Velocity tangential
        vel_t = vector_properties_.at(mpm::properties::Vector::Velocity)(dir_t,
                                                                         phase);
      } else {
        // General case, transform to local coordinate
        // Compute inverse rotation matrix
        const Eigen::Matrix<double, Tdim, Tdim> inverse_rotation_matrix =
            rotation_matrix_.inverse();
        // Transform to local coordinate
        Eigen::Matrix<double, Tdim, Tnphases> local_acceleration =
            inverse_rotation_matrix *
            vector_properties_.at(mpm::properties::Vector::Acceleration);
        Eigen::Matrix<double, Tdim, Tnphases> local_velocity =
            inverse_rotation_matrix *
            vector_properties_.at(mpm::properties::Vector::Velocity);
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
          vector_properties_.at(mpm::properties::Vector::Acceleration)(
              dir_t, phase) = acc_t;
        } else {
          // Local acceleration in terms of tangential and normal
          Eigen::Matrix<double, Tdim, Tnphases> acc;
          acc(dir_t, phase) = acc_t;
          acc(dir_n, phase) = acc_n;

          // General case, transform to global coordinate
          vector_properties_.at(mpm::properties::Vector::Acceleration)
              .col(phase) = rotation_matrix_ * acc.col(phase);
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
        acc = this->acceleration(phase);
        vel = this->velocity(phase);
      } else {
        // General case, transform to local coordinate
        // Compute inverse rotation matrix
        const Eigen::Matrix<double, Tdim, Tdim> inverse_rotation_matrix =
            rotation_matrix_.inverse();
        // Transform to local coordinate
        acc = inverse_rotation_matrix * this->acceleration(phase);
        vel = inverse_rotation_matrix * this->velocity(phase);
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
            acc(dir_t0) = -vel(dir_t0) / dt;
            acc(dir_t1) = -vel(dir_t1) / dt;
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
          vector_properties_.at(mpm::properties::Vector::Acceleration)
              .col(phase) = acc;
        } else {
          // General case, transform to global coordinate
          vector_properties_.at(mpm::properties::Vector::Acceleration)
              .col(phase) = rotation_matrix_ * acc;
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

//! Update nodal property at the nodes from particle
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_property(
    bool update, const std::string& property,
    const Eigen::MatrixXd& property_value, unsigned mat_id,
    unsigned nprops) noexcept {
  // Update/assign property
  std::lock_guard<std::mutex> guard(node_mutex_);
  property_handle_->update_property(property, prop_id_, mat_id, property_value,
                                    nprops);
}

//! Compute multimaterial change in momentum
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof,
               Tnphases>::compute_multimaterial_change_in_momentum() {
  // iterate over all materials in the material_ids set and update the change in
  // momentum
  std::lock_guard<std::mutex> guard(node_mutex_);
  for (auto mitr = material_ids_.begin(); mitr != material_ids_.end(); ++mitr) {
    const Eigen::Matrix<double, 1, 1> mass =
        property_handle_->property("masses", prop_id_, *mitr);
    const Eigen::Matrix<double, Tdim, 1> momentum =
        property_handle_->property("momenta", prop_id_, *mitr, Tdim);
    const Eigen::Matrix<double, Tdim, 1> change_in_momenta =
        vector_properties_.at(mpm::properties::Vector::Velocity) * mass -
        momentum;
    property_handle_->update_property("change_in_momenta", prop_id_, *mitr,
                                      change_in_momenta, Tdim);
  }
}

//! Compute multimaterial separation vector
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof,
               Tnphases>::compute_multimaterial_separation_vector() {
  // iterate over all materials in the material_ids set, update the
  // displacements and calculate the displacement of the center of mass for this
  // node
  std::lock_guard<std::mutex> guard(node_mutex_);
  for (auto mitr = material_ids_.begin(); mitr != material_ids_.end(); ++mitr) {
    const auto& material_displacement =
        property_handle_->property("displacements", prop_id_, *mitr, Tdim);
    const auto& material_mass =
        property_handle_->property("masses", prop_id_, *mitr);

    // displacement of the center of mass
    contact_displacement_ +=
        material_displacement /
        scalar_properties_.at(mpm::properties::Scalar::Mass)(0, 0);
    // assign nodal-multimaterial displacement by dividing it by this material's
    // mass
    property_handle_->assign_property(
        "displacements", prop_id_, *mitr,
        material_displacement / material_mass(0, 0), Tdim);
  }

  // iterate over all materials in the material_ids to compute the separation
  // vector
  for (auto mitr = material_ids_.begin(); mitr != material_ids_.end(); ++mitr) {
    const Eigen::Matrix<double, Tdim, 1> material_displacement =
        property_handle_->property("displacements", prop_id_, *mitr, Tdim);
    const Eigen::Matrix<double, 1, 1> material_mass =
        property_handle_->property("masses", prop_id_, *mitr);

    // Update the separation vector property
    const auto& separation_vector =
        (contact_displacement_ - material_displacement) *
        scalar_properties_.at(mpm::properties::Scalar::Mass)(0, 0) /
        (scalar_properties_.at(mpm::properties::Scalar::Mass)(0, 0) -
         material_mass(0, 0));
    property_handle_->update_property("separation_vectors", prop_id_, *mitr,
                                      separation_vector, Tdim);
  }
}

//! Compute multimaterial normal unit vector
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof,
               Tnphases>::compute_multimaterial_normal_unit_vector() {
  // Iterate over all materials in the material_ids set
  std::lock_guard<std::mutex> guard(node_mutex_);
  for (auto mitr = material_ids_.begin(); mitr != material_ids_.end(); ++mitr) {
    // calculte the normal unit vector
    VectorDim domain_gradient =
        property_handle_->property("domain_gradients", prop_id_, *mitr, Tdim);
    VectorDim normal_unit_vector = VectorDim::Zero();
    if (domain_gradient.norm() > std::numeric_limits<double>::epsilon())
      normal_unit_vector = domain_gradient.normalized();

    // assign nodal-multimaterial normal unit vector to property pool
    property_handle_->assign_property("normal_unit_vectors", prop_id_, *mitr,
                                      normal_unit_vector, Tdim);
  }
}

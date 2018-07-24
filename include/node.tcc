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

  this->initialise();
}

//! Initialise nodal properties
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::initialise() {
  mass_.setZero();
  volume_.setZero();
  external_force_.setZero();
  internal_force_.setZero();
  velocity_.setZero();
  momentum_.setZero();
  acceleration_.setZero();
  status_ = false;
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

//! Update external force (body force / traction force)
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::update_external_force(
    bool update, unsigned phase, const Eigen::VectorXd& force) {
  bool status = false;
  try {
    if (force.size() != external_force_.size()) {
      std::cout << external_force_.size() << "\t" << force.size() << "\n";
      throw std::runtime_error("Nodal force degrees of freedom don't match");
    }

    // Decide to update or assign
    double factor = 1.0;
    if (!update) factor = 0.;

    // Update/assign external force
    std::lock_guard<std::mutex> guard(node_mutex_);
    external_force_.col(phase) = external_force_.col(phase) * factor + force;
    status = true;
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
    status = false;
  }
  return status;
}

//! Update internal force (body force / traction force)
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::update_internal_force(
    bool update, unsigned phase, const Eigen::VectorXd& force) {
  bool status = false;
  try {
    if (force.size() != internal_force_.size()) {
      std::cout << internal_force_.size() << "\t" << force.size() << "\n";
      throw std::runtime_error("Nodal force degrees of freedom don't match");
    }

    // Decide to update or assign
    double factor = 1.0;
    if (!update) factor = 0.;

    // Update/assign internal force
    std::lock_guard<std::mutex> guard(node_mutex_);
    internal_force_.col(phase) = internal_force_.col(phase) * factor + force;
    status = true;
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
    status = false;
  }
  return status;
}

//! Assign nodal momentum
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::update_momentum(
    bool update, unsigned phase, const Eigen::VectorXd& momentum) {
  bool status = false;
  try {
    if (momentum.size() != momentum_.size()) {
      throw std::runtime_error("Nodal momentum degrees of freedom don't match");
    }

    // Decide to update or assign
    double factor = 1.0;
    if (!update) factor = 0.;

    // Update/assign momentum
    std::lock_guard<std::mutex> guard(node_mutex_);
    momentum_.col(phase) = momentum_.col(phase) * factor + momentum;
    status = true;
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
    status = false;
  }
  return status;
}

//! Compute velocity from momentum
//! velocity = momentum / mass
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::compute_velocity() {
  const double tolerance = std::numeric_limits<double>::lowest();
  for (unsigned phase = 0; phase < Tnphases; ++phase) {
    if (mass_(phase) > tolerance)
      velocity_.col(phase) = momentum_.col(phase) / mass_(phase);
  }
}

//! Update nodal acceleration
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::update_acceleration(
    bool update, unsigned phase, const Eigen::VectorXd& acceleration) {
  bool status = false;
  try {
    if (acceleration.size() != acceleration_.size()) {
      throw std::runtime_error(
          "Nodal acceleration degrees of freedom don't match");
    }

    // Decide to update or assign
    double factor = 1.0;
    if (!update) factor = 0.;

    //! Update/assign acceleration
    std::lock_guard<std::mutex> guard(node_mutex_);
    acceleration_.col(phase) = acceleration_.col(phase) * factor + acceleration;
    status = true;
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
    status = false;
  }
  return status;
}

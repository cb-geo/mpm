//! Constructor with id, coordinates and dof
//! \param[in] id Node id
//! \param[in] coord coordinates of the node
//! \param[in] dof Degrees of freedom
//! \tparam Tdim Dimension
//! \tparam Tdof Degrees of Freedom
//! \tparam Tnphases Number of phases
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
mpm::Node<Tdim, Tdof, Tnphases>::Node(
    Index id, const Eigen::Matrix<double, Tdim, 1>& coord)
    : NodeBase<Tdim>() {
  // Check if the dimension is between 1 & 3
  static_assert((Tdim >= 1 && Tdim <= 3), "Invalid global dimension");
  id_ = id;
  coordinates_ = coord;

  dof_ = Tdof;
  mass_.setZero();
  external_force_.setZero();
  internal_force_.setZero();
  velocity_.setZero();
  momentum_.setZero();
  acceleration_.setZero();
}

//! Initialise nodal properties
//! \tparam Tdim Dimension
//! \tparam Tdof Degrees of Freedom
//! \tparam Tnphases Number of phases
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::initialise() {
  external_force_.setZero();
  internal_force_.setZero();
  velocity_.setZero();
  momentum_.setZero();
  acceleration_.setZero();
}

//! Update mass
//! \tparam Tdim Dimension
//! \tparam Tdof Degrees of Freedom
//! \tparam Tnphases Number of phases
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_mass(
  bool update, unsigned nphase, double mass) {
    // Decide to update or assign
    double factor = 1.0;
    if (!update) factor = 0.;

    mass_(0, nphase) = mass_(0, nphase) * factor + mass; 
}

//! Update external force (body force / traction force)
//! \tparam Tdim Dimension
//! \tparam Tdof Degrees of Freedom
//! \tparam Tnphases Number of phases
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_external_force(
  bool update, unsigned nphase, const Eigen::VectorXd& force) {
  try {
    // Decide to update or assign
    double factor = 1.0;
    if (!update) factor = 0.;
    
    if (force.size() != external_force_.size()) {
      std::cout << external_force_.size() << "\t" << force.size() << "\n";
      throw std::runtime_error("Nodal force degrees of freedom don't match");
    }
    // Assign force
    external_force_.col(nphase) = external_force_.col(nphase) * factor + force;
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
  }
}

//! Update internal force (body force / traction force)
//! \tparam Tdim Dimension
//! \tparam Tdof Degrees of Freedom
//! \tparam Tnphases Number of phases
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_internal_force(
  bool update, unsigned nphase, const Eigen::VectorXd& force) {
  try {
    // Decide to update or assign
    double factor = 1.0;
    if (!update) factor = 0.;
    
    if (force.size() != internal_force_.size()) {
      std::cout << internal_force_.size() << "\t" << force.size() << "\n";
      throw std::runtime_error("Nodal force degrees of freedom don't match");
    }
    // Assign force
    internal_force_.col(nphase) = internal_force_.col(nphase) * factor + force;
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
  }
}

//! Assign nodal velocity
//! \tparam Tdim Dimension
//! \tparam Tdof Degrees of Freedom
//! \tparam Tnphases Number of phases
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::assign_velocity(
    unsigned nphase, const Eigen::VectorXd& velocity) {
  try {
    if (velocity.size() != velocity_.size()) {
      throw std::runtime_error("Nodal velocity degrees of freedom don't match");
    }
    // Assign velocity
    velocity_.col(nphase) = velocity;
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
  }
}

//! Assign nodal momentum
//! \tparam Tdim Dimension
//! \tparam Tdof Degrees of Freedom
//! \tparam Tnphases Number of phases
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::assign_momentum(
    unsigned nphase, const Eigen::VectorXd& momentum) {
  try {
    if (momentum.size() != momentum_.size()) {
      throw std::runtime_error("Nodal momentum degrees of freedom don't match");
    }
    // Assign momentum
    momentum_.col(nphase) = momentum;
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
  }
}

//! Assign nodal acceleration
//! \tparam Tdim Dimension
//! \tparam Tdof Degrees of Freedom
//! \tparam Tnphases Number of phases
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::assign_acceleration(
    unsigned nphase, const Eigen::VectorXd& acceleration) {
  try {
    if (acceleration.size() != acceleration_.size()) {
      throw std::runtime_error(
          "Nodal acceleration degrees of freedom don't match");
    }
    //! Assign acceleration
    acceleration_.col(nphase) = acceleration;
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
  }
}

//! Compute velocity from momentum
//! velocity = momentum / mass
//! \tparam Tdim Dimension
//! \tparam Tdof Degrees of Freedom
//! \tparam Tnphases Number of phases
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::compute_velocity() {
  const double tolerance = std::numeric_limits<double>::lowest();
  for (unsigned nphase = 0; nphase < Tnphases; ++nphase) {
    if (mass_(nphase) > tolerance)
      velocity_.col(nphase) = momentum_.col(nphase) / mass_(nphase);
  }
}

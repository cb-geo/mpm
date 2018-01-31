// Constructor with id, coordinates and dof
//! \param[in] id Node id
//! \param[in] coord coordinates of the node
//! \param[in] dof Degrees of freedom
//! \tparam Tdim Dimension
//! \tparam Tdof Degrees of Freedom
//! \tparam Tnphases Number of phases
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
mpm::Node<Tdim, Tdof, Tnphases>::Node(Index id, const Eigen::Matrix<double, Tdim, 1>& coord)
    : NodeBase<Tdim>(id, coord) {
  dof_ = Tdof;
  mass_.setZero();
  force_.setZero();
  velocity_.setZero();
  momentum_.setZero();
  acceleration_.setZero();
}

// Initialise nodal properties
//! \tparam Tdim Dimension
//! \tparam Tdof Degrees of Freedom
//! \tparam Tnphases Number of phases
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::initialise() {
  force_.setZero();
  velocity_.setZero();
  momentum_.setZero();
  acceleration_.setZero();
}

// Assign nodal force
//! \tparam Tdim Dimension
//! \tparam Tdof Degrees of Freedom
//! \tparam Tnphases Number of phases
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::assign_force(unsigned nphase, const Eigen::VectorXd& force) {
  try {
    if (force.size() != force_.size()) {
      std::cout << force_.size() << "\t" << force.size() << "\n";
      throw std::runtime_error("Nodal force degrees of freedom don't match");
    }
    // Assign force
    force_.col(nphase) = force;
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
  }
}

// Assign nodal velocity
//! \tparam Tdim Dimension
//! \tparam Tdof Degrees of Freedom
//! \tparam Tnphases Number of phases
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::assign_velocity(unsigned nphase, const Eigen::VectorXd& velocity) {
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

// Assign nodal momentum
//! \tparam Tdim Dimension
//! \tparam Tdof Degrees of Freedom
//! \tparam Tnphases Number of phases
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::assign_momentum(unsigned nphase, const Eigen::VectorXd& momentum) {
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

// Assign nodal acceleration
//! \tparam Tdim Dimension
//! \tparam Tdof Degrees of Freedom
//! \tparam Tnphases Number of phases
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::assign_acceleration(unsigned nphase, const Eigen::VectorXd& acceleration) {
  try {
    if (acceleration.size() != acceleration_.size()) {
      throw std::runtime_error(
          "Nodal acceleration degrees of freedom don't match");
    }
    // Assign acceleration
    acceleration_.col(nphase) = acceleration;
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
  }
}

// Constructor with id, coordinates and dof
//! \param[in] id Node id
//! \param[in] coord coordinates of the node
//! \param[in] dof Degrees of freedom
//! \tparam Tdim Dimension
template <unsigned Tdim>
mpm::Node<Tdim>::Node(Index id, const VectorDim& coord, unsigned dof)
    : NodeBase<Tdim>(id, coord), dof_{dof} {
  nphases_ = 1.;
  force_.resize(dof_);
  velocity_.resize(dof_);
  momentum_.resize(dof_);
  acceleration_.resize(dof_);
  this->initialise();
}

// Constructor with id, coordinates, dof and nphases
//! \param[in] id Node id
//! \param[in] coord coordinates of the node
//! \param[in] dof Degrees of freedom
//! \param[in] nphases Number of phases
//! \tparam Tdim Dimension
template <unsigned Tdim>
mpm::Node<Tdim>::Node(Index id, const VectorDim& coord, unsigned dof,
                      unsigned nphases)
    : NodeBase<Tdim>(id, coord), dof_{dof}, nphases_{nphases} {
  force_.resize(dof_ * nphases_);
  velocity_.resize(dof_ * nphases_);
  momentum_.resize(dof_ * nphases_);
  acceleration_.resize(dof_ * nphases_);
  this->initialise();
}

// Initialise nodal properties
//! \tparam Tdim Dimension
template <unsigned Tdim>
void mpm::Node<Tdim>::initialise() {
  force_.setZero();
  velocity_.setZero();
  momentum_.setZero();
  acceleration_.setZero();
}

// Assign nodal force
//! \tparam Tdim Dimension
template <unsigned Tdim>
void mpm::Node<Tdim>::assign_force(const Eigen::VectorXd& force) {
  try {
    if (force.size() != force_.size()) {
      throw std::runtime_error("Nodal force degrees of freedom don't match");
    }
    // Assign force
    force_ = force;
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
  }
}

// Assign nodal velocity
//! \tparam Tdim Dimension
template <unsigned Tdim>
void mpm::Node<Tdim>::assign_velocity(const Eigen::VectorXd& velocity) {
  try {
    if (velocity.size() != velocity_.size()) {
      throw std::runtime_error("Nodal velocity degrees of freedom don't match");
    }
    // Assign velocity
    velocity_ = velocity;
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
  }
}

// Assign nodal momentum
//! \tparam Tdim Dimension
template <unsigned Tdim>
void mpm::Node<Tdim>::assign_momentum(const Eigen::VectorXd& momentum) {
  try {
    if (momentum.size() != momentum_.size()) {
      throw std::runtime_error("Nodal momentum degrees of freedom don't match");
    }
    // Assign momentum
    momentum_ = momentum;
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
  }
}

// Assign nodal acceleration
//! \tparam Tdim Dimension
template <unsigned Tdim>
void mpm::Node<Tdim>::assign_acceleration(const Eigen::VectorXd& acceleration) {
  try {
    if (acceleration.size() != acceleration_.size()) {
      throw std::runtime_error(
          "Nodal acceleration degrees of freedom don't match");
    }
    // Assign acceleration
    acceleration_ = acceleration;
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
  }
}

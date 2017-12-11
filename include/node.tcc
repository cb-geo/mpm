// Constructor with id, coordinates and dof
//! \param[in] id Node id
//! \param[in] coord coordinates of the node
//! \param[in] dof Degrees of freedom
//! \tparam Tdim Dimension
template <unsigned Tdim>
mpm::Node<Tdim>::Node(Index id, const VectorDim& coord, unsigned dof)
    : NodeBase<Tdim>(id, coord), dof_{dof} {
  force_.resize(dof_);
  this->initialise();
}

// Initialise nodal properties
//! \tparam Tdim Dimension
template <unsigned Tdim>
void mpm::Node<Tdim>::initialise() {
  force_.setZero();
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

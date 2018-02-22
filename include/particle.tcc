//! Constructor with id and coordinates
//! \param[in] id Particle id
//! \param[in] coord coordinates of the particle
//! \tparam Tdim Dimension
//! \tparam Tnphases Number of phases
template <unsigned Tdim, unsigned Tnphases>
mpm::Particle<Tdim, Tnphases>::Particle(Index id, const VectorDim& coord)
    : mpm::ParticleBase<Tdim>(id, coord) {
  this->initialise();
}

//! Constructor with id, coordinates and status
//! \param[in] id Particle id
//! \param[in] coord coordinates of the particle
//! \param[in] status Particle status (active / inactive)
//! \tparam Tdim Dimension
//! \tparam Tnphases Number of phases
template <unsigned Tdim, unsigned Tnphases>
mpm::Particle<Tdim, Tnphases>::Particle(Index id, const VectorDim& coord,
                                        bool status)
    : mpm::ParticleBase<Tdim>(id, coord, status) {
  this->initialise();
}

// Assign a cell to particle
//! \param[in] cellptr Pointer to a cell
//! \tparam Tdim Dimension
//! \tparam Tnphases Number of phases
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::assign_cell(
    const std::shared_ptr<Cell<Tdim>>& cellptr) {
  cell_ = cellptr;
  cell_id_ = cellptr->id();
  return cell_->add_particle_id(this->id());
}

// Initialise particle properties
//! \tparam Tdim Dimension
//! \tparam Tnphases Number of phases
template <unsigned Tdim, unsigned Tnphases>
void mpm::Particle<Tdim, Tnphases>::initialise() {
  mass_.setZero();
  stress_.setZero();
  velocity_.setZero();
  momentum_.setZero();
  acceleration_.setZero();
}

// Assign particle stress
//! \tparam Tdim Dimension
//! \tparam Tnphases Number of phases
template <unsigned Tdim, unsigned Tnphases>
void mpm::Particle<Tdim, Tnphases>::assign_stress(
    unsigned nphase, const Eigen::VectorXd& stress) {
  try {
    if (stress.size() != stress_.size()) {
      std::cout << stress_.size() << "\t" << stress.size() << "\n";
      throw std::runtime_error(
          "Particle stress degrees of freedom don't match");
    }
    // Assign stress
    stress_.col(nphase) = stress;
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
  }
}

// Assign particle velocity
//! \tparam Tdim Dimension
//! \tparam Tnphases Number of phases
template <unsigned Tdim, unsigned Tnphases>
void mpm::Particle<Tdim, Tnphases>::assign_velocity(
    unsigned nphase, const Eigen::VectorXd& velocity) {
  try {
    if (velocity.size() != velocity_.size()) {
      throw std::runtime_error(
          "Particle velocity degrees of freedom don't match");
    }
    // Assign velocity
    velocity_.col(nphase) = velocity;
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
  }
}

// Assign particle momentum
//! \tparam Tdim Dimension
//! \tparam Tnphases Number of phases
template <unsigned Tdim, unsigned Tnphases>
void mpm::Particle<Tdim, Tnphases>::assign_momentum(
    unsigned nphase, const Eigen::VectorXd& momentum) {
  try {
    if (momentum.size() != momentum_.size()) {
      throw std::runtime_error(
          "Particle momentum degrees of freedom don't match");
    }
    // Assign momentum
    momentum_.col(nphase) = momentum;
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
  }
}

// Assign particle acceleration
//! \tparam Tdim Dimension
//! \tparam Tnphases Number of phases
template <unsigned Tdim, unsigned Tnphases>
void mpm::Particle<Tdim, Tnphases>::assign_acceleration(
    unsigned nphase, const Eigen::VectorXd& acceleration) {
  try {
    if (acceleration.size() != acceleration_.size()) {
      throw std::runtime_error(
          "Particle acceleration degrees of freedom don't match");
    }
    // Assign acceleration
    acceleration_.col(nphase) = acceleration;
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
  }
}



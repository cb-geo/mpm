//! Construct a particle with id and coordinates
template <unsigned Tdim, unsigned Tnphases>
mpm::Particle<Tdim, Tnphases>::Particle(Index id, const VectorDim& coord)
    : mpm::ParticleBase<Tdim>(id, coord) {
  this->initialise();
}

//! Construct a particle with id, coordinates and status
template <unsigned Tdim, unsigned Tnphases>
mpm::Particle<Tdim, Tnphases>::Particle(Index id, const VectorDim& coord,
                                        bool status)
    : mpm::ParticleBase<Tdim>(id, coord, status) {
  this->initialise();
}

// Initialise particle properties
template <unsigned Tdim, unsigned Tnphases>
void mpm::Particle<Tdim, Tnphases>::initialise() {
  mass_.setZero();
  stress_.setZero();
  velocity_.setZero();
  momentum_.setZero();
  acceleration_.setZero();
}

// Assign a cell to particle
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::assign_cell(
    std::shared_ptr<Cell<Tdim>> cellptr) {
  cell_ = cellptr;
  cell_id_ = cellptr->id();
  this->compute_reference_location();
  return cell_->add_particle_id(this->id());
}

// Compute reference location cell to particle
template <unsigned Tdim, unsigned Tnphases>
void mpm::Particle<Tdim, Tnphases>::compute_reference_location() {
  try {
    // Get reference location of a particle
    if (cell_->is_initialised()) {
      this->reference_location_ =
          cell_->local_coordinates_point(this->coordinates_);
    } else {
      throw std::runtime_error(
          "Cell is not initialised! "
          "Cannot compute local reference coordinates of the particle");
    }
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
  }
}

// Assign stress to the particle
template <unsigned Tdim, unsigned Tnphases>
void mpm::Particle<Tdim, Tnphases>::assign_stress(
    unsigned nphase, const Eigen::Matrix<double, 6, 1>& stress) {
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

// Assign velocity to the particle
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

// Assign momentum to the particle
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

//! Assign acceleration to the particle
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

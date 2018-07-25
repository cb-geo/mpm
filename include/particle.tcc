//! Construct a particle with id and coordinates
template <unsigned Tdim, unsigned Tnphases>
mpm::Particle<Tdim, Tnphases>::Particle(Index id, const VectorDim& coord)
    : mpm::ParticleBase<Tdim>(id, coord) {
  this->initialise();
  cell_ = nullptr;
  material_ = nullptr;
}

//! Construct a particle with id, coordinates and status
template <unsigned Tdim, unsigned Tnphases>
mpm::Particle<Tdim, Tnphases>::Particle(Index id, const VectorDim& coord,
                                        bool status)
    : mpm::ParticleBase<Tdim>(id, coord, status) {
  this->initialise();
  cell_ = nullptr;
  material_ = nullptr;
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

// Assign a material to particle
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::assign_material(
    const std::shared_ptr<Material>& material) {
  bool status = false;
  try {
    // Check if material is valid and properties are set
    if (material != nullptr && material->status()) {
      material_ = material;
      status = true;
    } else {
      throw std::runtime_error("Material is undefined!");
    }
  } catch (std::exception& exception) {
    std::cerr << __FILE__ << __LINE__ << exception.what() << '\n';
  }
  return status;
}

// Compute reference location cell to particle
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::compute_reference_location() {
  bool status = true;
  try {
    // Check if particle has a valid cell ptr
    if (cell_ != nullptr) {
#ifdef _MPM_ISOPARAMETRIC_
      // Get reference location of a particle with isoparametric transformation
      this->xi_ = cell_->transform_real_to_unit_cell(this->coordinates_);
#else
      // Get reference location of a particle on cartesian grid
      this->xi_ = cell_->local_coordinates_point(this->coordinates_);
#endif
    } else {
      throw std::runtime_error(
          "Cell is not initialised! "
          "Cannot compute local reference coordinates of the particle");
    }
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
    status = false;
  }
  return status;
}

// Compute shape functions and gradients
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::compute_shapefn() {
  bool status = true;
  try {
    // Check if particle has a valid cell ptr
    if (cell_ != nullptr) {
      // Compute local coordinates
      this->compute_reference_location();

      // Get shape function ptr of a cell
      const auto sfn = cell_->shapefn_ptr();

      // Compute shape function of the particle
      shapefn_ = sfn->shapefn(this->xi_);
      // Compute gradient shape function of the particle
      grad_shapefn_ = sfn->grad_shapefn(this->xi_);
      // Compute bmatrix of the particle
      bmatrix_ = sfn->bmatrix(this->xi_);
    } else {
      throw std::runtime_error(
          "Cell is not initialised! "
          "cannot compute shapefns for the particle");
    }
  } catch (std::exception& exception) {
    std::cerr << __FILE__ << __LINE__ << "\t" << exception.what() << '\n';
    status = false;
  }
  return status;
}

// Compute volume of particle
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::compute_volume() {
  bool status = true;
  try {
    // Check if particle has a valid cell ptr
    if (cell_ != nullptr) {
      // Volume of the cell / # of particles
      this->volume_ = cell_->volume() / cell_->nparticles();
    } else {
      throw std::runtime_error(
          "Cell is not initialised! "
          "cannot compute volume for the particle");
    }
  } catch (std::exception& exception) {
    std::cerr << __FILE__ << __LINE__ << "\t" << exception.what() << '\n';
    status = false;
  }
  return status;
}

// Compute mass of particle
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::compute_mass(unsigned phase) {
  bool status = true;
  try {
    // Check if particle volume is set and material ptr is valid
    if (volume_ != std::numeric_limits<double>::max() && material_ != nullptr) {
      // Mass = volume of particle * density
      this->mass_(phase) = volume_ * material_->property("density");
    } else {
      throw std::runtime_error(
          "Cell is not initialised! or material is invalid"
          "cannot compute mass for the particle");
    }
  } catch (std::exception& exception) {
    std::cerr << __FILE__ << __LINE__ << "\t" << exception.what() << '\n';
    status = false;
  }
  return status;
}

//! Map particle mass and momentum to nodes
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::map_mass_momentum_to_nodes(unsigned phase) {
  bool status = true;
  try {
    // Check if particle mass is set
    if (mass_(phase) != std::numeric_limits<double>::max()) {
      // Map particle mass and momentum to nodes
      this->cell_->map_mass_momentum_to_nodes(
          this->shapefn_, phase, mass_(phase), velocity_.col(phase));
    } else {
      throw std::runtime_error("Particle mass has not be computed");
    }
  } catch (std::exception& exception) {
    std::cerr << __FILE__ << __LINE__ << "\t" << exception.what() << '\n';
    status = false;
  }
  return status;
}

// Compute strain of the particle
template <unsigned Tdim, unsigned Tnphases>
void mpm::Particle<Tdim, Tnphases>::compute_strain(unsigned phase, double dt) {
  // Strain rate
  Eigen::VectorXd strain_rate = cell_->compute_strain_rate(bmatrix_, phase);
  // dstrain
  Eigen::Matrix<double, 6, 1> dstrain;
  dstrain.setZero();

  // TODO: Remove
  /*
  for (const auto& bm : bmatrix_) {
    std::cout << "Bmatrix: \n";
    for (unsigned i = 0; i < bm.rows(); ++i) {
      for (unsigned j = 0; j < bm.cols(); ++j) {
        std::cout << bm(i, j) << "\t";
      }
      std::cout << "\n";
    }
  }
  */

  switch (Tdim) {
    case (1): {
      dstrain(0) = strain_rate(0) * dt;
      break;
    }
    case (2): {
      dstrain(0) = strain_rate(0) * dt;
      dstrain(1) = strain_rate(1) * dt;
      dstrain(3) = strain_rate(2) * dt;
      break;
    }
    default: {
      dstrain = strain_rate * dt;
      break;
    }
  }

  // TODO: Remove
  /*
  for (unsigned i = 0; i < dstrain.size(); ++i)
    std::cout << "Strain ratae: " << dstrain(i) << "\n";
  */

  // Update strain
  strain_ += dstrain;
}

// Assign stress to the particle
template <unsigned Tdim, unsigned Tnphases>
void mpm::Particle<Tdim, Tnphases>::assign_stress(
    unsigned phase, const Eigen::Matrix<double, 6, 1>& stress) {
  // Assign stress
  stress_.col(phase) = stress;
}

// Assign velocity to the particle
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::assign_velocity(
    unsigned phase, const Eigen::VectorXd& velocity) {
  bool status = false;
  try {
    if (velocity.size() != velocity_.size()) {
      throw std::runtime_error(
          "Particle velocity degrees of freedom don't match");
    }
    // Assign velocity
    velocity_.col(phase) = velocity;
    status = true;
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
    status = false;
  }
  return status;
}

// Assign momentum to the particle
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::assign_momentum(
    unsigned phase, const Eigen::VectorXd& momentum) {
  bool status = false;
  try {
    if (momentum.size() != momentum_.size()) {
      throw std::runtime_error(
          "Particle momentum degrees of freedom don't match");
    }
    // Assign momentum
    momentum_.col(phase) = momentum;
    status = true;
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
    status = false;
  }
  return status;
}

//! Assign acceleration to the particle
template <unsigned Tdim, unsigned Tnphases>
bool mpm::Particle<Tdim, Tnphases>::assign_acceleration(
    unsigned phase, const Eigen::VectorXd& acceleration) {
  bool status = false;
  try {
    if (acceleration.size() != acceleration_.size()) {
      throw std::runtime_error(
          "Particle acceleration degrees of freedom don't match");
    }
    // Assign acceleration
    acceleration_.col(phase) = acceleration;
    status = true;
  } catch (std::exception& exception) {
    std::cerr << exception.what() << '\n';
    status = false;
  }
  return status;
}

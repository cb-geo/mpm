//! Construct a two phase particle with id and coordinates
template <unsigned Tdim>
mpm::FluidParticle<Tdim>::FluidParticle(Index id, const VectorDim& coord)
    : mpm::Particle<Tdim>(id, coord) {

  // Set initial private variables
  this->porosity_ = 1.;
  this->liquid_saturation_ = 1.;

  // Initialize vector data properties
  this->properties_["pressure"] = [&]() {
    Eigen::VectorXd vec_pressure = Eigen::VectorXd::Zero(3);
    vec_pressure[0] = this->pressure();
    // FIXME: This is to check free surface particles
    // TODO: To be removed somewhere
    vec_pressure[1] = this->free_surface();
    return vec_pressure;
  };

  // Logger
  std::string logger =
      "FluidParticle" + std::to_string(Tdim) + "d::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

// Compute mass of fluid particle
template <unsigned Tdim>
void mpm::FluidParticle<Tdim>::compute_mass() noexcept {
  // Check if particle volume is set and material ptr is valid
  assert(volume_ != std::numeric_limits<double>::max() && material_ != nullptr);
  // Mass = volume of particle * mass_density
  this->mass_density_ =
      material_->template property<double>(std::string("density")) *
      this->liquid_saturation_ * this->porosity_;
  this->mass_ = volume_ * mass_density_;
}

// Compute stress
template <unsigned Tdim>
void mpm::FluidParticle<Tdim>::compute_stress() noexcept {
  // Run particle compute stress
  mpm::Particle<Tdim>::compute_stress();

  // Calculate fluid turbulent stress
  this->stress_ += this->compute_turbulent_stress();
}

// Compute turbulent stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1>
    mpm::FluidParticle<Tdim>::compute_turbulent_stress() {
  // Compute turbulent stress depends on the model
  Eigen::Matrix<double, 6, 1> tstress;
  tstress.setZero();

  // TODO: To be implemented

  return tstress;
}

//! Map internal force
template <>
inline void mpm::FluidParticle<1>::map_internal_force() noexcept {
  // initialise a vector of total stress (deviatoric + turbulent - pressure)
  Eigen::Matrix<double, 6, 1> total_stress = this->stress_;
  total_stress(0) -= this->beta_ * state_variables_.at("pressure");
  // Compute nodal internal forces
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    // Compute force: -pstress * volume
    Eigen::Matrix<double, 1, 1> force;
    force[0] = dn_dx_(i, 0) * total_stress[0];
    force *= -1 * this->volume_;

    nodes_[i]->update_internal_force(true, mpm::ParticlePhase::SinglePhase,
                                     force);
  }
}

//! Map internal force
template <>
inline void mpm::FluidParticle<2>::map_internal_force() noexcept {
  // initialise a vector of total stress (deviatoric + turbulent - pressure)
  Eigen::Matrix<double, 6, 1> total_stress = this->stress_;
  total_stress(0) -= this->beta_ * state_variables_.at("pressure");
  total_stress(1) -= this->beta_ * state_variables_.at("pressure");

  // Compute nodal internal forces
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    // Compute force: -pstress * volume
    Eigen::Matrix<double, 2, 1> force;
    force[0] = dn_dx_(i, 0) * total_stress[0] + dn_dx_(i, 1) * total_stress[3];
    force[1] = dn_dx_(i, 1) * total_stress[1] + dn_dx_(i, 0) * total_stress[3];

    force *= -1. * this->volume_;

    nodes_[i]->update_internal_force(true, mpm::ParticlePhase::SinglePhase,
                                     force);
  }
}

//! Map internal force
template <>
inline void mpm::FluidParticle<3>::map_internal_force() noexcept {
  // initialise a vector of total stress (deviatoric + turbulent - pressure)
  Eigen::Matrix<double, 6, 1> total_stress = this->stress_;
  total_stress(0) -= this->beta_ * state_variables_.at("pressure");
  total_stress(1) -= this->beta_ * state_variables_.at("pressure");
  total_stress(2) -= this->beta_ * state_variables_.at("pressure");

  // Compute nodal internal forces
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    // Compute force: -pstress * volume
    Eigen::Matrix<double, 3, 1> force;
    force[0] = dn_dx_(i, 0) * total_stress[0] + dn_dx_(i, 1) * total_stress[3] +
               dn_dx_(i, 2) * total_stress[5];

    force[1] = dn_dx_(i, 1) * total_stress[1] + dn_dx_(i, 0) * total_stress[3] +
               dn_dx_(i, 2) * total_stress[4];

    force[2] = dn_dx_(i, 2) * total_stress[2] + dn_dx_(i, 1) * total_stress[4] +
               dn_dx_(i, 0) * total_stress[5];

    force *= -1. * this->volume_;

    nodes_[i]->update_internal_force(true, mpm::ParticlePhase::SinglePhase,
                                     force);
  }
}

// TODO: Check whether this is important
//! Assign particle pressure constraints
// template <unsigned Tdim>
// bool mpm::FluidParticle<Tdim>::assign_particle_pressure_constraint(
//     double pressure) {
//   bool status = true;
//   try {
//     this->pressure_constraint_ = pressure;
//   } catch (std::exception& exception) {
//     console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
//     status = false;
//   }
//   return status;
// }

//! Compute laplacian matrix element
template <unsigned Tdim>
bool mpm::FluidParticle<Tdim>::map_L_to_cell() {
  bool status = true;
  try {
    // Compute local matrix of Laplacian
    cell_->compute_L_element(dn_dx_, volume_, 1.0);
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Compute element F_s_element and F_m_element
template <unsigned Tdim>
bool mpm::FluidParticle<Tdim>::map_F_to_cell() {
  bool status = true;
  try {
    // Compute local matrix of F
    cell_->compute_F_element(
        shapefn_, dn_dx_,
        material_->template property<double>(std::string("density")) * volume_);
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Compute element K_cor_w_element_
template <unsigned Tdim>
bool mpm::FluidParticle<Tdim>::map_K_cor_to_cell() {
  bool status = true;
  try {
    cell_->compute_K_cor_element(shapefn_, dn_dx_, volume_);

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return status;
}

// Compute updated pore pressure of the particle based on nodal pressure
template <unsigned Tdim>
bool mpm::FluidParticle<Tdim>::compute_updated_pressure() {
  bool status = true;

  try {
    double pressure_increment = 0;
    for (unsigned i = 0; i < nodes_.size(); ++i) {
      pressure_increment += shapefn_(i) * nodes_[i]->pressure_increment();
    }
    // Get interpolated nodal pore pressure
    state_variables_.at("pressure") =
        state_variables_.at("pressure") * beta_ + pressure_increment;
    // Apply free surface
    if (this->free_surface()) state_variables_.at("pressure") = 0.0;
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// //! Assign particle permeability
// template <unsigned Tdim>
// bool mpm::FluidParticle<Tdim>::assign_permeability() {
//   bool status = true;
//   try {
//     // Check if material ptr is valid
//     if (material_ != nullptr) {
//       // Porosity parameter
//       const double k_p =
//           std::pow(this->porosity_, 3) / std::pow((1. - this->porosity_), 2);
//       switch (Tdim) {
//         case (1): {
//           permeability_(0) = material_->template property<double>("k_x");
//           c1_(0) = permeability_(0) / k_p;
//           break;
//         }
//         case (2): {
//           permeability_(0) = material_->template property<double>("k_x");
//           permeability_(1) = material_->template property<double>("k_y");
//           c1_(0) = permeability_(0) / k_p;
//           c1_(1) = permeability_(1) / k_p;
//           break;
//         }
//         default: {
//           permeability_(0) = material_->template property<double>("k_x");
//           permeability_(1) = material_->template property<double>("k_y");
//           permeability_(2) = material_->template property<double>("k_z");
//           c1_(0) = permeability_(0) / k_p;
//           c1_(1) = permeability_(1) / k_p;
//           c1_(2) = permeability_(2) / k_p;
//           break;
//         }
//       }
//     } else {
//       throw std::runtime_error("Material is invalid");
//     }
//   } catch (std::exception& exception) {
//     console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
//     status = false;
//   }
//   return status;
// }

// //! Update particle permeability
// template <unsigned Tdim>
// bool mpm::FluidParticle<Tdim>::update_permeability() {
//   bool status = true;
//   try {
//     // Porosity parameter
//     const double k_p =
//         std::pow(this->porosity_, 3) / std::pow((1. - this->porosity_), 2);

//     // Update permeability by KC equation
//     permeability_ = k_p * c1_;

//   } catch (std::exception& exception) {
//     console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
//     status = false;
//   }
//   return status;
// }

// //! Map drag force
// template <unsigned Tdim>
// bool mpm::FluidParticle<Tdim>::map_drag_force_coefficient() {
//   bool status = true;
//   try {
//     // Update permeability
//     this->update_permeability();
//     // Initialise drag force coefficient
//     VectorDim drag_force_coefficient;
//     drag_force_coefficient.setZero();

//     // Check if permeability coefficient is valid
//     for (unsigned i = 0; i < Tdim; ++i) {
//       if (permeability_(i) > 0.)
//         drag_force_coefficient(i) =
//             porosity_ * porosity_ * 9.81 *
//             // liquid_material_->template property<double>(
//             //     std::string("density")) /
//             material_->template property<double>(std::string("density")) /
//             permeability_(i);
//       else
//         throw std::runtime_error("Permeability coefficient is invalid");
//     }

//     // Map drag forces from particle to nodes
//     for (unsigned j = 0; j < nodes_.size(); ++j)
//       nodes_[j]->update_drag_force_coefficient(
//           true, drag_force_coefficient * this->volume_ * shapefn_(j));

//   } catch (std::exception& exception) {
//     console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
//     status = false;
//   }
//   return status;
// }
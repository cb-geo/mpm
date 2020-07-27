namespace mpm {
namespace twophaseparticle {

// Assign particle porosity
template <unsigned Tdim>
void assign_porosity(std::shared_ptr<mpm::ParticleBase<Tdim>> particle) {
  // Check if material ptr is valid
  if (particle->material() != nullptr) {
    double porosity = particle->material()->template property<double>(
        std::string("porosity"));

    // Check if the porosity value is valid
    if (porosity < 0. || porosity > 1.)
      throw std::runtime_error(
          "Particle porosity is negative or larger than one");

    // Update porosity
    particle->update_scalar_property(mpm::properties::Scalar::Porosity, false,
                                     porosity);
  } else {
    throw std::runtime_error("Material is invalid, could not assign porosity");
  }
}

//! Assign particle permeability
template <unsigned Tdim>
void assign_permeability(std::shared_ptr<mpm::ParticleBase<Tdim>> particle) {

  // Check if material ptr is valid
  if (particle->material() != nullptr) {
    // Porosity parameter k_p
    const double k_p =
        std::pow(particle->scalar_property(mpm::properties::Scalar::Porosity),
                 3) /
        std::pow(
            (1. - particle->scalar_property(mpm::properties::Scalar::Porosity)),
            2);
    // Initialise permeability vector
    Eigen::Matrix<double, Tdim, 1> permeability;
    permeability.setZero();

    // Different dimensions permeability
    switch (Tdim) {
      case (1): {
        permeability(0) =
            particle->material()->template property<double>("k_x") / k_p;
        break;
      }
      case (2): {
        permeability(0) =
            particle->material()->template property<double>("k_x") / k_p;
        permeability(1) =
            particle->material()->template property<double>("k_y") / k_p;

        break;
      }
      default: {
        permeability(0) =
            particle->material()->template property<double>("k_x") / k_p;
        permeability(1) =
            particle->material()->template property<double>("k_y") / k_p;
        permeability(2) =
            particle->material()->template property<double>("k_z") / k_p;
        break;
      }
    }
    // Assign permeability vector
    particle->update_vector_property(mpm::properties::Vector::Permeability,
                                     false, permeability);
  } else {
    throw std::runtime_error(
        "Material is invalid, could not assign permeability");
  }
}

// Compute particle mass for two phase
template <unsigned Tdim>
void compute_mass(std::shared_ptr<mpm::ParticleBase<Tdim>> particle) noexcept {
  // Check if particle volume is set and material ptr is valid
  assert(particle->volume() != std::numeric_limits<double>::max() &&
         particle->material() != nullptr);

  // Mass = volume of particle * mass_density
  // Solid density
  particle->update_scalar_property(
      mpm::properties::Scalar::MassDensity, false,
      particle->material()->template property<double>(std::string("density")));

  // Update particle solid mass
  particle->update_scalar_property(
      mpm::properties::Scalar::Mass, false,
      (1 - particle->scalar_property(mpm::properties::Scalar::Porosity)) *
          particle->volume() * particle->mass_density());

  // Update particle water mass
  particle->update_scalar_property(
      mpm::properties::Scalar::LiquidMass, false,
      particle->scalar_property(mpm::properties::Scalar::Porosity) *
          particle->volume() *
          particle->scalar_property(
              mpm::properties::Scalar::LiquidMassDensity));
}

//! Map liquid mass and momentum to nodes
template <unsigned Tdim>
void map_mass_momentum_to_nodes(
    std::shared_ptr<mpm::ParticleBase<Tdim>> particle) noexcept {
  // Check if particle mass is set
  assert(particle->mass() != std::numeric_limits<double>::max());

  // Map solid mass and momentum to nodes
  particle->map_scalar_property_to_nodes(mpm::properties::Scalar::Mass, true,
                                         mpm::ParticlePhase::Solid);
  particle->map_vector_property_to_nodes(
      mpm::properties::Vector::Momentum, true, mpm::ParticlePhase::Solid,
      particle->mass() * particle->velocity());

  // Check if particle mass is set
  assert(particle->scalar_property(mpm::properties::Scalar::LiquidMass) !=
         std::numeric_limits<double>::max());

  // Map liquid mass and momentum to nodes
  particle->map_scalar_property_to_nodes(
      mpm::properties::Scalar::LiquidMass, true, mpm::ParticlePhase::Liquid,
      particle->scalar_property(mpm::properties::Scalar::LiquidMass));
  particle->map_vector_property_to_nodes(
      mpm::properties::Vector::Momentum, true, mpm::ParticlePhase::Liquid,
      particle->scalar_property(mpm::properties::Scalar::LiquidMass) *
          particle->vector_property(mpm::properties::Vector::LiquidVelocity));
}

//! Map particle pore liquid pressure to nodes
template <unsigned Tdim>
void map_pore_pressure_to_nodes(
    std::shared_ptr<mpm::ParticleBase<Tdim>> particle) noexcept {
  // Check if particle mass is set
  assert(particle->scalar_property(mpm::properties::Scalar::LiquidMass) !=
         std::numeric_limits<double>::max());

  // Map particle liquid mass and pore pressure to nodes
  map_scalar_property_to_nodes(
      mpm::properties::Scalar::MassPressure, true, mpm::ParticlePhase::Liquid,
      particle->scalar_property(mpm::properties::Scalar::LiquidMass) *
          particle->scalar_property(mpm::properties::Scalar::PorePressure));
}

//! Map body force for both mixture and liquid
template <unsigned Tdim>
void map_body_force(std::shared_ptr<mpm::ParticleBase<Tdim>> particle,
                    const Eigen::Matrix<double, Tdim, 1>& pgravity) noexcept {
  //! Map body force for mixture
  particle->map_vector_property_to_nodes(
      mpm::properties::Vector::ExternalForce, true, mpm::ParticlePhase::Mixture,
      pgravity * (particle->mass() + particle->scalar_property(
                                         mpm::properties::Scalar::LiquidMass)));
  //! Map body force for liquid
  particle->map_vector_property_to_nodes(
      mpm::properties::Vector::ExternalForce, true, mpm::ParticlePhase::Liquid,
      pgravity *
          particle->scalar_property(mpm::properties::Scalar::LiquidMass));
}

//! Map drag force coefficient
template <unsigned Tdim>
void map_drag_force_coefficient(
    std::shared_ptr<mpm::ParticleBase<Tdim>> particle) {
  // Update permeability
  auto permeability =
      particle->vector_property(mpm::properties::Vector::Permeability) *
      std::pow(particle->scalar_property(mpm::properties::Scalar::Porosity),
               3) /
      std::pow(
          (1. - particle->scalar_property(mpm::properties::Scalar::Porosity)),
          2);
  // Initialise drag force coefficient
  Eigen::Matrix<double, Tdim, 1> drag_force_coefficient;
  drag_force_coefficient.setZero();

  // Check if permeability coefficient is valid
  for (unsigned i = 0; i < Tdim; ++i) {
    drag_force_coefficient(i) =
        std::pow(particle->scalar_property(mpm::properties::Scalar::Porosity),
                 2) *
        9.81 *
        particle->scalar_property(mpm::properties::Scalar::LiquidMassDensity) /
        permeability(i);
  }

  particle->map_vector_property_to_nodes(
      mpm::properties::Vector::DragForce, true, mpm::ParticlePhase::Solid,
      drag_force_coefficient * particle->volume());
}

// Update particle porosity
template <unsigned Tdim>
void update_porosity(std::shared_ptr<mpm::ParticleBase<Tdim>> particle,
                     double dt) {
  // Update particle porosity
  double porosity =
      1 - (1 - particle->scalar_property(mpm::properties::Scalar::Porosity)) /
              (1 + dt * particle->strain_rate().head(Tdim).sum());
  // Check if the value is valid
  if (porosity < 0) porosity = 1E-5;
  if (porosity < 1) porosity = 1 - 1E-5;

  // Assign new particle porosity
  particle->update_scalar_property(mpm::properties::Scalar::Porosity, false,
                                   porosity);
}

//! Initial pore pressure by water table
template <unsigned Tdim>
void initialise_pore_pressure_watertable(
    const unsigned dir_v, const unsigned dir_h,
    std::map<double, double>& refernece_points,
    std::shared_ptr<mpm::ParticleBase<Tdim>> particle) {
  // Initialise left boundary position (coordinate) and h0
  double left_boundary = std::numeric_limits<double>::lowest();
  double h0_left = 0.;
  // Initialise right boundary position (coordinate) and h0
  double right_boundary = std::numeric_limits<double>::max();
  double h0_right = 0.;
  // Position and h0 of particle (coordinate)
  const double position = particle->coordinates().at(dir_h);
  // Iterate over each refernece_points
  for (const auto& refernece_point : refernece_points) {
    // Find boundary
    if (refernece_point.first > left_boundary &&
        refernece_point.first <= position) {
      // Left boundary position and h0
      left_boundary = refernece_point.first;
      h0_left = refernece_point.second;
    } else if (refernece_point.first > position &&
               refernece_point.first <= right_boundary) {
      // Right boundary position and h0
      right_boundary = refernece_point.first;
      h0_right = refernece_point.second;
    }
  }
  // Check if the boundaries are assigned
  if (left_boundary != std::numeric_limits<double>::lowest()) {
    // Particle with left and right boundary
    if (right_boundary != std::numeric_limits<double>::max()) {
      particle->update_scalar_property(
          mpm::properties::Scalar::LiquidMassDensity, false,
          ((h0_right - h0_left) / (right_boundary - left_boundary) *
               (position - left_boundary) +
           h0_left - particle->coordinates_(dir_v)) *
              1000 * 9.81);
    } else
      // Particle with only left boundary
      particle->update_scalar_property(
          mpm::properties::Scalar::LiquidMassDensity, false,
          (h0_left - particle->coordinates_(dir_v)) * 1000 * 9.81);
  }
  // Particle with only right boundary
  else if (right_boundary != std::numeric_limits<double>::max())
    particle->update_scalar_property(
        mpm::properties::Scalar::LiquidMassDensity, false,
        (h0_right - particle->coordinates_(dir_v)) * 1000 * 9.81);

  else
    throw std::runtime_error(
        "Particle pore pressure can not be initialised by water table");
}

}  // namespace twophaseparticle
}  // namespace mpm
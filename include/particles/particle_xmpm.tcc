//! Construct a particle with id and coordinates
template <unsigned Tdim>
mpm::ParticleXMPM<Tdim>::ParticleXMPM(Index id, const VectorDim& coord)
    : mpm::Particle<Tdim>(id, coord) {

  this->initialise();
  // Logger
  std::string logger =
      "particlexmpm" + std::to_string(Tdim) + "d::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

//! Construct a particle with id, coordinates and status
template <unsigned Tdim>
mpm::ParticleXMPM<Tdim>::ParticleXMPM(Index id, const VectorDim& coord,
                                      bool status)
    : mpm::Particle<Tdim>(id, coord, status) {
  this->initialise();
  //! Logger
  std::string logger =
      "particlexmpm" + std::to_string(Tdim) + "d::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

// Initialise particle properties
template <unsigned Tdim>
void mpm::ParticleXMPM<Tdim>::initialise() {
  levelset_phi_ = 0.;
}

//! Map particle mass and momentum to nodes
template <unsigned Tdim>
void mpm::ParticleXMPM<Tdim>::map_mass_momentum_to_nodes() noexcept {
  // Check if particle mass is set
  assert(mass_ != std::numeric_limits<double>::max());

  // Map mass and momentum to nodes
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    nodes_[i]->update_mass(true, mpm::ParticlePhase::Solid,
                           mass_ * shapefn_[i]);
    nodes_[i]->update_momentum(true, mpm::ParticlePhase::Solid,
                               mass_ * shapefn_[i] * velocity_);
    if (nodes_[i]->discontinuity_enrich()) {
      // Unit 1x1 Eigen matrix to be used with scalar quantities
      Eigen::Matrix<double, 1, 1> nodal_mass;
      nodal_mass(0, 0) = sgn(levelset_phi_) * mass_ * shapefn_[i];
      // Map enriched mass and momentum to nodes
      nodes_[i]->update_discontinuity_property(true, "mass_enrich", nodal_mass,
                                               0, 1);
      nodes_[i]->update_discontinuity_property(true, "momenta_enrich",
                                               velocity_ * nodal_mass, 0, Tdim);
    }
  }
}

// Compute strain rate of the particle
template <>
inline Eigen::Matrix<double, 6, 1> mpm::ParticleXMPM<3>::compute_strain_rate(
    const Eigen::MatrixXd& dn_dx, unsigned phase) noexcept {
  // Define strain rate
  Eigen::Matrix<double, 6, 1> strain_rate = Eigen::Matrix<double, 6, 1>::Zero();
  const double tolerance = 1.E-16;
  Eigen::Vector3d vel;
  vel.setZero();
  // Compute corresponding nodal velocity
  for (unsigned i = 0; i < this->nodes_.size(); ++i) {
    if (nodes_[i]->discontinuity_enrich()) {
      double nodal_mass =
          nodes_[i]->mass(phase) +
          sgn(levelset_phi_) *
              nodes_[i]->discontinuity_property("mass_enrich", 1)(0, 0);
      if (nodal_mass < tolerance) continue;

      vel =
          (nodes_[i]->momentum(phase) +
           sgn(levelset_phi_) *
               nodes_[i]->discontinuity_property("momenta_enrich", 3).col(0)) /
          nodal_mass;
    } else {
      double nodal_mass = nodes_[i]->mass(phase);
      if (nodal_mass < tolerance) continue;
      vel = nodes_[i]->momentum(phase) / nodal_mass;
    }

    strain_rate[0] += dn_dx(i, 0) * vel[0];
    strain_rate[1] += dn_dx(i, 1) * vel[1];
    strain_rate[2] += dn_dx(i, 2) * vel[2];
    strain_rate[3] += dn_dx(i, 1) * vel[0] + dn_dx(i, 0) * vel[1];
    strain_rate[4] += dn_dx(i, 2) * vel[1] + dn_dx(i, 1) * vel[2];
    strain_rate[5] += dn_dx(i, 2) * vel[0] + dn_dx(i, 0) * vel[2];
  }

  for (unsigned i = 0; i < strain_rate.size(); ++i)
    if (std::fabs(strain_rate[i]) < 1.E-15) strain_rate[i] = 0.;
  return strain_rate;
}

// Compute strain of the particle
template <unsigned Tdim>
void mpm::ParticleXMPM<Tdim>::compute_strain(double dt) noexcept {
  // Assign strain rate
  strain_rate_ = this->compute_strain_rate(dn_dx_, mpm::ParticlePhase::Solid);
  // Update dstrain
  dstrain_ = strain_rate_ * dt;
  // Update strain
  strain_ += dstrain_;

  // Compute at centroid
  // Strain rate for reduced integration
  const Eigen::Matrix<double, 6, 1> strain_rate_centroid =
      this->compute_strain_rate(dn_dx_centroid_, mpm::ParticlePhase::Solid);

  // Assign volumetric strain at centroid
  dvolumetric_strain_ = dt * strain_rate_centroid.head(Tdim).sum();
  volumetric_strain_centroid_ += dvolumetric_strain_;
}

//! Map body force
template <unsigned Tdim>
void mpm::ParticleXMPM<Tdim>::map_body_force(
    const VectorDim& pgravity) noexcept {
  // Compute nodal body forces
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    nodes_[i]->update_external_force(true, mpm::ParticlePhase::Solid,
                                     (pgravity * mass_ * shapefn_(i)));
    if (nodes_[i]->discontinuity_enrich())
      nodes_[i]->update_discontinuity_property(
          true, "external_force_enrich",
          sgn(levelset_phi_) * pgravity * mass_ * shapefn_(i), 0, Tdim);
  }
}

//! Map internal force
template <>
inline void mpm::ParticleXMPM<1>::map_internal_force() noexcept {
  // Compute nodal internal forces
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    // Compute force: -pstress * volume
    Eigen::Matrix<double, 1, 1> force;
    force[0] = -1. * dn_dx_(i, 0) * volume_ * stress_[0];

    nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Solid, force);
    if (nodes_[i]->discontinuity_enrich())
      nodes_[i]->update_discontinuity_property(
          true, "internal_force_enrich", sgn(levelset_phi_) * force, 0, 1);
  }
}

//! Map internal force
template <>
inline void mpm::ParticleXMPM<2>::map_internal_force() noexcept {
  // Compute nodal internal forces
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    // Compute force: -pstress * volume
    Eigen::Matrix<double, 2, 1> force;
    force[0] = dn_dx_(i, 0) * stress_[0] + dn_dx_(i, 1) * stress_[3];
    force[1] = dn_dx_(i, 1) * stress_[1] + dn_dx_(i, 0) * stress_[3];

    force *= -1. * this->volume_;

    nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Solid, force);
    if (nodes_[i]->discontinuity_enrich())
      nodes_[i]->update_discontinuity_property(
          true, "internal_force_enrich", sgn(levelset_phi_) * force, 0, 2);
  }
}

//! Map internal force
template <>
inline void mpm::ParticleXMPM<3>::map_internal_force() noexcept {
  // Compute nodal internal forces
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    // Compute force: -pstress * volume
    Eigen::Matrix<double, 3, 1> force;
    force[0] = dn_dx_(i, 0) * stress_[0] + dn_dx_(i, 1) * stress_[3] +
               dn_dx_(i, 2) * stress_[5];

    force[1] = dn_dx_(i, 1) * stress_[1] + dn_dx_(i, 0) * stress_[3] +
               dn_dx_(i, 2) * stress_[4];

    force[2] = dn_dx_(i, 2) * stress_[2] + dn_dx_(i, 1) * stress_[4] +
               dn_dx_(i, 0) * stress_[5];

    force *= -1. * this->volume_;

    nodes_[i]->update_internal_force(true, mpm::ParticlePhase::Solid, force);
    if (nodes_[i]->discontinuity_enrich())
      nodes_[i]->update_discontinuity_property(
          true, "internal_force_enrich", sgn(levelset_phi_) * force, 0, 3);
  }
}

// Compute updated position of the particle
template <unsigned Tdim>
void mpm::ParticleXMPM<Tdim>::compute_updated_position(
    double dt, bool velocity_update) noexcept {
  // Check if particle has a valid cell ptr
  assert(cell_ != nullptr);
  // Get interpolated nodal velocity
  Eigen::Matrix<double, Tdim, 1> nodal_velocity =
      Eigen::Matrix<double, Tdim, 1>::Zero();
  const double tolerance = 1.E-16;
  unsigned int phase = mpm::ParticlePhase::Solid;
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    if (nodes_[i]->discontinuity_enrich()) {
      double nodal_mass =
          nodes_[i]->mass(phase) +
          sgn(levelset_phi_) *
              nodes_[i]->discontinuity_property("mass_enrich", 1)(0, 0);
      if (nodal_mass < tolerance) continue;

      nodal_velocity += shapefn_[i] *
                        (nodes_[i]->momentum(phase) +
                         sgn(levelset_phi_) * nodes_[i]->discontinuity_property(
                                                  "momenta_enrich", 3)) /
                        nodal_mass;
    } else {
      double nodal_mass = nodes_[i]->mass(phase);
      if (nodal_mass < tolerance) continue;
      nodal_velocity += shapefn_[i] * nodes_[i]->momentum(phase) / nodal_mass;
    }
  }
  // Acceleration update
  if (!velocity_update) {
    // Get interpolated nodal acceleration
    Eigen::Matrix<double, Tdim, 1> nodal_acceleration =
        Eigen::Matrix<double, Tdim, 1>::Zero();
    for (unsigned i = 0; i < nodes_.size(); ++i) {
      if (nodes_[i]->discontinuity_enrich()) {
        double nodal_mass =
            nodes_[i]->mass(phase) +
            sgn(levelset_phi_) *
                nodes_[i]->discontinuity_property("mass_enrich", 1)(0, 0);
        if (nodal_mass < tolerance) continue;

        auto force = nodes_[i]->internal_force(phase) +
                     sgn(levelset_phi_) * nodes_[i]->discontinuity_property(
                                              "internal_force_enrich", 3) +
                     nodes_[i]->external_force(phase) +
                     sgn(levelset_phi_) * nodes_[i]->discontinuity_property(
                                              "external_force_enrich", 3);

        nodal_acceleration += shapefn_[i] * force / nodal_mass;
      } else {
        double nodal_mass = nodes_[i]->mass(phase);
        if (nodal_mass < tolerance) continue;

        auto force =
            nodes_[i]->internal_force(phase) + nodes_[i]->external_force(phase);

        nodal_acceleration += shapefn_[i] * force / nodal_mass;
      }
    }

    // Update particle velocity from interpolated nodal acceleration
    this->velocity_ += nodal_acceleration * dt;
  }
  // Update particle velocity using interpolated nodal velocity
  else
    this->velocity_ = nodal_velocity;

  // New position  current position + velocity * dt
  this->coordinates_ += nodal_velocity * dt;
  // Update displacement (displacement is initialized from zero)
  this->displacement_ += nodal_velocity * dt;
}

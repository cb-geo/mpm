//! Construct a two phase particle with id and coordinates
template <unsigned Tdim>
mpm::TwoPhaseParticle<Tdim>::TwoPhaseParticle(Index id, const VectorDim& coord)
    : mpm::Particl<Tdim>(id, coord) {
  this->initialise_liquid_phase();

  // Set material pointer to null
  liquid_material_ = nullptr;
  // Logger
  std::string logger =
      "twophaseparticle" + std::to_string(Tdim) + "d::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

// Initialise liquid phase particle properties
template <unsigned Tdim>
void mpm::TwoPhaseParticle<Tdim>::initialise_liquid_phase() {
  liquid_mass_ = 0.;
  liquid_displacement_.setZero();
  liquid_velocity_.setZero();
  liquid_strain_rate_.setZero();
  liquid_strain_.setZero();
  set_traction_ = false;
  liquid_traction_.setZero();
  pore_pressure_ = 0.;
  liquid_saturation_ = 1.;

  this->liquid_properties_["velocities"] = [&]() { return liquid_velocity(); };
  this->liquid_properties_["pressure"] = [&]() {
    Eigen::VectorXd vec_pressure(1);
    vec_pressure << this->pore_pressure();
    return vec_pressure;
  };
}

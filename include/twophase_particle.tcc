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

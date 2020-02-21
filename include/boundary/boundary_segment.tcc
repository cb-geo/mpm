// constructor
template <unsigned Tdim>
mpm::BoundarySegment<Tdim>::BoundarySegment(
    mpm::Index id,
    const std::vector<std::shared_ptr<ParticleBase<Tdim>>>& pointptrs)
    : id_{id}, points_{pointptrs} {

  cells_.clear();

  //! Logger
  std::string logger = "boundary segment" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

//! Constructor
template <unsigned Tdim>
mpm::MPMExplicit<Tdim>::MPMExplicit(std::unique_ptr<IO>&& io)
    : mpm::MPM<Tdim>(std::move(io)) {
  //! Logger
  console_ = spdlog::get("MPMExplicit");
  try {
    analysis_ = io_->analysis();
    dt_ = analysis_["dt"].template get<double>();
    nsteps_ = analysis_["nsteps"].template get<mpm::Index>();
  } catch (std::domain_error& domain_error) {
    console_->error("Get analysis object: {}", domain_error.what());
  }
}

// Initialise
template <unsigned Tdim>
bool mpm::MPMExplicit<Tdim>::initialise() {
  return true;
}

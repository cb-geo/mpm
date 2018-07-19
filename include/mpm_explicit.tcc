//! Constructor
template <unsigned Tdim>
mpm::MPMExplicit<Tdim>::MPMExplicit(std::unique_ptr<IO>&& io)
    : mpm::MPM<Tdim>(std::move(io)) {
  //! Logger
  console_ = spdlog::get("MPMExplicit");
  try {
    analysis_ = io_->analysis();
    // Time-step size
    dt_ = analysis_["dt"].template get<double>();
    // Number of time steps
    nsteps_ = analysis_["nsteps"].template get<mpm::Index>();
    // Gravity
    if (analysis_.at("gravity").is_array() &&
        analysis_.at("gravity").size() == Tdim) {
      for (unsigned i = 0; i < analysis_.at("gravity").size(); ++i) {
        gravity_[i] = analysis_.at("gravity").at(i);
      }
    }
  } catch (std::domain_error& domain_error) {
    console_->error("Get analysis object: {}", domain_error.what());
  }
}

// Initialise
template <unsigned Tdim>
bool mpm::MPMExplicit<Tdim>::initialise() {
  bool status = false;
  try {
    auto mesh_props = io_->json_object("mesh");
    status = true;
  } catch (std::domain_error& domain_error) {
    console_->error("Get mesh object: {}", domain_error.what());
  }
  return status;
}

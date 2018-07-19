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
    // Get mesh properties
    auto mesh_props = io_->json_object("mesh");
    // Get Mesh reader from JSON object
    const std::string reader =
        mesh_props["mesh_reader"].template get<std::string>();
    // Create a mesh reader
    auto mesh_reader = Factory<mpm::ReadMesh<Tdim>>::instance()->create(reader);

    // Read particles
    auto particles = mesh_reader->read_particles(io_->file_name("particles"));

    for (const auto& particle : particles) {
      for (unsigned i = 0; i < particle.size(); ++i) {
        std::cout << particle(i) << "\t";
      }
      std::cout << "\n";
    }

    std::cout << "Number of particles " << particles.size() << "\n";

    status = true;
  } catch (std::domain_error& domain_error) {
    console_->error("Get mesh object: {}", domain_error.what());
  }
  return status;
}

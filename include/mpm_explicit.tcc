//! Constructor
template <unsigned Tdim>
mpm::MPMExplicit<Tdim>::MPMExplicit(std::unique_ptr<IO>&& io)
    : mpm::MPM<Tdim>(std::move(io)) {
  //! Logger
  console_ = spdlog::get("MPMExplicit");

  // Create a mesh with global id 0
  const mpm::Index id = 0;
  meshes_.emplace_back(std::make_unique<mpm::Mesh<Tdim>>(id));

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

    // Global Index
    mpm::Index gid = 0;
    // Node type
    const auto node_type = mesh_props["node_type"].template get<std::string>();
    // Create nodes from file
    meshes_.at(0)->create_nodes(
        gid,                                                    // global id
        node_type,                                              // node type
        mesh_reader->read_mesh_nodes(io_->file_name("mesh")));  // coordinates

    std::cout << "Number of nodes: " << meshes_.at(0)->nnodes() << "\n";

    // Shape function name
    const auto cell_type = mesh_props["cell_type"].template get<std::string>();
    // Shape function
    std::shared_ptr<mpm::ShapeFn<Tdim>> shapefn =
        Factory<mpm::ShapeFn<Tdim>>::instance()->create(cell_type);

    // Create cells from file
    meshes_.at(0)->create_cells(
        std::move(gid),  // global id
        shapefn,         // Shape function
        mesh_reader->read_mesh_cells(io_->file_name("mesh")));  // Node ids

    std::cout << "Number of cells " << meshes_.at(0)->ncells() << "\n";

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

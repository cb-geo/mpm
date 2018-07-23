//! Constructor
template <unsigned Tdim>
mpm::MPMExplicit<Tdim>::MPMExplicit(std::unique_ptr<IO>&& io)
    : mpm::MPM<Tdim>(std::move(io)) {
  //! Logger
  console_ = spdlog::get("MPMExplicit");

  // Create a mesh with global id 0
  const mpm::Index id = 0;
  meshes_.emplace_back(std::make_unique<mpm::Mesh<Tdim>>(id));

  // Empty all materials
  materials_.clear();

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
    console_->error(" {} {} Get analysis object: {}", __FILE__, __LINE__,
                    domain_error.what());
  }
}

// Initialise
template <unsigned Tdim>
bool mpm::MPMExplicit<Tdim>::initialise_mesh_particles() {
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

    // Shape function name
    const auto cell_type = mesh_props["cell_type"].template get<std::string>();
    // Shape function
    std::shared_ptr<mpm::ShapeFn<Tdim>> shapefn =
        Factory<mpm::ShapeFn<Tdim>>::instance()->create(cell_type);
    // Create cells from file
    meshes_.at(0)->create_cells(
        gid,      // global id
        shapefn,  // Shape function
        mesh_reader->read_mesh_cells(io_->file_name("mesh")));  // Node ids

    // Particle type
    const auto particle_type =
        mesh_props["particle_type"].template get<std::string>();
    // Create particles from file
    meshes_.at(0)->create_particles(gid,            // global id
                                    particle_type,  // particle type
                                    mesh_reader->read_particles(io_->file_name(
                                        "particles")));  // coordinates

    // Locate particles in cell
    auto unlocatable_particles = meshes_.at(0)->locate_particles_mesh();

    if (!unlocatable_particles.empty())
      throw std::runtime_error("Particle outside the mesh domain");

    status = true;
  } catch (std::exception& exception) {
    console_->error("{} {} Reading mesh and particles: {}", __FILE__, __LINE__,
                    exception.what());
  }
  return status;
}

template <unsigned Tdim>
bool mpm::MPMExplicit<Tdim>::initialise_materials() {
  bool status = false;
  materials_.clear();
  try {
    // Get materials properties
    auto materials = io_->json_object("materials");

    for (const auto material : materials) {
      // Create a new material from JSON object
      const std::string material_type =
          material["type"].template get<std::string>();

      unsigned material_id = material["id"].template get<unsigned>();

      // Create material
      auto mat = Factory<mpm::Material, unsigned>::instance()->create(
          material_type, std::move(material_id));

      // Add material to list
      materials_.emplace_back(mat);
    }
    status = true;
  } catch (std::exception& exception) {
    console_->error("{} {} Reading materials: {}", __FILE__, __LINE__,
                    exception.what());
  }
  return status;
}

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

// Initialise mesh and particles
template <unsigned Tdim>
bool mpm::MPMExplicit<Tdim>::initialise_mesh_particles() {
  bool status = true;
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
    bool node_status = meshes_.at(0)->create_nodes(
        gid,                                                    // global id
        node_type,                                              // node type
        mesh_reader->read_mesh_nodes(io_->file_name("mesh")));  // coordinates

    if (!node_status)
      throw std::runtime_error("Addition of nodes to mesh failed");

    // Shape function name
    const auto cell_type = mesh_props["cell_type"].template get<std::string>();
    // Shape function
    std::shared_ptr<mpm::ShapeFn<Tdim>> shapefn =
        Factory<mpm::ShapeFn<Tdim>>::instance()->create(cell_type);

    // Create cells from file
    bool cell_status = meshes_.at(0)->create_cells(
        gid,      // global id
        shapefn,  // Shape function
        mesh_reader->read_mesh_cells(io_->file_name("mesh")));  // Node ids

    if (!cell_status)
      throw std::runtime_error("Addition of cells to mesh failed");

    // Particle type
    const auto particle_type =
        mesh_props["particle_type"].template get<std::string>();
    // Create particles from file
    bool particle_status = meshes_.at(0)->create_particles(
        gid,            // global id
        particle_type,  // particle type
        mesh_reader->read_particles(
            io_->file_name("particles")));  // coordinates

    if (!particle_status)
      throw std::runtime_error("Addition of particles to mesh failed");

    // Locate particles in cell
    auto unlocatable_particles = meshes_.at(0)->locate_particles_mesh();

    if (!unlocatable_particles.empty())
      throw std::runtime_error("Particle outside the mesh domain");

  } catch (std::exception& exception) {
    console_->error("#{}: Reading mesh and particles: {}", __LINE__,
                    exception.what());
    status = false;
  }
  return status;
}

// Initialise materials
template <unsigned Tdim>
bool mpm::MPMExplicit<Tdim>::initialise_materials() {
  bool status = true;
  try {
    // Get materials properties
    auto materials = io_->json_object("materials");

    for (const auto material_props : materials) {
      // Get material type
      const std::string material_type =
          material_props["type"].template get<std::string>();

      // Get material id
      unsigned material_id = material_props["id"].template get<unsigned>();

      // Create a new material from JSON object
      auto mat = Factory<mpm::Material, unsigned>::instance()->create(
          material_type, std::move(material_id));

      // Initialise material properties
      mat->properties(material_props);

      // Add material to list
      auto result = materials_.insert(std::make_pair(mat->id(), mat));

      // If insert material failed
      if (!result.second) {
        status = false;
        throw std::runtime_error(
            "New material cannot be added, insertion failed");
      }
    }
  } catch (std::exception& exception) {
    console_->error("#{}: Reading materials: {}", __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! MPM Explicit solver
template <unsigned Tdim>
bool mpm::MPMExplicit<Tdim>::solve() {
  bool status = true;

  // Phase
  const unsigned phase = 0;
  // Initialise material
  bool mat_status = this->initialise_materials();
  if (!mat_status) status = false;

  // Initialise mesh and materials
  bool mesh_status = this->initialise_mesh_particles();
  if (!mesh_status) status = false;

  // Assign material to particles
  // Get mesh properties
  auto mesh_props = io_->json_object("mesh");
  // Material id
  const unsigned material_id =
      mesh_props["material_id"].template get<unsigned>();

  // Get material from list of materials
  auto material = materials_.at(material_id);

  // Iterate over each particle to assign material
  meshes_.at(0)->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::assign_material,
                std::placeholders::_1, material));

  for (mpm::Index steps = 0; steps < this->nsteps_; ++steps) {

    // Initialise nodes
    meshes_.at(0)->iterate_over_nodes(
        std::bind(&mpm::NodeBase<Tdim>::initialise, std::placeholders::_1));

    meshes_.at(0)->iterate_over_cells(
        std::bind(&mpm::Cell<Tdim>::activate_nodes, std::placeholders::_1));

    // Iterate over each particle to compute shapefn
    meshes_.at(0)->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::compute_shapefn, std::placeholders::_1));

    // Compute volume
    meshes_.at(0)->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::compute_volume, std::placeholders::_1));

    // Compute mass
    meshes_.at(0)->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::compute_mass, std::placeholders::_1, phase));
    // Assign mass and momentum to nodes
    meshes_.at(0)->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_mass_momentum_to_nodes,
                  std::placeholders::_1, phase));

    // Compute nodal velocity
    meshes_.at(0)->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_velocity,
                  std::placeholders::_1),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Iterate over each particle to calculate strain
    meshes_.at(0)->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::compute_strain,
                  std::placeholders::_1, phase, dt_));

    // Iterate over each particle to compute stress
    meshes_.at(0)->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::compute_stress,
                  std::placeholders::_1, phase));

    // Iterate over each particle to compute nodal body force
    meshes_.at(0)->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_body_force,
                  std::placeholders::_1, phase, this->gravity_));

    // Iterate over each particle to compute nodal internal force
    meshes_.at(0)->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_internal_force,
                  std::placeholders::_1, phase));

    // Iterate over active nodes to compute acceleratation and velocity
    meshes_.at(0)->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_acceleration_velocity,
                  std::placeholders::_1, phase, this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Iterate over each particle to compute updated position
    meshes_.at(0)->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::compute_updated_position,
                  std::placeholders::_1, phase, this->dt_));

    // Locate particles
    auto unlocatable_particles = meshes_.at(0)->locate_particles_mesh();

    if (!unlocatable_particles.empty())
      throw std::runtime_error("Particle outside the mesh domain");

    // Stats
    // TODO: Remove
    // Iterate over each particle stats
    meshes_.at(0)->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::stats, std::placeholders::_1));

    // TODO: Remove stats
    /*
    std::cout << "After: \n";
    meshes_.at(0)->iterate_over_nodes(
        std::bind(&mpm::NodeBase<Tdim>::stats, std::placeholders::_1));
    */
    meshes_.at(0)->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::stats, std::placeholders::_1),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));
  }
  return status;
}

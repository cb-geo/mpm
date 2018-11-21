//! Constructor
template <unsigned Tdim>
mpm::MPMExplicit<Tdim>::MPMExplicit(std::unique_ptr<IO>&& io)
    : mpm::MPM(std::move(io)) {
  //! Logger
  console_ = spdlog::get("MPMExplicit");

  // Create a mesh with global id 0
  const mpm::Index id = 0;
  // Set analysis step to start at 0
  step_ = 0;

  // Set mesh as isoparametric
  bool isoparametric = true;

  try {
    const auto mesh_props = io_->json_object("mesh");
    isoparametric = mesh_props["isoparametric"].template get<bool>();
  } catch (std::exception& exception) {
    console_->warn(
        "{} {} Isoparametric status of mesh: {}\n Setting mesh as "
        "isoparametric.",
        __FILE__, __LINE__, exception.what());
    isoparametric = true;
  }

  mesh_ = std::make_unique<mpm::Mesh<Tdim>>(id, isoparametric);

  // Empty all materials
  materials_.clear();

  try {
    analysis_ = io_->analysis();
    // Time-step size
    dt_ = analysis_["dt"].template get<double>();
    // Number of time steps
    nsteps_ = analysis_["nsteps"].template get<mpm::Index>();

    if (analysis_.at("gravity").is_array() &&
        analysis_.at("gravity").size() == gravity_.size()) {
      for (unsigned i = 0; i < gravity_.size(); ++i) {
        gravity_[i] = analysis_.at("gravity").at(i);
      }
    } else {
      throw std::runtime_error("Specified gravity dimension is invalid");
    }

    // Velocity update
    try {
      velocity_update_ = analysis_["velocity_update"].template get<bool>();
    } catch (std::exception& exception) {
      console_->warn(
          "{} #{}: Velocity update parameter is not specified, using default "
          "as false",
          __FILE__, __LINE__, exception.what());
      velocity_update_ = false;
    }

    post_process_ = io_->post_processing();
    // Output steps
    output_steps_ = post_process_["output_steps"].template get<mpm::Index>();

  } catch (std::domain_error& domain_error) {
    console_->error("{} {} Get analysis object: {}", __FILE__, __LINE__,
                    domain_error.what());
    abort();
  }

  // Default VTK attributes
  std::vector<std::string> vtk = {"velocities", "stresses", "strains"};
  try {
    if (post_process_.at("vtk").is_array() &&
        post_process_.at("vtk").size() > 0) {
      for (unsigned i = 0; i < post_process_.at("vtk").size(); ++i) {
        std::string attribute =
            post_process_["vtk"][i].template get<std::string>();
        if (std::find(vtk.begin(), vtk.end(), attribute) != vtk.end())
          vtk_attributes_.emplace_back(attribute);
        else
          throw std::runtime_error("Specificed VTK argument is incorrect");
      }
    } else {
      throw std::runtime_error(
          "Specificed VTK arguments are incorrect, using defaults");
    }
  } catch (std::exception& exception) {
    vtk_attributes_ = vtk;
    console_->warn("{} {}: {}", __FILE__, __LINE__, exception.what());
  }
}

// Initialise mesh
template <unsigned Tdim>
bool mpm::MPMExplicit<Tdim>::initialise_mesh() {
  // TODO: Fix phase
  const unsigned phase = 0;
  bool status = true;

  try {
    // Initialise MPI rank and size
    int mpi_rank = 0;
    int mpi_size = 1;

#ifdef USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    // Get number of MPI ranks
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif

    // Get mesh properties
    auto mesh_props = io_->json_object("mesh");
    // Get Mesh reader from JSON object
    const std::string reader =
        mesh_props["mesh_reader"].template get<std::string>();

    bool check_duplicates = true;
    try {
      check_duplicates = mesh_props["check_duplicates"].template get<bool>();
    } catch (std::exception& exception) {
      console_->warn(
          "{} #{}: Check duplicates, not specified setting default as true",
          __FILE__, __LINE__, exception.what());
      check_duplicates = true;
    }

    // Create a mesh reader
    auto mesh_reader = Factory<mpm::ReadMesh<Tdim>>::instance()->create(reader);

    auto nodes_begin = std::chrono::steady_clock::now();
    // Global Index
    mpm::Index gid = 0;
    // Node type
    const auto node_type = mesh_props["node_type"].template get<std::string>();
    // Create nodes from file
    bool node_status = mesh_->create_nodes(
        gid,                                                   // global id
        node_type,                                             // node type
        mesh_reader->read_mesh_nodes(io_->file_name("mesh")),  // coordinates
        check_duplicates);  // check duplicates

    if (!node_status)
      throw std::runtime_error("Addition of nodes to mesh failed");

    auto nodes_end = std::chrono::steady_clock::now();
    console_->info("Rank {} Read nodes: {} ms", mpi_rank,
                   std::chrono::duration_cast<std::chrono::milliseconds>(
                       nodes_end - nodes_begin)
                       .count());

    // Read and assign velocity constraints
    if (!io_->file_name("velocity_constraints").empty()) {
      bool velocity_constraints = mesh_->assign_velocity_constraints(
          mesh_reader->read_velocity_constraints(
              io_->file_name("velocity_constraints")));
      if (!velocity_constraints)
        throw std::runtime_error(
            "Velocity constraints are not properly assigned");
    }

    // Set nodal traction as false if file is empty
    if (io_->file_name("nodal_tractions").empty()) nodal_tractions_ = false;

    auto cells_begin = std::chrono::steady_clock::now();
    // Shape function name
    const auto cell_type = mesh_props["cell_type"].template get<std::string>();
    // Shape function
    std::shared_ptr<mpm::Element<Tdim>> element =
        Factory<mpm::Element<Tdim>>::instance()->create(cell_type);

    // Create cells from file
    bool cell_status = mesh_->create_cells(
        gid,                                                   // global id
        element,                                               // element tyep
        mesh_reader->read_mesh_cells(io_->file_name("mesh")),  // Node ids
        check_duplicates);  // Check duplicates

    if (!cell_status)
      throw std::runtime_error("Addition of cells to mesh failed");

    auto cells_end = std::chrono::steady_clock::now();
    console_->info("Rank {} Read cells: {} ms", mpi_rank,
                   std::chrono::duration_cast<std::chrono::milliseconds>(
                       cells_end - cells_begin)
                       .count());

  } catch (std::exception& exception) {
    console_->error("#{}: Reading mesh and particles: {}", __LINE__,
                    exception.what());
    status = false;
  }
  return status;
}

// Initialise particles
template <unsigned Tdim>
bool mpm::MPMExplicit<Tdim>::initialise_particles() {
  // TODO: Fix phase
  const unsigned phase = 0;
  bool status = true;

  try {
    // Initialise MPI rank and size
    int mpi_rank = 0;
    int mpi_size = 1;

#ifdef USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    // Get number of MPI ranks
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif

    // Get mesh properties
    auto mesh_props = io_->json_object("mesh");
    // Get Mesh reader from JSON object
    const std::string reader =
        mesh_props["mesh_reader"].template get<std::string>();

    bool check_duplicates = true;
    try {
      check_duplicates = mesh_props["check_duplicates"].template get<bool>();
    } catch (std::exception& exception) {
      console_->warn(
          "{} #{}: Check duplicates, not specified setting default as true",
          __FILE__, __LINE__, exception.what());
      check_duplicates = true;
    }

    // Create a mesh reader
    auto particle_reader =
        Factory<mpm::ReadMesh<Tdim>>::instance()->create(reader);

    auto particles_begin = std::chrono::steady_clock::now();

    // Get all particles
    std::vector<Eigen::Matrix<double, Tdim, 1>> all_particles;

    // Generate particles
    bool read_particles_file = false;
    try {
      unsigned nparticles_cell =
          mesh_props["generate_particles_cells"].template get<unsigned>();

      if (nparticles_cell > 0)
        all_particles = mesh_->generate_material_points(nparticles_cell);
      else
        throw std::runtime_error(
            "Specified # of particles per cell for generation is invalid!");

    } catch (std::exception& exception) {
      console_->warn("Generate particles is not set, reading particles file");
      read_particles_file = true;
    }

    // Read particles from file
    if (read_particles_file)
      all_particles =
          particle_reader->read_particles(io_->file_name("particles"));
    // Get all particle ids
    std::vector<mpm::Index> all_particles_ids(all_particles.size());
    std::iota(all_particles_ids.begin(), all_particles_ids.end(), 0);

    // Get local particles chunk
    std::vector<Eigen::Matrix<double, Tdim, 1>> particles;
    chunk_vector_quantities(all_particles, particles);

    // Get local particles ids chunks
    std::vector<mpm::Index> particles_ids;
    chunk_scalar_quantities(all_particles_ids, particles_ids);

    // Particle type
    const auto particle_type =
        mesh_props["particle_type"].template get<std::string>();

    // Create particles from file
    bool particle_status =
        mesh_->create_particles(particles_ids,      // global id
                                particle_type,      // particle type
                                particles,          // coordinates
                                check_duplicates);  // Check duplicates

    if (!particle_status)
      throw std::runtime_error("Addition of particles to mesh failed");

    auto particles_end = std::chrono::steady_clock::now();
    console_->info("Rank {} Read particles: {} ms", mpi_rank,
                   std::chrono::duration_cast<std::chrono::milliseconds>(
                       particles_end - particles_begin)
                       .count());
    try {
      // Read and assign particles cells
      if (!io_->file_name("particles_cells").empty()) {
        bool particles_cells =
            mesh_->assign_particles_cells(particle_reader->read_particles_cells(
                io_->file_name("particles_cells")));
        if (!particles_cells)
          throw std::runtime_error(
              "Cell ids are not properly assigned to particles");
      }
    } catch (std::exception& exception) {
      console_->error("{} #{}: Reading particles cells: {}", __FILE__, __LINE__,
                      exception.what());
    }

    auto particles_locate_begin = std::chrono::steady_clock::now();
    // Locate particles in cell
    auto unlocatable_particles = mesh_->locate_particles_mesh();

    if (!unlocatable_particles.empty())
      throw std::runtime_error("Particle outside the mesh domain");

    auto particles_locate_end = std::chrono::steady_clock::now();
    console_->info("Rank {} Locate particles: {} ms", mpi_rank,
                   std::chrono::duration_cast<std::chrono::milliseconds>(
                       particles_locate_end - particles_locate_begin)
                       .count());

    // Write particles and cells to file
    particle_reader->write_particles_cells(
        io_->output_file("particles-cells", ".txt", uuid_, 0, 0).string(),
        mesh_->particles_cells());

    auto particles_traction_begin = std::chrono::steady_clock::now();
    // Compute volume
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::compute_volume,
                  std::placeholders::_1, phase));

    // Read and assign particles volumes
    if (!io_->file_name("particles_volumes").empty()) {
      bool particles_volumes = mesh_->assign_particles_volumes(
          particle_reader->read_particles_volumes(
              io_->file_name("particles_volumes")));
      if (!particles_volumes)
        throw std::runtime_error("Particles volumes are not properly assigned");
    }

    // Read and assign particles tractions
    if (!io_->file_name("particles_tractions").empty()) {
      bool particles_tractions = mesh_->assign_particles_tractions(
          particle_reader->read_particles_tractions(
              io_->file_name("particles_tractions")));
      if (!particles_tractions)
        throw std::runtime_error(
            "Particles tractions are not properly assigned");
    }

    // Read and assign particles stresses
    if (!io_->file_name("particles_stresses").empty()) {

      // Get stresses of all particles
      const auto all_particles_stresses =
          particle_reader->read_particles_stresses(
              io_->file_name("particles_stresses"));
      // Chunked stresses
      std::vector<Eigen::Matrix<double, 6, 1>> particles_stresses;
      chunk_vector_quantities(all_particles_stresses, particles_stresses);

      // Read and assign particles stresses
      if (!mesh_->assign_particles_stresses(particles_stresses))
        throw std::runtime_error(
            "Particles stresses are not properly assigned");
    }

    auto particles_traction_end = std::chrono::steady_clock::now();
    console_->info("Rank {} Read particle traction and stresses: {} ms",
                   mpi_rank,
                   std::chrono::duration_cast<std::chrono::milliseconds>(
                       particles_traction_end - particles_traction_begin)
                       .count());

  } catch (std::exception& exception) {
    console_->error("#{}: Reading particles: {}", __LINE__, exception.what());
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
      auto material_id = material_props["id"].template get<unsigned>();

      // Create a new material from JSON object
      auto mat =
          Factory<mpm::Material<Tdim>, unsigned, const Json&>::instance()
              ->create(material_type, std::move(material_id), material_props);

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

//! Apply nodal tractions
template <unsigned Tdim>
bool mpm::MPMExplicit<Tdim>::apply_nodal_tractions() {
  bool status = true;
  try {
    // Read and assign nodes tractions
    if (!io_->file_name("nodal_tractions").empty()) {
      // Get mesh properties
      auto mesh_props = io_->json_object("mesh");
      // Get Mesh reader from JSON object
      const std::string reader =
          mesh_props["mesh_reader"].template get<std::string>();
      // Create a mesh reader
      auto node_reader =
          Factory<mpm::ReadMesh<Tdim>>::instance()->create(reader);

      bool nodal_tractions =
          mesh_->assign_nodal_tractions(node_reader->read_particles_tractions(
              io_->file_name("nodal_tractions")));
      if (!nodal_tractions)
        throw std::runtime_error("Nodal tractions are not properly assigned");
    } else
      nodal_tractions_ = false;
  } catch (std::exception& exception) {
    console_->error("#{}: Nodal traction: {}", __LINE__, exception.what());
    status = false;
    nodal_tractions_ = false;
  }
  return status;
}

//! Checkpoint resume
template <unsigned Tdim>
bool mpm::MPMExplicit<Tdim>::checkpoint_resume() {
  bool checkpoint = true;
  try {
    // TODO: Set phase
    const unsigned phase = 0;

    if (!analysis_["resume"]["resume"].template get<bool>())
      throw std::runtime_error("Resume analysis option is disabled!");

    // Get unique analysis id
    this->uuid_ = analysis_["resume"]["uuid"].template get<std::string>();
    // Get step
    this->step_ = analysis_["resume"]["step"].template get<mpm::Index>();

    // Input particle h5 file for resume
    std::string attribute = "particles";
    std::string extension = ".h5";

    auto particles_file =
        io_->output_file(attribute, extension, uuid_, step_, this->nsteps_)
            .string();
    // Load particle information from file
    mesh_->read_particles_hdf5(phase, particles_file);
    // Locate particles
    auto unlocatable_particles = mesh_->locate_particles_mesh();

    if (!unlocatable_particles.empty())
      throw std::runtime_error("Particle outside the mesh domain");

    // Increament step
    ++this->step_;

    console_->info("Checkpoint resume at step {} of {}", this->step_,
                   this->nsteps_);

  } catch (std::exception& exception) {
    console_->info("{} {} Resume failed, restarting analysis: {}", __FILE__,
                   __LINE__, exception.what());
    this->step_ = 0;
    checkpoint = false;
  }
  return checkpoint;
}

//! Write HDF5 files
template <unsigned Tdim>
void mpm::MPMExplicit<Tdim>::write_hdf5(mpm::Index step, mpm::Index max_steps) {
  // Write input geometry to vtk file
  std::string attribute = "particles";
  std::string extension = ".h5";

  auto particles_file =
      io_->output_file(attribute, extension, uuid_, step, max_steps).string();

  const unsigned phase = 0;
  mesh_->write_particles_hdf5(phase, particles_file);
}

#ifdef USE_VTK
//! Write VTK files
template <unsigned Tdim>
void mpm::MPMExplicit<Tdim>::write_vtk(mpm::Index step, mpm::Index max_steps) {

  // VTK PolyData writer
  auto vtk_writer = std::make_unique<VtkWriter>(mesh_->particle_coordinates());

  // Write mesh on step 0
  if (step == 0)
    vtk_writer->write_mesh(
        io_->output_file("mesh", ".vtp", uuid_, step, max_steps).string(),
        mesh_->nodal_coordinates(), mesh_->node_pairs());

  // Write input geometry to vtk file
  const std::string extension = ".vtp";
  const std::string attribute = "geometry";
  auto meshfile =
      io_->output_file(attribute, extension, uuid_, step, max_steps).string();
  vtk_writer->write_geometry(meshfile);

  // TODO fix phase
  unsigned phase = 0;

  for (const auto& attribute : vtk_attributes_) {
    // Write vector
    auto file =
        io_->output_file(attribute, extension, uuid_, step, max_steps).string();
    vtk_writer->write_vector_point_data(
        file, mesh_->particles_vector_data(attribute, phase), attribute);
  }
}
#endif

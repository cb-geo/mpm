//! Constructor
template <unsigned Tdim>
mpm::MPMBase<Tdim>::MPMBase(std::unique_ptr<IO>&& io)
    : mpm::MPM(std::move(io)) {
  //! Logger
  console_ = spdlog::get("MPMBase");

  // Create a mesh with global id 0
  const mpm::Index id = 0;

  // Set analysis step to start at 0
  step_ = 0;

  // Set mesh as isoparametric
  bool isoparametric = is_isoparametric();

  mesh_ = std::make_unique<mpm::Mesh<Tdim>>(id, isoparametric);

  // Empty all materials
  materials_.clear();

  try {
    analysis_ = io_->analysis();
    // Time-step size
    dt_ = analysis_["dt"].template get<double>();
    // Number of time steps
    nsteps_ = analysis_["nsteps"].template get<mpm::Index>();

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

    // Math functions
    try {
      if (analysis_["math_functions"].template get<bool>())
        this->initialise_math_functions();
    } catch (std::exception& exception) {
      console_->warn("{} #{}: No math functions are defined; set to default",
                     __FILE__, __LINE__, exception.what());
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
bool mpm::MPMBase<Tdim>::initialise_mesh() {
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

    // Read nodal euler angles and assign rotation matrices
    if (!io_->file_name("nodal_euler_angles").empty()) {
      bool rotation_matrices = mesh_->compute_nodal_rotation_matrices(
          mesh_reader->read_euler_angles(io_->file_name("nodal_euler_angles")));
      if (!rotation_matrices)
        throw std::runtime_error(
            "Euler angles are not properly assigned/computed");
    }

    // Read and assign velocity constraints
    if (!io_->file_name("velocity_constraints").empty()) {
      bool velocity_constraints = mesh_->assign_velocity_constraints(
          mesh_reader->read_velocity_constraints(
              io_->file_name("velocity_constraints")));
      if (!velocity_constraints)
        throw std::runtime_error(
            "Velocity constraints are not properly assigned");
    }

    // Read and assign friction constraints
    if (!io_->file_name("friction_constraints").empty()) {
      bool friction_constraints = mesh_->assign_friction_constraints(
          mesh_reader->read_friction_constraints(
              io_->file_name("friction_constraints")));
      if (!friction_constraints)
        throw std::runtime_error(
            "Friction constraints are not properly assigned");
    }

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

    // Compute cell neighbours
    mesh_->compute_cell_neighbours();

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
bool mpm::MPMBase<Tdim>::initialise_particles() {
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

    // Get particle properties
    auto particle_props = io_->json_object("particle");
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
        particle_props["particle_type"].template get<std::string>();

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

    // Read and assign particles velocity constraints
    if (!io_->file_name("particles_velocity_constraints").empty()) {
      bool particles_velocity_constraints =
          mesh_->assign_particles_velocity_constraints(
              particle_reader->read_velocity_constraints(
                  io_->file_name("particles_velocity_constraints")));
      if (!particles_velocity_constraints)
        throw std::runtime_error(
            "Particles velocity constraints are not properly assigned");
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

    // Read and assign particle sets
    if (!io_->file_name("entity_sets").empty()) {
      bool particle_sets = mesh_->create_particle_sets(
          (io_->entity_sets(io_->file_name("entity_sets"), "particle_sets")),
          check_duplicates);
    }
  } catch (std::exception& exception) {
    console_->error("#{}: Reading particles: {}", __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Initialise materials
template <unsigned Tdim>
bool mpm::MPMBase<Tdim>::initialise_materials() {
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

//! Apply properties to particles sets (e.g: material)
template <unsigned Tdim>
bool mpm::MPMBase<Tdim>::apply_properties_to_particles_sets() {
  bool status = false;
  // Set phase to zero
  unsigned phase = 0;
  // Assign material to particle sets
  try {
    // Get particle properties
    auto particle_props = io_->json_object("particle");
    // Get particle sets properties
    auto particle_sets = particle_props["particle_sets"];
    // Assign material to each particle sets
    for (const auto& psets : particle_sets) {
      // Get set material from list of materials
      auto set_material = materials_.at(psets["material_id"]);
      // Get sets ids
      std::vector<unsigned> sids = psets["set_id"];
      // Assign material to particles in the specific sets
      for (const auto& sitr : sids) {
        mesh_->iterate_over_particle_set(
            sitr, std::bind(&mpm::ParticleBase<Tdim>::assign_material,
                            std::placeholders::_1, phase, set_material));
      }
    }
    status = true;
  } catch (std::exception& exception) {
    console_->error("#{}: Particle sets material: {}", __LINE__,
                    exception.what());
  }
  return status;
}

//! Checkpoint resume
template <unsigned Tdim>
bool mpm::MPMBase<Tdim>::checkpoint_resume() {
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

    // Clear all particle ids
    mesh_->iterate_over_cells(
        std::bind(&mpm::Cell<Tdim>::clear_particle_ids, std::placeholders::_1));

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
void mpm::MPMBase<Tdim>::write_hdf5(mpm::Index step, mpm::Index max_steps) {
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
void mpm::MPMBase<Tdim>::write_vtk(mpm::Index step, mpm::Index max_steps) {

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

//! Return if a mesh is isoparametric
template <unsigned Tdim>
bool mpm::MPMBase<Tdim>::is_isoparametric() {
  bool isoparametric = true;

  try {
    const auto mesh_props = io_->json_object("mesh");
    isoparametric = mesh_props.at("isoparametric").template get<bool>();
  } catch (std::exception& exception) {
    console_->warn(
        "{} {} Isoparametric status of mesh: {}\n Setting mesh as "
        "isoparametric.",
        __FILE__, __LINE__, exception.what());
    isoparametric = true;
  }
  return isoparametric;
}

//! Initialise loads
template <unsigned Tdim>
bool mpm::MPMBase<Tdim>::initialise_loads() {
  bool status = true;
  try {
    auto loads = io_->json_object("external_loading_conditions");
    // Initialise gravity loading
    if (loads.at("gravity").is_array() &&
        loads.at("gravity").size() == gravity_.size()) {
      for (unsigned i = 0; i < gravity_.size(); ++i) {
        gravity_[i] = loads.at("gravity").at(i);
      }
    } else {
      throw std::runtime_error("Specified gravity dimension is invalid");
    }

    // Create a file reader
    const std::string reader =
        io_->json_object("mesh")["mesh_reader"].template get<std::string>();
    auto traction_reader =
        Factory<mpm::ReadMesh<Tdim>>::instance()->create(reader);

    // Read and assign particles surface tractions
    if (loads.find("particle_surface_traction") != loads.end()) {
      for (const auto& ptraction : loads["particle_surface_traction"]) {
        // Get the math function
        std::shared_ptr<FunctionBase> tfunction = nullptr;
        if(ptraction.find("math_function_id") != ptraction.end())
          tfunction = math_functions_.at(
              ptraction["math_function_id"].template get<unsigned>());
        // Read and assign particle surface tractions
        bool particles_tractions = mesh_->assign_particles_tractions(
            tfunction,
            traction_reader->read_tractions(
                io_->working_directory() +
                ptraction["input_file"].template get<std::string>()));
        if (!particles_tractions)
          throw std::runtime_error(
              "Particles tractions are not properly assigned");
      }
    } else
      console_->warn(
          "No particle surface traction is defined for the analysis");

    // Read and assign nodal concentrated forces
    if (loads.find("concentrated_nodal_forces") != loads.end()) {
      for (const auto& nforce : loads["concentrated_nodal_forces"]) {
        // Get the math function
        std::shared_ptr<FunctionBase> ffunction = nullptr;
        if (nforce.find("math_function_id") != nforce.end())
          ffunction = math_functions_.at(
              nforce["math_function_id"].template get<unsigned>());
        // Read and assign nodal concentrated forces
        bool nodal_force = mesh_->assign_nodal_concentrated_forces(
            ffunction, traction_reader->read_tractions(
                           io_->working_directory() +
                           nforce["input_file"].template get<std::string>()));
        if (!nodal_force)
          throw std::runtime_error(
              "Concentrated nodal forces are not properly assigned");
        set_node_concentrated_force_ = true;
      }
    } else
      console_->warn("No concentrated nodal force is defined for the analysis");
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Initialise math functions
template <unsigned Tdim>
bool mpm::MPMBase<Tdim>::initialise_math_functions() {
  bool status = true;
  try {
    // Get materials properties
    auto math_functions = io_->json_object("math_functions");
    for (const auto function_props : math_functions) {

      // Get math function id
      auto function_id = function_props["id"].template get<unsigned>();
      
      // Get function type
      const std::string function_type =
          function_props["type"].template get<std::string>();

      // Create a new function from JSON object
      auto function =
          Factory<mpm::FunctionBase, unsigned, const Json&>::instance()->create(
              function_type, std::move(function_id), function_props);

      // Add material to list
      auto insert_status =
          math_functions_.insert(std::make_pair(function->id(), function));

      // If insert material failed
      if (!insert_status.second) {
        status = false;
        throw std::runtime_error(
            "New math function cannot be added, insertion failed");
      }
    }
  } catch (std::exception& exception) {
    console_->error("#{}: Reading math functions: {}", __LINE__,
                    exception.what());
    status = false;
  }
  return status;
}

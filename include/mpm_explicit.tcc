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
  // Clear meshes
  meshes_.clear();
  meshes_.emplace_back(std::make_unique<mpm::Mesh<Tdim>>(id));

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

// Initialise mesh and particles
template <unsigned Tdim>
bool mpm::MPMExplicit<Tdim>::initialise_mesh_particles() {
  // TODO: Fix phase
  const unsigned phase = 0;
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

    // Read and assign velocity constraints
    if (!io_->file_name("velocity_constraints").empty()) {
      bool velocity_constraints = meshes_.at(0)->assign_velocity_constraints(
          mesh_reader->read_velocity_constraints(
              io_->file_name("velocity_constraints")));
      if (!velocity_constraints)
        throw std::runtime_error(
            "Velocity constraints are not properly assigned");
    }

    // Shape function name
    const auto cell_type = mesh_props["cell_type"].template get<std::string>();
    // Shape function
    std::shared_ptr<mpm::Element<Tdim>> element =
        Factory<mpm::Element<Tdim>>::instance()->create(cell_type);

    // Create cells from file
    bool cell_status = meshes_.at(0)->create_cells(
        gid,                                                    // global id
        element,                                                // element tyep
        mesh_reader->read_mesh_cells(io_->file_name("mesh")));  // Node ids

    if (!cell_status)
      throw std::runtime_error("Addition of cells to mesh failed");

#ifdef USE_MPI
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    // Get number of MPI ranks
    int mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    // Create MPI array type
    MPI_Datatype array_t;
    MPI_Type_vector(Tdim, 1, 1, MPI_DOUBLE, &array_t);
    MPI_Type_commit(&array_t);

    auto all_particles =
        mesh_reader->read_particles(io_->file_name("particles"));

    // Calculate chunk size to split router
    int chunk_size = all_particles.size() / mpi_size;
    MPI_Bcast(&chunk_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    std::vector<Eigen::Matrix<double, Tdim, 1>> particles(chunk_size);

    // Send particle chunks to different compute nodes
    MPI_Scatter(all_particles.data(), chunk_size, array_t, particles.data(),
                particles.size(), array_t, 0, MPI_COMM_WORLD);

    // Calculate the remaining chunk of od_pairs and add to rank 0
    int chunk_remainder = all_particles.size() % mpi_size;
    if (mpi_rank == 0) {
      particles.insert(particles.begin(), all_particles.end() - chunk_remainder,
                       all_particles.end());
    }

    // Particle type
    const auto particle_type =
        mesh_props["particle_type"].template get<std::string>();
    // Create particles from file
    bool particle_status =
        meshes_.at(0)->create_particles(gid,            // global id
                                        particle_type,  // particle type
                                        particles);     // coordinates

    if (!particle_status)
      throw std::runtime_error("Addition of particles to mesh failed");

    // Locate particles in cell
    auto unlocatable_particles = meshes_.at(0)->locate_particles_mesh();

    if (!unlocatable_particles.empty())
      throw std::runtime_error("Particle outside the mesh domain");
#else
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

#endif

    // Compute volume
    meshes_.at(0)->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::compute_volume,
                  std::placeholders::_1, phase));

    // Read and assign particles tractions
    if (!io_->file_name("particles_tractions").empty()) {
      bool particles_tractions = meshes_.at(0)->assign_particles_tractions(
          mesh_reader->read_particles_tractions(
              io_->file_name("particles_tractions")));
      if (!particles_tractions)
        throw std::runtime_error(
            "Particles tractions are not properly assigned");
    }

    // Read and assign particles stresses
    if (!io_->file_name("particles_stresses").empty()) {
      bool particles_stresses = meshes_.at(0)->assign_particles_stresses(
          mesh_reader->read_particles_stresses(
              io_->file_name("particles_stresses")));
      if (!particles_stresses)
        throw std::runtime_error(
            "Particles stresses are not properly assigned");
    }

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
      auto material_id = material_props["id"].template get<unsigned>();

      // Create a new material from JSON object
      auto mat = Factory<mpm::Material<Tdim>, unsigned>::instance()->create(
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
    meshes_.at(0)->read_particles_hdf5(phase, particles_file);
    // Locate particles
    auto unlocatable_particles = meshes_.at(0)->locate_particles_mesh();

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
  meshes_.at(0)->write_particles_hdf5(phase, particles_file);
}

#ifdef USE_VTK
//! Write VTK files
template <unsigned Tdim>
void mpm::MPMExplicit<Tdim>::write_vtk(mpm::Index step, mpm::Index max_steps) {
  const auto coordinates = meshes_.at(0)->particle_coordinates();
  // VTK PolyData writer
  auto vtk_writer = std::make_unique<VtkWriter>(coordinates);

  // Write input geometry to vtk file
  std::string extension = ".vtp";
  std::string attribute = "geometry";
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
        file, meshes_.at(0)->particles_vector_data(attribute, phase),
        attribute);
  }
}
#endif

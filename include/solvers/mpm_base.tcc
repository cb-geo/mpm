//! Constructor
template <unsigned Tdim>
mpm::MPMBase<Tdim>::MPMBase(const std::shared_ptr<IO>& io) : mpm::MPM(io) {
  //! Logger
  console_ = spdlog::get("MPMBase");

  // Create a mesh with global id 0
  const mpm::Index id = 0;

  // Set analysis step to start at 0
  step_ = 0;

  // Set mesh as isoparametric
  bool isoparametric = is_isoparametric();

  mesh_ = std::make_shared<mpm::Mesh<Tdim>>(id, isoparametric);

  // Create constraints
  constraints_ = std::make_shared<mpm::Constraints<Tdim>>(mesh_);

  // Empty all materials
  materials_.clear();

  try {
    analysis_ = io_->analysis();
    // Time-step size
    dt_ = analysis_["dt"].template get<double>();
    // Number of time steps
    nsteps_ = analysis_["nsteps"].template get<mpm::Index>();

    // nload balance
    if (analysis_.find("nload_balance_steps") != analysis_.end())
      nload_balance_steps_ =
          analysis_["nload_balance_steps"].template get<mpm::Index>();

    // Locate particles
    if (analysis_.find("locate_particles") != analysis_.end())
      locate_particles_ = analysis_["locate_particles"].template get<bool>();

    // Stress update method (USF/USL/MUSL)
    try {
      if (analysis_.find("stress_update") != analysis_.end())
        stress_update_ = mpm::stress_update.at(
            analysis_["stress_update"].template get<std::string>());
    } catch (std::exception& exception) {
      console_->warn(
          "{} #{}: {}. Stress update method is not specified, using USF as "
          "default\n",
          __FILE__, __LINE__, exception.what());
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

    // Damping
    try {
      if (analysis_.find("damping") != analysis_.end()) {
        if (!initialise_damping(analysis_.at("damping")))
          throw std::runtime_error("Damping parameters are not defined");
      }
    } catch (std::exception& exception) {
      console_->warn("{} #{}: Damping is not specified, using none as default",
                     __FILE__, __LINE__, exception.what());
    }

    // Math functions
    try {
      // Get materials properties
      auto math_functions = io_->json_object("math_functions");
      if (!math_functions.empty())
        this->initialise_math_functions(math_functions);
    } catch (std::exception& exception) {
      console_->warn("{} #{}: No math functions are defined", __FILE__,
                     __LINE__, exception.what());
    }

    post_process_ = io_->post_processing();
    // Output steps
    output_steps_ = post_process_["output_steps"].template get<mpm::Index>();

  } catch (std::domain_error& domain_error) {
    console_->error("{} {} Get analysis object: {}", __FILE__, __LINE__,
                    domain_error.what());
    abort();
  }

  // VTK state variables
  try {
    if (post_process_.at("vtk_statevars").is_array() &&
        post_process_.at("vtk_statevars").size() > 0) {
      for (unsigned i = 0; i < post_process_.at("vtk_statevars").size(); ++i) {
        std::string attribute =
            post_process_["vtk_statevars"][i].template get<std::string>();
        vtk_statevars_.emplace_back(attribute);
      }
    } else {
      throw std::runtime_error(
          "No VTK statevariable were specified, none will be generated");
    }
  } catch (std::exception& exception) {
    console_->warn("{} {}: {}", __FILE__, __LINE__, exception.what());
  }
}

// Initialise mesh
template <unsigned Tdim>
bool mpm::MPMBase<Tdim>::initialise_mesh() {
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
    const std::string io_type =
        mesh_props["io_type"].template get<std::string>();

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
    auto mesh_io = Factory<mpm::IOMesh<Tdim>>::instance()->create(io_type);

    auto nodes_begin = std::chrono::steady_clock::now();
    // Global Index
    mpm::Index gid = 0;
    // Node type
    const auto node_type = mesh_props["node_type"].template get<std::string>();

    // Mesh file
    std::string mesh_file =
        io_->file_name(mesh_props["mesh"].template get<std::string>());

    // Create nodes from file
    bool node_status =
        mesh_->create_nodes(gid,                                  // global id
                            node_type,                            // node type
                            mesh_io->read_mesh_nodes(mesh_file),  // coordinates
                            check_duplicates);                    // check dups

    if (!node_status) {
      status = false;
      throw std::runtime_error("Addition of nodes to mesh failed");
    }

    auto nodes_end = std::chrono::steady_clock::now();
    console_->info("Rank {} Read nodes: {} ms", mpi_rank,
                   std::chrono::duration_cast<std::chrono::milliseconds>(
                       nodes_end - nodes_begin)
                       .count());

    // Read and assign node sets
    this->node_entity_sets(mesh_props, check_duplicates);

    // Read nodal euler angles and assign rotation matrices
    this->node_euler_angles(mesh_props, mesh_io);

    // Read and assign velocity constraints
    this->nodal_velocity_constraints(mesh_props, mesh_io);

    // Read and assign friction constraints
    this->nodal_frictional_constraints(mesh_props, mesh_io);

    // Initialise cell
    auto cells_begin = std::chrono::steady_clock::now();
    // Shape function name
    const auto cell_type = mesh_props["cell_type"].template get<std::string>();
    // Shape function
    std::shared_ptr<mpm::Element<Tdim>> element =
        Factory<mpm::Element<Tdim>>::instance()->create(cell_type);

    // Create cells from file
    bool cell_status =
        mesh_->create_cells(gid,      // global id
                            element,  // element tyep
                            mesh_io->read_mesh_cells(mesh_file),  // Node ids
                            check_duplicates);                    // Check dups

    if (!cell_status) {
      status = false;
      throw std::runtime_error("Addition of cells to mesh failed");
    }

    // Compute cell neighbours
    mesh_->find_cell_neighbours();

    // Read and assign cell sets
    this->cell_entity_sets(mesh_props, check_duplicates);

    auto cells_end = std::chrono::steady_clock::now();
    console_->info("Rank {} Read cells: {} ms", mpi_rank,
                   std::chrono::duration_cast<std::chrono::milliseconds>(
                       cells_end - cells_begin)
                       .count());
  } catch (std::exception& exception) {
    console_->error("#{}: Reading mesh: {}", __LINE__, exception.what());
  }

  // Terminate if mesh creation failed
  if (!status) throw std::runtime_error("Initialisation of mesh failed");
  return status;
}

// Initialise particles
template <unsigned Tdim>
bool mpm::MPMBase<Tdim>::initialise_particles() {
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
    const std::string io_type =
        mesh_props["io_type"].template get<std::string>();

    bool check_duplicates = true;
    try {
      check_duplicates = mesh_props["check_duplicates"].template get<bool>();
    } catch (std::exception& exception) {
      console_->warn(
          "{} #{}: Check duplicates, not specified setting default as true",
          __FILE__, __LINE__, exception.what());
      check_duplicates = true;
    }

    auto particles_gen_begin = std::chrono::steady_clock::now();

    // Get particles properties
    auto json_particles = io_->json_object("particles");

    for (const auto& json_particle : json_particles) {
      // Generate particles
      bool gen_status =
          mesh_->generate_particles(io_, json_particle["generator"]);
      if (!gen_status) status = false;
    }

    auto particles_gen_end = std::chrono::steady_clock::now();
    console_->info("Rank {} Generate particles: {} ms", mpi_rank,
                   std::chrono::duration_cast<std::chrono::milliseconds>(
                       particles_gen_end - particles_gen_begin)
                       .count());

    auto particles_locate_begin = std::chrono::steady_clock::now();

    // Create a mesh reader
    auto particle_io = Factory<mpm::IOMesh<Tdim>>::instance()->create(io_type);

    // Read and assign particles cells
    this->particles_cells(mesh_props, particle_io);

    // Locate particles in cell
    auto unlocatable_particles = mesh_->locate_particles_mesh();

    if (!unlocatable_particles.empty())
      throw std::runtime_error("Particle outside the mesh domain");

    // Write particles and cells to file
    particle_io->write_particles_cells(
        io_->output_file("particles-cells", ".txt", uuid_, 0, 0).string(),
        mesh_->particles_cells());

    auto particles_locate_end = std::chrono::steady_clock::now();
    console_->info("Rank {} Locate particles: {} ms", mpi_rank,
                   std::chrono::duration_cast<std::chrono::milliseconds>(
                       particles_locate_end - particles_locate_begin)
                       .count());

    auto particles_volume_begin = std::chrono::steady_clock::now();
    // Compute volume
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::compute_volume, std::placeholders::_1));

    // Read and assign particles volumes
    this->particles_volumes(mesh_props, particle_io);

    // Read and assign particles stresses
    this->particles_stresses(mesh_props, particle_io);

    auto particles_volume_end = std::chrono::steady_clock::now();
    console_->info("Rank {} Read volume, velocity and stresses: {} ms",
                   mpi_rank,
                   std::chrono::duration_cast<std::chrono::milliseconds>(
                       particles_volume_end - particles_volume_begin)
                       .count());

    // Particle entity sets
    auto particles_sets_begin = std::chrono::steady_clock::now();
    this->particle_entity_sets(mesh_props, check_duplicates);
    auto particles_sets_end = std::chrono::steady_clock::now();

    // Read and assign particles velocity constraints
    this->particle_velocity_constraints(mesh_props, particle_io);

    console_->info("Rank {} Create particle sets: {} ms", mpi_rank,
                   std::chrono::duration_cast<std::chrono::milliseconds>(
                       particles_volume_end - particles_volume_begin)
                       .count());

    // Material id update using particle sets
    try {
      auto material_sets = io_->json_object("material_sets");
      if (!material_sets.empty()) {
        for (const auto& material_set : material_sets) {
          unsigned material_id =
              material_set["material_id"].template get<unsigned>();
          unsigned pset_id = material_set["pset_id"].template get<unsigned>();
          // Update material_id for particles in each pset
          mesh_->iterate_over_particle_set(
              pset_id,
              std::bind(&mpm::ParticleBase<Tdim>::assign_material,
                        std::placeholders::_1, materials_.at(material_id)));
        }
      }
    } catch (std::exception& exception) {
      console_->warn("{} #{}: Material sets are not specified", __FILE__, __LINE__,
                     exception.what());
    }

  } catch (std::exception& exception) {
    console_->error("#{}: MPM Base generating particles: {}", __LINE__,
                    exception.what());
    status = false;
  }
  if (!status) throw std::runtime_error("Initialisation of particles failed");
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
    // Copy materials to mesh
    mesh_->initialise_material_models(this->materials_);
  } catch (std::exception& exception) {
    console_->error("#{}: Reading materials: {}", __LINE__, exception.what());
    status = false;
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

    int mpi_rank = 0;
#ifdef USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif

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
  // Get active node pairs use true
  if (step % nload_balance_steps_ == 0)
    vtk_writer->write_mesh(
        io_->output_file("mesh", ".vtp", uuid_, step, max_steps).string(),
        mesh_->nodal_coordinates(), mesh_->node_pairs(true));

  // Write input geometry to vtk file
  const std::string extension = ".vtp";
  const std::string attribute = "geometry";
  auto meshfile =
      io_->output_file(attribute, extension, uuid_, step, max_steps).string();
  vtk_writer->write_geometry(meshfile);

  // MPI parallel vtk file
  int mpi_rank = 0;
  int mpi_size = 1;
  bool write_mpi_rank = false;

#ifdef USE_MPI
  // Get MPI rank
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  // Get number of MPI ranks
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif

  //! VTK vector variables
  std::vector<std::string> vtk_vector_data = {"displacements", "velocities"};

  // Write VTK attributes
  for (const auto& attribute : vtk_vector_data) {
    // Write vector
    auto file =
        io_->output_file(attribute, extension, uuid_, step, max_steps).string();
    vtk_writer->write_vector_point_data(
        file, mesh_->template particles_tensor_data<3>(attribute), attribute);

    // Write a parallel MPI VTK container file
#ifdef USE_MPI
    if (mpi_rank == 0 && mpi_size > 1) {
      auto parallel_file = io_->output_file(attribute, ".pvtp", uuid_, step,
                                            max_steps, write_mpi_rank)
                               .string();

      vtk_writer->write_parallel_vtk(parallel_file, attribute, mpi_size, step,
                                     max_steps);
    }
#endif
  }

  //! VTK tensor variables
  std::vector<std::string> vtk_tensor_data = {"stresses", "strains"};

  // Write VTK attributes
  for (const auto& attribute : vtk_tensor_data) {
    // Write vector
    auto file =
        io_->output_file(attribute, extension, uuid_, step, max_steps).string();
    vtk_writer->write_tensor_point_data(
        file, mesh_->template particles_tensor_data<6>(attribute), attribute);

    // Write a parallel MPI VTK container file
#ifdef USE_MPI
    if (mpi_rank == 0 && mpi_size > 1) {
      auto parallel_file = io_->output_file(attribute, ".pvtp", uuid_, step,
                                            max_steps, write_mpi_rank)
                               .string();

      vtk_writer->write_parallel_vtk(parallel_file, attribute, mpi_size, step,
                                     max_steps, 9);
    }
#endif
  }

  // VTK state variables
  for (const auto& attribute : vtk_statevars_) {
    // Write state variables
    auto file =
        io_->output_file(attribute, extension, uuid_, step, max_steps).string();
    vtk_writer->write_scalar_point_data(
        file, mesh_->particles_statevars_data(attribute), attribute);
    // Write a parallel MPI VTK container file
#ifdef USE_MPI
    if (mpi_rank == 0 && mpi_size > 1) {
      auto parallel_file = io_->output_file(attribute, ".pvtp", uuid_, step,
                                            max_steps, write_mpi_rank)
                               .string();
      unsigned ncomponents = 1;
      vtk_writer->write_parallel_vtk(parallel_file, attribute, mpi_size, step,
                                     max_steps, ncomponents);
    }
#endif
  }
}
#endif

#ifdef USE_PARTIO
//! Write Partio files
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::write_partio(mpm::Index step, mpm::Index max_steps) {

  // MPI parallel partio file
  int mpi_rank = 0;
  int mpi_size = 1;
#ifdef USE_MPI
  // Get MPI rank
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  // Get number of MPI ranks
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif

  // Get Partio file extensions
  const std::string extension = ".bgeo";
  const std::string attribute = "partio";
  // Create filename
  auto file =
      io_->output_file(attribute, extension, uuid_, step, max_steps).string();
  // Write partio file
  mpm::partio::write_particles(file, mesh_->particles_hdf5());
}
#endif  // USE_PARTIO

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
    const std::string io_type =
        io_->json_object("mesh")["io_type"].template get<std::string>();
    auto reader = Factory<mpm::IOMesh<Tdim>>::instance()->create(io_type);

    // Read and assign particles surface tractions
    if (loads.find("particle_surface_traction") != loads.end()) {
      for (const auto& ptraction : loads["particle_surface_traction"]) {
        // Get the math function
        std::shared_ptr<FunctionBase> tfunction = nullptr;
        // If a math function is defined set to function or use scalar
        if (ptraction.find("math_function_id") != ptraction.end())
          tfunction = math_functions_.at(
              ptraction.at("math_function_id").template get<unsigned>());
        // Set id
        int pset_id = ptraction.at("pset_id").template get<int>();
        // Direction
        unsigned dir = ptraction.at("dir").template get<unsigned>();
        // Traction
        double traction = ptraction.at("traction").template get<double>();

        // Create particle surface tractions
        bool particles_tractions = mesh_->create_particles_tractions(
            tfunction, pset_id, dir, traction);
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
        // Forces are specified in a file
        if (nforce.find("file") != nforce.end()) {
          std::string force_file =
              nforce.at("file").template get<std::string>();
          bool nodal_forces = mesh_->assign_nodal_concentrated_forces(
              reader->read_forces(io_->file_name(force_file)));
          if (!nodal_forces)
            throw std::runtime_error(
                "Nodal force file is invalid, forces are not properly "
                "assigned");
          set_node_concentrated_force_ = true;
        } else {
          // Get the math function
          std::shared_ptr<FunctionBase> ffunction = nullptr;
          if (nforce.find("math_function_id") != nforce.end())
            ffunction = math_functions_.at(
                nforce.at("math_function_id").template get<unsigned>());
          // Set id
          int nset_id = nforce.at("nset_id").template get<int>();
          // Direction
          unsigned dir = nforce.at("dir").template get<unsigned>();
          // Traction
          double force = nforce.at("force").template get<double>();

          // Read and assign nodal concentrated forces
          bool nodal_force = mesh_->assign_nodal_concentrated_forces(
              ffunction, nset_id, dir, force);
          if (!nodal_force)
            throw std::runtime_error(
                "Concentrated nodal forces are not properly assigned");
          set_node_concentrated_force_ = true;
        }
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
bool mpm::MPMBase<Tdim>::initialise_math_functions(const Json& math_functions) {
  bool status = true;
  try {
    // Get materials properties
    for (const auto& function_props : math_functions) {

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
            "Invalid properties for new math function, fn insertion failed");
      }
    }
  } catch (std::exception& exception) {
    console_->error("#{}: Reading math functions: {}", __LINE__,
                    exception.what());
    status = false;
  }
  return status;
}

//! Node entity sets
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::node_entity_sets(const Json& mesh_props,
                                          bool check_duplicates) {
  try {
    if (mesh_props.find("entity_sets") != mesh_props.end()) {
      std::string entity_sets =
          mesh_props["entity_sets"].template get<std::string>();
      if (!io_->file_name(entity_sets).empty()) {
        bool node_sets = mesh_->create_node_sets(
            (io_->entity_sets(io_->file_name(entity_sets), "node_sets")),
            check_duplicates);
        if (!node_sets)
          throw std::runtime_error("Node sets are not properly assigned");
      }
    } else
      throw std::runtime_error("Entity set JSON not found");
  } catch (std::exception& exception) {
    console_->warn("#{}: Entity sets are undefined {} ", __LINE__,
                   exception.what());
  }
}

//! Node Euler angles
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::node_euler_angles(
    const Json& mesh_props, const std::shared_ptr<mpm::IOMesh<Tdim>>& mesh_io) {
  try {
    if (mesh_props.find("boundary_conditions") != mesh_props.end() &&
        mesh_props["boundary_conditions"].find("nodal_euler_angles") !=
            mesh_props["boundary_conditions"].end()) {
      std::string euler_angles =
          mesh_props["boundary_conditions"]["nodal_euler_angles"]
              .template get<std::string>();
      if (!io_->file_name(euler_angles).empty()) {
        bool rotation_matrices = mesh_->compute_nodal_rotation_matrices(
            mesh_io->read_euler_angles(io_->file_name(euler_angles)));
        if (!rotation_matrices)
          throw std::runtime_error(
              "Euler angles are not properly assigned/computed");
      }
    } else
      throw std::runtime_error("Euler angles JSON not found");
  } catch (std::exception& exception) {
    console_->warn("#{}: Euler angles are undefined {} ", __LINE__,
                   exception.what());
  }
}

// Nodal velocity constraints
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::nodal_velocity_constraints(
    const Json& mesh_props, const std::shared_ptr<mpm::IOMesh<Tdim>>& mesh_io) {
  try {
    // Read and assign velocity constraints
    if (mesh_props.find("boundary_conditions") != mesh_props.end() &&
        mesh_props["boundary_conditions"].find("velocity_constraints") !=
            mesh_props["boundary_conditions"].end()) {
      // Iterate over velocity constraints
      for (const auto& constraints :
           mesh_props["boundary_conditions"]["velocity_constraints"]) {
        // Velocity constraints are specified in a file
        if (constraints.find("file") != constraints.end()) {
          std::string velocity_constraints_file =
              constraints.at("file").template get<std::string>();
          bool velocity_constraints =
              constraints_->assign_nodal_velocity_constraints(
                  mesh_io->read_velocity_constraints(
                      io_->file_name(velocity_constraints_file)));
          if (!velocity_constraints)
            throw std::runtime_error(
                "Velocity constraints are not properly assigned");

        } else {
          // Set id
          int nset_id = constraints.at("nset_id").template get<int>();
          // Direction
          unsigned dir = constraints.at("dir").template get<unsigned>();
          // Velocity
          double velocity = constraints.at("velocity").template get<double>();
          // Add velocity constraint to mesh
          auto velocity_constraint =
              std::make_shared<mpm::VelocityConstraint>(nset_id, dir, velocity);
          bool velocity_constraints =
              constraints_->assign_nodal_velocity_constraint(
                  nset_id, velocity_constraint);
          if (!velocity_constraints)
            throw std::runtime_error(
                "Nodal velocity constraint is not properly assigned");
        }
      }
    } else
      throw std::runtime_error("Velocity constraints JSON not found");
  } catch (std::exception& exception) {
    console_->warn("#{}: Velocity constraints are undefined {} ", __LINE__,
                   exception.what());
  }
}

// Nodal frictional constraints
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::nodal_frictional_constraints(
    const Json& mesh_props, const std::shared_ptr<mpm::IOMesh<Tdim>>& mesh_io) {
  try {
    // Read and assign friction constraints
    if (mesh_props.find("boundary_conditions") != mesh_props.end() &&
        mesh_props["boundary_conditions"].find("friction_constraints") !=
            mesh_props["boundary_conditions"].end()) {
      // Iterate over velocity constraints
      for (const auto& constraints :
           mesh_props["boundary_conditions"]["friction_constraints"]) {
        // Friction constraints are specified in a file
        if (constraints.find("file") != constraints.end()) {
          std::string friction_constraints_file =
              constraints.at("file").template get<std::string>();
          bool friction_constraints =
              constraints_->assign_nodal_friction_constraints(
                  mesh_io->read_friction_constraints(
                      io_->file_name(friction_constraints_file)));
          if (!friction_constraints)
            throw std::runtime_error(
                "Friction constraints are not properly assigned");

        } else {

          // Set id
          int nset_id = constraints.at("nset_id").template get<int>();
          // Direction
          unsigned dir = constraints.at("dir").template get<unsigned>();
          // Sign n
          int sign_n = constraints.at("sign_n").template get<int>();
          // Friction
          double friction = constraints.at("friction").template get<double>();
          // Add friction constraint to mesh
          auto friction_constraint = std::make_shared<mpm::FrictionConstraint>(
              nset_id, dir, sign_n, friction);
          bool friction_constraints =
              constraints_->assign_nodal_frictional_constraint(
                  nset_id, friction_constraint);
          if (!friction_constraints)
            throw std::runtime_error(
                "Nodal friction constraint is not properly assigned");
        }
      }
    } else
      throw std::runtime_error("Friction constraints JSON not found");

  } catch (std::exception& exception) {
    console_->warn("#{}: Friction conditions are undefined {} ", __LINE__,
                   exception.what());
  }
}

//! Cell entity sets
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::cell_entity_sets(const Json& mesh_props,
                                          bool check_duplicates) {
  try {
    if (mesh_props.find("entity_sets") != mesh_props.end()) {
      // Read and assign cell sets
      std::string entity_sets =
          mesh_props["entity_sets"].template get<std::string>();
      if (!io_->file_name(entity_sets).empty()) {
        bool cell_sets = mesh_->create_cell_sets(
            (io_->entity_sets(io_->file_name(entity_sets), "cell_sets")),
            check_duplicates);
        if (!cell_sets)
          throw std::runtime_error("Cell sets are not properly assigned");
      }
    } else
      throw std::runtime_error("Cell entity sets JSON not found");

  } catch (std::exception& exception) {
    console_->warn("#{}: Cell entity sets are undefined {} ", __LINE__,
                   exception.what());
  }
}

// Particles cells
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::particles_cells(
    const Json& mesh_props,
    const std::shared_ptr<mpm::IOMesh<Tdim>>& particle_io) {
  try {
    if (mesh_props.find("particle_cells") != mesh_props.end()) {
      std::string fparticles_cells =
          mesh_props["particle_cells"].template get<std::string>();

      if (!io_->file_name(fparticles_cells).empty()) {
        bool particles_cells =
            mesh_->assign_particles_cells(particle_io->read_particles_cells(
                io_->file_name(fparticles_cells)));
        if (!particles_cells)
          throw std::runtime_error(
              "Particle cells are not properly assigned to particles");
      }
    } else
      throw std::runtime_error("Particle cells JSON not found");

  } catch (std::exception& exception) {
    console_->warn("#{}: Particle cells are undefined {} ", __LINE__,
                   exception.what());
  }
}

// Particles volumes
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::particles_volumes(
    const Json& mesh_props,
    const std::shared_ptr<mpm::IOMesh<Tdim>>& particle_io) {
  try {
    if (mesh_props.find("particles_volumes") != mesh_props.end()) {
      std::string fparticles_volumes =
          mesh_props["particles_volumes"].template get<std::string>();
      if (!io_->file_name(fparticles_volumes).empty()) {
        bool particles_volumes =
            mesh_->assign_particles_volumes(particle_io->read_particles_volumes(
                io_->file_name(fparticles_volumes)));
        if (!particles_volumes)
          throw std::runtime_error(
              "Particles volumes are not properly assigned");
      }
    } else
      throw std::runtime_error("Particle volumes JSON not found");
  } catch (std::exception& exception) {
    console_->warn("#{}: Particle volumes are undefined {} ", __LINE__,
                   exception.what());
  }
}

// Particle velocity constraints
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::particle_velocity_constraints(
    const Json& mesh_props,
    const std::shared_ptr<mpm::IOMesh<Tdim>>& particle_io) {
  try {
    if (mesh_props.find("boundary_conditions") != mesh_props.end() &&
        mesh_props["boundary_conditions"].find(
            "particles_velocity_constraints") !=
            mesh_props["boundary_conditions"].end()) {

      // Iterate over velocity constraints
      for (const auto& constraints :
           mesh_props["boundary_conditions"]
                     ["particles_velocity_constraints"]) {

        // Set id
        int pset_id = constraints.at("pset_id").template get<int>();
        // Direction
        unsigned dir = constraints.at("dir").template get<unsigned>();
        // Velocity
        double velocity = constraints.at("velocity").template get<double>();
        // Add velocity constraint to mesh
        auto velocity_constraint =
            std::make_shared<mpm::VelocityConstraint>(pset_id, dir, velocity);
        mesh_->create_particle_velocity_constraint(pset_id,
                                                   velocity_constraint);
      }
    } else
      throw std::runtime_error("Particle velocity constraints JSON not found");
  } catch (std::exception& exception) {
    console_->warn("#{}: Particle velocity constraints are undefined {} ",
                   __LINE__, exception.what());
  }
}

// Particles stresses
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::particles_stresses(
    const Json& mesh_props,
    const std::shared_ptr<mpm::IOMesh<Tdim>>& particle_io) {
  try {
    if (mesh_props.find("particles_stresses") != mesh_props.end()) {
      std::string fparticles_stresses =
          mesh_props["particles_stresses"].template get<std::string>();
      if (!io_->file_name(fparticles_stresses).empty()) {

        // Get stresses of all particles
        const auto all_particles_stresses =
            particle_io->read_particles_stresses(
                io_->file_name(fparticles_stresses));

        // Read and assign particles stresses
        if (!mesh_->assign_particles_stresses(all_particles_stresses))
          throw std::runtime_error(
              "Particles stresses are not properly assigned");
      }
    } else
      throw std::runtime_error("Particle stresses JSON not found");

  } catch (std::exception& exception) {
    console_->warn("#{}: Particle stresses are undefined {} ", __LINE__,
                   exception.what());
  }
}

//! Particle entity sets
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::particle_entity_sets(const Json& mesh_props,
                                              bool check_duplicates) {
  // Read and assign particle sets
  try {
    if (mesh_props.find("entity_sets") != mesh_props.end()) {
      std::string entity_sets =
          mesh_props["entity_sets"].template get<std::string>();
      if (!io_->file_name(entity_sets).empty()) {
        bool particle_sets = mesh_->create_particle_sets(
            (io_->entity_sets(io_->file_name(entity_sets), "particle_sets")),
            check_duplicates);

        if (!particle_sets)
          throw std::runtime_error("Particle set creation failed");
      }
    } else
      throw std::runtime_error("Particle entity set JSON not found");

  } catch (std::exception& exception) {
    console_->warn("#{}: Particle sets are undefined {} ", __LINE__,
                   exception.what());
  }
}

// Initialise Damping
template <unsigned Tdim>
bool mpm::MPMBase<Tdim>::initialise_damping(const Json& damping_props) {

  // Read damping JSON object
  bool status = true;
  try {
    // Read damping type
    std::string type = damping_props.at("type").template get<std::string>();
    if (type == "Cundall") damping_type_ = mpm::Damping::Cundall;

    // Read damping factor
    damping_factor_ = damping_props.at("damping_factor").template get<double>();

  } catch (std::exception& exception) {
    console_->warn("#{}: Damping parameters are undefined {} ", __LINE__,
                   exception.what());
    status = false;
  }

  return status;
}

//! Domain decomposition
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::mpi_domain_decompose(bool initial_step) {
#ifdef USE_MPI
  // Initialise MPI rank and size
  int mpi_rank = 0;
  int mpi_size = 1;

  // Get MPI rank
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  // Get number of MPI ranks
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  if (mpi_size > 1 && mesh_->ncells() > 1) {

    // Initialize MPI
    MPI_Comm comm;
    MPI_Comm_dup(MPI_COMM_WORLD, &comm);

    auto mpi_domain_begin = std::chrono::steady_clock::now();
    console_->info("Rank {}, Domain decomposition started\n", mpi_rank);

    // Check if mesh has cells to partition
    if (mesh_->ncells() == 0)
      throw std::runtime_error("Container of cells is empty");

#ifdef USE_GRAPH_PARTITIONING
    // Create graph object if empty
    if (initial_step || graph_ == nullptr)
      graph_ = std::make_shared<Graph<Tdim>>(mesh_->cells());

    // Find number of particles in each cell across MPI ranks
    mesh_->find_nglobal_particles_cells();

    // Construct a weighted DAG
    graph_->construct_graph(mpi_size, mpi_rank);

    // Graph partitioning mode
    int mode = 4;  // FAST
    // Create graph partition
    graph_->create_partitions(&comm, mode);
    // Collect the partitions
    auto exchange_cells = graph_->collect_partitions(mpi_size, mpi_rank, &comm);

    // Identify shared nodes across MPI domains
    mesh_->find_domain_shared_nodes();
    // Identify ghost boundary cells
    mesh_->find_ghost_boundary_cells();

    // Delete all the particles which is not in local task parititon
    if (initial_step) mesh_->remove_all_nonrank_particles();
    // Transfer non-rank particles to appropriate cells
    else
      mesh_->transfer_nonrank_particles(exchange_cells);

#endif
    auto mpi_domain_end = std::chrono::steady_clock::now();
    console_->info("Rank {}, Domain decomposition: {} ms", mpi_rank,
                   std::chrono::duration_cast<std::chrono::milliseconds>(
                       mpi_domain_end - mpi_domain_begin)
                       .count());
  }
#endif  // MPI
}

//! MPM pressure smoothing
template <unsigned Tdim>
void mpm::MPMBase<Tdim>::pressure_smoothing(unsigned phase) {
  // Assign pressure to nodes
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::map_pressure_to_nodes,
                std::placeholders::_1, phase));

#ifdef USE_MPI
  int mpi_size = 1;

  // Get number of MPI ranks
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  // Run if there is more than a single MPI task
  if (mpi_size > 1) {
    // MPI all reduce nodal pressure
    mesh_->template nodal_halo_exchange<double, 1>(
        std::bind(&mpm::NodeBase<Tdim>::pressure, std::placeholders::_1, phase),
        std::bind(&mpm::NodeBase<Tdim>::assign_pressure, std::placeholders::_1,
                  phase, std::placeholders::_2));
  }
#endif

  // Smooth pressure over particles
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::compute_pressure_smoothing,
                std::placeholders::_1, phase));
}

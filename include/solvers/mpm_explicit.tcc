//! Constructor
template <unsigned Tdim>
mpm::MPMExplicit<Tdim>::MPMExplicit(const std::shared_ptr<IO>& io)
    : mpm::MPMBase<Tdim>(io) {
  //! Logger
  console_ = spdlog::get("MPMExplicit");
  //! Stress update
  if (this->stress_update_ == "usl")
    mpm_scheme_ = std::make_shared<mpm::MPMSchemeUSL<Tdim>>(mesh_, dt_);
  else
    mpm_scheme_ = std::make_shared<mpm::MPMSchemeUSF<Tdim>>(mesh_, dt_);

  //! Interface scheme
  if (this->interface_)
    contact_ = std::make_shared<mpm::ContactFriction<Tdim>>(mesh_);
  else
    contact_ = std::make_shared<mpm::Contact<Tdim>>(mesh_);

  //! Real-time monitoring
  real_time_monitoring_ = io_->analysis_bool("real_time_monitoring");
  if (io_->analysis_has_key("dashboard_url")) {
    dashboard_url_ = io_->analysis_string("dashboard_url");
  }
}

#include <curl/curl.h>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

//! MPM Explicit compute stress strain
template <unsigned Tdim>
void mpm::MPMExplicit<Tdim>::compute_stress_strain(unsigned phase) {
  // Iterate over each particle to calculate strain
  mesh_->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::compute_strain, std::placeholders::_1, dt_));

  // Iterate over each particle to update particle volume
  mesh_->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::update_volume, std::placeholders::_1));

  // Pressure smoothing
  if (pressure_smoothing_) this->pressure_smoothing(phase);

  // Iterate over each particle to compute stress
  mesh_->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::compute_stress, std::placeholders::_1));
}

//! Send real-time data to web dashboard
template <unsigned Tdim>
void mpm::MPMExplicit<Tdim>::send_real_time_data(unsigned step, double time) {
  if (!real_time_monitoring_) return;

  // Collect stress and strain data
  double avg_stress_xx = 0.0, avg_stress_yy = 0.0, avg_stress_zz = 0.0;
  double avg_strain_xx = 0.0, avg_strain_yy = 0.0, avg_strain_zz = 0.0;
  double max_stress_xx = -std::numeric_limits<double>::max();
  double min_stress_xx = std::numeric_limits<double>::max();
  size_t particle_count = 0;

  // Iterate over all particles to compute averages and extremes
  mesh_->iterate_over_particles([&](std::shared_ptr<mpm::ParticleBase<Tdim>> particle) {
    if (!particle->status()) return;

    auto stress = particle->stress();
    auto strain = particle->strain();

    avg_stress_xx += stress[0];
    avg_stress_yy += stress[1];
    avg_stress_zz += stress[2];
    avg_strain_xx += strain[0];
    avg_strain_yy += strain[1];
    avg_strain_zz += strain[2];

    max_stress_xx = std::max(max_stress_xx, stress[0]);
    min_stress_xx = std::min(min_stress_xx, stress[0]);

    particle_count++;
  });

  if (particle_count == 0) return;

  // Compute averages
  avg_stress_xx /= particle_count;
  avg_stress_yy /= particle_count;
  avg_stress_zz /= particle_count;
  avg_strain_xx /= particle_count;
  avg_strain_yy /= particle_count;
  avg_strain_zz /= particle_count;

  // Create JSON data
  json data = {
    {"step", step},
    {"total_steps", nsteps_},
    {"time", time},
    {"average_stress_xx", avg_stress_xx},
    {"average_stress_yy", avg_stress_yy},
    {"average_stress_zz", avg_stress_zz},
    {"average_strain_xx", avg_strain_xx},
    {"average_strain_yy", avg_strain_yy},
    {"average_strain_zz", avg_strain_zz},
    {"max_stress_xx", max_stress_xx},
    {"min_stress_xx", min_stress_xx}
  };

  // Convert JSON to string
  std::string json_str = data.dump();

  // Initialize CURL
  CURL *curl = curl_easy_init();
  if (curl) {
    // Set URL
    curl_easy_setopt(curl, CURLOPT_URL, dashboard_url_.c_str());

    // Set POST request
    curl_easy_setopt(curl, CURLOPT_POST, 1L);

    // Set POST data
    curl_easy_setopt(curl, CURLOPT_POSTFIELDS, json_str.c_str());

    // Set headers
    struct curl_slist *headers = NULL;
    headers = curl_slist_append(headers, "Content-Type: application/json");
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);

    // Perform request
    CURLcode res = curl_easy_perform(curl);
    if (res != CURLE_OK) {
      console_->warn("Failed to send real-time data: {}", curl_easy_strerror(res));
    }

    // Cleanup
    curl_slist_free_all(headers);
    curl_easy_cleanup(curl);
  }
}

//! MPM Explicit solver
template <unsigned Tdim>
bool mpm::MPMExplicit<Tdim>::solve() {
  bool status = true;

  console_->info("MPM analysis type {}", io_->analysis_type());

  // Initialise MPI rank and size
  int mpi_rank = 0;
  int mpi_size = 1;

#ifdef USE_MPI
  // Get MPI rank
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  // Get number of MPI ranks
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif

  // Phase
  const unsigned phase = 0;

  // Test if checkpoint resume is needed
  bool resume = false;
  if (analysis_.find("resume") != analysis_.end())
    resume = analysis_["resume"]["resume"].template get<bool>();

  // Enable repartitioning if resume is done with particles generated outside
  // the MPM code.
  bool repartition = false;
  if (analysis_.find("resume") != analysis_.end() &&
      analysis_["resume"].find("repartition") != analysis_["resume"].end())
    repartition = analysis_["resume"]["repartition"].template get<bool>();

  // Pressure smoothing
  pressure_smoothing_ = io_->analysis_bool("pressure_smoothing");

  // Interface
  interface_ = io_->analysis_bool("interface");

  // Initialise material
  this->initialise_materials();

  // Initialise mesh
  this->initialise_mesh();

  // Initialise particles
  if (!resume) this->initialise_particles();

  // Create nodal properties
  if (interface_) mesh_->create_nodal_properties();

  // Compute mass
  if (!resume)
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::compute_mass, std::placeholders::_1));

  bool initial_step = (resume == true) ? false : true;
  // Check point resume
  if (resume) {
    this->checkpoint_resume();
    if (repartition) {
      this->mpi_domain_decompose(initial_step);
    } else {
      mesh_->resume_domain_cell_ranks();
#ifdef USE_MPI
#ifdef USE_GRAPH_PARTITIONING
      MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif
    }
  } else {
    // Domain decompose
    this->mpi_domain_decompose(initial_step);
  }

  //! Particle entity sets and velocity constraints
  if (resume) {
    this->particle_entity_sets(false);
    this->particle_velocity_constraints();
  }

  // Initialise loading conditions
  this->initialise_loads();

  auto solver_begin = std::chrono::steady_clock::now();
  // Main loop
  for (; step_ < nsteps_; ++step_) {

    if (mpi_rank == 0) console_->info("Step: {} of {}.\n", step_, nsteps_);

#ifdef USE_MPI
#ifdef USE_GRAPH_PARTITIONING
    // Run load balancer at a specified frequency
    if (step_ % nload_balance_steps_ == 0 && step_ != 0)
      this->mpi_domain_decompose(false);
#endif
#endif

    // Inject particles
    mesh_->inject_particles(step_ * dt_);

    // Initialise nodes, cells and shape functions
    mpm_scheme_->initialise();

    // Initialise nodal properties and append material ids to node
    contact_->initialise();

    // Mass momentum and compute velocity at nodes
    mpm_scheme_->compute_nodal_kinematics(phase);

    // Map material properties to nodes
    contact_->compute_contact_forces();

    // Update stress first
    mpm_scheme_->precompute_stress_strain(phase, pressure_smoothing_);

    // Compute forces
    mpm_scheme_->compute_forces(gravity_, phase, step_,
                                set_node_concentrated_force_);

    // Particle kinematics
    mpm_scheme_->compute_particle_kinematics(velocity_update_, phase, "Cundall",
                                             damping_factor_);

    // Update Stress Last
    mpm_scheme_->postcompute_stress_strain(phase, pressure_smoothing_);

    // Send real-time data to dashboard
    send_real_time_data(step_, step_ * dt_);

    // Locate particles
    mpm_scheme_->locate_particles(this->locate_particles_);

#ifdef USE_MPI
#ifdef USE_GRAPH_PARTITIONING
    mesh_->transfer_halo_particles();
    MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif

    if (step_ % output_steps_ == 0) {
      // HDF5 outputs
      this->write_hdf5(this->step_, this->nsteps_);
#ifdef USE_VTK
      // VTK outputs
      this->write_vtk(this->step_, this->nsteps_);
#endif
#ifdef USE_PARTIO
      // Partio outputs
      this->write_partio(this->step_, this->nsteps_);
#endif
    }
  }
  auto solver_end = std::chrono::steady_clock::now();
  console_->info("Rank {}, Explicit {} solver duration: {} ms", mpi_rank,
                 mpm_scheme_->scheme(),
                 std::chrono::duration_cast<std::chrono::milliseconds>(
                     solver_end - solver_begin)
                     .count());

  return status;
}

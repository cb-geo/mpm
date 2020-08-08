#ifndef MPM_MPM_BASE_H_
#define MPM_MPM_BASE_H_

#include <numeric>

#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// MPI
#ifdef USE_MPI
#include "mpi.h"
#endif

#ifdef USE_GRAPH_PARTITIONING
#include "graph.h"
#endif

#include "constraints.h"
#include "mpm.h"
#include "particle.h"
#include "vector.h"

namespace mpm {

//! Stress update method
//! USF: Update Stress First
//! USL: Update Stress Last
//! MUSL: Modified Stress Last
enum class StressUpdate { USF, USL, MUSL };
extern std::map<std::string, StressUpdate> stress_update;

//! Damping type
//! None: No damping is specified
//! Cundall: Cundall damping
enum class Damping { None, Cundall };

//! MPMBase class
//! \brief A class that implements the fully base one phase mpm
//! \details A Base MPM class
//! \tparam Tdim Dimension
template <unsigned Tdim>
class MPMBase : public MPM {
 public:
  //! Default constructor
  MPMBase(const std::shared_ptr<IO>& io);

  //! Initialise mesh
  bool initialise_mesh() override;

  //! Initialise particles
  bool initialise_particles() override;

  //! Initialise materials
  bool initialise_materials() override;

  //! Initialise loading
  bool initialise_loads() override;

  //! Initialise math functions
  bool initialise_math_functions(const Json&) override;

  //! Solve
  bool solve() override { return true; }

  //! Checkpoint resume
  bool checkpoint_resume() override;

#ifdef USE_VTK
  //! Write VTK files
  void write_vtk(mpm::Index step, mpm::Index max_steps) override;
#endif

#ifdef USE_PARTIO
  //! Write PARTIO files
  void write_partio(mpm::Index step, mpm::Index max_steps) override;
#endif

  //! Write HDF5 files
  void write_hdf5(mpm::Index step, mpm::Index max_steps) override;

  //! Domain decomposition
  //! \param[in] initial_step Start of simulation or later steps
  void mpi_domain_decompose(bool initial_step = false) override;

  //! Pressure smoothing
  //! \param[in] phase Phase to smooth pressure
  void pressure_smoothing(unsigned phase);

 private:
  //! Return if a mesh will be isoparametric or not
  //! \retval isoparametric Status of mesh type
  bool is_isoparametric();

  //! Node entity sets
  //! \param[in] mesh_prop Mesh properties
  //! \param[in] check Check duplicates
  void node_entity_sets(const Json& mesh_prop, bool check);

  //! Node Euler angles
  //! \param[in] mesh_prop Mesh properties
  //! \param[in] mesh_io Mesh IO handle
  void node_euler_angles(const Json& mesh_prop,
                         const std::shared_ptr<mpm::IOMesh<Tdim>>& mesh_io);

  //! Nodal velocity constraints
  //! \param[in] mesh_prop Mesh properties
  //! \param[in] mesh_io Mesh IO handle
  void nodal_velocity_constraints(
      const Json& mesh_prop, const std::shared_ptr<mpm::IOMesh<Tdim>>& mesh_io);

  //! Nodal frictional constraints
  //! \param[in] mesh_prop Mesh properties
  //! \param[in] mesh_io Mesh IO handle
  void nodal_frictional_constraints(
      const Json& mesh_prop, const std::shared_ptr<mpm::IOMesh<Tdim>>& mesh_io);

  //! Cell entity sets
  //! \param[in] mesh_prop Mesh properties
  //! \param[in] check Check duplicates
  void cell_entity_sets(const Json& mesh_prop, bool check);

  //! Particles cells
  //! \param[in] mesh_prop Mesh properties
  //! \param[in] particle_io Particle IO handle
  void particles_cells(const Json& mesh_prop,
                       const std::shared_ptr<mpm::IOMesh<Tdim>>& particle_io);

  //! Particles volumes
  //! \param[in] mesh_prop Mesh properties
  //! \param[in] particle_io Particle IO handle
  void particles_volumes(const Json& mesh_prop,
                         const std::shared_ptr<mpm::IOMesh<Tdim>>& particle_io);

  //! Particle velocity constraints
  //! \param[in] mesh_prop Mesh properties
  //! \param[in] particle_io Particle IO handle
  void particle_velocity_constraints(
      const Json& mesh_prop,
      const std::shared_ptr<mpm::IOMesh<Tdim>>& particle_io);

  //! Particles stresses
  //! \param[in] mesh_prop Mesh properties
  //! \param[in] particle_io Particle IO handle
  void particles_stresses(
      const Json& mesh_prop,
      const std::shared_ptr<mpm::IOMesh<Tdim>>& particle_io);

  //! Particle entity sets
  //! \param[in] mesh_prop Mesh properties
  //! \param[in] check Check duplicates
  void particle_entity_sets(const Json& mesh_prop, bool check);

  //! Initialise damping
  //! \param[in] damping_props Damping properties
  bool initialise_damping(const Json& damping_props);

 protected:
  // Generate a unique id for the analysis
  using mpm::MPM::uuid_;
  //! Time step size
  using mpm::MPM::dt_;
  //! Current step
  using mpm::MPM::step_;
  //! Number of steps
  using mpm::MPM::nsteps_;
  //! Output steps
  using mpm::MPM::output_steps_;
  //! Load balancing steps
  using mpm::MPM::nload_balance_steps_;
  //! A unique ptr to IO object
  using mpm::MPM::io_;
  //! JSON analysis object
  using mpm::MPM::analysis_;
  //! JSON post-process object
  using mpm::MPM::post_process_;
  //! Logger
  using mpm::MPM::console_;

  //! Stress update method (default USF = 0, USL = 1, MUSL = 2)
  mpm::StressUpdate stress_update_{mpm::StressUpdate::USF};
  //! velocity update
  bool velocity_update_{false};
  //! Gravity
  Eigen::Matrix<double, Tdim, 1> gravity_;
  //! Mesh object
  std::shared_ptr<mpm::Mesh<Tdim>> mesh_;
  //! Constraints object
  std::shared_ptr<mpm::Constraints<Tdim>> constraints_;
  //! Materials
  std::map<unsigned, std::shared_ptr<mpm::Material<Tdim>>> materials_;
  //! Mathematical functions
  std::map<unsigned, std::shared_ptr<mpm::FunctionBase>> math_functions_;
  //! VTK state variables
  tsl::robin_map<unsigned, std::vector<std::string>> vtk_statevars_;
  //! Set node concentrated force
  bool set_node_concentrated_force_{false};
  //! Damping type
  mpm::Damping damping_type_{mpm::Damping::None};
  //! Damping factor
  double damping_factor_{0.};
  //! Locate particles
  bool locate_particles_{true};

#ifdef USE_GRAPH_PARTITIONING
  // graph pass the address of the container of cell
  std::shared_ptr<Graph<Tdim>> graph_{nullptr};
#endif
};  // MPMBase class
}  // namespace mpm

#include "mpm_base.tcc"

#endif  // MPM_MPM_BASE_H_

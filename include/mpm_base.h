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
#include "tbb/task_group.h"

#ifdef USE_GRAPH_PARTITIONING
#include "graph.h"
#endif

#include "container.h"
#include "friction_constraint.h"
#include "mpm.h"
#include "particle.h"
#include "velocity_constraint.h"

namespace mpm {

//! Qudrature rule
//! (1) StandardMPM considers material points with mass and volume as
//! integration points and quadrature weight is assumes as particle volume.
//! (2) Gauss quadrature computes the polynomial integration at fixed Gauss
//! points. This use is limited in MPM, however, it can be used when
//! the code is used for simple FEM analysis.
//! (3) Moving Gauss quadrature considers material points as moving
//! integration points with no mass/volume. The quadrature weight
//! is computed from Gauss qudrature rule.
//! (4) GaussMPM integration takes two integration methods for interior and
//! boundary elements. For interior elements, the integration is performaed at
//! fixed Gauss locations whereas for boundary elements, the integration is
//! performed at particle level.
enum class QuadratureRule { StandardMPM, Gauss, MovingGauss, GaussMPM };

//! Stress update method
//! USF: Update Stress First
//! USL: Update Stress Last
//! MUSL: Modified Stress Last
enum class StressUpdate { USF, USL, MUSL };
extern std::map<std::string, StressUpdate> stress_update;

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

  //! Initialise damping
  bool initialise_damping(const Json&) override;

  //! Solve
  bool solve() override { return true; }

  //! Checkpoint resume
  bool checkpoint_resume() override;

#ifdef USE_VTK
  //! Write VTK files
  void write_vtk(mpm::Index step, mpm::Index max_steps) override;
#endif

  //! Write HDF5 files
  void write_hdf5(mpm::Index step, mpm::Index max_steps) override;

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
  std::unique_ptr<mpm::Mesh<Tdim>> mesh_;
  //! Materials
  std::map<unsigned, std::shared_ptr<mpm::Material<Tdim>>> materials_;
  //! Mathematical functions
  std::map<unsigned, std::shared_ptr<mpm::FunctionBase>> math_functions_;
  //! VTK attributes
  std::vector<std::string> vtk_attributes_;
  //! Set node concentrated force
  bool set_node_concentrated_force_{false};
  //! Bool nodal tractions
  bool nodal_tractions_{true};
  //! Qudrature rule (default is standard mpm integration)
  mpm::QuadratureRule quadrature_rule_{QuadratureRule::StandardMPM};
  //! Polynomial order for Approximate Gauss quadrature rule
  unsigned quadrature_order_{2};
  // Level set methods
  bool ls_methods_{false};
  //! Damping type
  enum class Damping{ None, Cundall };
  Damping damping_type_{mpm::MPMBase<Tdim>::Damping::None};
  //! Damping factor
  double damping_factor_{0.};

#ifdef USE_GRAPH_PARTITIONING
  // graph pass the address of the container of cell
  std::shared_ptr<Graph<Tdim>> graph_{nullptr};
#endif
};  // MPMBase class
}  // namespace mpm

#include "mpm_base.tcc"

#endif  // MPM_MPM_BASE_H_

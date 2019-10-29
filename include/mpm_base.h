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

#include "container.h"
#include "generators/generator_factory.h"
#include "mpi_wrapper.h"
#include "mpm.h"
#include "particle.h"

namespace mpm {

//! Qudrature rule
//! (1) Standard MPM considers material points with mass and volume as
//! integration points and quadrature weight is assumes as particle volume.
//! (2) Gauss quadrature computes the polynomial integration at fixed Gauss
//! points. This use is limited in MPM, however, it can be used when
//! the code is used for simple FEM analysis.
//! (3) Moving Gauss quadrature considers material points as moving
//! integration points with no mass/volume. The quadrature weight
//! is computed from Gauss qudrature rule.
//! (4) Gauss_mpm integration takes two integration methods for interior and
//! boundary elements. For interior elements, the integration is performaed at
//! fixed Gauss locations whereas for boundary elements, the integration is
//! performed at particle level.
enum class QuadratureRule { Standard_mpm, Gauss, Moving_Gauss, Gauss_mpm };

//! Stress update method
//! usl: Update Stress Last
//! usf: Update Stress First
//! musl: Modified Stress Last
enum class StressUpdate { usf, usl, musl };

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

  //! Apply nodal tractions
  bool apply_nodal_tractions() override;

  //! Apply properties to particles sets (e.g: material)
  bool apply_properties_to_particles_sets() override;

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
  //! Stress update method (default usf = 0, usl = 1, musl = 2)
  mpm::StressUpdate stress_update_{mpm::StressUpdate::usf};
  //! velocity update
  bool velocity_update_{false};
  //! Gravity
  Eigen::Matrix<double, Tdim, 1> gravity_;
  //! Mesh object
  std::shared_ptr<mpm::Mesh<Tdim>> mesh_;
  //! Materials
  std::map<unsigned, std::shared_ptr<mpm::Material<Tdim>>> materials_;
  //! VTK attributes
  std::vector<std::string> vtk_attributes_;
  //! Bool nodal tractions
  bool nodal_tractions_{true};
  //! Qudrature rule (default is standard mpm integration)
  mpm::QuadratureRule quadrature_rule_{QuadratureRule::Standard_mpm};
  //! Polynomial order for Approximate Gauss quadrature rule
  unsigned quadrature_order_{2};
  // Level set methods
  bool ls_methods_{false};
};  // MPMBase class
}  // namespace mpm

#include "mpm_base.tcc"

#endif  // MPM_MPM_BASE_H_

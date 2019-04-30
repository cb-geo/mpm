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
#include "mpi_wrapper.h"
#include "mpm.h"
#include "particle.h"

namespace mpm {

//! MPMBase class
//! \brief A class that implements the fully base one phase mpm
//! \details A Base MPM class
//! \tparam Tdim Dimension
template <unsigned Tdim>
class MPMBase : public MPM {
 public:
  //! Default constructor
  MPMBase(std::unique_ptr<IO>&& io);

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

  //! velocity update
  bool velocity_update_{false};
  //! Gravity
  Eigen::Matrix<double, Tdim, 1> gravity_;
  //! Mesh object
  std::unique_ptr<mpm::Mesh<Tdim>> mesh_;
  //! Materials
  std::map<unsigned, std::shared_ptr<mpm::Material<Tdim>>> materials_;
  //! VTK attributes
  std::vector<std::string> vtk_attributes_;
  //! Bool nodal tractions
  bool nodal_tractions_{true};
};  // MPMBase class
}  // namespace mpm

#include "mpm_base.tcc"

#endif  // MPM_MPM_BASE_H_

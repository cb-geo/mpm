#ifndef MPM_MPM_EXPLICIT_H_
#define MPM_MPM_EXPLICIT_H_

#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

#include "container.h"
#include "mpm.h"
#include "particle.h"

namespace mpm {

//! MPMExplicit class
//! \brief A class that implements the fully explicit one phase mpm
//! \details A single-phase explicit MPM
//! \tparam Tdim Dimension
template <unsigned Tdim>
class MPMExplicit : public MPM {
 public:
  //! Default constructor
  MPMExplicit(std::unique_ptr<IO>&& io);

  //! Initialise mesh and particles
  bool initialise_mesh_particles() override;

  //! Initialise materials
  bool initialise_materials() override;

  //! Solve
  bool solve() override = 0;

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

  //! Gravity
  Eigen::Matrix<double, Tdim, 1> gravity_;
  //! Mesh object
  std::vector<std::unique_ptr<mpm::Mesh<Tdim>>> meshes_;
  //! Materials
  std::map<unsigned, std::shared_ptr<mpm::Material<Tdim>>> materials_;
  //! VTK attributes
  std::vector<std::string> vtk_attributes_;

};  // MPMExplicit class
}  // namespace mpm

#include "mpm_explicit.tcc"

#endif  // MPM_MPM_EXPLICIT_H_

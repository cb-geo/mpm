#ifndef MPM_MPM_SEMI_IMPLICIT_NAVIER_STOKES_H_
#define MPM_MPM_SEMI_IMPLICIT_NAVIER_STOKES_H_

#ifdef USE_GRAPH_PARTITIONING
#include "graph.h"
#endif

#include "mpm_base.h"

#include "assembler_base.h"
#include "cg_eigen.h"
#include "solver_base.h"

namespace mpm {

//! MPMSemiImplicit Navier Stokes class
//! \brief A class that implements the fractional step navier-stokes mpm
//! \tparam Tdim Dimension
template <unsigned Tdim>
class MPMSemiImplicitNavierStokes : public MPMBase<Tdim> {
 public:
  //! Default constructor
  MPMSemiImplicitNavierStokes(const std::shared_ptr<IO>& io);

  //! Return matrix assembler pointer
  std::shared_ptr<mpm::AssemblerBase<Tdim>> matrix_assembler() {
    return matrix_assembler_;
  }

  //! Solve
  bool solve() override;

  //! Class private functions
 private:
  //! Initialise matrix
  virtual bool initialise_matrix();

  //! Initialise matrix
  virtual bool reinitialise_matrix();

  //! Compute poisson equation
  bool compute_poisson_equation(std::string solver_type = "cg");

  //! Compute corrected velocity
  bool compute_correction_force();

  //! Class private variables
 private:
  // Generate a unique id for the analysis
  using mpm::MPMBase<Tdim>::uuid_;
  //! Time step size
  using mpm::MPMBase<Tdim>::dt_;
  //! Current step
  using mpm::MPMBase<Tdim>::step_;
  //! Number of steps
  using mpm::MPMBase<Tdim>::nsteps_;
  //! Output steps
  using mpm::MPMBase<Tdim>::output_steps_;
  //! A unique ptr to IO object
  using mpm::MPMBase<Tdim>::io_;
  //! JSON analysis object
  using mpm::MPMBase<Tdim>::analysis_;
  //! JSON post-process object
  using mpm::MPMBase<Tdim>::post_process_;
  //! Logger
  using mpm::MPMBase<Tdim>::console_;
  //! Stress update
  using mpm::MPMBase<Tdim>::stress_update_;
  //! velocity update
  using mpm::MPMBase<Tdim>::velocity_update_;
  //! Gravity
  using mpm::MPMBase<Tdim>::gravity_;
  //! Mesh object
  using mpm::MPMBase<Tdim>::mesh_;
  //! Materials
  using mpm::MPMBase<Tdim>::materials_;
  //! Pressure smoothing
  bool pressure_smoothing_{false};
  // Projection method paramter (beta)
  double beta_{1};
  //! Matrix assembler object
  std::shared_ptr<mpm::AssemblerBase<Tdim>> matrix_assembler_;
  //! Matrix solver object
  std::shared_ptr<mpm::SolverBase<Eigen::SparseMatrix<double>>> matrix_solver_;
  //! Volume tolerance for free surface
  double volume_tolerance_{0};

};  // MPMSemiImplicit class
}  // namespace mpm

#include "mpm_semi_implicit_navierstokes.tcc"

#endif  // MPM_MPM_SEMI_IMPLICIT_NAVIER_STOKES_H_

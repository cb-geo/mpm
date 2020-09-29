#ifndef MPM_MPM_EXPLICIT_TWOPHASE_H_
#define MPM_MPM_EXPLICIT_TWOPHASE_H_

#ifdef USE_GRAPH_PARTITIONING
#include "graph.h"
#endif

#include "solvers/mpm_base.h"

namespace mpm {

//! MPMExplicit class
//! \brief A class that implements the fully explicit one phase mpm
//! \details A single-phase explicit MPM
//! \tparam Tdim Dimension
template <unsigned Tdim>
class MPMExplicitTwoPhase : public MPMBase<Tdim> {
 public:
  //! Default constructor
  MPMExplicitTwoPhase(const std::shared_ptr<IO>& io);

  //! Solve
  bool solve() override;

  //! Compute stress strain
  void compute_stress_strain();

 protected:
  // Generate a unique id for the analysis
  using mpm::MPMBase<Tdim>::uuid_;
  //! Time step size
  using mpm::MPMBase<Tdim>::dt_;
  //! Current step
  using mpm::MPMBase<Tdim>::step_;
  //! Number of steps
  using mpm::MPMBase<Tdim>::nsteps_;
  //! Number of steps
  using mpm::MPMBase<Tdim>::nload_balance_steps_;
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

#ifdef USE_GRAPH_PARTITIONING
  //! Graph
  using mpm::MPMBase<Tdim>::graph_;
#endif

  //! velocity update
  using mpm::MPMBase<Tdim>::velocity_update_;
  //! Gravity
  using mpm::MPMBase<Tdim>::gravity_;
  //! Mesh object
  using mpm::MPMBase<Tdim>::mesh_;
  //! Materials
  using mpm::MPMBase<Tdim>::materials_;
  //! Node concentrated force
  using mpm::MPMBase<Tdim>::set_node_concentrated_force_;
  //! Damping type
  using mpm::MPMBase<Tdim>::damping_type_;
  //! Damping factor
  using mpm::MPMBase<Tdim>::damping_factor_;
  //! Locate particles
  using mpm::MPMBase<Tdim>::locate_particles_;

 private:
  //! Pressure smoothing
  bool pressure_smoothing_{false};
  //! Pore pressure smoothing
  bool pore_pressure_smoothing_{false};
  //! Compute free surface
  std::string free_surface_detection_;
  //! Volume tolerance for free surface
  double volume_tolerance_{0.};

};  // MPMExplicit class
}  // namespace mpm

#include "mpm_explicit_twophase.tcc"

#endif  // MPM_MPM_EXPLICIT_H_

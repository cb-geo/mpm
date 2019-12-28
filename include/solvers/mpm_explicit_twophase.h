#ifndef MPM_MPM_EXPLICIT_TWOPHASE_H_
#define MPM_MPM_EXPLICIT_TWOPHASE_H_

#include "mpm_base.h"

namespace mpm {

//! MPMExplicit class
//! \brief A class that implements the fully explicit one phase mpm
//! \details A single-phase explicit MPM
//! \tparam Tdim Dimension
template <unsigned Tdim>
class MPMExplicitTwoPhase : public MPMBase<Tdim> {
 public:
  //! Default constructor
  MPMExplicitTwoPhase(std::unique_ptr<IO>&& io);

  //! Solve
  bool solve() override;

 protected:
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
  //! VTK attributes
  using mpm::MPMBase<Tdim>::vtk_attributes_;

 private:
  //! Pressure smoothing
  bool pressure_smoothing_{false};

};  // MPMExplicit class
}  // namespace mpm

#include "mpm_explicit_twophase.tcc"

#endif  // MPM_MPM_EXPLICIT_H_

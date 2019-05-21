#ifndef MPM_MPM_EXPLICIT_2P_H_
#define MPM_MPM_EXPLICIT_2P_H_

#include "mpm_base.h"

namespace mpm {

//! MPMExplicit Two phases class
//! \brief A class that implements the fully explicit two phase mpm
//! \details A two-phase explicit MPM
//! \tparam Tdim Dimension
template <unsigned Tdim>
class MPMExplicit2P : public MPMBase<Tdim> {
 public:
  //! Default constructor
  MPMExplicit2P(std::unique_ptr<IO>&& io);

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
  //! Bool nodal tractions
  using mpm::MPMBase<Tdim>::nodal_tractions_;
  //! Number of phases
  using mpm::MPMBase<Tdim>::tnphases_;

 private:
  //! Boolean to switch between USL and USF
  bool usl_{false};
  //! Pressure smoothing
  bool pressure_smoothing_{false};

};  // MPMExplicit Two phases class
}  // namespace mpm

#include "mpm_explicit_2p.tcc"

#endif  // MPM_MPM_EXPLICIT_2P_H_

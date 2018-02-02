#ifndef MPM_PARTICLE_H_
#define MPM_PARTICLE_H_

#include "mpm.h"
#include "container.h"
#include "particle.h"


namespace mpm {

//! ExplicitOnePhaseMpm class
//! \Class that implements the fully explicit one phase mpm
//! \details ExplicitOnePhaseMpm class: dt
//! \tparam Tdim Dimension
template <unsigned Tdim>
class ExplicitOnePhaseMpm : public MPM<Tdim> {
 public:
  //! Default constructor
  ExplicitOnePhaseMpm() {};
  

  // container 1: container of particles
  // container 2: container of cells or nodes?

  // map particle mass to nodes
  // map particle momentum to nodes
  // compute nodal velocity

  // body force
  // traction force
  // internal force

  // solve equation of motion

  // update particles (stress , velocity, position)

 private:
  mpm::Container<mpm::ParticleBase<Tdim>> particles_;
  
}; // ExplicitOnePhaseMpm class
} // namespace mpm

#endif  // MPM_EXPLICITONEPHASE_MPM_H_

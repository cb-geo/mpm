#ifndef MPM_MPM_EXPLICIT_SINGLE_PHASE_H_
#define MPM_MPM_EXPLICIT_SINGLE_PHASE_H_

#include "container.h"
#include "mpm.h"
#include "particle.h"

namespace mpm {

//! MPMExplicitSinglePhase class
//! \brief A class that implements the fully explicit one phase mpm
//! \details ExplicitOnePhaseMpm class: dt
//! \tparam Tdim Dimension
template <unsigned Tdim>
class MPMExplicitSinglePhase : public MPM<Tdim> {
 public:
  //! Default constructor
  MPMExplicitSinglePhase(){};

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

};  // MPMExplicitSinglePhase class
}  // namespace mpm

#endif  // MPM_EXPLICITONEPHASE_MPM_H_

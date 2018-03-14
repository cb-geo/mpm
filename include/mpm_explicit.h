#ifndef MPM_MPM_EXPLICIT_H_
#define MPM_MPM_EXPLICIT_H_

#include "container.h"
#include "mpm.h"
#include "particle.h"

namespace mpm {

//! MPMExplicit class
//! \brief A class that implements the fully explicit one phase mpm
//! \details A single-phase explicit MPM
//! \tparam Tdim Dimension
template <unsigned Tdim>
class MPMExplicit : public MPM<Tdim> {
 public:
  //! Default constructor
  MPMExplicit() = default;

  // map particle mass to nodes
  // map particle momentum to nodes
  // compute nodal velocity

  // body force
  // traction force
  // internal force

  // solve equation of motion

  // update particles (stress , velocity, position)

};  // MPMExplicit class
}  // namespace mpm

#endif  // MPM_MPM_EXPLICIT_H_

#ifndef MPM_PARTICLE_H_
#define MPM_PARTICLE_H_

#include "mpm.h"
#include "cell.h"
#include "particle.h"

namespace mpm {

//! ExplicitOnePhaseMpm class
//! \Class that implements the fully explicit one phase mpm
//! \details ExplicitOnePhaseMpm class: dt
//! \tparam Tdim Dimension
template <unsigned Tdim>
class ExplicitOnePhaseMpm : public Mpm<Tdim> {

}; // ExplicitOnePhaseMpm class
} // namespace mpm

#endif  // MPM_EXPLICITONEPHASE_MPM_H_

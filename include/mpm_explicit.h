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

};  // MPMExplicit class
}  // namespace mpm

#endif  // MPM_MPM_EXPLICIT_H_

#ifndef MPM_MPM_EXPLICIT_USF_H_
#define MPM_MPM_EXPLICIT_USF_H_

#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

#include "container.h"
#include "mpm.h"
#include "mpm_explicit.h"
#include "particle.h"

namespace mpm {

//! MPMExplicitUsf class
//! \brief Explicit one phase mpm with USF
//! \details A single-phase explicit MPM with Update Stress Last
//! \tparam Tdim Dimension
template <unsigned Tdim>
class MPMExplicitUSF : public MPMExplicit<Tdim> {
 public:
  //! Constructor
  MPMExplicitUSF(std::unique_ptr<IO>&& io);

  //! Solve
  bool solve() override;

 protected:
  // Generate a unique id for the analysis
  using mpm::MPMExplicit<Tdim>::uuid_;
  //! Time step size
  using mpm::MPMExplicit<Tdim>::dt_;
  //! Current step
  using mpm::MPMExplicit<Tdim>::step_;
  //! Number of steps
  using mpm::MPMExplicit<Tdim>::nsteps_;
  //! Output steps
  using mpm::MPMExplicit<Tdim>::output_steps_;
  //! A unique ptr to IO object
  using mpm::MPMExplicit<Tdim>::io_;
  //! JSON analysis object
  using mpm::MPMExplicit<Tdim>::analysis_;
  //! JSON post-process object
  using mpm::MPMExplicit<Tdim>::post_process_;
  //! Logger
  using mpm::MPMExplicit<Tdim>::console_;

  //! Gravity
  using mpm::MPMExplicit<Tdim>::gravity_;
  //! Mesh object
  using mpm::MPMExplicit<Tdim>::mesh_;
  //! Materials
  using mpm::MPMExplicit<Tdim>::materials_;

};  // MPMExplicitUSF class
}  // namespace mpm

#include "mpm_explicit_usf.tcc"

#endif  // MPM_MPM_EXPLICIT_USF_H_

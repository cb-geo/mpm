#ifndef MPM_MPM_EXPLICIT_H_
#define MPM_MPM_EXPLICIT_H_

#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

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
  MPMExplicit(std::unique_ptr<IO>&& io);

  //! Initialise
  bool initialise();

 protected:
  // Generate a unique id for the analysis
  using mpm::MPM<Tdim>::uuid_;
  //! Time step size
  using mpm::MPM<Tdim>::dt_;
  //! Number of steps
  using mpm::MPM<Tdim>::nsteps_;
  //! Gravity
  using mpm::MPM<Tdim>::gravity_;
  //! Mesh object
  using mpm::MPM<Tdim>::meshes_;
  //! A unique ptr to IO object
  using mpm::MPM<Tdim>::io_;
  //! JSON analysis object
  using mpm::MPM<Tdim>::analysis_;
  //! Logger
  using mpm::MPM<Tdim>::console_;

};  // MPMExplicit class
}  // namespace mpm

#include "mpm_explicit.tcc"

#endif  // MPM_MPM_EXPLICIT_H_

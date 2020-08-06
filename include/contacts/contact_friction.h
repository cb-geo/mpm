#ifndef _MPM_CONTACT_FRICTION_
#define _MPM_CONTACT_FRICTION_

#include "contact.h"

namespace mpm {

//! ContactFriction class
//! \brief ContactFriction base class
//! \tparam Tdim Dimension
template <unsigned Tdim>
class ContactFriction : public Contact<Tdim> {
 public:
  //! Default constructor with mesh class
  ContactFriction(const std::shared_ptr<mpm::Mesh<Tdim>>& mesh);

  //! Intialize
  virtual inline void initialise() override;

  //! Compute contact forces
  virtual inline void compute_contact_forces() override;

 protected:
  //! Mesh object
  using mpm::Contact<Tdim>::mesh_;
  //! MPI Size
  using mpm::Contact<Tdim>::mpi_size_;
  //! MPI rank
  using mpm::Contact<Tdim>::mpi_rank_;
};  // Contactfriction class
}  // namespace mpm

#include "contact_friction.tcc"

#endif  // _MPM_CONTACT_FRICTION_

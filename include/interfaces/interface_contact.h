#ifndef _MPM_INTERFACE_CONTACT_
#define _MPM_INTERFACE_CONTACT_

#include "interface.h"

namespace mpm {

//! InterfaceContact class
//! \brief InterfaceContact base class
//! \tparam Tdim Dimension
template <unsigned Tdim>
class InterfaceContact : public Interface<Tdim> {
 public:
  //! Default constructor with mesh class
  InterfaceContact(const std::shared_ptr<mpm::Mesh<Tdim>>& mesh);

  //! Intialize
  virtual inline void initialise() override;

  //! Compute contact forces
  virtual inline void compute_contact_forces() override;

 protected:
  //! Mesh object
  using mpm::Interface<Tdim>::mesh_;
  //! MPI Size
  using mpm::Interface<Tdim>::mpi_size_;
  //! MPI rank
  using mpm::Interface<Tdim>::mpi_rank_;
};  // Interfacecontact class
}  // namespace mpm

#include "interface_contact.tcc"

#endif  // _MPM_INTERFACECONTACT_

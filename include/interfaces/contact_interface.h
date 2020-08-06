#ifndef _MPM_CONTACT_INTERFACE_
#define _MPM_CONTACT_INTERFACE_

#include "interface.h"

namespace mpm {

//! ContactInterface class
//! \brief Contactinterface base class
//! \tparam Tdim Dimension
template <unsigned Tdim>
class ContactInterface : public Interface<Tdim> {
 public:
  //! Default constructor with mesh class
  ContactInterface(const std::shared_ptr<mpm::Mesh<Tdim>>& mesh);

  //! Intialize
  virtual void initialise() override;

  //! Compute contact forces
  virtual void compute_contact_forces() override;

 protected:
  //! Mesh object
  using mpm::Interface<Tdim>::mesh_;
  //! MPI Size
  using mpm::Interface<Tdim>::mpi_size_;
  //! MPI rank
  using mpm::Interface<Tdim>::mpi_rank_;
};  // Contactinterface class
}  // namespace mpm

#include "contact_interface.tcc"

#endif  // _MPM_CONTACTINTERFACE_

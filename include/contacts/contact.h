#ifndef MPM_CONTACT_H_
#define MPM_CONTACT_H_

#include "mesh.h"

namespace mpm {

//! Contact class
//! \brief Contact base class
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Contact {
 public:
  //! Default constructor with mesh class
  Contact(const std::shared_ptr<mpm::Mesh<Tdim>>& mesh);

  //! Intialize
  virtual inline void initialise(){};

  //! Compute contact forces
  virtual inline void compute_contact_forces(){};

 protected:
  //! Mesh object
  std::shared_ptr<mpm::Mesh<Tdim>> mesh_;
};  // Contact class
}  // namespace mpm

#include "contact.tcc"

#endif  // MPM_CONTACT_H_

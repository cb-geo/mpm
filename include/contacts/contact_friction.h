#ifndef MPM_CONTACT_FRICTION_H_
#define MPM_CONTACT_FRICTION_H_

#include "contact.h"

namespace mpm {

//! ContactFriction class
//! \brief ContactFriction base class
//! \tparam Tdim Dimension
template <unsigned Tdim>
class ContactFriction : public Contact<Tdim> {
 public:
  //! Default constructor with mesh class
  ContactFriction(const std::shared_ptr<mpm::Mesh<Tdim>>& mesh, double friction,
                  std::string normal_type);

  //! Intialize
  virtual inline void initialise() override;

  //! Compute contact forces
  virtual inline void compute_contact_forces() override;

 protected:
  //! Mesh object
  using mpm::Contact<Tdim>::mesh_;
  //! Coefficient of friction
  double friction_{.0};
  //! Type of normal detection
  std::string normal_type_{"PVG"};
};  // Contactfriction class
}  // namespace mpm

#include "contact_friction.tcc"

#endif  // MPM_CONTACT_FRICTION_H_

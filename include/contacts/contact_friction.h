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
  ContactFriction(const std::shared_ptr<mpm::Mesh<Tdim>>& mesh, double friction);

  //! Intialize
  virtual inline void initialise() override;

  //! Compute nodal kinematics
  virtual inline void compute_nodal_kinematics() override;

  //! Compute contact forces
  //! \param[in] gravity Gravity vector
  //! \param[in] phase Index to indicate material phase
  //! \param[in] time Current time in the simulation
  //! \param[in] concentrated_nodal_forces Boolean for if a concentrated force
  //! is applied or not
  virtual inline void compute_contact_forces(
      const Eigen::Matrix<double, Tdim, 1>& gravity, unsigned phase,
      double time, bool concentrated_nodal_forces) override;

  //! Compute contact nodal kinematics
  //! \param[in] dt Timestep in analysis
  virtual inline void compute_contact_kinematics(double dt) override;

  //! Update particle position
  //! \param[in] dt Timestep in analysis
  //! \param[in] velocity_update Update particle velocity from nodal vel
  virtual inline void update_particles_contact(double dt,
                                               bool velocity_update) override;

 protected:
  //! Mesh object
  using mpm::Contact<Tdim>::mesh_;
  //! Coefficient of friction
  double friction_{.0};
};  // Contactfriction class
}  // namespace mpm

#include "contact_friction.tcc"

#endif  // MPM_CONTACT_FRICTION_H_

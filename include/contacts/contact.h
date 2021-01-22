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

  //! Compute nodal kinematics
  virtual inline void compute_nodal_kinematics(){};

  //! Compute contact forces
  //! \param[in] gravity Gravity vector
  //! \param[in] phase Index to indicate material phase
  //! \param[in] time Current time in the simulation
  //! \param[in] concentrated_nodal_forces Boolean for if a concentrated force
  //! is applied or not
  virtual inline void compute_contact_forces(
      const Eigen::Matrix<double, Tdim, 1>& gravity, unsigned phase,
      double time, bool concentrated_nodal_forces){};

  //! Compute contact nodal kinematics
  //! \param[in] dt Timestep in analysis
  virtual inline void compute_contact_kinematics(double dt){};

  //! Update particle position
  //! \param[in] dt Timestep in analysis
  //! \param[in] velocity_update Update particle velocity from nodal vel
  virtual inline void update_particles_contact(double dt, bool velocity_update){};

 protected:
  //! Mesh object
  std::shared_ptr<mpm::Mesh<Tdim>> mesh_;
};  // Contact class
}  // namespace mpm

#include "contact.tcc"

#endif  // MPM_CONTACT_H_

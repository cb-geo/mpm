#ifndef MPM_CONSTRAINTS_H_
#define MPM_CONSTRAINTS_H_

#include <memory>

#include "absorbing_constraint.h"
#include "friction_constraint.h"
#include "logger.h"
#include "mesh.h"
#include "velocity_constraint.h"

namespace mpm {

//! Constraints class to store velocity and frictional constraints
//! \brief Constraint class to store a constraint on mesh
template <unsigned Tdim>
class Constraints {
 public:
  // Constructor with mesh as input argument
  Constraints(std::shared_ptr<mpm::Mesh<Tdim>> mesh) {
    mesh_ = mesh;
    console_ =
        std::make_unique<spdlog::logger>("Constraints", mpm::stdout_sink);
  }

  //! Assign nodal velocity constraints
  //! \param[in] setid Node set id
  //! \param[in] velocity_constraints Velocity constraint at node, dir, velocity
  bool assign_nodal_velocity_constraint(
      int set_id, const std::shared_ptr<mpm::VelocityConstraint>& constraint);

  //! Assign velocity constraints to nodes
  //! \param[in] velocity_constraints Constraint at node, dir, and velocity
  bool assign_nodal_velocity_constraints(
      const std::vector<std::tuple<mpm::Index, unsigned, double>>&
          velocity_constraints);

  //! Assign nodal frictional constraints
  //! \param[in] setid Node set id
  //! \param[in] friction_constraints Constraint at node, dir, sign, friction
  bool assign_nodal_frictional_constraint(
      int nset_id,
      const std::shared_ptr<mpm::FrictionConstraint>& fconstraints);

  //! Assign friction constraints to nodes
  //! \param[in] friction_constraints Constraint at node, dir, sign, and
  //! friction
  bool assign_nodal_friction_constraints(
      const std::vector<std::tuple<mpm::Index, unsigned, int, double>>&
          friction_constraints);

  //! Assign nodal absorbing constraints
  //! \param[in] setid Node set id
  //! \param[in] absorbing_constraints Constraint at node, dir, delta, a, b
  bool assign_nodal_absorbing_constraint(
      int nset_id,
      const std::shared_ptr<mpm::AbsorbingConstraint>& aconstraints);

  //! Assign absorbing constraints to nodes
  //! \param[in] absorbing_constraints Constraint at node, dir, delta, a, and b
  bool assign_nodal_absorbing_constraints(
      const std::vector<std::tuple<mpm::Index, unsigned, double, double,
                                   double>>& absorbing_constraints);

 private:
  //! Mesh object
  std::shared_ptr<mpm::Mesh<Tdim>> mesh_;
  //! Logger
  std::unique_ptr<spdlog::logger> console_;
};
}  // namespace mpm

#include "constraints.tcc"

#endif  // MPM_CONSTRAINTS_H_

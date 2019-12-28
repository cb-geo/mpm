#ifndef MPM_VELOCITY_CONSTRAINT_H_
#define MPM_VELOCITY_CONSTRAINT_H_

namespace mpm {

//! VelocityConstraint class to store velocity constraint on a set
//! \brief VelocityConstraint class to store a constraint on a set
//! \details VelocityConstraint stores the constraint as a static value
class VelocityConstraint {
 public:
  // Constructor
  //! \param[setid] setid  set id
  //! \param[dir] dir Direction of constraint load
  //! \param[velocity] velocity Constraint  velocity
  VelocityConstraint(int setid, unsigned dir, double velocity)
      : setid_{setid}, dir_{dir}, velocity_{velocity} {};

  // Set id
  int setid() const { return setid_; }

  // Direction
  unsigned dir() const { return dir_; }

  // Return velocity
  double velocity() const { return velocity_; }

 private:
  // ID
  int setid_;
  // Direction
  unsigned dir_;
  // Velocity
  double velocity_;
};
}  // namespace mpm
#endif  // MPM_VELOCITY_CONSTRAINT_H_

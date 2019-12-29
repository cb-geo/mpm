#ifndef MPM_FRICTION_CONSTRAINT_H_
#define MPM_FRICTION_CONSTRAINT_H_

namespace mpm {

//! FrictionConstraint class to store friction constraint on a set
//! \brief FrictionConstraint class to store a constraint on a set
//! \details FrictionConstraint stores the constraint as a static value
class FrictionConstraint {
 public:
  // Constructor
  //! \param[in] setid  set id
  //! \param[in] dir Direction of constraint load
  //! \param[in] sign_n Sign of normal vector
  //! \param[in] friction Constraint  friction
  FrictionConstraint(int setid, unsigned dir, int sign_n, double friction)
      : setid_{setid}, dir_{dir}, sign_n_{sign_n}, friction_{friction} {};

  // Set id
  int setid() const { return setid_; }

  // Direction
  unsigned dir() const { return dir_; }

  // Sign of normal component
  int sign_n() const { return sign_n_; }

  // Return friction
  double friction() const { return friction_; }

 private:
  // ID
  int setid_;
  // Direction
  unsigned dir_;
  // Sign
  int sign_n_;
  // Friction
  double friction_;
};
}  // namespace mpm
#endif  // MPM_FRICTION_CONSTRAINT_H_

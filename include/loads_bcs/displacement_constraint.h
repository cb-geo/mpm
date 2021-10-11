#ifndef MPM_DISPLACEMENT_CONSTRAINT_H_
#define MPM_DISPLACEMENT_CONSTRAINT_H_

namespace mpm {

//! DisplacementConstraint class to store displacement constraint on a set
//! \brief DisplacementConstraint class to store a displacement constraint on a
//! set \details DisplacementConstraint stores the displacement constraint as a
//! static value
class DisplacementConstraint {
 public:
  // Constructor
  //! \param[in] setid  set id
  //! \param[in] dir Direction of displacement constraint load
  //! \param[in] velocity Constraint displacement
  DisplacementConstraint(int setid, unsigned dir, double displacement)
      : setid_{setid}, dir_{dir}, displacement_{displacement} {};

  // Set id
  int setid() const { return setid_; }

  // Direction
  unsigned dir() const { return dir_; }

  // Return velocity
  double displacement() const { return displacement_; }

 private:
  // ID
  int setid_;
  // Direction
  unsigned dir_;
  // Velocity
  double displacement_;
};
}  // namespace mpm
#endif  // MPM_DISPLACEMENT_CONSTRAINT_H_

#ifndef MPM_CONSTRAINT_H_
#define MPM_CONSTRAINT_H_

namespace mpm {

//! Constraint class to store constraint on a set
//! \brief Constraint class to store a constraint on a set
//! \details Constraint stores the constraint as a static value
class Constraint {
 public:
  // Constructor
  //! \param[setid] setid  set id
  //! \param[dir] dir Direction of constraint load
  //! \param[constraint] constraint  constraint
  Constraint(int setid, unsigned dir, double constraint)
      : setid_{setid}, dir_{dir}, constraint_{constraint} {};

  // Set id
  int setid() const { return setid_; }

  // Direction
  unsigned dir() const { return dir_; }

  // Return constraint
  double value() const { return constraint_; }

 private:
  // ID
  int setid_;
  // Direction
  unsigned dir_;
  // Constraint
  double constraint_;
};
}  // namespace mpm
#endif  // MPM_CONSTRAINT_H_

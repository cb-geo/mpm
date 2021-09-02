#ifndef MPM_PRESSURE_CONSTRAINT_H_
#define MPM_PRESSURE_CONSTRAINT_H_

namespace mpm {

//! PressureConstraint class to store pressure constraint on a set
//! \brief PressureConstraint class to store a constraint on a set
//! \details PressureConstraint stores the constraint as a static value
class PressureConstraint {
 public:
  // Constructor
  //! \param[in] setid  set id
  //! \param[in] pressure Constraint pressure
  PressureConstraint(int setid, unsigned phase, double pressure)
      : setid_{setid}, phase_{phase}, pressure_{pressure} {};

  // Set id
  int setid() const { return setid_; }

  // Return phase
  unsigned phase() const { return phase_; }

  // Return pressure
  double pressure() const { return pressure_; }

 private:
  // ID
  int setid_;
  // Phase
  unsigned phase_;
  // Velocity
  double pressure_;
};
}  // namespace mpm
#endif  // MPM_PRESSURE_CONSTRAINT_H_

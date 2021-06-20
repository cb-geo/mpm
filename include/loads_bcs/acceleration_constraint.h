#ifndef MPM_ACCELERATION_CONSTRAINT_H_
#define MPM_ACCELERATION_CONSTRAINT_H_

namespace mpm {

//! AccelerationConstraint class to store acceleration constraint on a set
//! \brief AccelerationConstraint class to store a constraint on a set
//! \details AccelerationConstraint stores the constraint as a static value
class AccelerationConstraint {
 public:
  // Constructor
  //! \param[in] setid  set id
  //! \param[in] acceleration_fn Math function if defined
  //! \param[in] dir Direction of constraint load
  //! \param[in] acceleration Constraint acceleration
  AccelerationConstraint(
      int setid, const std::shared_ptr<mpm::FunctionBase>& acceleration_fn,
      unsigned dir, double acceleration)
      : setid_{setid},
        acceleration_fn_{acceleration_fn},
        dir_{dir},
        acceleration_{acceleration} {};

  // Set id
  int setid() const { return setid_; }

  // Direction
  unsigned dir() const { return dir_; }

  // Return acceleration
  double acceleration(double current_time) const {
    // Static load when no math function is defined
    double scalar = (this->acceleration_fn_ != nullptr)
                        ? (this->acceleration_fn_)->value(current_time)
                        : 1.0;
    return acceleration_ * scalar;
  }

 private:
  // ID
  int setid_;
  // Math function
  std::shared_ptr<mpm::FunctionBase> acceleration_fn_;
  // Direction
  unsigned dir_;
  // Acceleration
  double acceleration_;
};
}  // namespace mpm
#endif  // MPM_ACCELERATION_CONSTRAINT_H_

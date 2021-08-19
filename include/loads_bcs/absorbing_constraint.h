#ifndef MPM_ABSORBING_CONSTRAINT_H_
#define MPM_ABSORBING_CONSTRAINT_H_

namespace mpm {

//! AbsorbingConstraint class to store friction constraint on a set
//! \brief AbsorbingConstraint class to store a constraint on a set
//! \details AbsorbingConstraint stores the constraint as a static value
class AbsorbingConstraint {
 public:
  // Constructor
  //! \param[in] setid  set id
  //! \param[in] dir Normal direction of boundary application
  //! \param[in] delta Virtual viscous layer thickness
  //! \param[in] h_min Cell height
  //! \param[in] a Equation constant
  //! \param[in] b Equation constant
  //! \param[in] position Nodal position along boundary
  AbsorbingConstraint(int setid, unsigned dir, double delta, double h_min,
                      double a = 1., double b = 1.,
                      std::string position = "null")
      : setid_{setid},
        dir_{dir},
        delta_{delta},
        h_min_{h_min},
        a_{a},
        b_{b},
        position_{position} {};

  // Set id
  int setid() const { return setid_; }

  // Direction
  unsigned dir() const { return dir_; }

  // Virtual viscous layer thickness
  double delta() const { return delta_; }

  // Cell height
  double h_min() const { return h_min_; }

  // Return a
  double a() const { return a_; }

  // Return b
  double b() const { return b_; }

  // Return position
  std::string position() const { return position_; }

 private:
  // ID
  int setid_;
  // Direction
  unsigned dir_;
  // Delta
  double delta_;
  // Cell height
  double h_min_;
  // a
  double a_;
  // b
  double b_;
  // Node position
  std::string position_;
};
}  // namespace mpm
#endif  // MPM_ABSORBING_CONSTRAINT_H_

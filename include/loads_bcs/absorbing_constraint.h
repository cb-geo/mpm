#ifndef MPM_ABSORBING_CONSTRAINT_H_
#define MPM_ABSORBING_CONSTRAINT_H_

#include "data_types.h"

namespace mpm {

//! AbsorbingConstraint class to store friction constraint on a set
//! \brief AbsorbingConstraint class to store a constraint on a set
//! \details AbsorbingConstraint stores the constraint as a static value
class AbsorbingConstraint {
 public:
  // Constructor
  //! \param[in] setid  set id
  //! \param[in] dir Direction of p-wave propagation in model
  //! \param[in] delta Virtual viscous layer thickness
  //! \param[in] h_min Characteristic length (cell height)
  //! \param[in] a Dimensionless dashpot weight factor, p-wave
  //! \param[in] b Dimensionless dashpot weight factor, s-wave
  //! \param[in] position Nodal position along boundary
  AbsorbingConstraint(int setid, unsigned dir, double delta, double h_min,
                      double a = 1., double b = 1.,
                      mpm::Position position = mpm::Position::None)
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
  mpm::Position position() const { return position_; }

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
  mpm::Position position_;
};
}  // namespace mpm
#endif  // MPM_ABSORBING_CONSTRAINT_H_

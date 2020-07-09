#ifndef MPM_PROPERTIES_H_
#define MPM_PROPERTIES_H_

namespace mpm {
namespace properties {
//! Scalar Properties
enum Scalar : unsigned int { Mass, Volume, MassDensity, Pressure };
//! Vector Properties
enum Vector : unsigned int {
  Displacement,
  Velocity,
  Acceleration,
  Momentum,
  ExternalForce,
  InternalForce
};
}  // namespace properties
}  // namespace mpm

#endif  // MPM_PROPERTIES_H_

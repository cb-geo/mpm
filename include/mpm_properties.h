#ifndef MPM_PROPERTIES_H_
#define MPM_PROPERTIES_H_

namespace mpm {
namespace properties {
//! Boolean Properties
enum Boolean : unsigned int { SetTraction, Friction, GenericBC };
//! Scalar Properties
enum Scalar : unsigned int {
  Mass,
  Volume,
  MassDensity,
  MassPressure,
  Pressure
};
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

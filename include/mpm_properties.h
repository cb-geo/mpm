#ifndef MPM_PROPERTIES_H_
#define MPM_PROPERTIES_H_

namespace mpm {
namespace properties {
//! Scalar Properties
enum Scalar : unsigned int {
  Mass = 0,
  Volume = 1,
  MassDensity = 2,
  MassPressure = 3,
  Pressure = 4
};
//! Vector Properties
enum Vector : unsigned int {
  Displacement = 0,
  Velocity = 1,
  Acceleration = 2,
  Momentum = 3,
  ExternalForce = 4,
  InternalForce = 5
};
}  // namespace properties
}  // namespace mpm

#endif  // MPM_PROPERTIES_H_

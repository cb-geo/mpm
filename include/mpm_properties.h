#ifndef MPM_PROPERTIES_H_
#define MPM_PROPERTIES_H_

namespace mpm {
namespace properties {
//! Scalar Properties
enum Scalar : unsigned int { Mass };
//! Vector Properties
enum Vector : unsigned int { Velocity, Displacement, Acceleration };
}  // namespace properties
}  // namespace mpm

#endif  // MPM_PROPERTIES_H_

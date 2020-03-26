#ifndef MPM_RADIAL_BASIS_FUNCTION_H_
#define MPM_RADIAL_BASIS_FUNCTION_H_

#include <vector>

#include "Eigen/Dense"

#include "logger.h"

namespace mpm {

// Namespace for radial basis function handling
namespace RadialBasisFunction {

//! Cubic Spline Radial Basis Function
//! Source: Monaghan, 1985; Monaghan, 1992
template <unsigned Tdim>
double cubic_spline(const double& smoothing_length,
                    const double& norm_distance) {

  // Assign multiplier depends on dimension
  double multiplier;
  if (Tdim == 2)
    multiplier = 15.0 / (7.0 * M_PI * std::pow(smoothing_length, 2));
  else if (Tdim == 3)
    multiplier = 3.0 / (2.0 * M_PI * std::pow(smoothing_length, 3));
  else
    throw std::runtime_error("Tdim is invalid");

  // Compute basis function
  double basis_function = multiplier;
  const double radius = norm_distance / smoothing_length;
  if (radius >= 0.0 && radius < 1.0)
    basis_function *=
        2.0 / 3.0 - std::pow(radius, 2) + 0.5 * std::pow(radius, 3);
  else if (radius >= 1.0 && radius < 2.0)
    basis_function *= 1.0 / 6.0 * std::pow((2.0 - radius), 3);
  else
    basis_function = 0.0;

  return basis_function;
}

//! Quintic Spline Radial Basis Function
//! Source: Liu, 2010
template <unsigned Tdim>
double quintic_spline(const double& smoothing_length,
                      const double& norm_distance) {

  // Assign multiplier depends on dimension
  double multiplier;
  if (Tdim == 2)
    multiplier = 1.0 / (478.0 * M_PI * std::pow(smoothing_length, 2));
  else if (Tdim == 3)
    multiplier = 3.0 / (359.0 * M_PI * std::pow(smoothing_length, 3));
  else
    throw std::runtime_error("Tdim is invalid");

  // Compute basis function
  double basis_function = multiplier;
  const double radius = norm_distance / smoothing_length;
  if (radius >= 0.0 && radius < 1.0)
    basis_function *=
        (std::pow(3.0 - radius, 5) - 6 * std::pow(2.0 - radius, 5) +
         15 * std::pow(1.0 - radius, 5));
  else if (radius >= 1.0 && radius < 2.0)
    basis_function *=
        (std::pow(3.0 - radius, 5) - 6 * std::pow(2.0 - radius, 5));
  else if (radius >= 2.0 && radius < 3.0)
    basis_function *= (std::pow(3.0 - radius, 5));
  else
    basis_function = 0.0;

  return basis_function;
}

//! Gaussian Kernel
//! Source: Liu, 2010
template <unsigned Tdim>
double gaussian_kernel(const double& smoothing_length,
                       const double& norm_distance) {

  // Assign multiplier depends on dimension
  double multiplier;
  if (Tdim == 2)
    multiplier = 1.0 / (M_PI * std::pow(smoothing_length, 2));
  else if (Tdim == 3)
    multiplier = 1.0 / std::pow((std::sqrt(M_PI) * smoothing_length), 3);
  else
    throw std::runtime_error("Tdim is invalid");

  // Compute basis function
  double basis_function = multiplier;
  const double radius = norm_distance / smoothing_length;
  if (radius >= 0.0 && radius < 3.0)
    basis_function *= std::exp(-std::pow(radius, 2));
  else
    basis_function = 0.0;

  return basis_function;
}

//! Super Gaussian Kernel
//! Source: Monaghan, 1992
template <unsigned Tdim>
double super_gaussian_kernel(const double& smoothing_length,
                             const double& norm_distance) {

  // Assign multiplier depends on dimension
  double multiplier;
  if (Tdim == 2)
    multiplier = 1.0 / std::pow((std::sqrt(M_PI) * smoothing_length), 2);
  else if (Tdim == 3)
    multiplier = 1.0 / std::pow((std::sqrt(M_PI) * smoothing_length), 3);
  else
    throw std::runtime_error("Tdim is invalid");

  // Compute basis function
  double basis_function = multiplier;
  const double radius = norm_distance / smoothing_length;
  if (radius >= 0.0 && radius < 3.0)
    basis_function *=
        std::exp(-std::pow(radius, 2)) * (Tdim / 2.0 + 1.0 - radius * radius);
  else
    basis_function = 0.0;

  return basis_function;
}

//! General Radial Basis Function Kernel call
template <unsigned Tdim>
double kernel(const double& smoothing_length, const double& norm_distance,
              const std::string type) {
  if (type == "cubic_spline") {
    return mpm::RadialBasisFunction::cubic_spline<Tdim>(smoothing_length,
                                                        norm_distance);
  } else if (type == "quintic_spline") {
    return mpm::RadialBasisFunction::quintic_spline<Tdim>(smoothing_length,
                                                          norm_distance);
  } else if (type == "gaussian") {
    return mpm::RadialBasisFunction::gaussian_kernel<Tdim>(smoothing_length,
                                                           norm_distance);
  } else if (type == "super_gaussian") {
    return mpm::RadialBasisFunction::super_gaussian_kernel<Tdim>(
        smoothing_length, norm_distance);
  } else {
    throw std::runtime_error(
        "RadialBasisFunction kernel type is invalid. Available types are: "
        "\"cubic_spline\", \"quintic_spline\", \"gaussian\", and, "
        "\"super_gaussian\".");
  }
}

}  // namespace RadialBasisFunction

}  // namespace mpm
#endif
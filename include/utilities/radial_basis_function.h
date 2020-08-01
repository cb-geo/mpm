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
double cubic_spline(const double smoothing_length, const double norm_distance) {

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
        (2.0 / 3.0 - std::pow(radius, 2) + 0.5 * std::pow(radius, 3));
  else if (radius >= 1.0 && radius < 2.0)
    basis_function *= (1.0 / 6.0 * std::pow((2.0 - radius), 3));
  else
    basis_function = 0.0;

  return basis_function;
}

//! Cubic Spline Radial Basis Function derivative
//! Source: Monaghan, 1985; Monaghan, 1992
template <unsigned Tdim>
double cubic_spline_derivative(const double smoothing_length,
                               const double norm_distance) {

  // Assign multiplier depends on dimension
  double multiplier;
  if (Tdim == 2)
    multiplier = 15.0 / (7.0 * M_PI * std::pow(smoothing_length, 2));
  else if (Tdim == 3)
    multiplier = 3.0 / (2.0 * M_PI * std::pow(smoothing_length, 3));
  else
    throw std::runtime_error("Tdim is invalid");

  // Compute basis function derivative
  double dw_dr = multiplier;
  const double radius = norm_distance / smoothing_length;
  if (radius >= 0.0 && radius < 1.0)
    dw_dr *= (-2.0 * radius + 1.5 * std::pow(radius, 2));
  else if (radius >= 1.0 && radius < 2.0)
    dw_dr *= -0.5 * std::pow((2.0 - radius), 2);
  else
    dw_dr = 0.0;

  return dw_dr;
}

//! Quintic Spline Radial Basis Function
//! Source: Liu, 2010
template <unsigned Tdim>
double quintic_spline(const double smoothing_length,
                      const double norm_distance) {

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
        (std::pow(3.0 - radius, 5) - 6.0 * std::pow(2.0 - radius, 5) +
         15.0 * std::pow(1.0 - radius, 5));
  else if (radius >= 1.0 && radius < 2.0)
    basis_function *=
        (std::pow(3.0 - radius, 5) - 6.0 * std::pow(2.0 - radius, 5));
  else if (radius >= 2.0 && radius < 3.0)
    basis_function *= (std::pow(3.0 - radius, 5));
  else
    basis_function = 0.0;

  return basis_function;
}

//! Quintic Spline Radial Basis Function derivative
//! Source: Liu, 2010
template <unsigned Tdim>
double quintic_spline_derivative(const double smoothing_length,
                                 const double norm_distance) {

  // Assign multiplier depends on dimension
  double multiplier;
  if (Tdim == 2)
    multiplier = 1.0 / (478.0 * M_PI * std::pow(smoothing_length, 2));
  else if (Tdim == 3)
    multiplier = 3.0 / (359.0 * M_PI * std::pow(smoothing_length, 3));
  else
    throw std::runtime_error("Tdim is invalid");

  // Compute basis function
  double dw_dr = multiplier;
  const double radius = norm_distance / smoothing_length;
  if (radius >= 0.0 && radius < 1.0)
    dw_dr *=
        (-5.0 * std::pow(3.0 - radius, 4) + 30. * std::pow(2.0 - radius, 4) -
         75. * std::pow(1.0 - radius, 4));
  else if (radius >= 1.0 && radius < 2.0)
    dw_dr *=
        (-5.0 * std::pow(3.0 - radius, 4) + 30. * std::pow(2.0 - radius, 4));
  else if (radius >= 2.0 && radius < 3.0)
    dw_dr *= (-5.0 * std::pow(3.0 - radius, 4));
  else
    dw_dr = 0.0;

  return dw_dr;
}

//! Gaussian Kernel
//! Source: Liu, 2010
template <unsigned Tdim>
double gaussian(const double smoothing_length, const double norm_distance) {

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

//! Gaussian Kernel derivative
//! Source: Liu, 2010
template <unsigned Tdim>
double gaussian_derivative(const double smoothing_length,
                           const double norm_distance) {

  // Assign multiplier depends on dimension
  double multiplier;
  if (Tdim == 2)
    multiplier = 1.0 / (M_PI * std::pow(smoothing_length, 2));
  else if (Tdim == 3)
    multiplier = 1.0 / std::pow((std::sqrt(M_PI) * smoothing_length), 3);
  else
    throw std::runtime_error("Tdim is invalid");

  // Compute basis function
  double dw_dr = multiplier;
  const double radius = norm_distance / smoothing_length;
  if (radius >= 0.0 && radius < 3.0)
    dw_dr *= -2.0 * radius * std::exp(-std::pow(radius, 2));
  else
    dw_dr = 0.0;

  return dw_dr;
}

//! Super Gaussian Kernel
//! Source: Monaghan, 1992
template <unsigned Tdim>
double super_gaussian(const double smoothing_length,
                      const double norm_distance) {

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
    basis_function *= radius * (2.0 * radius * radius - (double)Tdim - 4.0) *
                      std::exp(-std::pow(radius, 2));
  else
    basis_function = 0.0;

  return basis_function;
}

//! Super Gaussian Kernel derivative
//! Source: Monaghan, 1992
template <unsigned Tdim>
double super_gaussian_derivative(const double smoothing_length,
                                 const double norm_distance) {

  // Assign multiplier depends on dimension
  double multiplier;
  if (Tdim == 2)
    multiplier = 1.0 / std::pow((std::sqrt(M_PI) * smoothing_length), 2);
  else if (Tdim == 3)
    multiplier = 1.0 / std::pow((std::sqrt(M_PI) * smoothing_length), 3);
  else
    throw std::runtime_error("Tdim is invalid");

  // Compute basis function
  double dw_dr = multiplier;
  const double radius = norm_distance / smoothing_length;
  if (radius >= 0.0 && radius < 3.0)
    dw_dr *= std::exp(-std::pow(radius, 2)) *
             ((double)Tdim / 2.0 + 1.0 - radius * radius);
  else
    dw_dr = 0.0;

  return dw_dr;
}

//! General Radial Basis Function Kernel call
template <unsigned Tdim>
double kernel(const double smoothing_length, const double norm_distance,
              const std::string type = "cubic_spline") {
  if (type == "cubic_spline") {
    return mpm::RadialBasisFunction::cubic_spline<Tdim>(smoothing_length,
                                                        norm_distance);
  } else if (type == "quintic_spline") {
    return mpm::RadialBasisFunction::quintic_spline<Tdim>(smoothing_length,
                                                          norm_distance);
  } else if (type == "gaussian") {
    return mpm::RadialBasisFunction::gaussian<Tdim>(smoothing_length,
                                                    norm_distance);
  } else if (type == "super_gaussian") {
    return mpm::RadialBasisFunction::super_gaussian<Tdim>(smoothing_length,
                                                          norm_distance);
  } else {
    throw std::runtime_error(
        "RadialBasisFunction kernel type is invalid. Available types are: "
        "\"cubic_spline\", \"quintic_spline\", \"gaussian\", and, "
        "\"super_gaussian\".");
  }
}

//! General Radial Basis Function Kernel call
template <unsigned Tdim>
Eigen::Matrix<double, Tdim, 1> gradient(
    const double smoothing_length,
    const Eigen::Matrix<double, Tdim, 1>& relative_distance,
    const std::string type = "cubic_spline") {

  // Compute norm distance
  const double norm_distance = relative_distance.norm();
  double dw_dr;
  if (type == "cubic_spline") {
    dw_dr = mpm::RadialBasisFunction::cubic_spline_derivative<Tdim>(
        smoothing_length, norm_distance);
  } else if (type == "quintic_spline") {
    dw_dr = mpm::RadialBasisFunction::quintic_spline_derivative<Tdim>(
        smoothing_length, norm_distance);
  } else if (type == "gaussian") {
    dw_dr = mpm::RadialBasisFunction::gaussian_derivative<Tdim>(
        smoothing_length, norm_distance);
  } else if (type == "super_gaussian") {
    dw_dr = mpm::RadialBasisFunction::super_gaussian_derivative<Tdim>(
        smoothing_length, norm_distance);
  } else {
    throw std::runtime_error(
        "RadialBasisFunction gradient type is invalid. Available types are: "
        "\"cubic_spline\", \"quintic_spline\", \"gaussian\", and, "
        "\"super_gaussian\".");
  }

  // Gradient = dw_dr * r / ||r|| / h
  Eigen::Matrix<double, Tdim, 1> gradient = relative_distance;
  if (norm_distance > 1.e-12)
    gradient *= dw_dr / (norm_distance * smoothing_length);
  else
    gradient *= 0.0;

  return gradient;
}

}  // namespace RadialBasisFunction

}  // namespace mpm
#endif
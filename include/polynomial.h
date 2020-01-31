#ifndef MPM_POLYNOMIAL_H_
#define MPM_POLYNOMIAL_H_

#include <vector>

#include "Eigen/Dense"

#include "logger.h"

namespace mpm {

namespace Polynomial {

//! Base class that stores the definite integrals of monomials
//! over reference element
//! tparam Tdim Dimension
//! tparam Tporder Polynomial order
template <unsigned Tdim, unsigned Tporder, unsigned Tnterms>
struct IntegralBase {
  virtual Eigen::Matrix<double, Tnterms, 1> definite_integrals() const = 0;
  static const Eigen::Matrix<double, Tnterms, 1> Square_Definite_Integrals;
  static const Eigen::Matrix<double, Tnterms, 1> Tri_Definite_Integrals;
};

//! Evaluate monomials
//! \Details Computes the monomials of pth order polynomials for a given point
//! \Details 1st order polynomial= (1 + x) * (1 + y) * (1 + z)
//! \Details 2nd order polynomial= (1 + x + x^2) * (1 + y + y^2) * (1 + z + z^2)
template <unsigned Tdim>
inline Eigen::VectorXd evaluate_monomials(const unsigned porder,
                                          Eigen::Matrix<double, Tdim, 1> xi);

}  // namespace Polynomial
}  // namespace mpm

#endif

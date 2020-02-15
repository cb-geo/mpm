#ifndef MPM_POLYNOMIAL_INTERPOLATION_H_
#define MPM_POLYNOMIAL_INTERPOLATION_H_

#include <vector>

#include "Eigen/Dense"

#include "polynomial.h"

namespace mpm {

//! PolynomialInterpolation base class for polynomial interpolation
//! \brief Base class that interpolates a polynomial of given order
//! \tparam Tdim dimension
template <unsigned Tdim>
class PolynomialInterpolation {
 public:
  //! Constructor
  PolynomialInterpolation(){};

  //! destructor
  virtual ~PolynomialInterpolation(){};

  //! Delete copy constructor
  PolynomialInterpolation(const PolynomialInterpolation<Tdim>&) = delete;

  //! Delete assignement operator
  PolynomialInterpolation& operator=(const PolynomialInterpolation<Tdim>&) =
      delete;

};  // PolynomialInterpolation class
}  // namespace mpm

#endif  // MPM_POLYNOMIAL_INTERPOLATION_H_

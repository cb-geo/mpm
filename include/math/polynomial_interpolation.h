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

  //! Initialise interpolator
  //! \param[in] pcoord Coordinates of the point of interest
  //! \param[in] data_points Coordinates of the set of data points
  //! \param[in] spline_order Weight spline order
  //! \param[in] poly_order Order of polynomial to be interpolated
  //! \param[in] span Distance from the point of interest
  virtual void initialise(const VectorDim& pcoord,
                          const std::vector<VectorDim>& data_points,
                          unsigned spline_order, unsigned poly_order,
                          double span) = 0;

  //! Interpolate polynomial at a given point
  //! \param[in] data_values Known values of data points
  //! \retval Interpolated value at the point of interest
  virtual double interpolate_polynomial(
      const std::vector<double>& data_values) const = 0;

};  // PolynomialInterpolation class
}  // namespace mpm

#endif  // MPM_POLYNOMIAL_INTERPOLATION_H_

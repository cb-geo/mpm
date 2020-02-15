#ifndef MPM_MLS_POLYNOMIAL_INTERPOLATION_H_
#define MPM_MLS_POLYNOMIAL_INTERPOLATION_H_

#include <vector>

#include "Eigen/Dense"

#include "polynomial_interpolation.h"

namespace mpm {

//! MLSPolynomialInterpolation  class for polynomial interpolation
//! \brief Interpolates a polynomial using Moving Least Squares (MLS) method
//! \tparam Tdim dimension
//! \tparam Tnmonomials Number of monomials = (poly order + 1)^Tdim
template <unsigned Tdim, unsigned Tnmonomials>
class MLSPolyInterpolation : public PolynomialInterpolation<Tdim> {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Constructor
  //! \param[in] pcoord Coordinates of the point of interest
  //! \param[in] data_points Coordinates of the set of data points
  //! \param[in] spline_order Qeight spline order
  //! \param[in] poly_order Order of polynomial to be interpolated
  MLSPolyInterpolation(const VectorDim& pcoord,
                       const std::vector<VectorDim>& data_points,
                       unsigned spline_order, unsigned poly_order, double span);

  //! Destructor
  ~MLSPolyInterpolation() override{};

  //! Delete copy constructor
  MLSPolyInterpolation(const MLSPolyInterpolation<Tdim, Tnmonomials>&) = delete;

  //! Delete assignement operator
  MLSPolyInterpolation& operator=(
      const MLSPolyInterpolation<Tdim, Tnmonomials>&) = delete;

};  // MLSPolyInterpolation class
}  // namespace mpm

#include "mls_polyinterpolation.tcc"

#endif  // MPM_MLS_POLYNOMIAL_INTERPOLATION_H_

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
  MLSPolyInterpolation();

  //! Destructor
  ~MLSPolyInterpolation() override{};

  //! Delete copy constructor
  MLSPolyInterpolation(const MLSPolyInterpolation<Tdim, Tnmonomials>&) = delete;

  //! Delete assignement operator
  MLSPolyInterpolation& operator=(
      const MLSPolyInterpolation<Tdim, Tnmonomials>&) = delete;

  //! Initialise interpolator
  //! \param[in] pcoord Coordinates of the point of interest
  //! \param[in] data_points Coordinates of the set of data points
  //! \param[in] spline_order Weight spline order
  //! \param[in] poly_order Order of polynomial to be interpolated
  //! \param[in] span Distance from the point of interest
  void initialise(const VectorDim& pcoord,
                  const std::vector<VectorDim>& data_points,
                  unsigned spline_order, unsigned poly_order,
                  double span) override;

  //! Interpolate polynomial at a given point
  //! \param[in] data_values Known values associated with data points
  //! \retval Interpolated value at the point of interest
  double interpolate_polynomial(
      const std::vector<double>& data_values) const override;

 private:
  //! Initialise M matrix
  //! \param[in] datapoints Coordinates and values of data points
  //! \param[in] poly_order Polynomial order
  void initialise_M_matrix(const std::vector<VectorDim>& data_points,
                           unsigned poly_order);

  //! Initialise B matrix
  //! \param[in] datapoints Coordinates and values of data points
  //! \param[in] poly_order Polynomial order
  void initialise_B_vector(const std::vector<VectorDim>& data_points,
                           unsigned poly_order);

  //! Compute weights
  //! \param[in] pcoord Coordinates of the point of interest
  //! \param[in] data_points Coordinates of the set of data points
  //! \param[in] spline_order Qeight spline order
  void compute_weights(const VectorDim& pcoord,
                       const std::vector<VectorDim>& data_points,
                       unsigned spline_order, double span);

  //! Compute and return weight of a data point using cubic spline
  //! \param[in] distance Positive distance from the point of interest
  //! \retval Associated weight
  double cubic_spline_weight(const VectorDim& distance) const;

  //! Compute and return weight of a data point using quadratic spline
  //! \param[in] distance Positive distance from the point of interest
  //! \retval Associated weight
  double quadratic_spline_weight(const VectorDim& distance) const;

 private:
  //! M matrix
  Eigen::Matrix<double, Tnmonomials, Tnmonomials> M_matrix_;
  //! B matrix
  std::vector<Eigen::Matrix<double, Tnmonomials, 1>> B_vector_;
  //! Weights associated with each data point
  std::vector<double> weights_;
  //! Monomials at the point of interest
  Eigen::Matrix<double, Tnmonomials, 1> point_monomials_;
};  // MLSPolyInterpolation class
}  // namespace mpm

#include "mls_polyinterpolation.tcc"

#endif  // MPM_MLS_POLYNOMIAL_INTERPOLATION_H_

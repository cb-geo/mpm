// Constructor
template <unsigned Tdim, unsigned Tnmonomials>
mpm::MLSPolyInterpolation<Tdim, Tnmonomials>::MLSPolyInterpolation()
    : mpm::PolynomialInterpolation<Tdim>() {
  M_matrix_.setZero();
  B_vector_.clear();
  weights_.clear();
}

// Initialise interpolator
template <unsigned Tdim, unsigned Tnmonomials>
void mpm::MLSPolyInterpolation<Tdim, Tnmonomials>::initialise(
    const VectorDim& pcoord, const std::vector<VectorDim>& data_points,
    unsigned spline_order, unsigned poly_order, double span) {
  // Initialise monomial values of the polynomial at the given point
  this->point_monomials_ =
      mpm::Polynomial::evaluate_monomials<Tdim>(poly_order, pcoord);
  // Compute weights
  this->compute_weights(pcoord, data_points, spline_order, span);
  // Initialise M matrix for the given set of points
  this->initialise_M_matrix(data_points, poly_order);
  // Initialise B vector for the given set of points 
  this->initialise_B_vector(data_points, poly_order);
}

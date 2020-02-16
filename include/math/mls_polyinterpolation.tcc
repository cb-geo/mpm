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

  // Check spline order
  if (spline_order < 2 || spline_order > 3)
    throw std::runtime_error("MLS spline order is invalid");
  // Check coordinates list of data points
  if (data_points.empty())
    throw std::runtime_error("MLS Data points coordinates are empty");
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

// Compute associated weights of the data points set
template <unsigned Tdim, unsigned Tnmonomials>
void mpm::MLSPolyInterpolation<Tdim, Tnmonomials>::compute_weights(
    const VectorDim& point, const std::vector<VectorDim>& data_points,
    unsigned spline_order, double span) {

  weights_.clear();
  // Iterate over each data point and compute the associated weight
  for (const auto& pdata : data_points) {
    // Distance to data point
    VectorDim distance;
    for (unsigned dim = 0; dim < Tdim; ++dim)
      distance(dim) = std::fabs(pdata(dim) - point(dim)) / span;
    // Spline weight
    double weight = std::numeric_limits<double>::max();
    if (spline_order == 2) weight = this->quadratic_spline_weight(distance);
    if (spline_order == 3) weight = this->cubic_spline_weight(distance);
    weights_.emplace_back(weight);
  }
}

// Return cubic spline weight of the data point
template <unsigned Tdim, unsigned Tnmonomials>
double mpm::MLSPolyInterpolation<Tdim, Tnmonomials>::cubic_spline_weight(
    const Eigen::Matrix<double, Tdim, 1>& distance) const {

  Eigen::Matrix<double, Tdim, 1> weights;
  weights.setZero();
  for (unsigned dim = 0; dim < Tdim; ++dim) {
    if (distance(dim) <= 0.5)
      weights(dim) = 2. / 3. - (4. * distance(dim) * distance(dim)) +
                     (4. * distance(dim) * distance(dim) * distance(dim));
    else if (distance(dim) > 0.5 && distance(dim) < 1.0)
      weights(dim) = 4. / 3. - (4. * distance(dim)) +
                     (4. * distance(dim) * distance(dim)) -
                     (4. * distance(dim) * distance(dim) * distance(dim) / 3.);
    else
      weights(dim) = 0.;
  }
  return weights.prod();
}

// Return quadratic spline weight of the data point
template <unsigned Tdim, unsigned Tnmonomials>
double mpm::MLSPolyInterpolation<Tdim, Tnmonomials>::quadratic_spline_weight(
    const Eigen::Matrix<double, Tdim, 1>& distance) const {

  Eigen::Matrix<double, Tdim, 1> weights;
  weights.setZero();
  for (unsigned dim = 0; dim < Tdim; ++dim) {
    if (distance(dim) <= 0.5)
      weights(dim) = 3. / 4. - (distance(dim) * distance(dim));
    else if (distance(dim) > 0.5 && distance(dim) < 1.5)
      weights(dim) = 0.5 * (1.5 - distance(dim)) * (1.5 - distance(dim));
    else
      weights(dim) = 0.;
  }
  return weights.prod();
}

// Initialise M matrix
template <unsigned Tdim, unsigned Tnmonomials>
void mpm::MLSPolyInterpolation<Tdim, Tnmonomials>::initialise_M_matrix(
    const std::vector<VectorDim>& data_points, unsigned poly_order) {

  M_matrix_.setZero();
  unsigned p = 0;
  for (const auto& pdata_coord : data_points) {
    // Monomials
    Eigen::VectorXd monomials =
        mpm::Polynomial::evaluate_monomials<Tdim>(poly_order, pdata_coord);
    // Compute M matrix
    this->M_matrix_ += this->weights_[p] * monomials * monomials.transpose();
    ++p;
  }
  if (!std::fabs(M_matrix_.determinant()) > 0.)
    throw std::runtime_error("MLS M-matrix is non-invertible!");
}

// Initialise B vector
template <unsigned Tdim, unsigned Tnmonomials>
void mpm::MLSPolyInterpolation<Tdim, Tnmonomials>::initialise_B_vector(
    const std::vector<VectorDim>& data_points, unsigned poly_order) {

  B_vector_.clear();
  unsigned p = 0;
  for (const auto& pdata_coord : data_points) {
    // Monomials
    Eigen::VectorXd monomials =
        mpm::Polynomial::evaluate_monomials<Tdim>(poly_order, pdata_coord);
    // Compute M matrix
    B_vector_.emplace_back(weights_[p] * monomials);
    ++p;
  }
}

// Compute polynomial interpolation at the point of interest
template <unsigned Tdim, unsigned Tnmonomials>
double mpm::MLSPolyInterpolation<Tdim, Tnmonomials>::interpolate_polynomial(
    const std::vector<double>& data_values) const {
  double pvalue = std::numeric_limits<double>::max();

  if (data_values.size() != weights_.size() ||
      data_values.size() != B_vector_.size())
    throw std::runtime_error("MLS Data points are invalid");

  // B Matrix
  Eigen::Matrix<double, Tnmonomials, 1> B_matrix;
  B_matrix.setZero();
  for (unsigned i = 0; i < data_values.size(); ++i)
    B_matrix += data_values[i] * B_vector_[i];

  // Compute interpolated value at the point of interest
  pvalue = point_monomials_.transpose() * M_matrix_.inverse() * B_matrix;
  return pvalue;
}

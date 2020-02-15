template <unsigned Tdim>
inline Eigen::VectorXd mpm::Polynomial::evaluate_monomials(
    const unsigned porder, Eigen::Matrix<double, Tdim, 1> xi) {
  // number of terms in an univariate polynomial of order Tporder
  const unsigned unipoly_size = porder + 1;
  const unsigned nterms = pow((unipoly_size), Tdim);

  // multivariate polynomial is first written using univariate polynomials
  // Ex: 2nd order 2-dimensional multivariate polynomial
  // | 1   1   |
  // | x   y   |
  // | x^2 y^2 |
  Eigen::MatrixXd multi_polynomial(unipoly_size, Tdim);
  for (unsigned i = 0; i < unipoly_size; ++i) {
    for (unsigned j = 0; j < Tdim; ++j) multi_polynomial(i, j) = pow(xi(j), i);
  }

  // multivarate polynomial is now evaluated by multiplying the column vectors
  // Ex: monomials = (1 x x^2) * (1 y y^2)
  unsigned n = 0;
  Eigen::VectorXd monomials(nterms);

  // 2D evaluation of monomials
  if (Tdim == 2) {
    for (unsigned i = 0; i < unipoly_size; ++i) {
      for (unsigned j = 0; j < unipoly_size; ++j) {
        monomials(n) = multi_polynomial(i, 0) * multi_polynomial(j, 1);
        ++n;
      }
    }
  }

  // 3D evaluation of monomials
  n = 0;
  if (Tdim == 3) {
    for (unsigned i = 0; i < unipoly_size; ++i) {
      for (unsigned j = 0; j < unipoly_size; ++j) {
        for (unsigned k = 0; k < unipoly_size; ++k) {
          monomials(n) = multi_polynomial(i, 0) * multi_polynomial(j, 1) *
                         multi_polynomial(k, 2);
          ++n;
        }
      }
    }
  }
  return monomials;
}

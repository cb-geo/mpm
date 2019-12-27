//! Constructor with zero level set id and coefficient vector
template <unsigned Tdim>
mpm::LevelSet<Tdim>::LevelSet(unsigned id, const std::string& domain,
                              bool moving_statuss)
    : id_{id}, moving_{moving_statuss} {
  //! Logger
  std::string logger =
      "levelset" + std::to_string(Tdim) + "d::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);

  try {
    if (domain == "Integration_In")
      integration_domain_ = mpm::IntegrationDomain::Integration_In;
    else if (domain == "Integration_Out")
      integration_domain_ = mpm::IntegrationDomain::Integration_Out;
    else if (domain == "Boundary")
      integration_domain_ = mpm::IntegrationDomain::Boundary;
    else
      throw std::runtime_error(
          "Integration domain for level set is not properly specified, using "
          "default domain-in");
  } catch (std::exception& exception) {
    console_->warn("{} #{}: {}", __FILE__, __LINE__, exception.what());
  }
}

//! Evaluate signed distance function of a given point
template <unsigned Tdim>
double mpm::LevelSet<Tdim>::sign_distance(const VectorDim& point) const {
  // number of monomials
  const unsigned nterms = pow((poly_order_ + 1), Tdim);
  Eigen::VectorXd monomials =
      mpm::Polynomial::evaluate_monomials<Tdim>(poly_order_, point);

  double signed_distance = 0.0;
  for (unsigned i = 0; i < nterms; ++i)
    signed_distance += poly_coefficients_.at(i) * monomials(i);
  return signed_distance;
}

// TO DO
//! Evaluate normal vector to the zero level set at a given point
template <unsigned Tdim>
Eigen::Matrix<double, Tdim, 1> mpm::LevelSet<Tdim>::normal_vector(
    const VectorDim& point) const {
  VectorDim random;
  return random;
}

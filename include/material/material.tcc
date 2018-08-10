//! Get material property
template <unsigned Tdim>
double mpm::Material<Tdim>::property(const std::string& key) {
  double result = std::numeric_limits<double>::max();
  try {
    result = properties_[key].template get<double>();
  } catch (std::exception& except) {
    console_->error("Material parameter not found: {}", except.what());
  }
  return result;
}

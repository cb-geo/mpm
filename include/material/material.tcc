//! Get material property
template <unsigned Tdim>
double mpm::Material<Tdim>::property(const std::string& key) {
  double result = std::numeric_limits<double>::max();
  try {
    result = properties_[key].template get<double>();
  } catch (std::exception& except) {
    std::cerr << "Material parameter not found: " << except.what() << '\n';
  }
  return result;
}

//! Get material property
template <unsigned Tdim>
template <typename Ttype>
Ttype mpm::Material<Tdim>::property(const std::string& key) {
  try {
    return properties_[key].template get<Ttype>();
  } catch (std::exception& except) {
    console_->error("Property call to material parameter not found: {}",
                    except.what());
    throw std::runtime_error(
        "Property call to material parameter not found or invalid type");
  }
}

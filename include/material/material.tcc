//! Get material property
template <unsigned Tdim>
template <typename Targ>
Targ mpm::Material<Tdim>::property(const std::string& key) {
  try {
    return properties_[key].template get<Targ>();
  } catch (std::exception& except) {
    console_->error("Property call to material parameter not found: {}",
                    except.what());
    throw std::runtime_error(
        "Property call to material parameter not found or invalid type");
  }
}

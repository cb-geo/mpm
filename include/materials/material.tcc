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

//! Compute constitutive relations matrix
template <unsigned Tdim>
Eigen::Matrix<double, 6, 6> mpm::Material<Tdim>::compute_dmatrix(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr, mpm::dense_map* state_vars) {

  // Zero matrix
  return Eigen::Matrix<double, 6, 6>::Zero();
}
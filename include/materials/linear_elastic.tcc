//! Read material properties
template <unsigned Tdim>
mpm::LinearElastic<Tdim>::LinearElastic(unsigned id,
                                        const Json& material_properties)
    : Material<Tdim>(id, material_properties) {
  try {
    density_ = material_properties.at("density").template get<double>();
    youngs_modulus_ =
        material_properties.at("youngs_modulus").template get<double>();
    poisson_ratio_ =
        material_properties.at("poisson_ratio").template get<double>();

    // Calculate bulk modulus
    bulk_modulus_ = youngs_modulus_ / (3.0 * (1. - 2. * poisson_ratio_));

    // Special material properties
    if (material_properties.contains("earthquake")) {
      bool earthquake =
          material_properties.at("earthquake").template get<bool>();

      if (earthquake) {
        // Calculate constrained and shear modulus
        double constrained_modulus =
            youngs_modulus_ * (1. - poisson_ratio_) /
            ((1. + poisson_ratio_) * (1. - 2. * poisson_ratio_));
        double shear_modulus = youngs_modulus_ / (2.0 * (1. + poisson_ratio_));

        // Calculate wave velocities
        p_wave_velocity_ = sqrt(constrained_modulus / density_);
        s_wave_velocity_ = sqrt(shear_modulus / density_);
      }
    }

    properties_ = material_properties;
    // Set elastic tensor
    this->compute_elastic_tensor();
  } catch (Json::exception& except) {
    console_->error("Material parameter not set: {} {}\n", except.what(),
                    except.id);
  }
}

//! Return elastic tensor
template <unsigned Tdim>
bool mpm::LinearElastic<Tdim>::compute_elastic_tensor() {
  // Shear modulus
  const double G = youngs_modulus_ / (2.0 * (1. + poisson_ratio_));

  const double a1 = bulk_modulus_ + (4.0 / 3.0) * G;
  const double a2 = bulk_modulus_ - (2.0 / 3.0) * G;

  // clang-format off
  // compute elasticityTensor
  de_(0,0)=a1;    de_(0,1)=a2;    de_(0,2)=a2;    de_(0,3)=0;    de_(0,4)=0;    de_(0,5)=0;
  de_(1,0)=a2;    de_(1,1)=a1;    de_(1,2)=a2;    de_(1,3)=0;    de_(1,4)=0;    de_(1,5)=0;
  de_(2,0)=a2;    de_(2,1)=a2;    de_(2,2)=a1;    de_(2,3)=0;    de_(2,4)=0;    de_(2,5)=0;
  de_(3,0)= 0;    de_(3,1)= 0;    de_(3,2)= 0;    de_(3,3)=G;    de_(3,4)=0;    de_(3,5)=0;
  de_(4,0)= 0;    de_(4,1)= 0;    de_(4,2)= 0;    de_(4,3)=0;    de_(4,4)=G;    de_(4,5)=0;
  de_(5,0)= 0;    de_(5,1)= 0;    de_(5,2)= 0;    de_(5,3)=0;    de_(5,4)=0;    de_(5,5)=G;
  // clang-format on
  return true;
}

//! Compute stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::LinearElastic<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr, mpm::dense_map* state_vars) {
  const Vector6d dstress = this->de_ * dstrain;
  return (stress + dstress);
}

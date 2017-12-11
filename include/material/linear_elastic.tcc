//! Read material properties
//! \param[in] materail_properties Material properties
void mpm::LinearElastic::properties(const Json& materail_properties) {
  try {
    youngs_modulus_ = materail_properties["youngs_modulus"].template get<double>();
    poisson_ratio_ = materail_properties["poisson_ratio"].template get<double>();
  } catch (std::exception& except) {
    std::cerr << "Material parameter not set: " << except.what() << '\n';
  }
}

//! Return elastic tensor
//! \retval de_ Elastic tensor
mpm::Material::Matrix6x6 mpm::LinearElastic::elastic_tensor() {
  // Bulk and shear modulus
  const double K = youngs_modulus_ / (3.0 * (1. - 2. * poisson_ratio_));
  const double G = youngs_modulus_ / (2.0 * (1. + poisson_ratio_));

  const double a1 = K + (4.0 / 3.0) * G;
  const double a2 = K - (2.0 / 3.0) * G;

  // clang-format off
  // compute elasticityTensor
  de_(0,0)=a1;    de_(0,1)=a2;    de_(0,2)=a2;    de_(0,3)=0;    de_(0,4)=0;    de_(0,5)=0;
  de_(1,0)=a2;    de_(1,1)=a1;    de_(1,2)=a2;    de_(1,3)=0;    de_(1,4)=0;    de_(1,5)=0;
  de_(2,0)=a2;    de_(2,1)=a2;    de_(2,2)=a1;    de_(2,3)=0;    de_(2,4)=0;    de_(2,5)=0;
  de_(3,0)= 0;    de_(3,1)= 0;    de_(3,2)= 0;    de_(3,3)=G;    de_(3,4)=0;    de_(3,5)=0;
  de_(4,0)= 0;    de_(4,1)= 0;    de_(4,2)= 0;    de_(4,3)=0;    de_(4,4)=G;    de_(4,5)=0;
  de_(5,0)= 0;    de_(5,1)= 0;    de_(5,2)= 0;    de_(5,3)=0;    de_(5,4)=0;    de_(5,5)=G;
  // clang-format on

  return de_;
}

//! Compute stress
//! \param[in] strain Strain
//! \param[in] stress Stress
//! \retval updated_stress Updated value of stress
mpm::Material::Vector6d mpm::LinearElastic::compute_stress(
    const Vector6d& strain, const Vector6d& stress) {
  Vector6d dstress = de_ * strain;
  return stress + dstress;
}

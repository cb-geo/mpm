//! Constructor with id and material properties
template <unsigned Tdim>
mpm::NorSand<Tdim>::NorSand(unsigned id,
                                    const Json& material_properties)
    : Material<Tdim>(id, material_properties) {
  try {
    // General parameters
    // Density
    density_ = material_properties["density"].template get<double>();
    // Youngs modulus E
    youngs_modulus_ =
        material_properties.at("youngs_modulus").template get<double>();
    // Shear modulus constant An
    shear_modulus_constant_ =
        material_properties["shear_modulus_constant"].template get<double>();
    // Shear modulus exponent Gn
    shear_modulus_exponent_ =
        material_properties["shear_modulus_exponent"].template get<double>();
    // Reference pressure pref
    reference_pressure_ =
        material_properties["reference_pressure"].template get<double>();
    // Poisson ratio
    poisson_ratio_ =
        material_properties["poisson_ratio"].template get<double>();
    // Critical state coefficient M
    M_ = material_properties["M"].template get<double>();
    // Volumetric coupling (dilatancy) parameter N
    N_ = material_properties["N"].template get<double>();
    // Minimum void ratio
    e_min_ =
        material_properties["e_min"].template get<double>();
    // Maximum void ratio
    e_max_ =
        material_properties["e_max"].template get<double>();
    // Crushing pressure
    crushing_pressure_ = material_properties["crushing_pressure"].template get<double>();
    // Dilatancy coefficient chi
    chi_ = material_properties["chi"].template get<double>();
    // Hardening modulus
    hardening_modulus_ = material_properties["hardening_modulus"].template get<double>();
    // Initial void ratio
    void_ratio_initial_ = material_properties["void_ratio_initial"].template get<double>();
    
    // Properties
    properties_ = material_properties;

  } catch (std::exception& except) {
    console_->error("Material parameter not set: {}\n", except.what());
  }
}

//! Initialise state variables
template <unsigned Tdim>
mpm::dense_map mpm::NorSand<Tdim>::initialise_state_variables() {
  mpm::dense_map state_vars = { // Current void ratio
                                {"void_ratio", void_ratio_initial_}
                              };

  return state_vars;
}

//! Compute elastic tensor
template <unsigned Tdim>
bool mpm::NorSand<Tdim>::compute_elastic_tensor() {
  // Shear modulus
  const double G = shear_modulus_;
  const double a1 = bulk_modulus_ + (4.0 / 3.0) * G;
  const double a2 = bulk_modulus_ - (2.0 / 3.0) * G;
  // compute elastic stiffness matrix
  // clang-format off
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
Eigen::Matrix<double, 6, 1> mpm::NorSand<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr, mpm::dense_map* state_vars) {

  // Compression positive
  Vector6d stress_neg = -1 * stress;

  // Compute the mean pressure (must not be zero, compression positive)
  double mean_p = abs(stress_neg(0) + stress_neg(1) + stress_neg(2)) / 3.;
  // Compute deviatoric q (must not be zero)
  double deviatoric_q = sqrt(0.5 * (pow((stress_neg(0) - stress_neg(1)), 2) + pow((stress_neg(1) - stress_neg(2)), 2) + pow((stress_neg(2) - stress_neg(0)), 2) +
                             6 * (pow(stress_neg(3), 2) + pow(stress_neg(4), 2) + pow(stress_neg(5), 2))));

  // Compute pressure image
  double p_image = mean_p * pow(1/(1 - N_) - ((N_ - 1)/N_) * deviatoric_q / M_ / mean_p, ((N_ - 1)/N_));
  // Compute void ratio image
  double e_image = e_max_ - (e_max_ - e_min_) / log(crushing_pressure_ / p_image);

  // Shear modulus
  // shear_modulus_ = shear_modulus_constant_ * pow(mean_p / reference_pressure_, shear_modulus_exponent_);
  shear_modulus_ = youngs_modulus_ / (2.0 * (1. + poisson_ratio_));
  // Bulk modulus
  bulk_modulus_ = shear_modulus_ * (2.0 * (1 + poisson_ratio_)) / (3.0 * (1. - 2. * poisson_ratio_));
  // Set elastic tensor
  this->compute_elastic_tensor();

  // Get void ratio and update state variables
  double dvolumetric_strain = dstrain(0) + dstrain(1) + dstrain(2); 
  double void_ratio = (*state_vars).at("void_ratio") - (1 + (*state_vars).at("void_ratio")) * dvolumetric_strain;  
  (*state_vars).at("void_ratio") = void_ratio;

  // Compute psi
  double psi_image = void_ratio - e_image;

  // Estimate dilatancy at peak
  double D_min = chi_ * psi_image;
  // Estimate maximum image pressure
  double p_image_max = mean_p * pow((1 + D_min * N_ / M_), ((N_ - 1) / N_));

  // Compute derivatives
  double dF_dp = -1. * M_ / N_ * (1 + (N_ - 1) / (1 - N_)) * pow((mean_p / p_image), (N_ / (1 - N_)));
  Vector6d dp_dsigma = Vector6d::Zero();
  dp_dsigma(0) = 1./3.;
  dp_dsigma(1) = 1./3.;
  dp_dsigma(2) = 1./3.;
  double dF_dq = 1.;
  Vector6d dq_dsigma = Vector6d::Zero();
  dq_dsigma(0) = 3./2./deviatoric_q * (stress_neg(0) - mean_p);
  dq_dsigma(1) = 3./2./deviatoric_q * (stress_neg(1) - mean_p);
  dq_dsigma(2) = 3./2./deviatoric_q * (stress_neg(2) - mean_p);  
  dq_dsigma(3) = 3./2./deviatoric_q * stress_neg(3);
  dq_dsigma(4) = 3./2./deviatoric_q * stress_neg(4);
  dq_dsigma(5) = 3./2./deviatoric_q * stress_neg(5);

  Vector6d dF_dsigma = dF_dp * dp_dsigma + dF_dq * dq_dsigma;

  double dF_dpi = M_ / N_ * (N_ - 1)  / (1 - N_) * pow((mean_p / p_image), (1 / (1 - N_)));
  double dpi_depsd = hardening_modulus_ * (p_image_max - p_image);

  double dF_dsigma_ii = dF_dsigma(0) + dF_dsigma(1) + dF_dsigma(2);
  double dF_dsigma_deviatoric = sqrt(2/3) * sqrt( pow(dF_dsigma(0) - dF_dsigma_ii/3, 2) + 
                                                  pow(dF_dsigma(1) - dF_dsigma_ii/3, 2) + 
                                                  pow(dF_dsigma(2) - dF_dsigma_ii/3, 2) + 
                                                  2 * pow(dF_dsigma(3), 2) + 
                                                  2 * pow(dF_dsigma(4), 2) + 
                                                  2 * pow(dF_dsigma(5), 2) );
  
  // Construct Dp matrix
  // Matrix6x6 Dp = de_ * (dF_dsigma.transpose() * de_ * dF_dsigma) / 
  //                (dF_dsigma.transpose() * de_ * dF_dsigma - dF_dpi * dpi_depsd * dF_dsigma_deviatoric);
  Matrix6x6 Dp = de_ * (de_ * dF_dsigma * dF_dsigma.transpose()) / 
                 (dF_dsigma.transpose() * de_ * dF_dsigma - dF_dpi * dpi_depsd * dF_dsigma_deviatoric);

  // Compute D matrix used in stress update
  Matrix6x6 D_matrix = de_ - Dp;

  // Update stress
  Vector6d updated_stress = stress + D_matrix * dstrain;

  // std::cout << "Current stress: " << stress << '\n';
  // std::cout << "Updated stress: " << updated_stress << '\n';

  return updated_stress;
}
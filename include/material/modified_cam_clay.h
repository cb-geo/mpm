#ifndef MPM_MATERIAL_MODIFIED_CAM_CLAY_H_
#define MPM_MATERIAL_MODIFIED_CAM_CLAY_H_

#include <limits>

#include <cmath>

#include "Eigen/Dense"

#include "material.h"

namespace mpm {

//! ModifiedCamClay class
//! \brief Modified Cam Clay material model
//! \tparam Tdim Dimension
template <unsigned Tdim>
class ModifiedCamClay : public Material<Tdim> {
 public:
  //! Define a vector of 6 dof
  using Vector6d = Eigen::Matrix<double, 6, 1>;
  //! Define a Matrix of 6 x 6
  using Matrix6x6 = Eigen::Matrix<double, 6, 6>;

  //! Failure state
  enum FailureState { Elastic = 0, Yield = 1 };

  //! Constructor with id and material properties
  //! \param[in] material_properties Material properties
  ModifiedCamClay(unsigned id, const Json& material_properties);

  //! Destructor
  ~ModifiedCamClay() override{};

  //! Delete copy constructor
  ModifiedCamClay(const ModifiedCamClay&) = delete;

  //! Delete assignement operator
  ModifiedCamClay& operator=(const ModifiedCamClay&) = delete;

  //! Initialise history variables
  //! \retval state_vars State variables with history
  mpm::dense_map initialise_state_variables() override;

  //! Compute stress
  //! \param[in] stress Stress
  //! \param[in] dstrain Strain
  //! \param[in] particle Constant point to particle base
  //! \param[in] state_vars History-dependent state variables
  //! \retval updated_stress Updated value of stress
  Vector6d compute_stress(const Vector6d& stress, const Vector6d& dstrain,
                          const ParticleBase<Tdim>* ptr,
                          mpm::dense_map* state_vars) override;

  //! Compute stress invariants (j3, q, theta, and epsilon)
  //! \param[in] stress Stress
  //! \param[in] state_vars History-dependent state variables
  //! \retval status of computation of stress invariants
  bool compute_stress_invariants(const Vector6d& stress,
                                 mpm::dense_map* state_vars);

  //! Compute deviatoric stress tensor
  //! \param[in] stress Stress
  //! \param[in] state_vars History-dependent state variables
  //! \retval Deviatoric stress tensor
  Eigen::Matrix<double, 6, 1> compute_deviatoric_stress_tensor(
      const Vector6d& stress, mpm::dense_map* state_vars);

  //! Compute yield function and yield state
  //! \param[in] state_vars History-dependent state variables
  //! \retval yield_type Yield type (elastic or yield)
  FailureState compute_yield_state(mpm::dense_map* state_vars);

  //! Compute bonding parameters
  //! \param[in] chi Degredation
  //! \param[in] state_vars History-dependent state variables
  void compute_bonding_parameters(const double chi, mpm::dense_map* state_vars);

  //! Compute subloading parameters
  //! \param[in] subloading_r Subloading ratio
  //! \param[in] state_vars History-dependent state variables
  void compute_subloading_parameters(const double subloading_r,
                                     mpm::dense_map* state_vars);

  //! Compute dF/dmul
  //! \param[in] state_vars History-dependent state variables
  //! \param[in] df_dmul dF / ddelta_phi
  void compute_df_dmul(const mpm::dense_map* state_vars, double* df_dmul);

  //! Compute dF/dSigma
  //! \param[in] state_vars History-dependent state variables
  //! \param[in] stress Stress
  //! \param[in] df_dsigma dF/dSigma
  void compute_df_dsigma(const mpm::dense_map* state_vars,
                         const Vector6d& stress, Vector6d* df_dsigma);

  //! Compute G and dG/dpc
  //! \param[in] state_vars History-dependent state variables
  //! \param[in] pc_n Preconsolidation pressure of last step
  //! \param[in] p_trial Volumetric trial stress
  //! \param[in] g_function G
  //! \param[in] dg_dpc dG / dpc
  void compute_dg_dpc(const mpm::dense_map* state_vars, const double pc_n,
                      const double p_trial, double* g_function, double* dg_dpc);

  //! Compute elastic tensor
  //! \param[in] state_vars History-dependent state variables
  //! \retval status of computation
  bool compute_elastic_tensor(mpm::dense_map* state_vars);

  //! Compute plastic tensor
  //! \param[in] stress Stress
  //! \param[in] state_vars History-dependent state variables
  //! \retval status of computation
  bool compute_plastic_tensor(const Vector6d& stress,
                              mpm::dense_map* state_vars);

 protected:
  //! material id
  using Material<Tdim>::id_;
  //! Material properties
  using Material<Tdim>::properties_;
  //! Logger
  using Material<Tdim>::console_;

 private:
  //! Elastic stiffness matrix
  Matrix6x6 de_;
  //! Plastic stiffness matrix
  Matrix6x6 dp_;
  //! General parameters
  //! Density
  double density_{std::numeric_limits<double>::max()};
  //! Youngs modulus
  double youngs_modulus_{std::numeric_limits<double>::max()};
  //! Poisson ratio
  double poisson_ratio_{std::numeric_limits<double>::max()};
  //! Modified Cam Clay parameters
  //! Reference mean pressure
  double p_ref_{std::numeric_limits<double>::max()};
  //! Reference void ratio
  double e_ref_{std::numeric_limits<double>::max()};
  //! Initial void ratio
  double e0_{std::numeric_limits<double>::max()};
  //! Initial preconsolidation pressure
  double pc0_{std::numeric_limits<double>::max()};
  //! OCR
  double ocr_{std::numeric_limits<double>::max()};
  //! M (or Mtc in "Three invariants type")
  double m_{std::numeric_limits<double>::max()};
  //! Lambda
  double lambda_{std::numeric_limits<double>::max()};
  //! Kappa
  double kappa_{std::numeric_limits<double>::max()};
  //! Three invariants
  bool three_invariants_{false};
  //! Subloading surface properties
  //! Subloading status
  bool subloading_{false};
  //! Material constant controling plastic deformation
  double subloading_u_{std::numeric_limits<double>::epsilon()};
  //! Bonding properties
  //! Bonding status
  bool bonding_{false};
  //! Material constants a
  double mc_a_{std::numeric_limits<double>::epsilon()};
  //! Material constants b
  double mc_b_{std::numeric_limits<double>::epsilon()};
  //! Material constants c
  double mc_c_{std::numeric_limits<double>::epsilon()};
  //! Material constants d
  double mc_d_{std::numeric_limits<double>::epsilon()};
  //! Degradation
  double m_degradation_{std::numeric_limits<double>::epsilon()};
  // Increment in shear modulus
  double m_shear_ = {std::numeric_limits<double>::epsilon()};
  //! Hydrate saturation
  double s_h_{std::numeric_limits<double>::epsilon()};
};  // ModifiedCamClay class
}  // namespace mpm

#include "modified_cam_clay.tcc"

#endif  // MPM_MATERIAL_MODIFIED_CAM_CLAY_H_

#ifndef MPM_MATERIAL_NORSAND_H_
#define MPM_MATERIAL_NORSAND_H_

#include <cmath>

#include "Eigen/Dense"

#include "material.h"

namespace mpm {

//! NorSand class
//! \brief Mohr Coulomb material model
//! \details Mohr Coulomb material model with softening
//! \tparam Tdim Dimension
template <unsigned Tdim>
class NorSand : public Material<Tdim> {
 public:
  //! Define a vector of 6 dof
  using Vector6d = Eigen::Matrix<double, 6, 1>;
  //! Define a Matrix of 6 x 6
  using Matrix6x6 = Eigen::Matrix<double, 6, 6>;

  //! Failure state
  enum FailureState { Elastic = 0, Yield = 1 };

  //! Constructor with id and material properties
  //! \param[in] material_properties Material properties
  NorSand(unsigned id, const Json& material_properties);

  //! Destructor
  ~NorSand() override{};

  //! Delete copy constructor
  NorSand(const NorSand&) = delete;

  //! Delete assignement operator
  NorSand& operator=(const NorSand&) = delete;

  //! Initialise history variables
  //! \retval state_vars State variables with history
  mpm::dense_map initialise_state_variables() override;

  //! Thermodynamic pressure
  //! \param[in] volumetric_strain dVolumetric_strain
  //! \retval pressure Pressure for volumetric strain
  double thermodynamic_pressure(double volumetric_strain) const override {
    return 0;
  };

  //! Compute stress
  //! \param[in] stress Stress
  //! \param[in] dstrain Strain
  //! \param[in] particle Constant point to particle base
  //! \param[in] state_vars History-dependent state variables
  //! \retval updated_stress Updated value of stress
  Vector6d compute_stress(const Vector6d& stress, const Vector6d& dstrain,
                          const ParticleBase<Tdim>* ptr,
                          mpm::dense_map* state_vars) override;

  //! Compute stress invariants (p, q, etc)
  //! \param[in] stress Stress
  //! \param[in] state_vars History-dependent state variables
  //! \retval status of computation of stress invariants
  bool compute_stress_invariants(const Vector6d& stress,
                                 mpm::dense_map* state_vars);

  //! Compute p_cohesion and p_dilation
  //! \param[in] state_vars History-dependent state variables
  //! \retval status of computation of stress invariants
  bool compute_p_bond(mpm::dense_map* state_vars);

  //! Compute state variables (void ratio, p_image, e_image, etc)
  //! \param[in] stress Stress
  //! \param[in] state_vars History-dependent state variables
  //! \param[in] yield_type Yild type (elastic or yield)
  //! \retval status of computation of stress invariants
  bool compute_state_variables(const Vector6d& stress, const Vector6d& dstrain,
                               mpm::dense_map* state_vars,
                               FailureState yield_type);

  //! Compute yield function and yield state
  //! \param[in] state_vars History-dependent state variables
  //! \param[in] stress Stress
  //! \retval yield_type Yield type (elastic or yield)
  FailureState compute_yield_state(double* yield_function,
                                   const Vector6d& stress,
                                   const mpm::dense_map* state_vars);

 protected:
  //! material id
  using Material<Tdim>::id_;
  //! Material properties
  using Material<Tdim>::properties_;
  //! Logger
  using Material<Tdim>::console_;

 private:
  //! Compute elastic tensor
  bool compute_elastic_tensor();

  //! Compute plastic tensor
  bool compute_plastic_tensor(const Vector6d& stress,
                              mpm::dense_map* state_vars);

  //! Elastic stiffness matrix
  Matrix6x6 de_;
  //! Plastic stiffness matrix
  Matrix6x6 dp_;
  //! Density
  double density_{std::numeric_limits<double>::max()};
  //! Youngs modulus
  double youngs_modulus_{std::numeric_limits<double>::max()};
  //! Shear modulusc constant A
  double shear_modulus_constant_{std::numeric_limits<double>::max()};
  //! Shear modulus exponent Gn
  double shear_modulus_exponent_{std::numeric_limits<double>::max()};
  //! Bulk modulus
  double bulk_modulus_{std::numeric_limits<double>::max()};
  //! Shear modulus
  double shear_modulus_{std::numeric_limits<double>::max()};
  //! Reference pressure pref
  double reference_pressure_{std::numeric_limits<double>::max()};
  //! Poisson ratio
  double poisson_ratio_{std::numeric_limits<double>::max()};
  //! Critical state friction angle
  double friction_cs_{std::numeric_limits<double>::max()};
  //! Critical state coefficient M in triaxial compression
  double Mtc_{std::numeric_limits<double>::max()};
  //! Critical state coefficient M in triaxial extension
  double Mte_{std::numeric_limits<double>::max()};
  //! Volumetric coupling (dilatancy) parameter N
  double N_{std::numeric_limits<double>::max()};
  //! Minimum void ratio
  double e_min_{std::numeric_limits<double>::max()};
  //! Maximum void ratio
  double e_max_{std::numeric_limits<double>::max()};
  //! Crushing pressure
  double crushing_pressure_{std::numeric_limits<double>::max()};
  //! Dilatancy coefficient
  double chi_{std::numeric_limits<double>::max()};
  //! Hardening modulus
  double hardening_modulus_{std::numeric_limits<double>::max()};
  //! Initial void ratio
  double void_ratio_initial_{std::numeric_limits<double>::max()};
  //! Initial image pressure
  double p_image_initial_{std::numeric_limits<double>::max()};
  //! Flag for bonded model
  bool bond_model_{false};
  //! Initial p_cohesion
  double p_cohesion_initial_{std::numeric_limits<double>::max()};
  //! Initial p_dilation
  double p_dilation_initial_{std::numeric_limits<double>::max()};
  //! Cohesion degradation parameter m upon shearing
  double m_cohesion_{std::numeric_limits<double>::max()};
  //! Dilation degradation parameter m upon shearing
  double m_dilation_{std::numeric_limits<double>::max()};

};  // NorSand class
}  // namespace mpm

#include "norsand.tcc"

#endif  // MPM_MATERIAL_NORSAND_H_

#ifndef MPM_MATERIAL_MOHR_COULOMB_H_
#define MPM_MATERIAL_MOHR_COULOMB_H_

#include <limits>

#include "Eigen/Dense"

#include "material.h"

namespace mpm {

//! MohrCoulomb class
//! \brief Linear Elastic material model
//! \details MohrCoulomb class stresses and strains
//! \tparam Tdim Dimension
template <unsigned Tdim>
class MohrCoulomb : public Material<Tdim> {
 public:
  //! Define a vector of 6 dof
  using Vector6d = Eigen::Matrix<double, 6, 1>;
  //! Define a Matrix of 6 x 6
  using Matrix6x6 = Eigen::Matrix<double, 6, 6>;

  //! Constructor with id and material properties
  //! \param[in] material_properties Material properties
  MohrCoulomb(unsigned id, const Json& material_properties);

  //! Destructor
  ~MohrCoulomb() override{};

  //! Delete copy constructor
  MohrCoulomb(const MohrCoulomb&) = delete;

  //! Delete assignement operator
  MohrCoulomb& operator=(const MohrCoulomb&) = delete;

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
  //! \retval updated_stress Updated value of stress
  Vector6d compute_stress(const Vector6d& stress, const Vector6d& dstrain,
                          const ParticleBase<Tdim>* ptr) override;

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
  //! Compute j2, j3, rho, theta
  bool compute_rho_theta(const Vector6d stress, double& _j2, double& _j3,
			 double& _rho, double& _theta);
  //! Check the yield type(tension/shear) and return the value of two yield functions
  int check_yield(Eigen::Matrix<double, 2, 1>& _yield_function,
       const double _epsilon, const double _rho, const double _theta);
  //! Compute dF/dSigma and dP/dSigma
  bool compute_df_dp(const int _yield_type, const double _j2, const double _j3,
         const double _rho, const double _theta, const Vector6d stress,
         Vector6d& _dF_dSigma, Vector6d& _dP_dSigma,
         const double _epds, double& _softening,
         const ParticleBase<Tdim>* ptr);

 private:
  //! Elastic stiffness matrix
  Matrix6x6 de_;
  //! Density
  double density_{std::numeric_limits<double>::max()};
  //! Youngs modulus
  double youngs_modulus_{std::numeric_limits<double>::max()};
  //! Bulk modulus
  double bulk_modulus_{std::numeric_limits<double>::max()};
  //! Shear modulus
  double shear_modulus_{std::numeric_limits<double>::max()};
  //! Poisson ratio
  double poisson_ratio_{std::numeric_limits<double>::max()};
  //! Friction angle phi
  double friction_angle_{std::numeric_limits<double>::max()};
  //! Dilation angle psi
  double dilation_angle_{std::numeric_limits<double>::max()};
  //! Cohesion
  double cohesion_{std::numeric_limits<double>::max()};
  //! Residual friction angle phi
  double residual_friction_angle_{std::numeric_limits<double>::max()};
  //! Residual dilation angle psi
  double residual_dilation_angle_{std::numeric_limits<double>::max()};
  //! Residual cohesion
  double residual_cohesion_{std::numeric_limits<double>::max()};
  //! Peak plastic deviatoric strain
  double peak_EPDS_{std::numeric_limits<double>::max()};
  //! Critical plastic deviatoric strain
  double crit_EPDS_{std::numeric_limits<double>::max()};
  //! Tension cutoff
  double tension_cutoff_{std::numeric_limits<double>::max()};
  //! Porosity
  double porosity_{std::numeric_limits<double>::max()};
  //! Permeability
  double permeability_{std::numeric_limits<double>::max()};
  //! Friction angle phi
  double phi_{std::numeric_limits<double>::max()};
  //! Dilation angle psi
  double psi_{std::numeric_limits<double>::max()};
  //! Cohesion
  double c_{std::numeric_limits<double>::max()};
private:
  //! tracke the equivalent_plastic_deviatoric_strain
  double epds_last = 0.;
  //! value of PI
  const double PI = std::atan(1.0) * 4.;

};  // MohrCoulomb class
}  // namespace mpm

#include "MohrCoulomb.tcc"

#endif  // MPM_MATERIAL_MOHR_COULOMB_H_

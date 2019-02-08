#ifndef MPM_MATERIAL_LINEAR_ELASTIC_H_
#define MPM_MATERIAL_LINEAR_ELASTIC_H_

#include <limits>

#include "Eigen/Dense"

#include "material.h"

namespace mpm {

//! LinearElastic class
//! \brief Linear Elastic material model
//! \details LinearElastic class stresses and strains
//! \tparam Tdim Dimension
template <unsigned Tdim>
class LinearElastic : public Material<Tdim> {
 public:
  //! Define a vector of 6 dof
  using Vector6d = Eigen::Matrix<double, 6, 1>;
  //! Define a Matrix of 6 x 6
  using Matrix6x6 = Eigen::Matrix<double, 6, 6>;

  //! Constructor with id
  //! \param[in] material_properties Material properties
  LinearElastic(unsigned id, const Json& material_properties);

  //! Destructor
  ~LinearElastic() override{};

  //! Delete copy constructor
  LinearElastic(const LinearElastic&) = delete;

  //! Delete assignement operator
  LinearElastic& operator=(const LinearElastic&) = delete;

  //! Initialise history variables
  //! \param[in] state_vars State variables with history
  bool initialise_state_variables(std::map<std::string, double>* state_vars) {
    return true;
  }

  //! Thermodynamic pressure
  //! \param[in] volumetric_strain dVolumetric_strain
  //! \retval pressure Pressure for volumetric strain
  double thermodynamic_pressure(double volumetric_strain) const override;

  //! Compute stress
  //! \param[in] stress Stress
  //! \param[in] dstrain Strain
  //! \param[in] particle Constant point to particle base
  //! \param[in] state_vars History-dependent state variables
  //! \retval updated_stress Updated value of stress
  Vector6d compute_stress(const Vector6d& stress, const Vector6d& dstrain,
                          const ParticleBase<Tdim>* ptr,
                          std::map<std::string, double>* state_vars) override;

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

 private:
  //! Elastic stiffness matrix
  Matrix6x6 de_;
  //! Density
  double density_{std::numeric_limits<double>::max()};
  //! Youngs modulus
  double youngs_modulus_{std::numeric_limits<double>::max()};
  //! Poisson ratio
  double poisson_ratio_{std::numeric_limits<double>::max()};
  //! Bulk modulus
  double bulk_modulus_{std::numeric_limits<double>::max()};
};  // LinearElastic class
}  // namespace mpm

#include "linear_elastic.tcc"

#endif  // MPM_MATERIAL_LINEAR_ELASTIC_H_

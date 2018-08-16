#ifndef MPM_MATERIAL_BINGHAM_H_
#define MPM_MATERIAL_BINGHAM_H_

#include <iostream>
#include <limits>

#include "Eigen/Dense"

#include "material.h"

namespace mpm {

//! Bingham class
//! \brief Bingham fluid material model
//! \details Bingham class stresses and strains
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Bingham : public Material<Tdim> {
 public:
  //! Define a vector of 6 dof
  using Vector6d = Eigen::Matrix<double, 6, 1>;
  //! Define a Matrix of 6 x 6
  using Matrix6x6 = Eigen::Matrix<double, 6, 6>;

  //! Constructor with id
  Bingham(unsigned id) : Material<Tdim>(id){};

  //! Destructor
  virtual ~Bingham() override{};

  //! Delete copy constructor
  Bingham(const Bingham&) = delete;

  //! Delete assignement operator
  Bingham& operator=(const Bingham&) = delete;

  //! Read material properties
  //! \param[in] material_properties Material properties
  void properties(const Json& material_properties) override;

  //! Compute elastic tensor
  //! \retval de_ Elastic tensor
  Matrix6x6 elastic_tensor() override;

  //! Compute stress
  //! \param[in] stress Stress
  //! \param[in] dstrain Strain
  //! \retval updated_stress Updated value of stress
  Vector6d compute_stress(const Vector6d& stress,
                          const Vector6d& dstrain) override;
  //! Compute stress
  //! \param[in] stress Stress
  //! \param[in] dstrain Strain
  //! \param[in] particle Constant point to particle base
  //! \retval updated_stress Updated value of stress
  Vector6d compute_stress(const Vector6d& stress, const Vector6d& dstrain,
                          const ParticleBase<Tdim>* ptr) override;

  //! Check if this material needs a particle handle
  bool property_handle() const override { return true; }

 protected:
  //! material id
  using Material<Tdim>::id_;
  //! material status
  using Material<Tdim>::status_;
  //! Material properties
  using Material<Tdim>::properties_;
  //! Logger
  using Material<Tdim>::console_;

 private:
  //! Dirac delta function in Voigt notation
  Eigen::Matrix<double, 6, 1> dirac_delta() const;

  //! Density
  double density_{std::numeric_limits<double>::max()};
  //! Youngs modulus
  double youngs_modulus_{std::numeric_limits<double>::max()};
  //! Poisson ratio
  double poisson_ratio_{std::numeric_limits<double>::max()};
  //! Tau0 - shear yield stress in unit of [Pa]
  double tau0_{std::numeric_limits<double>::max()};
  //! mu - constant plastic viscosity [N s / m^2 or kg / m / s]
  double mu_{std::numeric_limits<double>::max()};
  //! Critical yielding shear rate
  double critical_shear_rate_{std::numeric_limits<double>::max()};

};  // Bingham class
}  // namespace mpm

#include "bingham.tcc"

#endif  // MPM_MATERIAL_BINGHAM_H_

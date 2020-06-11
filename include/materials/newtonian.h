#ifndef MPM_MATERIAL_NEWTONIAN_H_
#define MPM_MATERIAL_NEWTONIAN_H_

#include <limits>

#include "Eigen/Dense"

#include "material.h"

namespace mpm {

//! Newtonian class
//! \brief Newtonian fluid material model
//! \details Newtonian class stresses and strains
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Newtonian : public Material<Tdim> {
 public:
  //! Define a vector of 6 dof
  using Vector6d = Eigen::Matrix<double, 6, 1>;
  //! Define a Matrix of 6 x 6
  using Matrix6x6 = Eigen::Matrix<double, 6, 6>;

  //! Constructor with id and material properties
  //! \param[in] id Material ID
  //! \param[in] material_properties Material properties
  Newtonian(unsigned id, const Json& material_properties);

  //! Destructor
  ~Newtonian() override{};

  //! Delete copy constructor
  Newtonian(const Newtonian&) = delete;

  //! Delete assignement operator
  Newtonian& operator=(const Newtonian&) = delete;

  //! Initialise history variables
  //! \retval state_vars State variables with history
  mpm::dense_map initialise_state_variables() override;

  //! State variables
  std::vector<std::string> state_variables() const override;

  //! Compute stress
  //! \param[in] stress Stress
  //! \param[in] dstrain Strain
  //! \param[in] particle Constant point to particle base
  //! \param[in] state_vars History-dependent state variables
  //! \retval updated_stress Updated value of stress
  Vector6d compute_stress(const Vector6d& stress, const Vector6d& dstrain,
                          const Particle<Tdim>* ptr,
                          mpm::dense_map* state_vars) override;

 protected:
  //! material id
  using Material<Tdim>::id_;
  //! Material properties
  using Material<Tdim>::properties_;
  //! Logger
  using Material<Tdim>::console_;

 private:
  //! Thermodynamic pressure
  //! \param[in] volumetric_strain dVolumetric_strain
  //! \retval pressure Pressure for volumetric strain
  double thermodynamic_pressure(double volumetric_strain) const;

  //! Density
  double density_{std::numeric_limits<double>::max()};
  //! Fluid Bulk modulus
  double bulk_modulus_{std::numeric_limits<double>::max()};
  //! Fluid Dynamic viscosity
  double dynamic_viscosity_{std::numeric_limits<double>::max()};
  //! Compressibility multiplier
  double compressibility_multiplier_{1.0};

};  // Newtonian class
}  // namespace mpm

#include "newtonian.tcc"

#endif  // MPM_MATERIAL_NEWTONIAN_H_

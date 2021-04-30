#ifndef MPM_MATERIAL_RIGID_H_
#define MPM_MATERIAL_RIGID_H_

#include <limits>

#include "Eigen/Dense"

#include "material.h"

namespace mpm {

//! Rigid class
//! \brief Rigid material model
//! \details Rigid class stresses and strains
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Rigid : public Material<Tdim> {
 public:
  //! Define a vector of 6 dof
  using Vector6d = Eigen::Matrix<double, 6, 1>;
  //! Define a Matrix of 6 x 6
  using Matrix6x6 = Eigen::Matrix<double, 6, 6>;

  //! Constructor with id
  //! \param[in] material_properties Material properties
  Rigid(unsigned id, const Json& material_properties);

  //! Destructor
  ~Rigid() override{};

  //! Delete copy constructor
  Rigid(const Rigid&) = delete;

  //! Delete assignement operator
  Rigid& operator=(const Rigid&) = delete;

  //! Initialise history variables
  //! \retval state_vars State variables with history
  mpm::dense_map initialise_state_variables() override {
    mpm::dense_map state_vars;
    return state_vars;
  }

  //! State variables
  std::vector<std::string> state_variables() const override;

  //! Compute stress
  //! \param[in] stress Stress
  //! \param[in] dstrain Strain
  //! \param[in] particle Constant point to particle base
  //! \param[in] state_vars History-dependent state variables
  //! \retval updated_stress Updated value of stress
  Vector6d compute_stress(const Vector6d& stress, const Vector6d& dstrain,
                          const ParticleBase<Tdim>* ptr,
                          mpm::dense_map* state_vars) override;

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
  //! Density
  double density_{std::numeric_limits<double>::max()};
};  // Rigid class
}  // namespace mpm

#include "rigid.tcc"

#endif  // MPM_MATERIAL_RIGID_H_

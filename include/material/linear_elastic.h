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
  LinearElastic(unsigned id) : Material<Tdim>(id) {
    //! Logger
    std::string logger = "material::" + std::to_string(id);
    console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
  };

  //! Destructor
  ~LinearElastic() override{};

  //! Delete copy constructor
  LinearElastic(const LinearElastic&) = delete;

  //! Delete assignement operator
  LinearElastic& operator=(const LinearElastic&) = delete;

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
  bool property_handle() const override { return false; }


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
  //! Elastic stiffness matrix
  Matrix6x6 de_;
  //! Density
  double density_{std::numeric_limits<double>::max()};
  //! Youngs modulus
  double youngs_modulus_{std::numeric_limits<double>::max()};
  //! Poisson ratio
  double poisson_ratio_{std::numeric_limits<double>::max()};
};  // LinearElastic class
}  // namespace mpm

#include "linear_elastic.tcc"

#endif  // MPM_MATERIAL_LINEAR_ELASTIC_H_

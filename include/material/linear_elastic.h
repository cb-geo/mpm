#ifndef MPM_MATERIAL_LINEAR_ELASTIC_H_
#define MPM_MATERIAL_LINEAR_ELASTIC_H_

#include <iostream>
#include <limits>

#include "Eigen/Dense"

#include "material.h"

namespace mpm {

//! LinearElastic class
//! \brief Linear Elastic material model
//! \details LinearElastic class stresses and strains
class LinearElastic : public Material {
 public:
  //! Define a vector of 6 dof
  using Vector6d = Eigen::Matrix<double, 6, 1>;
  //! Define a Matrix of 6 x 6
  using Matrix6x6 = Eigen::Matrix<double, 6, 6>;

  //! Constructor with id
  LinearElastic(unsigned id) : Material(id){};

  //! Destructor
  virtual ~LinearElastic(){};

  //! Delete copy constructor
  LinearElastic(const LinearElastic&) = delete;

  //! Delete assignement operator
  LinearElastic& operator=(const LinearElastic&) = delete;

  //! Read material properties
  //! \param[in] materail_properties Material properties
  void properties(const Json& material_properties );

  //! Get material property
  //! \param[in] key Material properties key
  //! \retval result Value of material property
  double property(const std::string& key);

  //! Compute elastic tensor
  //! \retval de_ Elastic tensor
  Matrix6x6 elastic_tensor();

  //! Compute stress
  //! \param[in] strain Strain
  //! \param[in] stress Stress
  //! \retval updated_stress Updated value of stress
  void compute_stress(Vector6d& stress, const Vector6d& strain);

 protected:
  //! material id
  using Material::id_;
  //! material status
  using Material::status_;

 private:
  //! Elastic stiffness matrix
  Matrix6x6 de_;
  //! Material properties
  Json properties_;
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

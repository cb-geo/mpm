#ifndef MPM_MATERIAL_BINGHAM_H_
#define MPM_MATERIAL_BINGHAM_H_

#include <limits>

#include "Eigen/Dense"

#include "material.h"

namespace mpm {

//! Bingham class
//! \brief Bingham fluid material model
//! \details Bingham class stresses and strains
class Bingham : public Material {
 public:
  //! Define a vector of 6 dof
  using Vector6d = Eigen::Matrix<double, 6, 1>;
  //! Define a Matrix of 6 x 6
  using Matrix6x6 = Eigen::Matrix<double, 6, 6>;

  //! Constructor with id
  Bingham(unsigned id) : Material(id){};

  //! Destructor
  virtual ~Bingham(){};

  //! Delete copy constructor
  Bingham(const Bingham&) = delete;

  //! Delete assignement operator
  Bingham& operator=(const Bingham&) = delete;

  //! Read material properties
  //! \param[in] materail_properties Material properties
  void properties(const Json& material_properties );

  //! Compute stress
  //! \param[in] strain Strain
  //! \param[in] stress Stress
  //! \retval updated_stress Updated value of stress
  void compute_stress(Vector6d& stress, const Vector6d& dstrain, const Eigen::Matrix3d& strain_rate);

 protected:
  //! material id
  using Material::id_;
  //! material status
  using Material::status_;
  //! Material properties
  using Material::properties_;
 private:
  //! Elastic stiffness matrix
  Matrix6x6 de_;
  //! Density
  double density_{std::numeric_limits<double>::max()};
  //! Youngs modulus
  double youngs_modulus_{std::numeric_limits<double>::max()};
  //! Poisson ratio
  double poisson_ratio_{std::numeric_limits<double>::max()};
  //! Tau0 - shear threshold 
  double tau0_{std::numeric_limits<double>::max()};
  //! mu parameter 
  double mu_{std::numeric_limits<double>::max()};
  //! Strain cutoff 
  double strain_cutoff_{std::numeric_limits<double>::max()};

};  // Bingham class
}  // namespace mpm

#include "bingham.tcc"

#endif  // MPM_MATERIAL_BINGHAM_H_

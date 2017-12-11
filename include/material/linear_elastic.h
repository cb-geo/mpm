#ifndef MPM_LINEAR_ELASTIC_H_
#define MPM_LINEAR_ELASTIC_H_

#include <limits>

#include "Eigen/Dense"

#include "material.h"

namespace mpm {

// LinearElastic class
//! \brief Linear Elastic material model
//! \details LinearElastic class stresses and strains
class LinearElastic : public Material {
 public:
  //! Define a vector of 6 dof
  using Vector6d = Eigen::Matrix<double, 6, 1>;
  //! Define a Matrix of 6 x 6
  using Matrix6x6 = Eigen::Matrix<double, 6, 6>;

  // Constructor with id
  LinearElastic(unsigned id) : Material(id){};

  //! Destructor
  virtual ~LinearElastic(){};

  //! Delete copy constructor
  LinearElastic(const LinearElastic&) = delete;

  //! Delete assignement operator
  LinearElastic& operator=(const LinearElastic&) = delete;

  //! Return id of the linearelastic
  unsigned id() const { return id_; }

  //! Compute elastic tensor
  Matrix6x6 elastic_tensor() {
    Matrix6x6 De;
    return De;
  }

  //! Compute stress
  void compute_stress(const Vector6d& strain, const Vector6d& stress) {}

 protected:
  //! linearelastic id
  using Material::id_;
};  // LinearElastic class
}  // mpm namespace

#endif  // MPM_LINEAR_ELASTIC_H_

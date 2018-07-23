#ifndef MPM_MATERIAL_MATERIAL_H_
#define MPM_MATERIAL_MATERIAL_H_

#include <limits>

#include "Eigen/Dense"
#include "json.hpp"

#include "factory.h"

// JSON
using Json = nlohmann::json;

namespace mpm {

//! Material base class
//! \brief Base class that stores the information about materials
//! \details Material class stresses and strains
class Material {
 public:
  //! Define a vector of 6 dof
  using Vector6d = Eigen::Matrix<double, 6, 1>;
  //! Define a Matrix of 6 x 6
  using Matrix6x6 = Eigen::Matrix<double, 6, 6>;

  // Constructor with id
  //! \param[in] id Material id
  Material(unsigned id) : id_{id} {}

  //! Destructor
  virtual ~Material(){};

  //! Delete copy constructor
  Material(const Material&) = delete;

  //! Delete assignement operator
  Material& operator=(const Material&) = delete;

  //! Return id of the material
  unsigned id() const { return id_; }

  //! Return status as true, when properties are assigned
  bool status() const { return status_; }
  
  //! Read material properties
  virtual void properties(const Json&) = 0;

  //! Compute elastic tensor
  virtual Matrix6x6 elastic_tensor() = 0;

  //! Compute stress
  virtual void compute_stress(Vector6d& stress, const Vector6d& strain) = 0;

 protected:
  //! material id
  unsigned id_{std::numeric_limits<unsigned>::max()};
  //! status
  bool status_{false};
};  // Material class
}  // namespace mpm

#endif  // MPM_MATERIAL_MATERIAL_H_

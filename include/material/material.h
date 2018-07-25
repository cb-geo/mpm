#ifndef MPM_MATERIAL_MATERIAL_H_
#define MPM_MATERIAL_MATERIAL_H_

#include <iostream>
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
  //! \param[in] materail_properties Material properties
  virtual void properties(const Json&) = 0;

  //! Get material property
  //! \param[in] key Material properties key
  //! \retval result Value of material property
  double property(const std::string& key) {
    double result = std::numeric_limits<double>::max();
    try {
      result = properties_[key].template get<double>();
    } catch (std::exception& except) {
      std::cerr << "Material parameter not found: " << except.what() << '\n';
    }
    return result;
  }

  //! Compute elastic tensor
  //! \retval de_ Elastic tensor
  virtual Matrix6x6 elastic_tensor() = 0;

  //! Compute stress
  //! \param[in] strain Strain
  //! \param[in] stress Stress
  //! \retval updated_stress Updated value of stress
  virtual void compute_stress(Vector6d& stress, const Vector6d& strain) = 0;

  virtual void compute_stress(Vector6d& stress, const Vector6d& dstrain, const Eigen::Matrix3d& strain_rate) = 0;


 protected:
  //! material id
  unsigned id_{std::numeric_limits<unsigned>::max()};
  //! status
  bool status_{false};
  //! Material properties
  Json properties_;
};  // Material class
}  // namespace mpm

#endif  // MPM_MATERIAL_MATERIAL_H_

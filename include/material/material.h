#ifndef MPM_MATERIAL_MATERIAL_H_
#define MPM_MATERIAL_MATERIAL_H_

#include <iostream>
#include <limits>

#include "Eigen/Dense"
#include "json.hpp"

#include "factory.h"
#include "particle.h"
#include "particle_base.h"

// JSON
using Json = nlohmann::json;

namespace mpm {

// Forward declaration of ParticleBase
template <unsigned Tdim>
class ParticleBase;

//! Material base class
//! \brief Base class that stores the information about materials
//! \details Material class stresses and strains
//! \tparam Tdim Dimension
template <unsigned Tdim>
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
  double property(const std::string& key);

  //! Compute elastic tensor
  //! \retval de_ Elastic tensor
  virtual Matrix6x6 elastic_tensor() = 0;

  //! Compute stress
  //! \param[in] strain Strain
  //! \param[in] stress Stress
  //! \retval updated_stress Updated value of stress
  virtual Vector6d compute_stress(const Vector6d& stress,
                                  const Vector6d& strain) = 0;

  //! Compute stress
  //! \param[in] strain Strain
  //! \param[in] stress Stress
  //! \param[in] particle Constant point to particle base
  //! \retval updated_stress Updated value of stress
  virtual Vector6d compute_stress(const Vector6d& stress,
                                  const Vector6d& strain,
                                  const ParticleBase<Tdim>* ptr) = 0;

  // Test function
  virtual void testfn(const ParticleBase<Tdim>* ptr) = 0;

 protected:
  //! material id
  unsigned id_{std::numeric_limits<unsigned>::max()};
  //! status
  bool status_{false};
  //! Material properties
  Json properties_;
};  // Material class
}  // namespace mpm

#include "material.tcc"

#endif  // MPM_MATERIAL_MATERIAL_H_

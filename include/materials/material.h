#ifndef MPM_MATERIAL_MATERIAL_H_
#define MPM_MATERIAL_MATERIAL_H_

#include <limits>

#include "Eigen/Dense"
#include "json.hpp"

#include "factory.h"
#include "logger.h"
#include "map.h"
#include "particle.h"
#include "particle_base.h"
#include "material_utility.h"

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
  Material(unsigned id, const Json& material_properties) : id_{id} {
    //! Logger
    std::string logger = "material::" + std::to_string(id);
    console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
  }

  //! Destructor
  virtual ~Material(){};

  //! Delete copy constructor
  Material(const Material&) = delete;

  //! Delete assignement operator
  Material& operator=(const Material&) = delete;

  //! Return id of the material
  unsigned id() const { return id_; }

  //! Get material property
  //! \tparam Ttype Return type for proerpty
  //! \param[in] key Material property key
  //! \retval result Value of material property
  template <typename Ttype>
  Ttype property(const std::string& key);

  //! Initialise history variables
  virtual mpm::dense_map initialise_state_variables() = 0;

  //! Compute stress
  //! \param[in] stress Stress
  //! \param[in] dstrain Strain
  //! \param[in] particle Constant point to particle base
  //! \param[in] state_vars History-dependent state variables
  //! \retval updated_stress Updated value of stress
  virtual Vector6d compute_stress(const Vector6d& stress,
                                  const Vector6d& dstrain,
                                  const ParticleBase<Tdim>* ptr,
                                  mpm::dense_map* state_vars) = 0;

 protected:
  //! material id
  unsigned id_{std::numeric_limits<unsigned>::max()};
  //! Material properties
  Json properties_;
  //! Logger
  std::unique_ptr<spdlog::logger> console_;
};  // Material class
}  // namespace mpm

#include "material.tcc"

#endif  // MPM_MATERIAL_MATERIAL_H_

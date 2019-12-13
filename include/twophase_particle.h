#ifndef MPM_TWOPHASE_PARTICLE_H_
#define MPM_TWOPHASE_PARTICLE_H_

#include <array>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "logger.h"
#include "particle.h"

namespace mpm {

//! TwoPhaseParticle class
//! \brief Class that stores the information about second-phase (water)
//! particles
//! \tparam Tdim Dimension
template <unsigned Tdim>
class TwoPhaseParticle : public mpm::Particle<Tdim> {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Construct a particle with id and coordinates
  //! \param[in] id Particle id
  //! \param[in] coord Coordinates of the particles
  TwoPhaseParticle(Index id, const VectorDim& coord);

  //! Destructor
  ~TwoPhaseParticle() override{};

  //! Delete copy constructor
  TwoPhaseParticle(const TwoPhaseParticle<Tdim>&) = delete;

  //! Delete assignment operator
  TwoPhaseParticle& operator=(const TwoPhaseParticle<Tdim>&) = delete;

  //! Initialise liquid phase
  void initialise_liquid_phase() override;

  //! Assign material
  //! \param[in] material Pointer to a material
  bool assign_liquid_material(
      const std::shared_ptr<Material<Tdim>>& material) override;

  //! Compute liquid mass
  bool compute_liquid_mass() override;

  //! Return liquid mass
  double liquid_mass() const override { return liquid_mass_; }

  //! Assign liquid mass and momentum to nodes
  bool map_liquid_mass_momentum_to_nodes() override;

  //! Assign pore pressure to nodes
  bool map_pore_pressure_to_nodes() override;

  //! Compute pore pressure somoothening by interpolating nodal pressure
  bool compute_pore_pressure_smoothening() override;

  //! Map liquid body force
  void map_liquid_body_force() override;

  //! Map liquid internal force
  void map_liquid_internal_force() override;

  //! Map liquid traction force
  void map_liquid_traction_force() override;

  //! Map liquid drag force
  void map_liquid_drag_force() override;

  //! Compute pore pressure
  void compute_pore_pressure() override;

  //! Map two phase mixture internal force
  void map_twophase_mixture_internal_force() override;

  //! Assign liquid traction
  //! \param[in] direction Index corresponding to the direction of traction
  //! \param[in] traction Particle traction in specified direction
  //! \retval status Assignment status
  bool assign_liquid_traction(unsigned direction, double traction) override;

 private:
  //! Shape functions
  using Particle<Tdim>::shapefn_;

  //! Liquid mass
  double liquid_mass_;
  //! Degree of saturation
  double liquid_saturation{1.0};
  //! Liquid velocity
  Eigen::Matrix<double, Tdim, 1> liquid_velocity_;
  //! Pore pressure
  double pore_pressure_;
  //! Traction
  Eigen::Matrix<double, Tdim, 1> liquid_traction_;
  //! Material
  std::shared_ptr<Material<Tdim>> liquid_material_;
  //! Map of vector properties
  std::map<std::string, std::function<Eigen::VectorXd()>> liquid_properties_;
  //! Logger
  std::unique_ptr<spdlog::logger> console_;
};  // TwoPhaseParticle class
}  // namespace mpm

#include "twophase_particle.tcc"

#endif  // MPM_TWOPHASE_PARTICLE_H__

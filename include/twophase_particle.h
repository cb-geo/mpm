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

  //! Assign saturation degree
  bool assign_saturation_degree() override;

  //! Assign pore pressure
  //! \param[in] pressure Pore liquid pressure
  void assign_pore_pressure(const double& pressure) override {
    this->pore_pressure_ = pressure;
  }

  //! Assign liquid traction
  //! \param[in] direction Index corresponding to the direction of traction
  //! \param[in] traction Particle traction in specified direction
  //! \retval status Assignment status
  bool assign_liquid_traction(unsigned direction, double traction) override;

  //! Assign mixture traction
  //! \param[in] direction Index corresponding to the direction of traction
  //! \param[in] traction Particle traction in specified direction
  //! \retval status Assignment status
  bool assign_mixture_traction(unsigned direction, double traction) override;

  //! Compute liquid mass
  bool compute_liquid_mass() override;

  //! Return liquid mass
  //! \retval liquid mass Liquid phase mass
  double liquid_mass() const override { return liquid_mass_; }

  //! Assign liquid mass and momentum to nodes
  bool map_liquid_mass_momentum_to_nodes() override;

  //! Assign pore pressure to nodes
  bool map_pore_pressure_to_nodes() override;

  //! Compute pore pressure somoothening by interpolating nodal pressure
  bool compute_pore_pressure_smoothing() override;

  //! Map liquid body force
  //! \param[in] pgravity Gravity of a particle
  void map_liquid_body_force(const VectorDim& pgravity) override;

  //! Map liquid phase traction force
  void map_liquid_traction_force() override;

  //! Map liquid internal force
  bool map_liquid_internal_force() override;

  //! Compute pore pressure
  //! \param[in] dt Time step size
  bool compute_pore_pressure(double dt) override;

  //! Return liquid pore pressure
  //! \retval pore pressure Pore liquid pressure
  double pore_pressure() const override { return pore_pressure_; }

  //! Map two phase mixture body force
  //! \param[in] mixture Identification for Mixture
  //! \param[in] pgravity Gravity of the particle
  void map_mixture_body_force(unsigned mixture,
                              const VectorDim& pgravity) override;

  //! Map two phase mixture traction force
  //! \param[in] mixture Identification for Mixture
  void map_mixture_traction_force(unsigned mixture) override;

  //! Map two phase mixture internal force
  //! \param[in] mixture Identification for Mixture
  bool map_mixture_internal_force(unsigned mixture) override;

  //! Map drag force coefficient
  bool map_drag_force_coefficient() override;

  //! Compute updated velocity of the particle using nodal acceleration
  //! \param[in] dt Analysis time step
  //! \retval status Compute status
  bool compute_updated_liquid_kinematics(double dt) override;

  //! Compute updated velocity of the particle based on nodal velocity
  //! \param[in] dt Analysis time step
  //! \retval status Compute status
  bool compute_updated_liquid_velocity(double dt) override;

  //! Assign particle liquid phase velocity constraints
  //! Directions can take values between 0 and Dim
  //! \param[in] dir Direction of particle velocity constraint
  //! \param[in] velocity Applied particle liquid phase velocity constraint
  //! \retval status Assignment status
  bool assign_particle_liquid_velocity_constraint(unsigned dir,
                                                  double velocity) override;

  //! Apply particle liquid phase velocity constraints
  void apply_particle_liquid_velocity_constraints() override;

  //! Return vector data of particle liquid phase
  //! \param[in] property Property string
  //! \retval vecdata Vector data of particle liquid phase property
  Eigen::VectorXd liquid_vector_data(const std::string& property) override;

  //! Return velocity of the particle liquid phase
  //! \retval liquid velocity Liquid phase velocity
  VectorDim liquid_velocity() const override { return liquid_velocity_; }

 private:
  //! Shape functions
  using Particle<Tdim>::shapefn_;
  //! B matrix
  using Particle<Tdim>::bmatrix_;
  //! Effective stress of soil skeleton
  using Particle<Tdim>::stress_;
  //! Soil skeleton strain rate
  using Particle<Tdim>::strain_rate_;
  //! Soil skeleton material
  using Particle<Tdim>::material_;
  //! Particle total volume
  using ParticleBase<Tdim>::volume_;
  //! Porosity
  using ParticleBase<Tdim>::porosity_;
  //! Cell
  using ParticleBase<Tdim>::cell_;
  //! Material
  std::shared_ptr<Material<Tdim>> liquid_material_;
  //! Liquid material id
  unsigned liquid_material_id_{std::numeric_limits<unsigned>::max()};
  //! Liquid mass
  double liquid_mass_;
  //! Liquid mass density (bulk density = liquid mass / total volume)
  double liquid_mass_density_;
  //! Degree of saturation
  double liquid_saturation_{1.0};
  //! Displacement
  Eigen::Matrix<double, Tdim, 1> liquid_displacement_;
  //! Liquid velocity
  Eigen::Matrix<double, Tdim, 1> liquid_velocity_;
  //! Particle liquid phase velocity constraints
  std::map<unsigned, double> liquid_velocity_constraints_;
  //! Pore pressure
  double pore_pressure_;
  //! Set liquid phase traction
  bool set_liquid_traction_;
  //! Set mixture traction
  bool set_mixture_traction_;
  //! Traction for liquid phase
  Eigen::Matrix<double, Tdim, 1> liquid_traction_;
  //! Traction for mixture (soil skeleton + pore liquid)
  Eigen::Matrix<double, Tdim, 1> mixture_traction_;
  //! Liquid strain rate
  Eigen::Matrix<double, 6, 1> liquid_strain_rate_;  // delete if not needed
  //! Liquid strain rate
  Eigen::Matrix<double, 6, 1> liquid_strain_;  // delete if not needed
  //! Map of vector properties
  std::map<std::string, std::function<Eigen::VectorXd()>> liquid_properties_;
  //! Logger
  std::unique_ptr<spdlog::logger> console_;
};  // TwoPhaseParticle class
}  // namespace mpm

#include "twophase_particle.tcc"

#endif  // MPM_TWOPHASE_PARTICLE_H__

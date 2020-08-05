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

  //! Construct a twophase particle with id and coordinates
  //! \param[in] id Particle id
  //! \param[in] coord Coordinates of the particles
  TwoPhaseParticle(Index id, const VectorDim& coord);

  //! Construct a twophase particle with id, coordinates and status
  //! \param[in] id Particle id
  //! \param[in] coord coordinates of the particle
  //! \param[in] status Particle status (active / inactive)
  TwoPhaseParticle(Index id, const VectorDim& coord, bool status);

  //! Destructor
  ~TwoPhaseParticle() override{};

  //! Delete copy constructor
  TwoPhaseParticle(const TwoPhaseParticle<Tdim>&) = delete;

  //! Delete assignment operator
  TwoPhaseParticle& operator=(const TwoPhaseParticle<Tdim>&) = delete;

  //! Initialise particle from HDF5 data
  //! \param[in] particle HDF5 data of particle
  //! \retval status Status of reading HDF5 particle
  bool initialise_particle(const HDF5Particle& particle) override;

  //! Initialise particle HDF5 data and material
  //! \param[in] particle HDF5 data of particle
  //! \param[in] material Material associated with the particle
  //! \retval status Status of reading HDF5 particle
  bool initialise_particle(
      const HDF5Particle& particle,
      const std::shared_ptr<Material<Tdim>>& material) override;

  //! Initialise liquid phase
  void initialise() override;

  //! Retrun particle data as HDF5
  //! \retval particle HDF5 data of the particle
  HDF5Particle hdf5() const override;

  //! Assign saturation degree
  bool assign_saturation_degree() override;

  //! Compute both solid and liquid mass
  void compute_mass() noexcept override;

  //! Map particle mass and momentum to nodes (both solid and liquid)
  void map_mass_momentum_to_nodes() noexcept override;

  //! Map body force
  //! \param[in] pgravity Gravity of a particle
  void map_body_force(const VectorDim& pgravity) noexcept override;

  //! Map traction force
  void map_traction_force() noexcept override;

  //! Map internal force
  inline void map_internal_force() noexcept override;

  //! Compute updated position of the particle and kinematics of both solid and
  //! liquid phase \param[in] dt Analysis time step \param[in] velocity_update
  //! Update particle velocity from nodal vel when true
  void compute_updated_position(double dt,
                                bool velocity_update = false) noexcept override;

  //! Assign velocity to the particle liquid phase
  //! \param[in] velocity A vector of particle liquid phase velocity
  //! \retval status Assignment status
  bool assign_liquid_velocity(const VectorDim& velocity) override;

  //! Compute pore pressure
  //! \param[in] dt Time step size
  void compute_pore_pressure(double dt) noexcept override;

  //! Map drag force coefficient
  bool map_drag_force_coefficient() override;

  //! Assign particles initial pore pressure by watertable
  //! \param[in] dir_v Vertical direction (Gravity direction) of the watertable
  //! \param[in] dir_h Horizontal direction of the watertable
  //! \param[in] gravity Gravity vector
  //! \param[in] reference_points
  //! (Horizontal coordinate of borehole + height of 0 pore pressure)
  bool initialise_pore_pressure_watertable(
      const unsigned dir_v, const unsigned dir_h, VectorDim& gravity,
      std::map<double, double>& reference_points);

  //! Update porosity
  //! \param[in] dt Analysis time step
  void update_porosity(double dt) override;

  //! Assign particle permeability
  //! \retval status Assignment status
  bool assign_permeability() override;

  //! Assign porosity
  bool assign_porosity() override;

  //! Map particle pressure to nodes
  bool map_pressure_to_nodes(unsigned phase = mpm::ParticlePhase::Solid,
                             double dt = 0, Index step = 0) noexcept override;

  //! Apply particle velocity constraints
  //! \param[in] dir Direction of particle velocity constraint
  //! \param[in] velocity Applied particle velocity constraint
  void apply_particle_velocity_constraints(unsigned dir,
                                           double velocity) override;

  //! Assign traction to the particle
  //! \param[in] direction Index corresponding to the direction of traction
  //! \param[in] traction Particle traction in specified direction
  //! \retval status Assignment status
  bool assign_traction(unsigned direction, double traction) override;

  //! Return velocity of the particle liquid phase
  //! \retval liquid velocity Liquid phase velocity
  VectorDim liquid_velocity() const override { return liquid_velocity_; }

  //! Return liquid mass
  //! \retval liquid mass Liquid phase mass
  double liquid_mass() const override { return liquid_mass_; }

  //! Reture porosity
  //! \retval porosity Porosity
  double porosity() const override { return porosity_; }

 private:
  //! Assign liquid mass and momentum to nodes
  virtual void map_liquid_mass_momentum_to_nodes() noexcept;

  //! Map two phase mixture body force
  //! \param[in] mixture Identification for Mixture
  //! \param[in] pgravity Gravity of the particle
  virtual void map_mixture_body_force(unsigned mixture,
                                      const VectorDim& pgravity) noexcept;

  //! Map liquid body force
  //! \param[in] pgravity Gravity of a particle
  virtual void map_liquid_body_force(const VectorDim& pgravity) noexcept;

  //! Map two phase mixture traction force
  virtual void map_mixture_traction_force() noexcept;

  //! Map two phase liquid traction force
  virtual void map_liquid_traction_force() noexcept;

  //! Map liquid internal force
  virtual void map_liquid_internal_force() noexcept;

  //! Map two phase mixture internal force
  //! \param[in] mixture Identification for Mixture
  virtual void map_mixture_internal_force(unsigned mixture) noexcept;

  //! Compute updated velocity of the particle based on nodal velocity
  //! \param[in] dt Analysis time step
  //! \retval status Compute status
  virtual void compute_updated_liquid_velocity(double dt,
                                               bool velocity_update) noexcept;

 protected:
  //! coordinates
  using ParticleBase<Tdim>::coordinates_;
  //! Cell
  using ParticleBase<Tdim>::cell_;
  //! Nodes
  using ParticleBase<Tdim>::nodes_;
  //! State variables
  using ParticleBase<Tdim>::state_variables_;
  //! Shape functions
  using Particle<Tdim>::shapefn_;
  //! Effective stress of soil skeleton
  using Particle<Tdim>::stress_;
  //! Soil skeleton strain rate
  using Particle<Tdim>::strain_rate_;
  //! Soil skeleton material
  using Particle<Tdim>::material_;
  //! Material id
  using ParticleBase<Tdim>::material_id_;
  //! Particle total volume
  using Particle<Tdim>::volume_;
  //! Particle mass for solid phase
  using Particle<Tdim>::mass_;
  //! Particle mass density
  using Particle<Tdim>::mass_density_;
  //! dN/dX
  using Particle<Tdim>::dn_dx_;
  //! dN/dX at cell centroid
  using Particle<Tdim>::dn_dx_centroid_;
  //! Set traction
  using Particle<Tdim>::set_traction_;
  //! Surface Traction (given as a stress; force/area)
  using Particle<Tdim>::traction_;
  //! Size of particle
  using Particle<Tdim>::size_;
  //! Particle velocity constraints
  using Particle<Tdim>::particle_velocity_constraints_;

  //! Liquid mass
  double liquid_mass_;
  //! Liquid mass density (bulk density = liquid mass / total volume)
  double liquid_mass_density_;
  //! Degree of saturation
  double liquid_saturation_{1.0};
  //! Material point porosity (volume of voids / total volume)
  double porosity_{0.0};
  //! Set liquid traction
  bool set_liquid_traction_{false};
  //! Liquid traction
  Eigen::Matrix<double, Tdim, 1> liquid_traction_;
  //! Liquid velocity
  Eigen::Matrix<double, Tdim, 1> liquid_velocity_;
  //! Pore pressure constraint
  double pore_pressure_constraint_{std::numeric_limits<unsigned>::max()};
  //! Permeability parameter c1 (k = k_p * c1)
  VectorDim permeability_;
  //! Logger
  std::unique_ptr<spdlog::logger> console_;

};  // TwoPhaseParticle class
}  // namespace mpm

#include "twophase_particle.tcc"

#endif  // MPM_TWOPHASE_PARTICLE_H__

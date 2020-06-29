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
//! \brief Class that stores the information about twophase particles
//! Derive from particle
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
  void initialise_liquid_phase() override;

  //! Retrun particle data as HDF5
  //! \retval particle HDF5 data of the particle
  HDF5Particle hdf5() const override;

  //! Assign material
  //! \param[in] material Pointer to a material
  bool assign_liquid_material(
      const std::shared_ptr<Material<Tdim>>& material) override;

  //! Assign porosity
  bool assign_porosity() override;

  //! Assign pore pressure
  //! \param[in] pressure Pore liquid pressure
  void assign_pore_pressure(double pressure) override {
    this->pore_pressure_ = pressure;
  }

  //! Assign liquid traction
  //! \param[in] direction Index corresponding to the direction of traction
  //! \param[in] traction Particle traction in specified direction
  //! \retval status Assignment status
  bool assign_liquid_traction(unsigned direction, double traction) override;

  //! Return liquid phase traction
  VectorDim liquid_traction() const { return liquid_traction_; };

  //! Compute both solid and liquid mass
  void compute_mass() noexcept override;

  //! Map particle mass and momentum to nodes (both solid and liquid)
  void map_mass_momentum_to_nodes() noexcept override;

  //! Map body force
  //! \param[in] pgravity Gravity of a particle
  void map_body_force(const VectorDim& pgravity) noexcept override;

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

  //! Return velocity of the particle liquid phase
  //! \retval liquid velocity Liquid phase velocity
  VectorDim liquid_velocity() const override { return liquid_velocity_; }

  //! Return liquid strain of the particle
  Eigen::Matrix<double, 6, 1> liquid_strain() const override {
    return liquid_strain_;
  }

  //! Return liquid mass
  //! \retval liquid mass Liquid phase mass
  double liquid_mass() const override { return liquid_mass_; }

  //! Assign pore pressure to nodes
  void map_pore_pressure_to_nodes() noexcept override;

  //! Compute pore pressure somoothening by interpolating nodal pressure
  bool compute_pore_pressure_smoothing() noexcept override;

  //! Compute pore pressure
  //! \param[in] dt Time step size
  void compute_pore_pressure(double dt) noexcept override;

  //! Return pore pressure
  //! \retval pore_pressure Pore pressure
  double pore_pressure() const override { return pore_pressure_; }

  //! Return excessive pore pressure
  //! \retval excessive pore pressure Excessive pore pressure
  double excessive_pore_pressure() const override {
    return excessive_pore_pressure_;
  }

  //! Return free surface status
  //! \retval free surface Free surface status
  bool free_surface() const override { return free_surface_; }

  //! Map drag force coefficient
  bool map_drag_force_coefficient() override;

  //! Update particle permeability
  //! \retval status Update status
  VectorDim update_permeability() override;

  //! Update porosity
  //! \param[in] dt Analysis time step
  bool update_porosity(double dt) override;

  //! Assign particle liquid phase velocity constraints
  //! Directions can take values between 0 and Dim
  //! \param[in] dir Direction of particle velocity constraint
  //! \param[in] velocity Applied particle liquid phase velocity constraint
  //! \retval status Assignment status
  bool assign_particle_liquid_velocity_constraint(unsigned dir,
                                                  double velocity) override;

  //! Apply particle liquid phase velocity constraints
  void apply_particle_liquid_velocity_constraints() override;

  //! Assign particle pressure constraints
  //! \retval status Assignment status
  bool assign_particle_pore_pressure_constraint(double pressure) override;

  //! Assign particles initial pore pressure by watertable
  bool initialise_pore_pressure_watertable(
      const unsigned dir_v, const unsigned dir_h,
      std::map<double, double>& refernece_points);

 private:
  //! Assign particle permeability
  //! \retval status Assignment status
  virtual bool assign_permeability();

  //! Compute liquid mass
  virtual void compute_liquid_mass() noexcept;

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

  //! Map liquid phase traction force
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
  //! Shape functions
  using Particle<Tdim>::shapefn_;
  //! Effective stress of soil skeleton
  using Particle<Tdim>::stress_;
  //! Soil skeleton strain rate
  using Particle<Tdim>::strain_rate_;
  //! Soil skeleton material
  using Particle<Tdim>::material_;
  //! Particle total volume
  using Particle<Tdim>::volume_;
  //! Particle mass density
  using Particle<Tdim>::mass_density_;
  //! dN/dX
  using Particle<Tdim>::dn_dx_;
  //! dN/dX at cell centroid
  using Particle<Tdim>::dn_dx_centroid_;
  //! Solid mass
  using Particle<Tdim>::mass_;
  //! Set traction
  using Particle<Tdim>::set_traction_;
  //! Material
  std::shared_ptr<Material<Tdim>> liquid_material_;
  //! Liquid material id
  unsigned liquid_material_id_{std::numeric_limits<unsigned>::max()};
  //! Liquid mass density (bulk density = liquid mass / total volume)
  double liquid_mass_density_{0.};
  //! Liquid mass
  double liquid_mass_{0.};
  //! Porosity
  double porosity_{0.};
  //! Permeability c1 (k = k_p * c1_)
  Eigen::Matrix<double, Tdim, 1> permeability_c1_;
  //! Liquid velocity
  Eigen::Matrix<double, Tdim, 1> liquid_velocity_;
  //! Particle liquid phase velocity constraints
  std::map<unsigned, double> liquid_velocity_constraints_;
  //! Pore pressure
  double pore_pressure_{0.};
  //! Excessive pore pressure
  double excessive_pore_pressure_{0.};
  //! Free surface
  bool free_surface_{false};
  //! Pore pressure constraint
  double pore_pressure_constraint_{std::numeric_limits<double>::max()};
  //! Set liquid phase traction
  bool set_liquid_traction_{false};
  //! Traction for liquid phase
  Eigen::Matrix<double, Tdim, 1> liquid_traction_;
  //! Liquid strain rate
  Eigen::Matrix<double, 6, 1> liquid_strain_rate_;
  //! Liquid strain rate
  Eigen::Matrix<double, 6, 1> liquid_strain_;
  //! Logger
  std::unique_ptr<spdlog::logger> console_;
};  // TwoPhaseParticle class
}  // namespace mpm

#include "twophase_particle.tcc"

#endif  // MPM_TWOPHASE_PARTICLE_H__

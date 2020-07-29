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

  //! Return mixture traction
  VectorDim mixture_traction() const { return mixture_traction_; };

  //! Compute both solid and liquid mass
  void compute_mass() noexcept override;

  //! Map particle mass and momentum to nodes (both solid and liquid)
  void map_mass_momentum_to_nodes() noexcept override;

  //! Map body force
  //! \param[in] pgravity Gravity of a particle
  void map_body_force(const VectorDim& pgravity) noexcept override;

  //! Assign traction to the particle
  //! \param[in] direction Index corresponding to the direction of traction
  //! \param[in] traction Particle traction in specified direction
  //! \retval status Assignment status
  bool assign_traction(unsigned direction, double traction) override;

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

  //! Compute pore pressure somoothening by interpolating nodal pressure
  bool compute_pore_pressure_smoothing() noexcept override;

  //! Compute pore pressure
  //! \param[in] dt Time step size
  void compute_pore_pressure(double dt) noexcept override;

  //! Return liquid pore pressure
  //! \retval pore pressure Pore liquid pressure
  double pore_pressure() const override { return pore_pressure_; }

  //! Map drag force coefficient
  bool map_drag_force_coefficient() override;

  //! Update particle permeability
  //! \retval status Update status
  bool update_permeability() override;

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

  //! Initial pore pressure
  //! \param[in] pore pressure Initial pore pressure
  void initial_pore_pressure(double pore_pressure) override {
    this->pore_pressure_ = pore_pressure;
  }

  //! Assign particles initial pore pressure by watertable
  bool initialise_pore_pressure_watertable(
      const unsigned dir_v, const unsigned dir_h,
      std::map<double, double>& refernece_points);

  //! Update porosity
  //! \param[in] dt Analysis time step
  bool update_porosity(double dt) override;

  //! Assign particle permeability
  //! \retval status Assignment status
  bool assign_permeability() override;

  //! Assign porosity
  bool assign_porosity() override;

  //! Map particle pressure to nodes
  bool map_pressure_to_nodes(
      unsigned phase = mpm::ParticlePhase::Solid) noexcept override;

  //! Compute pressure smoothing of the particle based on nodal pressure
  //! $$\hat{p}_p = \sum_{i = 1}^{n_n} N_i(x_p) p_i$$
  bool compute_pressure_smoothing(
      unsigned phase = mpm::ParticlePhase::Solid) noexcept override;

 private:
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

  //! Assign mixture traction
  //! \param[in] direction Index corresponding to the direction of traction
  //! \param[in] traction Particle traction in specified direction
  //! \retval status Assignment status
  virtual bool assign_mixture_traction(unsigned direction, double traction);

  //! Map liquid phase traction force
  virtual void map_liquid_traction_force() noexcept;

  //! Map two phase mixture traction force
  //! \param[in] mixture Identification for Mixture
  virtual void map_mixture_traction_force(unsigned mixture) noexcept;

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
  //! Particle mass density
  using Particle<Tdim>::mass_density_;
  //! Particle mass for solid phase
  using Particle<Tdim>::mass_;
  //! dN/dX
  using Particle<Tdim>::dn_dx_;
  //! dN/dX at cell centroid
  using Particle<Tdim>::dn_dx_centroid_;

  //! Liquid mass
  double liquid_mass_;
  //! Liquid mass density (bulk density = liquid mass / total volume)
  double liquid_mass_density_;
  //! Degree of saturation
  double liquid_saturation_{1.0};
  //! Material point porosity (volume of voids / total volume)
  double porosity_{0.0};
  //! Liquid velocity
  Eigen::Matrix<double, Tdim, 1> liquid_velocity_;
  //! Particle liquid phase velocity constraints
  std::map<unsigned, double> liquid_velocity_constraints_;
  //! Pore pressure constraint
  double pore_pressure_constraint_{std::numeric_limits<unsigned>::max()};
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
  //! Logger
  std::unique_ptr<spdlog::logger> console_;
  //! Permeability
  VectorDim permeability_;
  //! Permeability parameter
  VectorDim c1_;
  //! reference pore pressure
  double reference_pore_pressure_{0};
};  // TwoPhaseParticle class
}  // namespace mpm

#include "twophase_particle.tcc"

#endif  // MPM_TWOPHASE_PARTICLE_H__

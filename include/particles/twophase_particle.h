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

  //! Retrun particle data as HDF5
  //! \retval particle HDF5 data of the particle
  HDF5Particle hdf5() const override;

  //! Initialise properties
  void initialise() override;

  //! Map internal force
  inline void map_internal_force() noexcept override;

  //! Compute updated position of the particle and kinematics of both solid and
  //! liquid phase \param[in] dt Analysis time step \param[in] velocity_update
  //! Update particle velocity from nodal vel when true
  void compute_updated_position(double dt,
                                bool velocity_update = false) noexcept override;
  //-------------------------------------------------------------------------------
  //! TODO
  //! Assign liquid traction
  //! \param[in] direction Index corresponding to the direction of traction
  //! \param[in] traction Particle traction in specified direction
  //! \retval status Assignment status
  bool assign_liquid_traction(unsigned direction, double traction) override;

  //! Compute pore pressure
  //! \param[in] dt Time step size
  void compute_pore_pressure(double dt) noexcept override;

  //! Compute pore pressure somoothening by interpolating nodal pressure
  bool compute_pore_pressure_smoothing() noexcept override;

  //! Assign particle liquid phase velocity constraints
  //! Directions can take values between 0 and Dim
  //! \param[in] dir Direction of particle velocity constraint
  //! \param[in] velocity Applied particle liquid phase velocity constraint
  //! \retval status Assignment status
  bool assign_particle_liquid_velocity_constraint(unsigned dir,
                                                  double velocity) override;

  //! Assign particle pressure constraints
  //! \retval status Assignment status
  bool assign_particle_pore_pressure_constraint(double pressure) override;
  //-------------------------------------------------------------------------------

 private:
  //! Compute updated velocity of the particle based on nodal velocity
  //! \param[in] dt Analysis time step
  //! \retval status Compute status
  void compute_updated_liquid_velocity(double dt,
                                       bool velocity_update) noexcept;

  //! Apply particle liquid phase velocity constraints
  void apply_particle_liquid_velocity_constraints();

  //! Map liquid phase traction force
  void map_liquid_traction_force() noexcept;

  //! Return mass of the particle
  double liquid_mass() const {
    return this->scalar_property(mpm::properties::Scalar::LiquidMass);
  }

  //! Return porosity of the particle
  double porosity() const {
    return this->scalar_property(mpm::properties::Scalar::Porosity);
  }

  //! Return liquid strain of the particle
  Eigen::Matrix<double, 6, 1> liquid_strain() const { return liquid_strain_; }

  //! Return liquid velocity of the particle
  VectorDim liquid_velocity() const {
    return this->vector_property(mpm::properties::Vector::LiquidVelocity);
  }

  //! Return pore pressure of the particle
  double pore_pressure() const {
    return this->scalar_property(mpm::properties::Scalar::PorePressure);
  }

  //! Return free surface status
  bool free_surface() const { return free_surface_; }

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
  //! dN/dX
  using Particle<Tdim>::dn_dx_;
  //! dN/dX at cell centroid
  using Particle<Tdim>::dn_dx_centroid_;
  //! Scalar properties
  using ParticleBase<Tdim>::scalar_properties_;
  //! Vector properties
  using ParticleBase<Tdim>::vector_properties_;
  //! Particle liquid phase velocity constraints
  std::map<unsigned, double> liquid_velocity_constraints_;
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
#include "twophase_particle_functions.tcc"

#endif  // MPM_TWOPHASE_PARTICLE_H__

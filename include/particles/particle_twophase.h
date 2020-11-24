#ifndef MPM_PARTICLE_TWOPHASE_H_
#define MPM_PARTICLE_TWOPHASE_H_

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

  //! Initialise particle from POD data
  //! \param[in] particle POD data of particle
  //! \retval status Status of reading POD particle
  bool initialise_particle(PODParticle& particle) override;

  //! Initialise particle POD data and material
  //! \param[in] particle POD data of particle
  //! \param[in] materials Material associated with the particle arranged in a
  //! vector
  //! \retval status Status of reading POD particle
  bool initialise_particle(
      PODParticle& particle,
      const std::vector<std::shared_ptr<Material<Tdim>>>& materials) override;

  //! Initialise particle liquid phase on top of the regular solid phase
  void initialise() override;

  //! Return particle data as POD
  //! \retval particle POD of the particle
  std::shared_ptr<void> pod() const override;

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
      const unsigned dir_v, const unsigned dir_h, const VectorDim& gravity,
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
  bool map_pressure_to_nodes(
      unsigned phase = mpm::ParticlePhase::Solid) noexcept override;

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

  //! Type of particle
  std::string type() const override {
    return (Tdim == 2) ? "P2D2PHASE" : "P3D2PHASE";
  }

  //! Serialize
  //! \retval buffer Serialized buffer data
  std::vector<uint8_t> serialize() override;

  //! Deserialize
  //! \param[in] buffer Serialized buffer data
  //! \param[in] material Particle material pointers
  void deserialize(
      const std::vector<uint8_t>& buffer,
      std::vector<std::shared_ptr<mpm::Material<Tdim>>>& materials) override;

  //! ----------------------------------------------------------------
  //! Semi-Implicit integration functions based on Chorin's Projection
  //! ----------------------------------------------------------------

  //! Assigning beta parameter to particle
  //! \param[in] pressure parameter determining type of projection
  void assign_projection_parameter(double parameter) override {
    this->projection_param_ = parameter;
  };

  //! Map drag matrix to cell assuming linear-darcy drag force
  bool map_drag_matrix_to_cell() override;

  //! Update pressure after solving poisson equation
  bool compute_updated_pressure() override;

  //! Map laplacian element matrix to cell (used in poisson equation LHS)
  bool map_laplacian_to_cell() override;

  //! Map poisson rhs element matrix to cell (used in poisson equation RHS)
  bool map_poisson_right_to_cell() override;

  //! Map correction matrix element matrix to cell (used to correct velocity)
  bool map_correction_matrix_to_cell() override;

 protected:
  //! Compute pack size
  //! \retval pack size of serialized object
  int compute_pack_size() const override;

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
  virtual void map_mixture_internal_force() noexcept;

  //! Compute updated velocity of the particle based on nodal velocity
  //! \param[in] dt Analysis time step
  //! \retval status Compute status
  virtual void compute_updated_liquid_velocity(double dt,
                                               bool velocity_update) noexcept;

 protected:
  //! particle id
  using ParticleBase<Tdim>::id_;
  //! coordinates
  using ParticleBase<Tdim>::coordinates_;
  //! Status
  using ParticleBase<Tdim>::status_;
  //! Cell
  using ParticleBase<Tdim>::cell_;
  //! Cell id
  using ParticleBase<Tdim>::cell_id_;
  //! Nodes
  using ParticleBase<Tdim>::nodes_;
  //! State variables
  using ParticleBase<Tdim>::state_variables_;
  //! Shape functions
  using Particle<Tdim>::shapefn_;
  //! dN/dX
  using Particle<Tdim>::dn_dx_;
  //! dN/dX at cell centroid
  using Particle<Tdim>::dn_dx_centroid_;
  //! Size of particle in natural coordinates
  using Particle<Tdim>::natural_size_;
  //! Materials
  using Particle<Tdim>::material_;
  //! Material ids
  using ParticleBase<Tdim>::material_id_;
  //! Particle mass for solid phase
  using Particle<Tdim>::mass_;
  //! Particle total volume
  using Particle<Tdim>::volume_;
  //! Particle mass density
  using Particle<Tdim>::mass_density_;
  //! Displacement
  using Particle<Tdim>::displacement_;
  //! Velocity
  using Particle<Tdim>::velocity_;
  //! Effective stress of soil skeleton
  using Particle<Tdim>::stress_;
  //! Solid skeleton strains
  using Particle<Tdim>::strain_;
  //! Volumetric strain at centroid
  using Particle<Tdim>::volumetric_strain_centroid_;
  //! Soil skeleton strain rate
  using Particle<Tdim>::strain_rate_;
  //! Set traction
  using Particle<Tdim>::set_traction_;
  //! Surface Traction (given as a stress; force/area)
  using Particle<Tdim>::traction_;
  //! Size of particle
  using Particle<Tdim>::size_;
  //! Particle velocity constraints
  using Particle<Tdim>::particle_velocity_constraints_;
  //! Size of particle
  using Particle<Tdim>::pack_size_;

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
  //! Projection parameter for semi-implicit update
  double projection_param_{1.0};
  //! Pore pressure constraint
  double pore_pressure_constraint_{std::numeric_limits<unsigned>::max()};
  //! Permeability parameter c1 (k = k_p * c1)
  VectorDim permeability_;
  //! Logger
  std::unique_ptr<spdlog::logger> console_;

};  // TwoPhaseParticle class
}  // namespace mpm

#include "particle_twophase.tcc"

#endif  // MPM_PARTICLE_TWOPHASE_H__

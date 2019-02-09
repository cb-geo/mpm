#ifndef MPM_PARTICLE_H_
#define MPM_PARTICLE_H_

#include <array>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "cell.h"
#include "logger.h"
#include "particle_base.h"

namespace mpm {

//! Global index type for the particle
using Index = unsigned long long;

//! Particle class
//! \brief Base class that stores the information about particles
//! \details Particle class: id_ and coordinates.
//! \tparam Tdim Dimension
//! \tparam Tnphases Number of phases
template <unsigned Tdim, unsigned Tnphases>
class Particle : public ParticleBase<Tdim> {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Define DOFs
  static const unsigned Tdof = (Tdim == 1) ? 1 : 3 * (Tdim - 1);

  //! Construct a particle with id and coordinates
  //! \param[in] id Particle id
  //! \param[in] coord coordinates of the particle
  Particle(Index id, const VectorDim& coord);

  //! Construct a particle with id, coordinates and status
  //! \param[in] id Particle id
  //! \param[in] coord coordinates of the particle
  //! \param[in] status Particle status (active / inactive)
  Particle(Index id, const VectorDim& coord, bool status);

  //! Destructor
  ~Particle() override{};

  //! Delete copy constructor
  Particle(const Particle<Tdim, Tnphases>&) = delete;

  //! Delete assignment operator
  Particle& operator=(const Particle<Tdim, Tnphases>&) = delete;

  //! Initialise particle from HDF5 data
  //! \param[in] particle HDF5 data of particle
  //! \retval status Status of reading HDF5 particle
  bool initialise_particle(const HDF5Particle& particle) override;

  //! Initialise properties
  void initialise() override;

  //! Compute reference coordinates in a cell
  bool compute_reference_location() override;

  //! Return reference location
  VectorDim reference_location() const override { return xi_; }

  //! Assign a cell to particle
  //! If point is in new cell, assign new cell and remove particle id from old
  //! cell. If point can't be found in the new cell, check if particle is still
  //! valid in the old cell, if it is leave it as is. If not, set cell as null
  //! \param[in] cellptr Pointer to a cell
  bool assign_cell(const std::shared_ptr<Cell<Tdim>>& cellptr) override;

  //! Assign cell id
  //! \param[in] id Cell id
  bool assign_cell_id(Index id) override;

  //! Return cell id
  Index cell_id() const override { return cell_id_; }

  //! Return cell ptr status
  bool cell_ptr() const override { return cell_ != nullptr; }

  //! Remove cell associated with the particle
  void remove_cell() override;

  //! Compute shape functions of a particle, based on local coordinates
  bool compute_shapefn() override;

  //! Assign volume
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] volume Volume of particle for the phase
  bool assign_volume(unsigned phase, double volume) override;

  //! Return volume
  //! \param[in] phase Index corresponding to the phase
  double volume(unsigned phase) const override { return volume_(phase); }

  //! Return size of particle in natural coordinates
  VectorDim natural_size() const override { return natural_size_; }

  //! Compute volume as cell volume / nparticles
  //! \param[in] phase Index corresponding to the phase
  bool compute_volume(unsigned phase) override;

  //! Update volume based on centre volumetric strain rate
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] dt Analysis time step
  bool update_volume_strainrate(unsigned phase, double dt) override;

  //! Compute mass as volume * density
  //! \param[in] phase Index corresponding to the phase
  bool compute_mass(unsigned phase) override;

  //! Map particle mass and momentum to nodes
  //! \param[in] phase Index corresponding to the phase
  bool map_mass_momentum_to_nodes(unsigned phase) override;

  //! Assign nodal mass to particles
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] mass Mass from the particles in a cell
  //! \retval status Assignment status
  void assign_mass(unsigned phase, double mass) override {
    mass_(phase) = mass;
  }

  //! Return mass of the particles
  //! \param[in] phase Index corresponding to the phase
  double mass(unsigned phase) const override { return mass_(phase); }

  //! Assign material
  //! \param[in] material Pointer to a material
  bool assign_material(
      const std::shared_ptr<Material<Tdim>>& material) override;

  //! Return pressure of the particles
  //! \param[in] phase Index corresponding to the phase
  double pressure(unsigned phase) const override { return pressure_(phase); }

  //! Compute strain
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] dt Analysis time step
  void compute_strain(unsigned phase, double dt) override;

  //! Return strain of the particle
  //! \param[in] phase Index corresponding to the phase
  Eigen::Matrix<double, 6, 1> strain(unsigned phase) const override {
    return strain_.col(phase);
  }

  //! Return strain rate of the particle
  //! \param[in] phase Index corresponding to the phase
  Eigen::Matrix<double, 6, 1> strain_rate(unsigned phase) const override {
    return strain_rate_.col(phase);
  };

  //! Return volumetric strain of centroid
  //! \param[in] phase Index corresponding to the phase
  //! \retval volumetric strain at centroid
  double volumetric_strain_centroid(unsigned phase) const override {
    return volumetric_strain_centroid_(phase);
  }

  //! Initial stress
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] stress Initial sress corresponding to the phase
  void initial_stress(unsigned phase,
                      const Eigen::Matrix<double, 6, 1>& stress) override {
    this->stress_.col(phase) = stress;
  }

  //! Compute stress
  bool compute_stress(unsigned phase) override;

  //! Return stress of the particle
  //! \param[in] phase Index corresponding to the phase
  Eigen::Matrix<double, 6, 1> stress(unsigned phase) const override {
    return stress_.col(phase);
  }

  //! Map body force
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] pgravity Gravity of a particle
  void map_body_force(unsigned phase, const VectorDim& pgravity) override;

  //! Map internal force
  //! \param[in] phase Index corresponding to the phase
  bool map_internal_force(unsigned phase) override;

  //! Assign velocity to the particle
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] velocity A vector of particle velocity
  //! \retval status Assignment status
  bool assign_velocity(unsigned phase, const VectorDim& velocity) override;

  //! Return velocity of the particle
  //! \param[in] phase Index corresponding to the phase
  VectorDim velocity(unsigned phase) const override {
    return velocity_.col(phase);
  }

  //! Assign traction to the particle
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] direction Index corresponding to the direction of traction
  //! \param[in] traction Particle traction in specified direction
  //! \retval status Assignment status
  bool assign_traction(unsigned phase, unsigned direction,
                       double traction) override;

  //! Return traction of the particle
  //! \param[in] phase Index corresponding to the phase
  VectorDim traction(unsigned phase) const override {
    return traction_.col(phase);
  }

  //! Map traction force
  //! \param[in] phase Index corresponding to the phase
  void map_traction_force(unsigned phase) override;

  //! Compute updated position of the particle
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] dt Analysis time step
  bool compute_updated_position(unsigned phase, double dt) override;

  //! Compute updated position of the particle based on nodal velocity
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] dt Analysis time step
  bool compute_updated_position_velocity(unsigned phase, double dt) override;

  //! Return a state variable
  //! \param[in] var State variable
  //! \retval Quantity of the state history variable
  double state_variable(const std::string& var) const override {
    return state_variables_.at(var);
  }

 private:
  //! Update pressure of the particles
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] dvolumetric_strain dvolumetric strain in a cell
  bool update_pressure(unsigned phase, double dvolumetric_strain);

  //! particle id
  using ParticleBase<Tdim>::id_;
  //! coordinates
  using ParticleBase<Tdim>::coordinates_;
  //! Reference coordinates (in a cell)
  using ParticleBase<Tdim>::xi_;
  //! Cell
  using ParticleBase<Tdim>::cell_;
  //! Cell id
  using ParticleBase<Tdim>::cell_id_;
  //! Status
  using ParticleBase<Tdim>::status_;
  //! Material
  using ParticleBase<Tdim>::material_;
  //! State variables
  using ParticleBase<Tdim>::state_variables_;
  //! Mass
  Eigen::Matrix<double, 1, Tnphases> mass_;
  //! Volume
  Eigen::Matrix<double, 1, Tnphases> volume_;
  //! Size of particle
  Eigen::Matrix<double, 1, Tdim> size_;
  //! Size of particle in natural coordinates
  Eigen::Matrix<double, 1, Tdim> natural_size_;
  //! Pressure
  Eigen::Matrix<double, 1, Tnphases> pressure_;
  //! Stresses
  Eigen::Matrix<double, 6, Tnphases> stress_;
  //! Strains
  Eigen::Matrix<double, 6, Tnphases> strain_;
  //! Volumetric strain at centroid
  Eigen::Matrix<double, Tnphases, 1> volumetric_strain_centroid_;
  //! Strain rate
  Eigen::Matrix<double, 6, Tnphases> strain_rate_;
  //! dstrains
  Eigen::Matrix<double, 6, Tnphases> dstrain_;
  //! Velocity
  Eigen::Matrix<double, Tdim, Tnphases> velocity_;
  //! Set traction
  bool set_traction_{false};
  //! Traction
  Eigen::Matrix<double, Tdim, Tnphases> traction_;
  //! Shape functions
  Eigen::VectorXd shapefn_;
  //! B-Matrix
  std::vector<Eigen::MatrixXd> bmatrix_;
  //! Logger
  std::unique_ptr<spdlog::logger> console_;
};  // Particle class
}  // namespace mpm

#include "particle.tcc"

#endif  // MPM_PARTICLE_H__

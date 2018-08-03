#ifndef MPM_PARTICLE_H_
#define MPM_PARTICLE_H_

#include <array>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "cell.h"
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

  //! Delete assignement operator
  Particle& operator=(const Particle<Tdim, Tnphases>&) = delete;

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

  //! Return cell id
  Index cell_id() const override { return cell_id_; }

  //! Remove cell associated with the particle
  void remove_cell() override;

  //! Compute shape functions of a particle, based on local coordinates
  bool compute_shapefn() override;

  //! Assign volume
  void assign_volume(double volume) override { volume_ = volume; }

  //! Return volume
  double volume() const override { return volume_; }

  //! Compute volume as cell volume / nparticles
  bool compute_volume() override;

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
    mass_(0, phase) = mass;
  }

  //! Return mass of the particlesx
  //! \param[in] phase Index corresponding to the phase
  double mass(unsigned phase) const override { return mass_(0, phase); }

  //! Assign material
  //! \param[in] material Pointer to a material
  bool assign_material(
      const std::shared_ptr<Material<Tdim>>& material) override;

  //! Compute strain
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] dt Analysis time step
  void compute_strain(unsigned phase, double dt) override;

  //! Return strain of the particle
  //! \param[in] phase Index corresponding to the phase
  Eigen::Matrix<double, 6, 1> strain(unsigned phase) const override {
    return strain_.col(phase);
  }

  //! Compute stress
  bool compute_stress(unsigned phase) override;

  //! Assign stress to the particle
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] stress A vector of particle stress
  void assign_stress(unsigned phase,
                     const Eigen::Matrix<double, 6, 1>& stress) override;

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
  bool assign_velocity(unsigned phase,
                       const Eigen::VectorXd& velocity) override;

  //! Return velocity of the particle
  //! \param[in] phase Index corresponding to the phase
  Eigen::VectorXd velocity(unsigned phase) const override {
    return velocity_.col(phase);
  }

  //! Assign momentum to the particle
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] momentum A vector of particle momentum
  //! \retval status Assignment status
  bool assign_momentum(unsigned phase,
                       const Eigen::VectorXd& momentum) override;

  //! Return momentum of the particle
  //! \param[in] phase Index corresponding to the phase
  Eigen::VectorXd momentum(unsigned phase) const override {
    return momentum_.col(phase);
  }

  //! Assign acceleration to the particle
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] acceleration A vector of particle acceleration
  //! \retval status Assignment status
  bool assign_acceleration(unsigned phase,
                           const Eigen::VectorXd& acceleration) override;

  //! Return acceleration of the particle
  //! \param[in] phase Index corresponding to the phase
  Eigen::VectorXd acceleration(unsigned phase) const override {
    return acceleration_.col(phase);
  }

  //! Compute updated position of the particle
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] dt Analysis time step
  bool compute_updated_position(unsigned phase, double dt) override;

  //! TODO: Remove
  void stats() override {
    std::string out = "Particle: " + std::to_string(id_) + "\t" +
                      std::to_string(status_) + " position: ";
    std::string val = "";
    for (unsigned i = 0; i < coordinates_.size(); ++i)
      val += std::to_string(coordinates_(i, 0)) + "\t";
    out += val + "\n";

    std::cout << out;
  }

 private:
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
  //! Volume
  using ParticleBase<Tdim>::volume_;
  //! Material
  using ParticleBase<Tdim>::material_;
  //! Mass
  Eigen::Matrix<double, 1, Tnphases> mass_;
  //! Stresses
  Eigen::Matrix<double, 6, Tnphases> stress_;
  //! Strains
  Eigen::Matrix<double, 6, Tnphases> strain_;
  //! dstrains
  Eigen::Matrix<double, 6, Tnphases> dstrain_;
  //! Velocity
  Eigen::Matrix<double, Tdim, Tnphases> velocity_;
  //! Momentum
  Eigen::Matrix<double, Tdim, Tnphases> momentum_;
  //! Acceleration
  Eigen::Matrix<double, Tdim, Tnphases> acceleration_;
  //! Shape functions
  Eigen::VectorXd shapefn_;
  //! Gradient of shape functions
  Eigen::MatrixXd grad_shapefn_;
  //! B-Matrix
  std::vector<Eigen::MatrixXd> bmatrix_;

};  // Particle class
}  // namespace mpm

#include "particle.tcc"

#endif  // MPM_PARTICLE_H__

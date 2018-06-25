#ifndef MPM_PARTICLE_H_
#define MPM_PARTICLE_H_

#include <array>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "cell.h"
#include "particle_base.h"
#include "serialize.h"

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
  virtual ~Particle(){};

  //! Delete copy constructor
  Particle(const Particle<Tdim, Tnphases>&) = delete;

  //! Delete assignement operator
  Particle& operator=(const Particle<Tdim, Tnphases>&) = delete;

  //! Initialise properties
  void initialise();

  //! Compute reference coordinates in a cell
  void compute_reference_location();

  //! Return reference location
  VectorDim reference_location() const { return reference_location_; }

  // Assign a cell to particle
  //! \param[in] cellptr Pointer to a cell
  bool assign_cell(std::shared_ptr<Cell<Tdim>> cellptr);

  //! Return cell id
  Index cell_id() const { return cell_id_; }

  //! Assign nodal mass to particles
  //! \param[in] nphase Index corresponding to the phase
  //! \param[in] mass Mass from the particles in a cell
  //! \retval status Assignment status
  void assign_mass(unsigned nphase, double mass) { mass_(0, nphase) = mass; }

  //! Return mass of the particlesx
  //! \param[in] nphase Index corresponding to the phase
  double mass(unsigned nphase) const { return mass_(0, nphase); }

  //! Assign stress to the particle
  //! \param[in] nphase Index corresponding to the phase
  //! \param[in] stress A vector of particle stress
  void assign_stress(unsigned nphase,
                     const Eigen::Matrix<double, 6, 1>& stress);

  //! Return stress of the particle
  //! \param[in] nphase Index corresponding to the phase
  Eigen::Matrix<double, 6, 1> stress(unsigned nphase) const {
    return stress_.col(nphase);
  }

  //! Assign velocity to the particle
  //! \param[in] nphase Index corresponding to the phase
  //! \param[in] velocity A vector of particle velocity
  //! \retval status Assignment status
  bool assign_velocity(unsigned nphase, const Eigen::VectorXd& velocity);

  //! Return velocity of the particle
  //! \param[in] nphase Index corresponding to the phase
  Eigen::VectorXd velocity(unsigned nphase) const {
    return velocity_.col(nphase);
  }

  //! Assign momentum to the particle
  //! \param[in] nphase Index corresponding to the phase
  //! \param[in] momentum A vector of particle momentum
  //! \retval status Assignment status
  bool assign_momentum(unsigned nphase, const Eigen::VectorXd& momentum);

  //! Return momentum of the particle
  //! \param[in] nphase Index corresponding to the phase
  Eigen::VectorXd momentum(unsigned nphase) const {
    return momentum_.col(nphase);
  }

  //! Assign acceleration to the particle
  //! \param[in] nphase Index corresponding to the phase
  //! \param[in] acceleration A vector of particle acceleration
  //! \retval status Assignment status
  bool assign_acceleration(unsigned nphase,
                           const Eigen::VectorXd& acceleration);

  //! Return acceleration of the particle
  //! \param[in] nphase Index corresponding to the phase
  Eigen::VectorXd acceleration(unsigned nphase) const {
    return acceleration_.col(nphase);
  }

  friend class boost::serialization::access;
  //! Serialization / desierailization of particle
  //! \tparam Archive Boost archive object
  //! \param[in] ar Boost archive type
  //! \param[in] version Version numbering for the class object
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& id_;
    ar& coordinates_;
    std::cout << "Derived\n";
  }

 private:
  //! particle id
  using ParticleBase<Tdim>::id_;
  //! coordinates
  using ParticleBase<Tdim>::coordinates_;
  //! Reference coordinates (in a cell)
  Eigen::Matrix<double, Tdim, 1> reference_location_;
  //! Cell
  std::shared_ptr<Cell<Tdim>> cell_;
  //! Cell id
  using ParticleBase<Tdim>::cell_id_;
  //! Status
  using ParticleBase<Tdim>::status_;
  //! Mass
  Eigen::Matrix<double, 1, Tnphases> mass_;
  //! Stresses
  Eigen::Matrix<double, 6, Tnphases> stress_;
  //! Velocity
  Eigen::Matrix<double, Tdim, Tnphases> velocity_;
  //! Momentum
  Eigen::Matrix<double, Tdim, Tnphases> momentum_;
  //! Acceleration
  Eigen::Matrix<double, Tdim, Tnphases> acceleration_;

};  // Particle class
}  // namespace mpm

#include "particle.tcc"

#endif  // MPM_PARTICLE_H__

#ifndef MPM_PARTICLEBASE_H_
#define MPM_PARTICLEBASE_H_

#include <array>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "cell.h"

namespace mpm {

//! Global index type for the particleBase
using Index = unsigned long long;

//! ParticleBase class
//! \brief Base class that stores the information about particleBases
//! \details ParticleBase class: id_ and coordinates.
//! \tparam Tdim Dimension
template <unsigned Tdim>
class ParticleBase {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Constructor with id and coordinates
  //! \param[in] id Particle id
  //! \param[in] coord coordinates of the particle
  ParticleBase(Index id, const VectorDim& coord);

  //! Constructor with id, coordinates and status
  //! \param[in] id Particle id
  //! \param[in] coord coordinates of the particle
  //! \param[in] status Particle status (active / inactive)
  ParticleBase(Index id, const VectorDim& coord, bool status);

  //! Destructor
  virtual ~ParticleBase(){};

  //! Delete copy constructor
  ParticleBase(const ParticleBase<Tdim>&) = delete;

  //! Delete assignement operator
  ParticleBase& operator=(const ParticleBase<Tdim>&) = delete;

  //! Return id of the particleBase
  Index id() const { return id_; }

  //! Assign coordinates
  //! \param[in] coord Assign coord as coordinates of the particleBase
  void assign_coordinates(const VectorDim& coord) { coordinates_ = coord; }

  //! Return coordinates
  //! \retval coordinates_ return coordinates of the particleBase
  VectorDim coordinates() const { return coordinates_; }

  //! Compute reference coordinates in a cell
  virtual void compute_reference_location() = 0;

  //! Return reference location
  virtual VectorDim reference_location() const = 0;

  //! Assign cell
  virtual bool assign_cell(std::shared_ptr<Cell<Tdim>> cellptr) = 0;

  //! Return cell id
  virtual Index cell_id() const = 0;

  //! Compute shape functions
  virtual bool compute_shapefn() = 0;

  //! Assign status
  void assign_status(bool status) { status_ = status; }

  //! Status
  bool status() const { return status_; }

  //! Initialise properties
  virtual void initialise() = 0;

  //! Assign mass
  virtual void assign_mass(unsigned nphase, double mass) = 0;

  //! Return mass
  virtual double mass(unsigned nphase) const = 0;

  //! Assign stress
  virtual void assign_stress(unsigned nphase,
                             const Eigen::Matrix<double, 6, 1>& stress) = 0;

  //! Return stress
  virtual Eigen::Matrix<double, 6, 1> stress(unsigned nphase) const = 0;

  //! Assign velocity
  virtual bool assign_velocity(unsigned nphase,
                               const Eigen::VectorXd& velocity) = 0;

  //! Return velocity
  virtual Eigen::VectorXd velocity(unsigned nphase) const = 0;

  //! Assign momentum
  virtual bool assign_momentum(unsigned nphase,
                               const Eigen::VectorXd& momentum) = 0;

  //! Return momentum
  virtual Eigen::VectorXd momentum(unsigned nphase) const = 0;

  //! Assign acceleration
  virtual bool assign_acceleration(unsigned nphase,
                                   const Eigen::VectorXd& acceleration) = 0;

  //! Return acceleration
  virtual Eigen::VectorXd acceleration(unsigned nphase) const = 0;

 protected:
  //! particleBase id
  Index id_{std::numeric_limits<Index>::max()};
  //! coordinates
  VectorDim coordinates_;
  //! Cell id
  Index cell_id_{std::numeric_limits<Index>::max()};
  //! Status
  bool status_{true};

};  // ParticleBase class
}  // namespace mpm

#include "particle_base.tcc"

#endif  // MPM_PARTICLEBASE_H__

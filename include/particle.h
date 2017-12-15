#ifndef MPM_PARTICLE_H_
#define MPM_PARTICLE_H_

#include <array>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "cell.h"
#include "serialize.h"

namespace mpm {

//! Global index type for the particle
using Index = unsigned long long;

//! Particle class
//! \brief Base class that stores the information about particles
//! \details Particle class: id_ and coordinates.
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Particle {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Constructor with id and coordinates
  Particle(Index id, const VectorDim& coord);

  //! Constructor with id, coordinates and status
  Particle(Index id, const VectorDim& coord, bool status);

  //! Destructor
  virtual ~Particle(){};

  //! Delete copy constructor
  Particle(const Particle<Tdim>&) = delete;

  //! Delete assignement operator
  Particle& operator=(const Particle<Tdim>&) = delete;

  //! Return id of the particle
  Index id() const { return id_; }

  //! Assign coordinates
  //! \param[in] coord Assign coord as coordinates of the particle
  void coordinates(const VectorDim& coord) { coordinates_ = coord; }

  //! Return coordinates
  //! \retval coordinates_ return coordinates of the particle
  VectorDim coordinates() const { return coordinates_; }

  //! Assign cell
  bool assign_cell(const std::shared_ptr<Cell<Tdim>>& cellptr);

  //! Assign status
  void assign_status(bool status) { status_ = status; }

  //! Status
  bool status() const { return status_; }

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& id_;
    ar& coordinates_;
  }

 private:
  //! particle id
  Index id_{std::numeric_limits<Index>::max()};

  //! coordinates
  VectorDim coordinates_;

  //! Cell
  std::shared_ptr<Cell<Tdim>> cell_;

  //! Status
  bool status_{true};

  //! Serialize
  //! \tparam Archive Boost Archive
  //! \param[in] ar Archive
  //! \param[in] version Version of class
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) const {
    // note, version is always the latest when saving
    ar& id_;
    ar& coordinates_;
  }
};  // Particle class
}  // mpm namespace

#include "particle.tcc"

#endif  // MPM_PARTICLE_H__

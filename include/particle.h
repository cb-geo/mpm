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

  //! Constructor with id and coordinates
  Particle(Index id, const VectorDim& coord);

  //! Constructor with id, coordinates and status
  Particle(Index id, const VectorDim& coord, bool status);

  //! Destructor
  virtual ~Particle(){};

  //! Delete copy constructor
  Particle(const Particle<Tdim, Tnphases>&) = delete;

  //! Delete assignement operator
  Particle& operator=(const Particle<Tdim, Tnphases>&) = delete;

  //! Assign cell
  bool assign_cell(const std::shared_ptr<Cell<Tdim>>& cellptr);

  friend class boost::serialization::access;
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

  //! Cell
  std::shared_ptr<Cell<Tdim>> cell_;

  //! Status
  using ParticleBase<Tdim>::status_;

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
    std::cout << "Derived\n";
  }
};  // Particle class
}  // namespace mpm

#include "particle.tcc"

#endif  // MPM_PARTICLE_H__

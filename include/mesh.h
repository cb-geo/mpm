#ifndef MPM_MESH_H_
#define MPM_MESH_H_

#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "Eigen/Dense"

#include "cell.h"
#include "container.h"
#include "particle.h"

namespace mpm {

//! Mesh class
//! \brief Base class that stores the information about meshes
//! \details Mesh class: id_ and coordinates.
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Mesh {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Constructor with id
  Mesh(unsigned id);

  //! Destructor
  ~Mesh(){};

  //! Delete copy constructor
  Mesh(const Mesh<Tdim>&) = delete;

  //! Delete assignement operator
  Mesh& operator=(const Mesh<Tdim>&) = delete;

  //! Return id of the mesh
  unsigned id() const { return id_; }

  //! Add neighbouring mesh
  bool add_neighbour(unsigned id, const std::shared_ptr<Mesh<Tdim>>& neighbour);

  //! Number of neighbours
  unsigned nneighbours() const { return neighbour_meshes_.size(); }

  //! Add particle
  bool add_particle(const std::shared_ptr<mpm::Particle<Tdim>>& particle);

  //! Add particle
  bool remove_particle(const std::shared_ptr<mpm::Particle<Tdim>>& particle);

  //! Number of particles
  mpm::Index nparticles() const { return particles_.size(); }

  //! Active mesh (if a particle is present)
  bool status() const { return particles_.size(); }

 protected:
  //! mesh id
  unsigned id_{std::numeric_limits<unsigned>::max()};

  //! Container of mesh neighbours
  Handler<Mesh<Tdim>> neighbour_meshes_;

  //! Container of particles
  Container<Particle<Tdim>> particles_;
};  // Mesh class
}  // mpm namespace

#include "mesh.tcc"

#endif  // MPM_MESH_H_

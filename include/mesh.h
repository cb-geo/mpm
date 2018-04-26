#ifndef MPM_MESH_H_
#define MPM_MESH_H_

#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "Eigen/Dense"
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>

#include "cell.h"
#include "container.h"
#include "node.h"
#include "particle.h"
#include "particle_base.h"

namespace mpm {

//! Mesh class
//! \brief Base class that stores the information about meshes
//! \details Mesh class which stores the particles, nodes, cells and its
//! neighbours \tparam Tdim Dimension
template <unsigned Tdim>
class Mesh {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  // Construct a mesh with a global unique id
  //! \param[in] id Global mesh id
  Mesh(unsigned id);

  //! Default destructor
  ~Mesh() = default;

  //! Delete copy constructor
  Mesh(const Mesh<Tdim>&) = delete;

  //! Delete assignement operator
  Mesh& operator=(const Mesh<Tdim>&) = delete;

  //! Return id of the mesh
  unsigned id() const { return id_; }

  //! Add a neighbour mesh, using the local id for the new mesh and a mesh
  //! pointer
  //! \param[in] local_id local id of the mesh
  //! \param[in] neighbour A shared pointer to the neighbouring mesh
  //! \retval insertion_status Return the successful addition of a node
  bool add_neighbour(unsigned local_id,
                     const std::shared_ptr<Mesh<Tdim>>& neighbour);

  //! Return the number of neighbouring meshes
  unsigned nneighbours() const { return neighbour_meshes_.size(); }

  //! Add a particle to the mesh
  //! \param[in] particle A shared pointer to particle
  //! \retval insertion_status Return the successful addition of a particle
  bool add_particle(const std::shared_ptr<mpm::ParticleBase<Tdim>>& particle);

  //! Remove a particle from the mesh
  //! \param[in] particle A shared pointer to particle
  //! \retval insertion_status Return the successful addition of a particle
  bool remove_particle(
      const std::shared_ptr<mpm::ParticleBase<Tdim>>& particle);

  //! Number of particles in the mesh
  mpm::Index nparticles() const { return particles_.size(); }

  //! Add a node to the mesh
  //! \param[in] node A shared pointer to node
  //! \retval insertion_status Return the successful addition of a node
  bool add_node(const std::shared_ptr<mpm::NodeBase<Tdim>>& node);

  //! Remove a node from the mesh
  //! \param[in] node A shared pointer to node
  //! \retval insertion_status Return the successful addition of a node
  bool remove_node(const std::shared_ptr<mpm::NodeBase<Tdim>>& node);

  //! Return the number of nodes
  mpm::Index nnodes() const { return nodes_.size(); }

  //! Iterate over nodes
  template <typename Toper>
  Toper iterate_over_nodes(Toper oper);

  //! Add a cell from the mesh
  //! \param[in] cell A shared pointer to cell
  //! \retval insertion_status Return the successful addition of a cell
  bool add_cell(const std::shared_ptr<mpm::Cell<Tdim>>& cell);

  //! Remove a cell from the mesh
  //! \param[in] cell A shared pointer to cell
  //! \retval insertion_status Return the successful addition of a cell
  bool remove_cell(const std::shared_ptr<mpm::Cell<Tdim>>& cell);

  //! Number of cells in the mesh
  mpm::Index ncells() const { return cells_.size(); }

  //! Return status of the mesh. A mesh is active, if at least one particle is
  //! present
  bool status() const { return particles_.size(); }

 protected:
  //! mesh id
  unsigned id_{std::numeric_limits<unsigned>::max()};

  //! Container of mesh neighbours
  Handler<Mesh<Tdim>> neighbour_meshes_;

  //! Container of particles
  Container<ParticleBase<Tdim>> particles_;

  //! Container of nodes
  Container<NodeBase<Tdim>> nodes_;

  //! Container of cells
  Container<Cell<Tdim>> cells_;
};  // Mesh class
}  // namespace mpm

#include "mesh.tcc"

#endif  // MPM_MESH_H_

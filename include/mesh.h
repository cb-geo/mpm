#ifndef MPM_MESH_H_
#define MPM_MESH_H_

#include <algorithm>
#include <array>
#include <limits>
#include <memory>
#include <numeric>
#include <vector>

// Eigen
#include "Eigen/Dense"
// MPI
#ifdef USE_MPI
#include "mpi.h"
#endif
// TBB
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>

#include "cell.h"
#include "container.h"
#include "factory.h"
#include "hdf5.h"
#include "logger.h"
#include "material/material.h"
#include "node.h"
#include "particle.h"
#include "particle_base.h"

#include "map.h"

namespace mpm {

//! Mesh class
//! \brief Base class that stores the information about meshes
//! \details Mesh class which stores the particles, nodes, cells and neighbours
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Mesh {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  // Construct a mesh with a global unique id
  //! \param[in] id Global mesh id
  //! \param[in] isoparametric Mesh is isoparametric
  Mesh(unsigned id, bool isoparametric = true);

  //! Default destructor
  ~Mesh() = default;

  //! Delete copy constructor
  Mesh(const Mesh<Tdim>&) = delete;

  //! Delete assignement operator
  Mesh& operator=(const Mesh<Tdim>&) = delete;

  //! Return id of the mesh
  unsigned id() const { return id_; }

  //! Return if a mesh is isoparametric
  bool is_isoparametric() const { return isoparametric_; }

  //! Create nodes from coordinates
  //! \param[in] gnid Global node id
  //! \param[in] node_type Node type
  //! \param[in] coordinates Nodal coordinates
  //! \param[in] check_duplicates Parameter to check duplicates
  //! \retval status Create node status
  bool create_nodes(mpm::Index gnid, const std::string& node_type,
                    const std::vector<VectorDim>& coordinates,
                    bool check_duplicates = true);

  //! Add a node to the mesh
  //! \param[in] node A shared pointer to node
  //! \param[in] check_duplicates Parameter to check duplicates
  //! \retval insertion_status Return the successful addition of a node
  bool add_node(const std::shared_ptr<mpm::NodeBase<Tdim>>& node,
                bool check_duplicates = true);

  //! Remove a node from the mesh
  //! \param[in] node A shared pointer to node
  //! \retval insertion_status Return the successful addition of a node
  bool remove_node(const std::shared_ptr<mpm::NodeBase<Tdim>>& node);

  //! Return the number of nodes
  mpm::Index nnodes() const { return nodes_.size(); }

  //! Iterate over nodes
  //! \tparam Toper Callable object typically a baseclass functor
  template <typename Toper>
  void iterate_over_nodes(Toper oper);

  //! Iterate over nodes with predicate
  //! \tparam Toper Callable object typically a baseclass functor
  //! \tparam Tpred Predicate
  template <typename Toper, typename Tpred>
  void iterate_over_nodes_predicate(Toper oper, Tpred pred);

  //! Create a list of active nodes in mesh
  void find_active_nodes();

  //! Iterate over active nodes
  //! \tparam Toper Callable object typically a baseclass functor
  template <typename Toper>
  void iterate_over_active_nodes(Toper oper);

#ifdef USE_MPI
  //! All reduce over nodal scalar property
  //! \tparam Tgetfunctor Functor for getter
  //! \tparam Tsetfunctor Functor for setter
  //! \param[in] getter Getter function
  template <typename Tgetfunctor, typename Tsetfunctor>
  void allreduce_nodal_scalar_property(Tgetfunctor getter, Tsetfunctor setter);
#endif

#ifdef USE_MPI
  //! All reduce over nodal vector property
  //! \tparam Tgetfunctor Functor for getter
  //! \tparam Tsetfunctor Functor for setter
  //! \param[in] getter Getter function
  template <typename Tgetfunctor, typename Tsetfunctor>
  void allreduce_nodal_vector_property(Tgetfunctor getter, Tsetfunctor setter);
#endif

  //! Create cells from list of nodes
  //! \param[in] gcid Global cell id
  //! \param[in] element Element type
  //! \param[in] cells Node ids of cells
  //! \param[in] check_duplicates Parameter to check duplicates
  //! \retval status Create cells status
  bool create_cells(mpm::Index gnid,
                    const std::shared_ptr<mpm::Element<Tdim>>& element,
                    const std::vector<std::vector<mpm::Index>>& cells,
                    bool check_duplicates = true);

  //! Add a cell from the mesh
  //! \param[in] cell A shared pointer to cell
  //! \param[in] check_duplicates Parameter to check duplicates
  //! \retval insertion_status Return the successful addition of a cell
  bool add_cell(const std::shared_ptr<mpm::Cell<Tdim>>& cell,
                bool check_duplicates = true);

  //! Remove a cell from the mesh
  //! \param[in] cell A shared pointer to cell
  //! \retval insertion_status Return the successful addition of a cell
  bool remove_cell(const std::shared_ptr<mpm::Cell<Tdim>>& cell);

  //! Number of cells in the mesh
  mpm::Index ncells() const { return cells_.size(); }

  //! Iterate over cells
  //! \tparam Toper Callable object typically a baseclass functor
  template <typename Toper>
  void iterate_over_cells(Toper oper);

  //! Create particles from coordinates
  //! \param[in] gpids Global particle ids
  //! \param[in] particle_type Particle type
  //! \param[in] coordinates Nodal coordinates
  //! \param[in] check_duplicates Parameter to check duplicates
  //! \retval status Create particle status
  bool create_particles(const std::vector<mpm::Index>& gpids,
                        const std::string& particle_type,
                        const std::vector<VectorDim>& coordinates,
                        bool check_duplicates = true);

  //! Add a particle to the mesh
  //! \param[in] particle A shared pointer to particle
  //! \param[in] checks Parameter to check duplicates and addition
  //! \retval insertion_status Return the successful addition of a particle
  bool add_particle(const std::shared_ptr<mpm::ParticleBase<Tdim>>& particle,
                    bool checks = true);

  //! Remove a particle from the mesh
  //! \param[in] particle A shared pointer to particle
  //! \retval insertion_status Return the successful addition of a particle
  bool remove_particle(
      const std::shared_ptr<mpm::ParticleBase<Tdim>>& particle);

  //! Number of particles in the mesh
  mpm::Index nparticles() const { return particles_.size(); }

  //! Locate particles in a cell
  //! Iterate over all cells in a mesh to find the cell in which particles
  //! are located.
  //! \retval particles Particles which cannot be located in the mesh
  std::vector<std::shared_ptr<mpm::ParticleBase<Tdim>>> locate_particles_mesh();

  //! Iterate over particles
  //! \tparam Toper Callable object typically a baseclass functor
  template <typename Toper>
  void iterate_over_particles(Toper oper);

  //! Return coordinates of particles
  std::vector<Eigen::Matrix<double, 3, 1>> particle_coordinates();

  //! Return particles vector data
  //! \param[in] attribute Name of the vector data attribute
  //! \param[in] phase Index corresponding to the phase
  std::vector<Eigen::Matrix<double, 3, 1>> particles_vector_data(
      const std::string& attribute, unsigned phase);

  //! Assign velocity constraints to nodes
  //! \param[in] velocity_constraints Constraint at node, dir, and velocity
  bool assign_velocity_constraints(
      const std::vector<std::tuple<mpm::Index, unsigned, double>>&
          velocity_constraints);

  //! Assign velocity constraints to cells
  //! \param[in] velocity_constraints Constraint at cell id, face id, dir, and
  //! velocity
  bool assign_cell_velocity_constraints(
      const std::vector<std::tuple<mpm::Index, unsigned, unsigned, double>>&
          velocity_constraints);

  //! Assign particles volumes
  //! \param[in] particle_volumes Volume at dir on particle
  bool assign_particles_volumes(
      const std::vector<std::tuple<mpm::Index, double>>& particle_volumes);

  //! Assign particles tractions
  //! \param[in] particle_tractions Traction at dir on particle
  bool assign_particles_tractions(
      const std::vector<std::tuple<mpm::Index, unsigned, double>>&
          particle_tractions);

  //! Assign nodal traction force
  //! \param[in] nodal_tractions Traction at dir on nodes
  bool assign_nodal_tractions(
      const std::vector<std::tuple<mpm::Index, unsigned, double>>&
          nodal_tractions);

  //! Assign particles stresses
  //! \param[in] particle_stresses Initial stresses of particle
  bool assign_particles_stresses(
      const std::vector<Eigen::Matrix<double, 6, 1>>& particle_stresses);

  //! Assign particles cells
  //! \param[in] particles_cells Particles and cells
  bool assign_particles_cells(
      const std::vector<std::array<mpm::Index, 2>>& particles_cells);

  //! Return particles cells
  //! \retval particles_cells Particles and cells
  std::vector<std::array<mpm::Index, 2>> particles_cells() const;

  //! Return status of the mesh. A mesh is active, if at least one particle is
  //! present
  bool status() const { return particles_.size(); }

  //! Generate points
  //! \param[in] nquadratures Number of points per direction in cell
  //! \retval point Material point coordinates
  std::vector<VectorDim> generate_material_points(unsigned nquadratures = 1);

  //! Add a neighbour mesh, using the local id for the new mesh and a mesh
  //! pointer
  //! \param[in] local_id local id of the mesh
  //! \param[in] neighbour A shared pointer to the neighbouring mesh
  //! \retval insertion_status Return the successful addition of a node
  bool add_neighbour(unsigned local_id,
                     const std::shared_ptr<Mesh<Tdim>>& neighbour);

  //! Return the number of neighbouring meshes
  unsigned nneighbours() const { return neighbour_meshes_.size(); }

  //! Write HDF5 particles
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] filename Name of HDF5 file to write particles data
  //! \retval status Status of writing HDF5 output
  bool write_particles_hdf5(unsigned phase, const std::string& filename);

  //! Read HDF5 particles
  //! \param[in] phase Index corresponding to the phase
  //! \param[in] filename Name of HDF5 file to write particles data
  //! \retval status Status of reading HDF5 output
  bool read_particles_hdf5(unsigned phase, const std::string& filename);

  //! Return nodal coordinates
  std::vector<Eigen::Matrix<double, 3, 1>> nodal_coordinates() const;

  //! Return node pairs
  std::vector<std::array<mpm::Index, 2>> node_pairs() const;

  //! Return map of particle sets
  std::map<mpm::Index, std::vector<mpm::Index>> particle_sets() const {
    return this->particle_sets_;
  };

  //! Create map of particle sets
  void create_particle_sets(
      const std::map<mpm::Index, std::vector<mpm::Index>>& particle_sets) {
    this->particle_sets_ = particle_sets;
  };

  //! Return map of particles for fast retrieval
  Map<ParticleBase<Tdim>> map_particles() { return this->map_particles_; };

 private:
  // Locate a particle in mesh cells
  bool locate_particle_cells(
      const std::shared_ptr<mpm::ParticleBase<Tdim>>& particle);

 private:
  //! mesh id
  unsigned id_{std::numeric_limits<unsigned>::max()};
  //! Isoparametric mesh
  bool isoparametric_{true};
  //! Container of mesh neighbours
  Map<Mesh<Tdim>> neighbour_meshes_;
  //! Container of particles
  Container<ParticleBase<Tdim>> particles_;
  //! Map of particles for fast retrieval
  Map<ParticleBase<Tdim>> map_particles_;
  //! Map of particle sets
  std::map<mpm::Index, std::vector<mpm::Index>> particle_sets_;
  //! Container of nodes
  Container<NodeBase<Tdim>> nodes_;
  //! Container of active nodes
  Container<NodeBase<Tdim>> active_nodes_;
  //! Map of node sets
  std::map<mpm::Index, std::vector<mpm::Index>> node_sets_;
  //! Map of nodes for fast retrieval
  Map<NodeBase<Tdim>> map_nodes_;
  //! Map of cells for fast retrieval
  Map<Cell<Tdim>> map_cells_;
  //! Container of cells
  Container<Cell<Tdim>> cells_;
  //! Logger
  std::unique_ptr<spdlog::logger> console_;
};  // Mesh class
}  // namespace mpm

#include "mesh.tcc"

#endif  // MPM_MESH_H_

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
#include <Eigen/Sparse>
// MPI
#ifdef USE_MPI
#include "mpi.h"
#endif
// OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif
// TSL Maps
#include <tsl/robin_map.h>
// JSON
#include "json.hpp"
using Json = nlohmann::json;

#include "cell.h"
#include "factory.h"
#include "friction_constraint.h"
#include "function_base.h"
#include "generators/injection.h"
#include "geometry.h"
#include "hdf5_particle.h"
#include "io.h"
#include "io_mesh.h"
#include "logger.h"
#include "material.h"
#include "mpi_datatypes.h"
#include "nodal_properties.h"
#include "node.h"
#include "particle.h"
#include "particle_base.h"
#include "radial_basis_function.h"
#include "traction.h"
#include "vector.h"
#include "velocity_constraint.h"

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

  //! Return container of nodes
  mpm::Vector<NodeBase<Tdim>> nodes() { return nodes_; }

  //! Return the number of nodes in rank
  mpm::Index nnodes_rank();

  //! Iterate over nodes
  //! \tparam Toper Callable object typically a baseclass functor
  template <typename Toper>
  void iterate_over_nodes(Toper oper);

  //! Iterate over nodes with predicate
  //! \tparam Toper Callable object typically a baseclass functor
  //! \tparam Tpred Predicate
  template <typename Toper, typename Tpred>
  void iterate_over_nodes_predicate(Toper oper, Tpred pred);

  //! Return a vector of nodes
  //! \param[in] set_id Set of id of nodes (-1 for all nodes)
  Vector<NodeBase<Tdim>> nodes(unsigned set_id) const {
    return (set_id == -1) ? this->nodes_ : node_sets_.at(set_id);
  }

  //! Return a nodal shared_ptr
  std::shared_ptr<NodeBase<Tdim>> node(unsigned node_id) {
    return map_nodes_[node_id];
  }

  //! Create a list of active nodes in mesh
  void find_active_nodes();

  //! Iterate over active nodes
  //! \tparam Toper Callable object typically a baseclass functor
  template <typename Toper>
  void iterate_over_active_nodes(Toper oper);

#ifdef USE_MPI
  //! All reduce over nodal property
  //! \tparam Ttype Type of property to accumulate
  //! \tparam Tnparam Size of individual property
  //! \tparam Tgetfunctor Functor for getter
  //! \tparam Tsetfunctor Functor for setter
  //! \param[in] getter Getter function
  template <typename Ttype, unsigned Tnparam, typename Tgetfunctor,
            typename Tsetfunctor>
  void nodal_halo_exchange(Tgetfunctor getter, Tsetfunctor setter);
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

  //! Number of cells in mesh rank
  mpm::Index ncells_rank(bool active_cells = false);

  //! Iterate over cells
  //! \tparam Toper Callable object typically a baseclass functor
  template <typename Toper>
  void iterate_over_cells(Toper oper);

  //! Find cell neighbours
  void find_cell_neighbours();

  //! Find global nparticles across MPI ranks / cell
  void find_nglobal_particles_cells();

  //! Create particles from coordinates
  //! \param[in] particle_type Particle type
  //! \param[in] coordinates Nodal coordinates
  //! \param[in] material_id ID of the material
  //! \param[in] pset_id Set ID of the particles
  //! \param[in] check_duplicates Parameter to check duplicates
  //! \retval status Create particle status
  bool create_particles(const std::string& particle_type,
                        const std::vector<VectorDim>& coordinates,
                        const std::vector<unsigned>& material_ids,
                        unsigned pset_id, bool check_duplicates = true);

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

  //! Remove a particle by id
  bool remove_particle_by_id(mpm::Index id);

  //! Remove a particle from the mesh
  //! \param[in] pids Vector of particle ids
  void remove_particles(const std::vector<mpm::Index>& pids);

  //! Remove all particles in a cell in nonlocal rank
  void remove_all_nonrank_particles();

  //! Transfer halo particles to different ranks
  void transfer_halo_particles();

  //! Transfer particles to different ranks in nonlocal rank cells
  //! \param[in] exchange_cells Vector of cell ids that needs exchange
  void transfer_nonrank_particles(
      const std::vector<mpm::Index>& exchange_cells);

  //! Find shared nodes across MPI domains in the mesh
  void find_domain_shared_nodes();

  //! Find number of domain shared nodes in local rank
  mpm::Index nshared_nodes() const { return domain_shared_nodes_.size(); }

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

  //! Iterate over particles with predicate
  //! \tparam Toper Callable object typically a baseclass functor
  //! \tparam Tpred Predicate
  template <typename Toper, typename Tpred>
  void iterate_over_particles_predicate(Toper oper, Tpred pred);

  //! Iterate over particle set
  //! \tparam Toper Callable object typically a baseclass functor
  //! \param[in] set_id particle set id
  template <typename Toper>
  void iterate_over_particle_set(int set_id, Toper oper);

  //! Return coordinates of particles
  std::vector<Eigen::Matrix<double, 3, 1>> particle_coordinates();

  //! Return particles tensor data
  //! \param[in] attribute Name of the tensor data attribute
  template <unsigned Tsize>
  std::vector<Eigen::Matrix<double, Tsize, 1>> particles_tensor_data(
      const std::string& attribute);

  //! Return particles scalar data
  //! \param[in] attribute Name of the scalar data attribute
  std::vector<double> particles_statevars_data(const std::string& attribute);

  //! Compute and assign rotation matrix to nodes
  //! \param[in] euler_angles Map of node number and respective euler_angles
  bool compute_nodal_rotation_matrices(
      const std::map<mpm::Index, Eigen::Matrix<double, Tdim, 1>>& euler_angles);

  //! Assign particles volumes
  //! \param[in] particle_volumes Volume at dir on particle
  bool assign_particles_volumes(
      const std::vector<std::tuple<mpm::Index, double>>& particle_volumes);

  //! Create particles tractions
  //! \param[in] mfunction Math function if defined
  //! \param[in] setid Particle set id
  //! \param[in] dir Direction of traction load
  //! \param[in] traction Particle traction
  bool create_particles_tractions(
      const std::shared_ptr<FunctionBase>& mfunction, int set_id, unsigned dir,
      double traction);

  //! Apply traction to particles
  //! \param[in] current_time Current time
  void apply_traction_on_particles(double current_time);

  //! Create particle velocity constraints tractions
  //! \param[in] setid Node set id
  bool create_particle_velocity_constraint(
      int set_id, const std::shared_ptr<mpm::VelocityConstraint>& constraint);

  //! Apply particles velocity constraints
  void apply_particle_velocity_constraints();

  //! Assign nodal concentrated force
  //! \param[in] nodal_forces Force at dir on nodes
  bool assign_nodal_concentrated_forces(
      const std::vector<std::tuple<mpm::Index, unsigned, double>>&
          nodal_forces);

  //! Assign nodal concentrated force
  //! \param[in] mfunction Math function if defined
  //! \param[in] setid Node set id
  //! \param[in] dir Direction of force
  //! \param[in] node_forces Concentrated force at dir on nodes
  bool assign_nodal_concentrated_forces(
      const std::shared_ptr<FunctionBase>& mfunction, int set_id, unsigned dir,
      double force);

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
  //! \param[in] particle_type Particle type
  //! \param[in] material_id ID of the material
  //! \param[in] cset_id Set ID of the cell [-1 for all cells]
  //! \param[in] pset_id Set ID of the particles
  //! \retval point Material point coordinates
  bool generate_material_points(unsigned nquadratures,
                                const std::string& particle_type,
                                const std::vector<unsigned>& material_ids,
                                int cset_id, unsigned pset_id);

  //! Initialise material models
  //! \param[in] materials Material models
  void initialise_material_models(
      const std::map<unsigned, std::shared_ptr<mpm::Material<Tdim>>>&
          materials) {
    materials_ = materials;
  }

  //! Find particle neighbours
  //! \param[in] cell of interest
  void find_particle_neighbours();

  //! Find particle neighbours
  //! \param[in] cell of interest
  void find_particle_neighbours(const std::shared_ptr<mpm::Cell<Tdim>>& cell);

  //! Add a neighbour mesh, using the local id for the new mesh and a mesh
  //! pointer
  //! \param[in] local_id local id of the mesh
  //! \param[in] neighbour A shared pointer to the neighbouring mesh
  //! \retval insertion_status Return the successful addition of a node
  bool add_neighbour(unsigned local_id,
                     const std::shared_ptr<Mesh<Tdim>>& neighbour);

  //! Return the number of neighbouring meshes
  unsigned nneighbours() const { return neighbour_meshes_.size(); }

  //! Find ghost boundary cells
  void find_ghost_boundary_cells();

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

  //! Return HDF5 particles
  //! \retval particles_hdf5 Vector of HDF5 particles
  std::vector<mpm::HDF5Particle> particles_hdf5() const;

  //! Return nodal coordinates
  std::vector<Eigen::Matrix<double, 3, 1>> nodal_coordinates() const;

  //! Return node pairs
  std::vector<std::array<mpm::Index, 2>> node_pairs(bool active = false) const;

  //! Create map of vector of particles in sets
  //! \param[in] map of particles ids in sets
  //! \param[in] check_duplicates Parameter to check duplicates
  //! \retval status Status of  create particle sets
  bool create_particle_sets(
      const tsl::robin_map<mpm::Index, std::vector<mpm::Index>>& particle_sets,
      bool check_duplicates);

  //! Create map of vector of nodes in sets
  //! \param[in] map of nodes ids in sets
  //! \param[in] check_duplicates Parameter to check duplicates
  //! \retval status Status of  create node sets
  bool create_node_sets(
      const tsl::robin_map<mpm::Index, std::vector<mpm::Index>>& node_sets,
      bool check_duplicates);

  //! Create map of vector of cells in sets
  //! \param[in] map of cells ids in sets
  //! \param[in] check_duplicates Parameter to check duplicates
  //! \retval status Status of  create cell sets
  bool create_cell_sets(
      const tsl::robin_map<mpm::Index, std::vector<mpm::Index>>& cell_sets,
      bool check_duplicates);

  //! Get the vector of cell
  mpm::Vector<Cell<Tdim>> cells();

  //! Return particle cell ids
  std::map<mpm::Index, mpm::Index>* particles_cell_ids();

  //! Return nghost cells
  unsigned nghost_cells() const { return ghost_cells_.size(); }

  //! Return nlocal ghost cells
  unsigned nlocal_ghost_cells() const { return local_ghost_cells_.size(); }

  //! Generate particles
  //! \param[in] io IO object handle
  //! \param[in] generator Point generator object
  bool generate_particles(const std::shared_ptr<mpm::IO>& io,
                          const Json& generator);

  //! Inject particles
  void inject_particles(double current_time);

  //! Compute free surface
  //! \param[in] method Type of method to use
  //! \param[in] volume_tolerance for volume_fraction approach
  //! \retval status Status of compute_free_surface
  bool compute_free_surface(
      std::string method,
      double volume_tolerance = std::numeric_limits<unsigned>::epsilon());

  //! Compute free surface by density method
  //! \details Using simple approach of volume fraction approach as (Kularathna
  //! & Soga, 2017) and density ratio comparison (Hamad, 2015). This method is
  //! fast, but less accurate.
  //! \param[in] volume_tolerance for volume_fraction approach
  //! \retval status Status of compute_free_surface
  bool compute_free_surface_by_density(
      double volume_tolerance = std::numeric_limits<unsigned>::epsilon());

  //! Compute free surface by geometry method
  //! \details Using a more expensive approach using neighbouring particles and
  //! current geometry. This method combine multiple checks in order to simplify
  //! and fasten the process: (1) Volume fraction approach as (Kularathna & Soga
  //! 2017), (2) Density comparison approach as (Hamad, 2015), and (3) Geometry
  //! based approach as (Marrone et al. 2010)
  //! \param[in] volume_tolerance for volume_fraction approach
  //! \retval status Status of compute_free_surface
  bool compute_free_surface_by_geometry(
      double volume_tolerance = std::numeric_limits<unsigned>::epsilon());

  //! Get free surface node set
  std::set<mpm::Index> free_surface_nodes();

  //! Get free surface cell set
  std::set<mpm::Index> free_surface_cells();

  //! Get free surface particle set
  std::set<mpm::Index> free_surface_particles();

  //! Create a list of active nodes in mesh and assign active node id
  //! (rank-wise)
  unsigned assign_active_nodes_id();

  //! Return container of active nodes
  mpm::Vector<NodeBase<Tdim>> active_nodes() { return active_nodes_; }

  //! Return global node indices
  std::vector<Eigen::VectorXi> global_node_indices() const;

  //! Compute correction force in the node
  bool compute_nodal_correction_force(
      const Eigen::SparseMatrix<double>& correction_matrix,
      const Eigen::VectorXd& pressure_increment, double dt);

  // Create the nodal properties' map
  void create_nodal_properties();

  // Initialise the nodal properties' map
  void initialise_nodal_properties();

 private:
  // Read particles from file
  //! \param[in] pset_id Set ID of the particles
  bool read_particles_file(const std::shared_ptr<mpm::IO>& io,
                           const Json& generator, unsigned pset_id);

  // Locate a particle in mesh cells
  bool locate_particle_cells(
      const std::shared_ptr<mpm::ParticleBase<Tdim>>& particle);

 private:
  //! mesh id
  unsigned id_{std::numeric_limits<unsigned>::max()};
  //! Isoparametric mesh
  bool isoparametric_{true};
  //! Vector of mesh neighbours
  Map<Mesh<Tdim>> neighbour_meshes_;
  //! Vector of particles
  Vector<ParticleBase<Tdim>> particles_;
  //! Vector of particles ids and cell ids
  std::map<mpm::Index, mpm::Index> particles_cell_ids_;
  //! Vector of particle sets
  tsl::robin_map<unsigned, std::vector<mpm::Index>> particle_sets_;
  //! Map of particles for fast retrieval
  Map<ParticleBase<Tdim>> map_particles_;
  //! Vector of nodes
  Vector<NodeBase<Tdim>> nodes_;
  //! Vector of domain shared nodes
  Vector<NodeBase<Tdim>> domain_shared_nodes_;
  //! Boundary nodes
  Vector<NodeBase<Tdim>> boundary_nodes_;
  //! Vector of node sets
  tsl::robin_map<unsigned, Vector<NodeBase<Tdim>>> node_sets_;
  //! Vector of active nodes
  Vector<NodeBase<Tdim>> active_nodes_;
  //! Map of nodes for fast retrieval
  Map<NodeBase<Tdim>> map_nodes_;
  //! Map of cells for fast retrieval
  Map<Cell<Tdim>> map_cells_;
  //! Vector of cells
  Vector<Cell<Tdim>> cells_;
  //! Vector of ghost cells sharing the current MPI rank
  Vector<Cell<Tdim>> ghost_cells_;
  //! Vector of local ghost cells
  Vector<Cell<Tdim>> local_ghost_cells_;
  //! Vector of cell sets
  tsl::robin_map<unsigned, Vector<Cell<Tdim>>> cell_sets_;
  //! Map of ghost cells to the neighbours ranks
  std::map<unsigned, std::vector<unsigned>> ghost_cells_neighbour_ranks_;
  //! Faces and cells
  std::multimap<std::vector<mpm::Index>, mpm::Index> faces_cells_;
  //! Materials
  std::map<unsigned, std::shared_ptr<mpm::Material<Tdim>>> materials_;
  //! Loading (Particle tractions)
  std::vector<std::shared_ptr<mpm::Traction>> particle_tractions_;
  //! Particle velocity constraints
  std::vector<std::shared_ptr<mpm::VelocityConstraint>>
      particle_velocity_constraints_;
  //! Vector of generators for particle injections
  std::vector<mpm::Injection> particle_injections_;
  //! Nodal property pool
  std::shared_ptr<mpm::NodalProperties> nodal_properties_{nullptr};
  //! Logger
  std::unique_ptr<spdlog::logger> console_;
  //! Maximum number of halo nodes
  unsigned nhalo_nodes_{0};
  //! Maximum number of halo nodes
  unsigned ncomms_{0};
};  // Mesh class
}  // namespace mpm

#include "mesh.tcc"

#endif  // MPM_MESH_H_

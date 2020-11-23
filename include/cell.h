#ifndef MPM_CELL_H_
#define MPM_CELL_H_

#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <vector>

#include "Eigen/Dense"
#include "Eigen/LU"

#include "affine_transform.h"
#include "element.h"
#include "geometry.h"
#include "logger.h"
#include "map.h"
#include "node_base.h"
#include "quadrature.h"

namespace mpm {

//! Cell class
//! \brief Base class that stores the information about cells
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Cell {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Define DOFs
  static const unsigned Tdof = (Tdim == 1) ? 1 : 3 * (Tdim - 1);

  //! Constructor with id, number of nodes and element
  //! \param[in] id Global cell id
  //! \param[in] nnodes Number of nodes per cell
  //! \param[in] elementptr Pointer to an element type
  //! \param[in] isoparametric Cell type is isoparametric
  Cell(Index id, unsigned nnodes,
       const std::shared_ptr<const Element<Tdim>>& elementptr,
       bool isoparametric = true);

  //! Default destructor
  ~Cell() = default;

  //! Delete copy constructor
  Cell(const Cell<Tdim>&) = delete;

  //! Delete assignement operator
  Cell& operator=(const Cell<Tdim>&) = delete;

  //! Return id of the cell
  Index id() const { return id_; }

  //! Initialise cell properties
  bool initialise();

  //! Return the initialisation status of cells
  //! \retval initialisation_status Cell has nodes, element type and volumes
  bool is_initialised() const;

  //! Assign quadrature
  void assign_quadrature(unsigned nquadratures);

  //! Generate points
  std::vector<Eigen::Matrix<double, Tdim, 1>> generate_points();

  //! Return the number of particles
  unsigned nparticles() const { return particles_.size(); }

  //! Assign global nparticles
  //! \param[in] nparticles Number of global particles of cell
  void nglobal_particles(unsigned nparticles) {
    nglobal_particles_ = nparticles;
  }

  //! nglobal particles
  //! \retval nglobal_particles_ Number of global particles of cell
  unsigned nglobal_particles() const { return nglobal_particles_; }

  //! Return the status of a cell: active (if a particle is present)
  bool status() const { return particles_.size(); }

  //! Return particles_
  std::vector<Index> particles() const { return particles_; }

  //! Number of nodes
  unsigned nnodes() const { return nodes_.size(); }

  //! Return nodes of the cell
  std::vector<std::shared_ptr<mpm::NodeBase<Tdim>>> nodes() const {
    return nodes_;
  }

  //! Return nodes id in a cell
  std::set<mpm::Index> nodes_id() const {
    std::set<mpm::Index> nodes_id_lists;
    for (const auto& node : nodes_) nodes_id_lists.insert(node->id());
    return nodes_id_lists;
  }

  //! Side node pair ids
  std::vector<std::array<mpm::Index, 2>> side_node_pairs() const;

  //! Activate nodes if particle is present
  void activate_nodes();

  //! Return a pointer to element type of a cell
  std::shared_ptr<const Element<Tdim>> element_ptr() { return element_; }

  //! Return the number of shape functions, returns zero if the element type is
  //! not set.
  unsigned nfunctions() const {
    return (this->element_ != nullptr ? this->element_->nfunctions() : 0);
  };

  //! Add a node pointer to cell
  //! \param[in] local_id local id of the node
  //! \param[in] node A shared pointer to the node
  //! \retval insertion_status Return the successful addition of a node
  bool add_node(unsigned local_id, const std::shared_ptr<NodeBase<Tdim>>& node);

  //! Add a neighbour cell
  //! \param[in] neighbour_id id of the neighbouring cell
  //! \retval insertion_status Return the successful addition of a node
  bool add_neighbour(mpm::Index neighbour_id);

  //! Number of neighbours
  unsigned nneighbours() const { return neighbours_.size(); }

  //! Return neighbour ids
  std::set<mpm::Index> neighbours() const { return neighbours_; }

  //! Add an id of a particle in the cell
  //! \param[in] id Global id of a particle
  //! \retval status Return the successful addition of a particle id
  bool add_particle_id(Index id);

  //! Remove a particle id from the cell (moved to a different cell / killed)
  //! \param[in] id Global id of a particle
  void remove_particle_id(Index id);

  //! Clear all particle ids in the cell
  void clear_particle_ids() { particles_.clear(); }

  //! Compute the volume of the cell
  void compute_volume();

  //! Return the volume of the cell
  double volume() const { return volume_; }

  //! Compute the centroid of the cell
  void compute_centroid();

  //! Return the centroid of the cell
  Eigen::Matrix<double, Tdim, 1> centroid() const { return centroid_; }

  //! Return the dN/dx at the centroid of the cell
  Eigen::MatrixXd dn_dx_centroid() const { return dn_dx_centroid_; }

  //! Compute mean length of cell
  void compute_mean_length();

  //! Return the mean_length
  double mean_length() const { return mean_length_; }

  //! Return nodal coordinates
  Eigen::MatrixXd nodal_coordinates() const { return nodal_coordinates_; }

  //! Check if a point is in a cartesian cell by checking the domain ranges
  //! \param[in] point Coordinates of point
  inline bool point_in_cartesian_cell(
      const Eigen::Matrix<double, Tdim, 1>& point);

  //! Check if a point is in a isoparametric cell
  //! Use an affine transformation and NR to check if a transformed point is in
  //! a unit cell. This is useful for points on the surface, where
  //! volume calculations are tricky. The transformed point should be between -1
  //! and 1 in a unit cell
  //! \param[in] point Coordinates of point
  //! \param[in|out] xi Local coordinates of point
  //! \retval status Return if a point is in cell or not
  inline bool is_point_in_cell(const Eigen::Matrix<double, Tdim, 1>& point,
                               Eigen::Matrix<double, Tdim, 1>* xi);

  //! Return the local coordinates of a point in a cell
  //! \param[in] point Coordinates of a point
  //! \retval xi Local coordinates of a point
  inline Eigen::Matrix<double, Tdim, 1> local_coordinates_point(
      const Eigen::Matrix<double, Tdim, 1>& point);

  //! Return the local coordinates of a point in a unit cell
  //! Using newton iteration / affine transformation
  //! \param[in] point Coordinates of a point
  //! \retval xi Local coordinates of a point
  inline Eigen::Matrix<double, Tdim, 1> transform_real_to_unit_cell(
      const Eigen::Matrix<double, Tdim, 1>& point);

  //! Assign MPI rank to nodes
  void assign_mpi_rank_to_nodes();

  //! Compute normal vector
  void compute_normals();

  //! Return sorted face node ids
  std::vector<std::vector<mpm::Index>> sorted_face_node_ids();

  //! Assign ranks
  //! \param[in] Rank of cell
  void rank(unsigned mpi_rank);

  //! Return rank
  unsigned rank() const;

  //! Return previous mpi rank
  unsigned previous_mpirank() const;

  /**
   * \defgroup MultiPhase Functions dealing with multi-phase MPM
   */
  /**@{*/

  //! Assign solving status
  //! \ingroup MultiPhase
  //! \param[in] status Cell solving status for parallel free-surface detection
  void assign_solving_status(bool status) { solving_status_ = status; }

  //! Return solving status of a cell: active (if a particle is present in all
  //! MPI rank)
  //! \ingroup MultiPhase
  bool solving_status() const { return solving_status_; }

  //! Assign free surface
  //! \ingroup MultiPhase
  //! \param[in] free_surface boolean indicating free surface cell
  void assign_free_surface(bool free_surface) { free_surface_ = free_surface; };

  //! Return free surface bool
  //! \ingroup MultiPhase
  //! \retval free_surface_ indicating free surface cell
  bool free_surface() const { return free_surface_; };

  //! Assign volume traction to node
  //! \ingroup MultiPhase
  //! \param[in] volume_fraction cell volume fraction
  void assign_volume_fraction(double volume_fraction) {
    volume_fraction_ = volume_fraction;
  };

  //! Return cell volume fraction
  //! \ingroup MultiPhase
  //! \retval volume_fraction_ cell volume fraction
  double volume_fraction() const { return volume_fraction_; };

  //! Map cell volume to the nodes
  //! \ingroup MultiPhase
  //! \param[in] phase to map volume
  void map_cell_volume_to_nodes(unsigned phase);

  //! Initialize local elemental matrices
  //! \ingroup MultiPhase
  bool initialise_element_matrix();

  //! Return local node indices
  //! \ingroup MultiPhase
  Eigen::VectorXi local_node_indices();

  //! Return local laplacian
  //! \ingroup MultiPhase
  const Eigen::MatrixXd& laplacian_matrix() { return laplacian_matrix_; };

  //! Compute local laplacian matrix (Used in poisson equation)
  //! \ingroup MultiPhase
  //! \param[in] grad_shapefn shape function gradient
  //! \param[in] pvolume volume weight
  //! \param[in] multiplier multiplier
  void compute_local_laplacian(const Eigen::MatrixXd& grad_shapefn,
                               double pvolume,
                               double multiplier = 1.0) noexcept;

  //! Return local laplacian RHS matrix
  //! \ingroup MultiPhase
  const Eigen::MatrixXd& poisson_right_matrix() {
    return poisson_right_matrix_;
  };

  //! Compute local poisson RHS matrix (Used in poisson equation)
  //! \ingroup MultiPhase
  //! \param[in] shapefn shape function
  //! \param[in] grad_shapefn shape function gradient
  //! \param[in] pvolume volume weight
  void compute_local_poisson_right(const Eigen::VectorXd& shapefn,
                                   const Eigen::MatrixXd& grad_shapefn,
                                   double pvolume,
                                   double multiplier = 1.0) noexcept;

  //! Return local correction matrix
  //! \ingroup MultiPhase
  const Eigen::MatrixXd& correction_matrix() { return correction_matrix_; };

  //! Compute local correction matrix (Used to correct velocity)
  //! \ingroup MultiPhase
  void compute_local_correction_matrix(const Eigen::VectorXd& shapefn,
                                       const Eigen::MatrixXd& grad_shapefn,
                                       double pvolume) noexcept;

  /**@}*/

 private:
  //! Approximately check if a point is in a cell
  //! \param[in] point Coordinates of point
  bool approx_point_in_cell(const Eigen::Matrix<double, Tdim, 1>& point);

 private:
  //! Mutex
  std::mutex cell_mutex_;
  //! cell id
  Index id_{std::numeric_limits<Index>::max()};
  //! MPI Rank
  unsigned rank_{0};
  //! Previous MPI Rank
  unsigned previous_mpirank_{0};
  //! Isoparametric
  bool isoparametric_{true};
  //! Number of nodes
  unsigned nnodes_{0};
  //! Volume
  double volume_{std::numeric_limits<double>::lowest()};
  //! Centroid
  VectorDim centroid_;
  //! mean_length of cell
  double mean_length_{std::numeric_limits<double>::max()};
  //! particles ids in cell
  std::vector<Index> particles_;
  //! Number of global nparticles
  unsigned nglobal_particles_{0};
  //! Container of node pointers (local id, node pointer)
  std::vector<std::shared_ptr<NodeBase<Tdim>>> nodes_;
  //! Nodal coordinates
  Eigen::MatrixXd nodal_coordinates_;
  //! Container of cell neighbour ids
  std::set<mpm::Index> neighbours_;
  //! Shape function
  std::shared_ptr<const Element<Tdim>> element_{nullptr};
  //! Quadrature
  std::shared_ptr<Quadrature<Tdim>> quadrature_{nullptr};
  //! dN/dx
  Eigen::MatrixXd dn_dx_centroid_;
  //! Velocity constraints
  //! key: face_id, value: pair of direction [0/1/2] and velocity value
  std::map<unsigned, std::vector<std::pair<unsigned, double>>>
      velocity_constraints_;
  //! Normal of face
  //! first-> face_id, second->vector of the normal
  std::map<unsigned, Eigen::VectorXd> face_normals_;

  /**
   * \defgroup MultiPhaseVariables Variables for multi-phase MPM
   * @{
   */
  //! Solving status
  bool solving_status_{false};
  //! Free surface bool
  bool free_surface_{false};
  //! Volume fraction
  double volume_fraction_{0.0};
  //! Local laplacian matrix
  Eigen::MatrixXd laplacian_matrix_;
  //! Local poisson RHS matrix
  Eigen::MatrixXd poisson_right_matrix_;
  //! Local correction RHS matrix
  Eigen::MatrixXd correction_matrix_;
  /**@}*/

  //! Logger
  std::unique_ptr<spdlog::logger> console_;
};  // Cell class
}  // namespace mpm

#include "cell.tcc"
#include "cell_multiphase.tcc"

#endif  // MPM_CELL_H_
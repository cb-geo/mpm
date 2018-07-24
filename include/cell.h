#ifndef MPM_CELL_H_
#define MPM_CELL_H_

#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "Eigen/Dense"
#include "Eigen/LU"

#include "affine_transform.h"
#include "handler.h"
#include "node_base.h"
#include "shapefn.h"

namespace mpm {

//! Global index type for the cell
using Index = unsigned long long;

//! Cell class
//! \brief Base class that stores the information about cells
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Cell {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Define DOF for stresses
  static const unsigned Tdof = (Tdim == 2) ? 3 : 6;

  //! Constructor with id, number of nodes and shapefn
  //! \param[in] id Global cell id
  //! \param[in] nnodes Number of nodes per cell
  //! \param[in] shapefnptr Pointer to a shape function
  Cell(Index id, unsigned nnodes,
       const std::shared_ptr<ShapeFn<Tdim>>& shapefnptr);

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
  //! \retval initialisation_status Cell has nodes, shape functions and volumes
  bool is_initialised() const;

  //! Return the number of particles
  unsigned nparticles() const { return particles_.size(); }

  //! Return the status of a cell: active (if a particle is present)
  bool status() const { return particles_.size(); }

  //! Number of nodes
  unsigned nnodes() const { return nodes_.size(); }

  //! Assign shape function
  //! \param[in] shapefnptr Pointer to a shape function
  bool shapefn(const std::shared_ptr<ShapeFn<Tdim>>& shapefnptr);

  //! Return a pointer to shape function of a cell
  std::shared_ptr<ShapeFn<Tdim>> shapefn_ptr() { return shapefn_; }

  //! Return the number of shape functions, returns zero if the shapefunction is
  //! not set.
  unsigned nfunctions() const {
    return (this->shapefn_ != nullptr ? this->shapefn_->nfunctions() : 0);
  };

  //! Add a node pointer to cell
  //! \param[in] local_id local id of the node
  //! \param[in] node A shared pointer to the node
  //! \retval insertion_status Return the successful addition of a node
  bool add_node(unsigned local_id, std::shared_ptr<NodeBase<Tdim>> node);

  //! Add a neighbour cell
  //! \param[in] local_id local id of the neighbouring cell
  //! \param[in] neighbour A shared pointer to the neighbouring cell
  //! \retval insertion_status Return the successful addition of a node
  bool add_neighbour(unsigned local_id,
                     const std::shared_ptr<Cell<Tdim>>& neighbour);

  //! Number of neighbours
  unsigned nneighbours() const { return neighbour_cells_.size(); }

  //! Add an id of a particle in the cell
  //! \param[in] id Global id of a particle
  //! \retval status Return the successful addition of a particle id
  bool add_particle_id(Index id);

  //! Remove a particle id from the cell (moved to a different cell / killed)
  //! \param[in] id Global id of a particle
  void remove_particle_id(Index id);

  //! Compute the volume of the cell
  void compute_volume();

  //! Return the volume of the cell
  double volume() const { return volume_; }

  //! Compute the centroid of the cell
  void compute_centroid();

  //! Return the centroid of the cell
  Eigen::Matrix<double, Tdim, 1> centroid() const { return centroid_; }

  //! Compute mean length of cell
  void compute_mean_length();

  //! Return the mean_length
  double mean_length() const { return mean_length_; }

  //! Check if a point is in a 2D cell
  //! Cell is broken into sub-triangles with point as one of the
  //! vertex The sum of the sub-volume should be equal to the volume of the cell
  //! for a point to be in the cell
  bool point_in_cell(const Eigen::Matrix<double, 2, 1>& point);

  //! Check if a point is in a 3D cell
  //! Cell is broken into sub-tetrahedrons with point as one of the
  //! vertex The sum of the sub-volume should be equal to the volume of the cell
  //! for a point to be in the cell
  bool point_in_cell(const Eigen::Matrix<double, 3, 1>& point);

  //! Return the local coordinates of a point in a 1D cell
  //! \param[in] point Coordinates of a point
  //! \retval xi Local coordinates of a point
  Eigen::Matrix<double, 1, 1> local_coordinates_point(
      const Eigen::Matrix<double, 1, 1>& point);

  //! Return the local coordinates of a point in a 2D cell
  //! \param[in] point Coordinates of a point
  //! \retval xi Local coordinates of a point
  Eigen::Matrix<double, 2, 1> local_coordinates_point(
      const Eigen::Matrix<double, 2, 1>& point);

  //! Return the local coordinates of a point in a 3D cell
  //! \param[in] point Coordinates of a point
  //! \retval xi Local coordinates of a point
  Eigen::Matrix<double, 3, 1> local_coordinates_point(
      const Eigen::Matrix<double, 3, 1>& point);

  //! Return the local coordinates of a point in a unit cell
  //! Using newton iteration
  //! \param[in] point Coordinates of a point
  //! \retval xi Local coordinates of a point
  Eigen::Matrix<double, 1, 1> transform_real_to_unit_cell(
      const Eigen::Matrix<double, 1, 1>& point);

  //! Return the local coordinates of a point in a unit cell
  //! Using newton iteration
  //! \param[in] point Coordinates of a point
  //! \retval xi Local coordinates of a point
  Eigen::Matrix<double, 2, 1> transform_real_to_unit_cell(
      const Eigen::Matrix<double, 2, 1>& point);

  //! Return the local coordinates of a point in a unit cell
  //! Using newton iteration
  //! \param[in] point Coordinates of a point
  //! \retval xi Local coordinates of a point
  Eigen::Matrix<double, 3, 1> transform_real_to_unit_cell(
      const Eigen::Matrix<double, 3, 1>& point);

  //! Map particle mass to nodes
  //! \param[in] shapefn Shapefns at local coordinates of particle
  //! \param[in] phase Phase associate to the particle
  //! \param[in] pmass mass of particle
  void map_particle_mass_to_nodes(const Eigen::VectorXd& xi, unsigned phase,
                                  double pmass);

  //! Map particle volume to nodes
  //! \param[in] xi local coordinates of particle
  //! \param[in] phase Phase associate to the particle
  //! \param[in] pvolume volume of particle
  void map_particle_volume_to_nodes(const VectorDim& xi, unsigned phase,
                                    double pvolume);

  //! Compute the nodal momentum with particle mass & velocity for a phase
  //! \param[in] xi local coordinates of particle
  //! \param[in] phase Phase associate to the particle
  //! \param[in] pmass mass of a particle
  //! \param[in] velocity velocity of a particle
  void compute_nodal_momentum(const VectorDim& xi, unsigned phase, double pmass,
                              const Eigen::VectorXd& pvelocity);

  //! Compute the nodal body force of a cell from particle mass and gravity
  //! \param[in] xi local coordinates of particle
  //! \param[in] phase Phase associate to the particle
  //! \param[in] pmass Mass of a particle
  //! \param[in] pgravity Gravity of a particle
  void compute_nodal_body_force(const VectorDim& xi, unsigned phase,
                                double pmass, const VectorDim& pgravity);

  //! Compute the noal internal force  of a cell from particle stress and volume
  //! \param[in] xi local coordinates of particle
  //! \param[in] phase Phase associate to the particle
  //! \param[in] pvolume Volume of particle
  //! \param[in] pstress Stress of particle
  void compute_nodal_internal_force(unsigned phase, double pvolume,
                                    const VectorDim& xi,
                                    const Eigen::Matrix<double, 6, 1>& pstress);

  //! Return velocity at given location by interpolating from nodes
  //! \param[in] xi local coordinates of particle
  //! \param[in] phase Phase associate to the particle
  //! \retval velocity Interpolated velocity at xi
  Eigen::VectorXd interpolate_nodal_velocity(const VectorDim& xi,
                                             unsigned phase);

  //! Return acceleration at given location by interpolating from nodes
  //! \param[in] xi local coordinates of particle
  //! \param[in] phase Phase associate to the particle
  //! \retval acceleration Interpolated acceleration at xi
  Eigen::VectorXd interpolate_nodal_acceleration(const VectorDim& xi,
                                                 unsigned phase);

 protected:
  //! cell id
  Index id_{std::numeric_limits<Index>::max()};

  //! Number of nodes
  unsigned nnodes_{0};

  //! Volume
  double volume_{std::numeric_limits<double>::max()};

  //! Centroid
  VectorDim centroid_;

  //! mean_length of cell
  double mean_length_{std::numeric_limits<double>::max()};

  //! particles ids in cell
  std::vector<Index> particles_;

  //! Container of node pointers (local id, node pointer)
  Handler<NodeBase<Tdim>> nodes_;

  //! Container of cell neighbours
  Handler<Cell<Tdim>> neighbour_cells_;

  //! Shape function
  std::shared_ptr<ShapeFn<Tdim>> shapefn_{nullptr};
};  // Cell class
}  // namespace mpm

#include "cell.tcc"

#endif  // MPM_CELL_H_

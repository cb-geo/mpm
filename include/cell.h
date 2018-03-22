#ifndef MPM_CELL_H_
#define MPM_CELL_H_

#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "Eigen/Dense"
#include "Eigen/LU"

#include "handler.h"
#include "node_base.h"
#include "shapefn.h"

namespace mpm {

//! Global index type for the cell
using Index = unsigned long long;

//! Cell class
//! \brief Base class that stores the information about cells
//! \details Cell class: id_ and coordinates.
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Cell {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Define DOF for stresses
  static const unsigned Tdof = (Tdim == 2) ? 3 : 6;

  //! Constructor with id and number of nodes
  Cell(Index id, unsigned nnodes);

  //! Constructor with id, number of nodes and shapefn
  Cell(Index id, unsigned nnodes,
       const std::shared_ptr<ShapeFn<Tdim>>& shapefnptr);

  //! Destructor
  ~Cell(){};

  //! Delete copy constructor
  Cell(const Cell<Tdim>&) = delete;

  //! Delete assignement operator
  Cell& operator=(const Cell<Tdim>&) = delete;

  //! Return id of the cell
  Index id() const { return id_; }

  //! Number of nodes
  unsigned nnodes() const { return nodes_.size(); }

  //! Assign shape function
  bool shapefn(const std::shared_ptr<ShapeFn<Tdim>>& shapefnptr);

  //! Number of shape functions
  unsigned nfunctions() const {
    return (this->shapefn_ != nullptr ? this->shapefn_->nfunctions() : 0);
  };

  //! Add node to cell
  bool add_node(unsigned local_id, const std::shared_ptr<NodeBase<Tdim>>& node);

  //! Add neighbouring cell
  bool add_neighbour(unsigned id, const std::shared_ptr<Cell<Tdim>>& neighbour);

  //! Number of neighbours
  unsigned nneighbours() const { return neighbour_cells_.size(); }

  //! Add particle id
  bool add_particle_id(Index id);

  //! Remove particle id
  void remove_particle_id(Index id);

  //! Active cell (if a particle is present)
  bool status() const { return particles_.size(); }

  //! Compute volume
  void compute_volume();

  //! Return volume
  double volume() const { return volume_; }

  //! Point in cell 2D
  bool point_in_cell(const Eigen::Matrix<double, 2, 1>& point);

  //! Point in cell 3D
  bool point_in_cell(const Eigen::Matrix<double, 3, 1>& point);

  //! Return the local coordinates of a point in a 2D cell
  Eigen::Matrix<double, 2, 1> local_coordinates_point(
      const Eigen::Matrix<double, 2, 1>& point);

  //! Return the local coordinates of a point in a 3D cell
  Eigen::Matrix<double, 3, 1> local_coordinates_point(
      const Eigen::Matrix<double, 3, 1>& point);

  //! Map particle mass to nodes
  void map_particle_mass_to_nodes(const VectorDim& xi, unsigned nphase,
                                  double pmass);

  //! Assign particle momentum to nodes
  void map_momentum_to_nodes(const VectorDim& xi, unsigned nphase, double pmass,
                             const Eigen::VectorXd& pvelocity);

  //! Map body force to nodes
  void map_body_force_to_nodes(const VectorDim& xi, unsigned nphase,
                               double pmass, const VectorDim& pgravity);

  //! Map internal force to nodes
  void map_internal_force_to_nodes(unsigned nphase, double pvolume,
                                   const VectorDim& xi,
                                   const Eigen::VectorXd& pstress);

  //! Return velocity at given location by interpolating from nodes
  Eigen::VectorXd interpolate_nodal_velocity(const VectorDim& xi,
                                             unsigned nphase);

  //! Return acceleration at given location by interpolating from nodes
  Eigen::VectorXd interpolate_nodal_acceleration(const VectorDim& xi,
                                                 unsigned nphase);

 protected:
  //! cell id
  Index id_{std::numeric_limits<Index>::max()};

  //! Number of nodes
  unsigned nnodes_{0};

  //! Volume
  double volume_{std::numeric_limits<double>::max()};
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

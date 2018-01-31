#ifndef MPM_CELL_H_
#define MPM_CELL_H_

#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "Eigen/Dense"

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

  //! Constructor with id and coordinates
  Cell(Index id, unsigned nnodes);

  //! Constructor with id, coordinates and shapefn
  Cell(Index id, unsigned nnodes,
       const std::shared_ptr<ShapeFn<Tdim>>& shapefnptr);

  //! Destructor
  virtual ~Cell(){};

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

  //! Add particle id
  void remove_particle_id(Index id);

  //! Active cell (if a particle is present)
  bool status() const { return particles_.size(); }

 protected:
  //! cell id
  Index id_{std::numeric_limits<Index>::max()};

  //! Number of nodes
  unsigned nnodes_{0};

  //! particles ids in cell
  std::vector<Index> particles_;

  //! Container of node pointers (local id, node pointer)
  Handler<NodeBase<Tdim>> nodes_;

  //! Container of cell neighbours
  Handler<Cell<Tdim>> neighbour_cells_;

  //! Shape function
  std::shared_ptr<ShapeFn<Tdim>> shapefn_{nullptr};
};  // Cell class
}  // mpm namespace

#include "cell.tcc"

#endif  // MPM_CELL_H_

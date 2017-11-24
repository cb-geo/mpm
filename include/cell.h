#ifndef MPM_CELL_H_
#define MPM_CELL_H_

#include <array>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "Eigen/Dense"

#include "handler.h"
#include "node.h"

namespace mpm {
  
// Global index type for the cell
using Index = unsigned long long;

// Cell class
//! \brief Base class that stores the information about cells
//! \details Cell class: id_ and coordinates.
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Cell {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;
  
  // Constructor with id and coordinates
  //! \param[in] id Global cell id
  //! \param[in] nnodes Number of nodes per cell
  Cell(Index id, unsigned nnodes) : id_{id}, nnodes_{nnodes} {
    // Check if the dimension is between 1 & 3
    static_assert((Tdim >= 1 && Tdim <= 3), "Invalid global dimension");
  };

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

  //! Add node to cell
  bool add_node(unsigned local_id, const std::shared_ptr<Node<Tdim>>& node);

  //! Add neighbouring cell
  void add_neighbour(const std::shared_ptr<Cell<Tdim>>& neighbour) {
    neighbour_cells_.emplace_back(neighbour);
  }

  //! Number of neighbours
  unsigned nneighbours() const { return neighbour_cells_.size(); }

 private:
  //! cell id
  Index id_ { std::numeric_limits<Index>::max() };

  //! Number of nodes
  unsigned nnodes_{0};

  //! Container of node pointers (local id, node pointer)
  Handler<Node<Tdim>> nodes_;

  //! Container of cell neighbours
  std::vector<std::weak_ptr<Cell<Tdim>>> neighbour_cells_;

}; // Cell class
} // mpm namespace

#include "cell.tcc"

#endif  // MPM_CELL_H_

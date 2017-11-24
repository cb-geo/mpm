#ifndef MPM_CELL_H_
#define MPM_CELL_H_

#include <array>
#include <iostream>
#include <limits>
#include <vector>

#include "Eigen/Dense"

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
  Cell(Index id)
      : id_{id} {
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

 protected:
  //! cell id
  Index id_;
}; // Cell class
} // mpm namespace
#endif  // MPM_CELL_H_

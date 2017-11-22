#ifndef MPM_NODE_H_
#define MPM_NODE_H_

#include <array>
#include <iostream>
#include <limits>
#include <vector>

#include "Eigen/Dense"

namespace mpm {
  
// Global index type for the node
using Index = unsigned long long;

// Node class
//! \brief Base class that stores the information about nodes
//! \details Node class: id_ and coordinates.
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Node {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;
  
  // Constructor with id and coordinates
  //! \param[in] id Node id
  //! \param[in] coord coordinates of the node
  Node(Index id, const VectorDim& coord)
      : id_{id} {
    // Check if the dimension is between 1 & 3
    static_assert((Tdim >= 1 && Tdim <= 3), "Invalid global dimension");
    coordinates_ = coord;
  };

  //! Destructor
  virtual ~Node(){};

  //! Return id of the node
  Index id() const { return id_; }

  //! Assign coordinates
  //! \param[in] coord Assign coord as coordinates of the node
  void coordinates(const VectorDim& coord) {
    coordinates_ = coord;
  }

  //! Return coordinates
  //! \param[out] coordinates_ return coordinates of the node
  VectorDim coordinates() const { return coordinates_; }

 private:
  //! Copy constructor
  Node(const Node<Tdim>&);

  //! Assignement operator
  Node& operator=(const Node<Tdim>&);

 protected:
  //! node id
  Index id_;

  //! nodal coordinates
  VectorDim coordinates_;
}; // Node class
} // mpm namespace
#endif  // MPM_NODE_H_

#ifndef MPM_NODE_H_
#define MPM_NODE_H_

#include <array>
#include <iostream>
#include <limits>
#include <vector>

#include "Eigen/Dense"

#include "node_base.h"

namespace mpm {

// Global index type for the node
using Index = unsigned long long;

// Node class
//! \brief Base class that stores the information about nodes
//! \details Node class: id_ and coordinates.
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Node : public NodeBase<Tdim> {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  // Constructor with id and coordinates
  //! \param[in] id Node id
  //! \param[in] coord coordinates of the node
  Node(Index id, const VectorDim& coord) : NodeBase<Tdim>(id, coord){};

  //! Destructor
  virtual ~Node(){};

  //! Delete copy constructor
  Node(const Node<Tdim>&) = delete;

  //! Delete assignement operator
  Node& operator=(const Node<Tdim>&) = delete;

 protected:
  //! node id
  using NodeBase<Tdim>::id_;

  //! nodal coordinates
  using NodeBase<Tdim>::coordinates_;
};  // Node class
}  // mpm namespace
#endif  // MPM_NODE_H_

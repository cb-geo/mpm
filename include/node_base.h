#ifndef MPM_NODE_BASE_H_
#define MPM_NODE_BASE_H_

#include <array>
#include <iostream>
#include <limits>
#include <vector>

#include "Eigen/Dense"

namespace mpm {

// Global index type for the nodebase
using Index = unsigned long long;

// NodeBase class
//! \brief Base class that stores the information about nodes
//! \details NodeBase class: id_ and coordinates.
//! \tparam Tdim Dimension
template <unsigned Tdim>
class NodeBase {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  // Constructor with id and coordinates
  //! \param[in] id Node id
  //! \param[in] coord coordinates of the nodebase
  NodeBase(Index id, const VectorDim& coord) : id_{id} {
    // Check if the dimension is between 1 & 3
    static_assert((Tdim >= 1 && Tdim <= 3), "Invalid global dimension");
    coordinates_ = coord;
  };

  //! Destructor
  virtual ~NodeBase(){};

  //! Delete copy constructor
  NodeBase(const NodeBase<Tdim>&) = delete;

  //! Delete assignement operator
  NodeBase& operator=(const NodeBase<Tdim>&) = delete;

  //! Return id of the nodebase
  Index id() const { return id_; }

  //! Assign coordinates
  //! \param[in] coord Assign coord as coordinates of the nodebase
  void coordinates(const VectorDim& coord) { coordinates_ = coord; }

  //! Return coordinates
  //! \param[out] coordinates_ return coordinates of the nodebase
  VectorDim coordinates() const { return coordinates_; }

 protected:
  //! nodebase id
  Index id_{std::numeric_limits<Index>::max()};

  //! nodal coordinates
  VectorDim coordinates_;
};  // NodeBase class
}  // mpm namespace
#endif  // MPM_NODE_BASE_H_

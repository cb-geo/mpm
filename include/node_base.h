#ifndef MPM_NODE_BASE_H_
#define MPM_NODE_BASE_H_

#include <array>
#include <iostream>
#include <limits>
#include <vector>

namespace mpm {
  
// Global index type for the node
using Index = long long;

// Node Base class
//! \brief Base class that stores the information about nodes
//! \details NodeBase class: id_ and coordinates.
//! \tparam Tdim Dimension
template <unsigned Tdim>
class NodeBase {
 public:
  // Constructor with id and coordinates
  //! \param[in] id Node id
  //! \param[in] coord coordinates of the node
  NodeBase(const Index& id, const std::array<double, Tdim>& coord)
      : id_{id} {
    // Check if the dimension is between 1 & 3
    static_assert((Tdim >= 1 && Tdim <= 3), "Invalid global dimension");
    coordinates_ = coord;
  };

  //! Destructor
  virtual ~NodeBase(){};

  //! Return id of the node
  long long id() const { return id_; }

  //! Assign coordinates
  //! \param[in] coord Assign coord as coordinates of the node
  void coordinates(const std::array<double, Tdim>& coord) {
    coordinates_ = coord;
  }

  //! Return coordinates
  //! \param[out] coordinates_ return coordinates of the node
  std::array<double, Tdim> coordinates() const { return coordinates_; }

  //! Info
  void info() {
    std::cout << "Node id: " << id_ << ", coordinates: ";
    for (const auto& coord : coordinates_) std::cout << coord << ", ";
    std::cout << std::endl;
  }

 private:
  //! Copy constructor
  NodeBase(const NodeBase<Tdim>&);

  //! Assignement operator
  NodeBase& operator=(const NodeBase<Tdim>&);

 protected:
  //! node id
  Index id_;

  //! nodal coordinates
  std::array<double, Tdim> coordinates_;
}; // NodeBase class
} // mpm namespace
#endif  // MPM_NODE_BASE_H_

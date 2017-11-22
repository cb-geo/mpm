#ifndef MPM_NODE_HANDLER_H_
#define MPM_NODE_HANDLER_H_

#include <unordered_map>

#include "node.h"

namespace mpm {

// Global index type for the node
using Index = long long;

// Node handler class
//! \brief A class that offers a container and iterators for nodes
//! \tparam Tdim Dimension
template <unsigned Tdim>
class NodeHandler {
 public:
  //! Default constructor
  NodeHandler() = default;
  
  //! Insert node
  bool insert_node (const std::shared_ptr<Node<Tdim>>& node);

  //! Return number of nodes
  std::size_t nnodes() const { return nodes_.size(); }
  
  //! Return begin iterator of nodes
  typename std::unordered_map<Index, std::shared_ptr<Node<Tdim>>>::const_iterator
      nodes_begin() {
    return nodes_.begin();
  }

  //! Return end iterator of nodes
  typename std::unordered_map<Index, std::shared_ptr<Node<Tdim>>>::const_iterator
      nodes_end() {
    return nodes_.end();
  }

 private:
  // Unordered map of index and pointer to nodes
  std::unordered_map<Index, std::shared_ptr<Node<Tdim>>> nodes_;
};  // Node handler class

#include "node_handler.tcc"
  
}  // mpm namespace
#endif  // MPM_NODE_HANDLER_H_

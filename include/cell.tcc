//! Add a node pointer
//! \param[in] local_id local id of the node
//! \param[in] ptr A shared pointer
//! \retval insertion_status Return the successful addition of a node
//! \tparam Tdim Dimension
template<unsigned Tdim>
bool mpm::Cell<Tdim>::add_node(unsigned local_id, const std::shared_ptr<mpm::Node<Tdim>>& node_ptr) {
  bool insertion_status = false;
  try {
    // If number of node ptrs in a cell is less than the maximum number of nodes per cell
    // The local id should be between 0 and maximum number of nodes
    if (nodes_.size() < this->nnodes_ &&
        (local_id >= 0 && local_id < this->nnodes_)) {
      insertion_status = nodes_.insert(local_id, node_ptr);
    } else {
      throw std::runtime_error(
          "Number nodes in a cell exceeds the maximum allowed per cell");
    }
  } catch (std::exception& exception) {
    std::cerr << exception.what() << "\n";
  }
  return insertion_status;
}

//! Insert node
//! \param[in] node A shared pointer to a node
//! \tparam Tdim Dimension
template <unsigned Tdim>
bool mpm::NodeHandler<Tdim>::insert_node(const std::shared_ptr<Node<Tdim>>& node) {
  bool insertion_status = nodes_.insert(std::make_pair(node->id(), node)).second;
  return insertion_status;
}

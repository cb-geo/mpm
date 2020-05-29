//! Construct a semi-implicit eigen matrix assembler
template <unsigned Tdim>
mpm::AssemblerParallelSemiImplicitNavierStokes<
    Tdim>::AssemblerParallelSemiImplicitNavierStokes()
    : mpm::AssemblerEigenSemiImplicitNavierStokes<Tdim>() {
  //! Logger
  console_ =
      spdlog::stdout_color_mt("AssemblerParallelSemiImplicitNavierStokes");
}

//! Assign global node indices
template <unsigned Tdim>
bool mpm::AssemblerParallelSemiImplicitNavierStokes<
    Tdim>::assign_global_node_indices(unsigned nactive_node,
                                      unsigned nglobal_active_node) {
  bool status = true;
  try {
    // Total number of active node
    active_dof_ = nactive_node;
    global_active_dof_ = nglobal_active_node;

    global_node_indices_ = mesh_->global_node_indices();

    // Initialise mapping vector
    rank_global_mapper_.resize(active_dof_);

    // Nodes container
    const auto& nodes = mesh_->active_nodes();
    for (int counter = 0; counter < nodes.size(); counter++) {
      // Assign get nodal global index
      rank_global_mapper_[counter] = nodes[counter]->global_active_id();
    }

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

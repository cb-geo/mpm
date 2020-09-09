//! Constructor of contact with mesh
template <unsigned Tdim>
mpm::ContactFriction<Tdim>::ContactFriction(
    const std::shared_ptr<mpm::Mesh<Tdim>>& mesh)
    : mpm::Contact<Tdim>(mesh) {}

//! Initialize nodal properties
template <unsigned Tdim>
inline void mpm::ContactFriction<Tdim>::initialise() {
  // Initialise nodal properties
  mesh_->initialise_nodal_properties();

  // Append material ids to nodes
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::append_material_id_to_nodes,
                std::placeholders::_1));
}

//! Compute contact forces
template <unsigned Tdim>
inline void mpm::ContactFriction<Tdim>::compute_contact_forces(
    const Eigen::Matrix<double, Tdim, 1>& gravity, unsigned phase, double time,
    bool concentrated_nodal_forces) {

  // Map multimaterial properties from particles to nodes
  mesh_->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::map_multimaterial_mass_momentum_to_nodes,
      std::placeholders::_1));

  // Map multimaterial displacements from particles to nodes
  mesh_->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::map_multimaterial_displacements_to_nodes,
      std::placeholders::_1));

  // Map multimaterial domain gradients from particles to nodes
  mesh_->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Tdim>::map_multimaterial_domain_gradients_to_nodes,
      std::placeholders::_1));

  // Compute multimaterial change in momentum
  mesh_->iterate_over_nodes(
      std::bind(&mpm::NodeBase<Tdim>::compute_multimaterial_change_in_momentum,
                std::placeholders::_1));

  // Compute multimaterial separation vector
  mesh_->iterate_over_nodes(
      std::bind(&mpm::NodeBase<Tdim>::compute_multimaterial_separation_vector,
                std::placeholders::_1));

  // Compute multimaterial normal unit vector
  mesh_->iterate_over_nodes(
      std::bind(&mpm::NodeBase<Tdim>::compute_multimaterial_normal_unit_vector,
                std::placeholders::_1));

  // Map multimaterial body force
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::map_multimaterial_body_force,
                std::placeholders::_1, gravity));

  // Apply particle traction and map to multimaterial nodes
  mesh_->apply_multimaterial_traction_on_particles();

  // Iterate over each node to add concentrated node force to multimaterial
  // external force
  if (concentrated_nodal_forces)
    mesh_->iterate_over_nodes(
        std::bind(&mpm::NodeBase<Tdim>::apply_multimaterial_concentrated_force,
                  std::placeholders::_1, phase, time));

  // Iterate over each particle to compute multimaterial nodal internal force
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::map_multimaterial_internal_force,
                std::placeholders::_1));
}

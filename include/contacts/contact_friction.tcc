//! Constructor of contact with mesh
template <unsigned Tdim>
mpm::ContactFriction<Tdim>::ContactFriction(
    const std::shared_ptr<mpm::Mesh<Tdim>>& mesh, double friction,
    std::string normal_type)
    : mpm::Contact<Tdim>(mesh) {
  friction_ = friction;
  normal_type_ = normal_type;
}

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
inline void mpm::ContactFriction<Tdim>::compute_contact_forces() {

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
                std::placeholders::_1, normal_type_));
}

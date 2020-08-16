//! Initialise shared pointer to nodal properties pool for discontinuity
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::initialise_discontinuity_property_handle(
    unsigned prop_id,
    std::shared_ptr<mpm::NodalProperties> property_handle) noexcept {
  // the property handle and the property id is set in the node
  this->property_handle_ = property_handle;
  this->discontinuity_prop_id_ = prop_id;
}

//! Update nodal property at the nodes from particle for discontinuity
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::update_discontinuity_property(
    bool update, const std::string& property,
    const Eigen::MatrixXd& property_value, unsigned discontinuity_id,
    unsigned nprops) noexcept {
  // Update/assign property
  node_mutex_.lock();
  property_handle_->update_property(property, discontinuity_prop_id_,
                                    discontinuity_id, property_value, nprops);
  node_mutex_.unlock();
}

// Return data in the nodal properties map at a specific index
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
Eigen::MatrixXd mpm::Node<Tdim, Tdof, Tnphases>::discontinuity_property(
    const std::string& property, unsigned nprops) noexcept {
  // Const pointer to location of property: node_id * nprops x mat_id
  auto property_value =
      property_handle_->property(property, discontinuity_prop_id_, 0, nprops);
  ;
  // mpm::MapProperty property_handle(position, nprops);
  return property_value;
}

//! Compute acceleration and velocity with cundall damping factor
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
bool mpm::Node<Tdim, Tdof, Tnphases>::intergrate_momentum_discontinuity(
    unsigned phase, double dt) noexcept {
  momentum_.col(phase) =
      momentum_.col(phase) +
      (internal_force_.col(phase) + external_force_.col(phase)) * dt;
  if (discontinuity_enrich_) {
    property_handle_->update_property(
        "momenta_enrich", discontinuity_prop_id_, 0,
        (property_handle_->property("internal_force_enrich",
                                    discontinuity_prop_id_, 0, Tdim) +
         property_handle_->property("external_force_enrich",
                                    discontinuity_prop_id_, 0, Tdim)) *
            dt,
        Tdim);
  }
  // Apply velocity constraints, which also sets acceleration to 0,
  // when velocity is set.
  this->apply_velocity_constraints_discontinuity();

  //need to be done
  Eigen::Matrix<double, 3, 1> normal{0.44721359474414313, 0,
                                    0.89442719147920724};
  property_handle_->assign_property("normal_unit_vectors_discontinuity",
                                    discontinuity_prop_id_, 0, normal, Tdim);

  this->self_contact_discontinuity(dt);

  this->apply_velocity_constraints_discontinuity();

  return true;
}
//! Apply velocity constraints
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof,
               Tnphases>::apply_velocity_constraints_discontinuity() {
  // Set velocity constraint
  for (const auto& constraint : this->velocity_constraints_) {
    // Direction value in the constraint (0, Dim * Nphases)
    const unsigned dir = constraint.first;
    // Direction: dir % Tdim (modulus)
    const auto direction = static_cast<unsigned>(dir % Tdim);
    // Phase: Integer value of division (dir / Tdim)
    const auto phase = static_cast<unsigned>(dir / Tdim);

    if (!generic_boundary_constraints_) {
      // Velocity constraints are applied on Cartesian boundaries
      // this->velocity_(direction, phase) = constraint.second;
      // need to do for one direction

      this->momentum_(direction, phase) = this->mass(phase) * constraint.second;
      property_handle_->assign_property(
          "momenta_enrich", discontinuity_prop_id_ * Tdim + direction, 0,
          property_handle_->property("mass_enrich", discontinuity_prop_id_, 0,
                                     1) *
              constraint.second,
          1);
      // Set acceleration to 0 in direction of velocity constraint
      // this->acceleration_(direction, phase) = 0.;
      this->internal_force_(direction, phase) = 0;
      this->external_force_(direction, phase) = 0;

      Eigen::Matrix<double, 1, 1> momentum;
      momentum.setZero();
      property_handle_->assign_property(
          "internal_force_enrich", discontinuity_prop_id_ * Tdim + direction, 0,
          momentum, 1);
      property_handle_->assign_property(
          "external_force_enrich", discontinuity_prop_id_ * Tdim + direction, 0,
          momentum, 1);
    } else {  // need to do
      // Velocity constraints on general boundaries
      // Compute inverse rotation matrix
      const Eigen::Matrix<double, Tdim, Tdim> inverse_rotation_matrix =
          rotation_matrix_.inverse();
      // Transform to local coordinate
      Eigen::Matrix<double, Tdim, Tnphases> local_velocity =
          inverse_rotation_matrix * this->velocity_;
      Eigen::Matrix<double, Tdim, Tnphases> local_acceleration =
          inverse_rotation_matrix * this->acceleration_;
      // Apply boundary condition in local coordinate
      local_velocity(direction, phase) = constraint.second;
      local_acceleration(direction, phase) = 0.;
      // Transform back to global coordinate
      this->velocity_ = rotation_matrix_ * local_velocity;
      this->acceleration_ = rotation_matrix_ * local_acceleration;
    }
  }
}
//! Apply velocity constraints
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::self_contact_discontinuity(
    double dt) noexcept {

  if (!discontinuity_enrich_) return;

  unsigned phase = 0;
  const double tolerance = 1.0E-15;

  Eigen::Matrix<double, 1, 1> mass_enrich =
      property_handle_->property("mass_enrich", discontinuity_prop_id_, 0, 1);
  Eigen::Matrix<double, Tdim, 1> momenta_enrich = property_handle_->property(
      "momenta_enrich", discontinuity_prop_id_, 0, Tdim);
  Eigen::Matrix<double, Tdim, 1> normal_vector = property_handle_->property(
      "normal_unit_vectors_discontinuity", discontinuity_prop_id_, 0, Tdim);

  // mass for different sides
  auto mass_positive = mass_.col(phase) + mass_enrich;
  auto mass_negative = mass_.col(phase) - mass_enrich;

  if (mass_positive(phase) < tolerance || mass_negative(phase) < tolerance)
    return;

  // velocity for different sides
  auto velocity_positive =
      (momentum_.col(phase) + momenta_enrich) / mass_positive(phase);
  auto velocity_negative =
      (momentum_.col(phase) - momenta_enrich) / mass_negative(phase);

  // relative normal velocity
  if ((velocity_positive - velocity_negative).col(phase).dot(normal_vector) >=
      0)
    return;

  // the contact momentum for sticking contact
  auto momentum_contact = (mass_enrich(phase) * momentum_.col(phase) -
                           mass_(phase) * momenta_enrich) /
                          mass_(phase);
  auto force_contact = momentum_contact / dt;

  //! friction_coef < 0: move together without slide
  // need to be done
  double friction_coef = 0;

  if (friction_coef < 0) {
    property_handle_->update_property("momenta_enrich", discontinuity_prop_id_,
                                      0, momentum_contact.col(phase), Tdim);
    property_handle_->update_property("external_force_enrich",
                                      discontinuity_prop_id_, 0,
                                      force_contact.col(phase), Tdim);
  } else {
    double momentum_contact_norm =
        momentum_contact.col(phase).dot(normal_vector);
    double force_contact_norm = momentum_contact_norm / dt;

    // the maximum frictional contact force
    double max_frictional_force = friction_coef * abs(force_contact_norm);

    auto momentum_tangential =
        momentum_contact.col(phase) - momentum_contact_norm * normal_vector;
    auto force_tangential = momentum_tangential / dt;

    double force_tangential_value = force_tangential.norm();

    double frictional_force = force_tangential_value < max_frictional_force
                                  ? force_tangential_value
                                  : max_frictional_force;

    //! adjust the momentum and force
    property_handle_->update_property(
        "momenta_enrich", discontinuity_prop_id_, 0,
        momentum_contact_norm * normal_vector +
            frictional_force * force_tangential.col(phase).normalized() * dt,
        Tdim);
    property_handle_->update_property(
        "external_force_enrich", discontinuity_prop_id_, 0,
        force_contact_norm * normal_vector +
            frictional_force * force_tangential.col(phase).normalized(),
        Tdim);
  }
}

  //! Compute normal direction of each enrich node
  //! Apply velocity constraints
template <unsigned Tdim, unsigned Tdof, unsigned Tnphases>
void mpm::Node<Tdim, Tdof, Tnphases>::compute_normal_vector() noexcept {
  
    }
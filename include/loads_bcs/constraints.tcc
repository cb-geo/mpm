//! Assign nodal velocity constraints
template <unsigned Tdim>
bool mpm::Constraints<Tdim>::assign_nodal_velocity_constraint(
    int set_id, const std::shared_ptr<mpm::VelocityConstraint>& vconstraint) {
  bool status = true;
  try {
    int set_id = vconstraint->setid();
    auto nset = mesh_->nodes(set_id);
    if (nset.size() == 0)
      throw std::runtime_error(
          "Node set is empty for assignment of velocity constraints");

    unsigned dir = vconstraint->dir();
    double velocity = vconstraint->velocity();
    for (auto nitr = nset.cbegin(); nitr != nset.cend(); ++nitr) {
      if (!(*nitr)->assign_velocity_constraint(dir, velocity))
        throw std::runtime_error(
            "Failed to initialise velocity constraint at node");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assign velocity constraints to nodes
template <unsigned Tdim>
bool mpm::Constraints<Tdim>::assign_nodal_velocity_constraints(
    const std::vector<std::tuple<mpm::Index, unsigned, double>>&
        velocity_constraints) {
  bool status = true;
  try {
    for (const auto& velocity_constraint : velocity_constraints) {
      // Node id
      mpm::Index nid = std::get<0>(velocity_constraint);
      // Direction
      unsigned dir = std::get<1>(velocity_constraint);
      // Velocity
      double velocity = std::get<2>(velocity_constraint);

      // Apply constraint
      if (!mesh_->node(nid)->assign_velocity_constraint(dir, velocity))
        throw std::runtime_error(
            "Nodal velocity constraints assignment failed");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assign friction constraints to nodes
template <unsigned Tdim>
bool mpm::Constraints<Tdim>::assign_nodal_frictional_constraint(
    int nset_id, const std::shared_ptr<mpm::FrictionConstraint>& fconstraint) {
  bool status = true;
  try {
    int set_id = fconstraint->setid();
    auto nset = mesh_->nodes(set_id);
    if (nset.size() == 0)
      throw std::runtime_error(
          "Node set is empty for assignment of velocity constraints");
    unsigned dir = fconstraint->dir();
    int nsign_n = fconstraint->sign_n();
    double friction = fconstraint->friction();
    for (auto nitr = nset.cbegin(); nitr != nset.cend(); ++nitr) {
      if (!(*nitr)->assign_friction_constraint(dir, nsign_n, friction))
        throw std::runtime_error(
            "Failed to initialise velocity constraint at node");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assign friction constraints to nodes
template <unsigned Tdim>
bool mpm::Constraints<Tdim>::assign_nodal_friction_constraints(
    const std::vector<std::tuple<mpm::Index, unsigned, int, double>>&
        friction_constraints) {
  bool status = true;
  try {
    for (const auto& friction_constraint : friction_constraints) {
      // Node id
      mpm::Index nid = std::get<0>(friction_constraint);
      // Direction
      unsigned dir = std::get<1>(friction_constraint);
      // Sign
      int sign = std::get<2>(friction_constraint);
      // Friction
      double friction = std::get<3>(friction_constraint);

      // Apply constraint
      if (!mesh_->node(nid)->assign_friction_constraint(dir, sign, friction))
        throw std::runtime_error(
            "Nodal friction constraints assignment failed");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Apply absorbing constraints to nodes
template <unsigned Tdim>
bool mpm::Constraints<Tdim>::assign_nodal_absorbing_constraint(
    int nset_id,
    const std::shared_ptr<mpm::AbsorbingConstraint>& absorbing_constraint) {
  bool status = true;
  try {
    int set_id = absorbing_constraint->setid();
    auto nset = mesh_->nodes(set_id);
    if (nset.size() == 0)
      throw std::runtime_error(
          "Node set is empty for application of absorbing constraints");
    unsigned dir = absorbing_constraint->dir();
    double delta = absorbing_constraint->delta();
    double h_min = absorbing_constraint->h_min();
    double a = absorbing_constraint->a();
    double b = absorbing_constraint->b();
    assert(delta > h_min / (2 * a));
    assert(delta > h_min / (2 * b));
    for (auto nitr = nset.cbegin(); nitr != nset.cend(); ++nitr) {
      if (!(*nitr)->apply_absorbing_constraint(dir, delta, h_min, a, b))
        throw std::runtime_error(
            "Failed to apply absorbing constraint at node");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assign absorbing constraints to nodes
template <unsigned Tdim>
bool mpm::Constraints<Tdim>::assign_nodal_absorbing_constraints(
    const std::vector<std::tuple<mpm::Index, unsigned, double, double, double,
                                 double>>& absorbing_constraints) {
  bool status = true;
  try {
    for (const auto& absorbing_constraint : absorbing_constraints) {
      // Node id
      mpm::Index nid = std::get<0>(absorbing_constraint);
      // Direction
      unsigned dir = std::get<1>(absorbing_constraint);
      // Delta
      double delta = std::get<2>(absorbing_constraint);
      // h_min
      double h_min = std::get<3>(absorbing_constraint);
      // a
      double a = std::get<4>(absorbing_constraint);
      // b
      double b = std::get<5>(absorbing_constraint);
      assert(delta > h_min / (2 * a));
      assert(delta > h_min / (2 * b));
      // Apply constraint
      if (!mesh_->node(nid)->apply_absorbing_constraint(dir, delta, h_min, a,
                                                        b))
        throw std::runtime_error(
            "Nodal absorbing constraints assignment failed");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}
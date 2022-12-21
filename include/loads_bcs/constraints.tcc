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

//! Assign nodal acceleration constraints
template <unsigned Tdim>
bool mpm::Constraints<Tdim>::assign_nodal_acceleration_constraint(
    int set_id,
    const std::shared_ptr<mpm::AccelerationConstraint>& constraint) {
  bool status = true;
  try {
    int set_id = constraint->setid();
    auto nset = mesh_->nodes(set_id);
    if (nset.size() == 0)
      throw std::runtime_error(
          "Node set is empty for assignment of acceleration constraints");

    unsigned dir = constraint->dir();
    double acceleration = constraint->acceleration(0);
    for (auto nitr = nset.cbegin(); nitr != nset.cend(); ++nitr) {
      if (!(*nitr)->assign_acceleration_constraint(dir, acceleration))
        throw std::runtime_error(
            "Failed to initialise acceleration constraint at node");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assign acceleration constraints to nodes
template <unsigned Tdim>
bool mpm::Constraints<Tdim>::assign_nodal_acceleration_constraints(
    const std::vector<std::tuple<mpm::Index, unsigned, double>>&
        acceleration_constraints) {
  bool status = true;
  try {
    for (const auto& acceleration_constraint : acceleration_constraints) {
      // Node id
      mpm::Index nid = std::get<0>(acceleration_constraint);
      // Direction
      unsigned dir = std::get<1>(acceleration_constraint);
      // Acceleration
      double acceleration = std::get<2>(acceleration_constraint);

      // Apply constraint
      if (!mesh_->node(nid)->assign_acceleration_constraint(dir, acceleration))
        throw std::runtime_error(
            "Nodal acceleration constraints assignment failed");
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
          "Node set is empty for assignment of friction constraints");
    unsigned dir = fconstraint->dir();
    int nsign_n = fconstraint->sign_n();
    double friction = fconstraint->friction();
    for (auto nitr = nset.cbegin(); nitr != nset.cend(); ++nitr) {
      if (!(*nitr)->assign_friction_constraint(dir, nsign_n, friction))
        throw std::runtime_error(
            "Failed to initialise friction constraint at node");
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

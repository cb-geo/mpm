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

//! Assign nodal pressure constraints
template <unsigned Tdim>
bool mpm::Constraints<Tdim>::assign_nodal_pressure_constraint(
    const std::shared_ptr<FunctionBase>& mfunction,
    const std::shared_ptr<mpm::PressureConstraint>& pconstraint) {
  bool status = true;
  try {
    int set_id = pconstraint->setid();
    auto nset = mesh_->nodes(set_id);
    if (nset.size() != 0) {
      unsigned phase = pconstraint->phase();
      double pressure = pconstraint->pressure();
      tbb::parallel_for(
          tbb::blocked_range<int>(size_t(0), size_t(nset.size())),
          [&](const tbb::blocked_range<int>& range) {
            for (int i = range.begin(); i != range.end(); ++i) {
              status = nset[i]->assign_pressure_constraint(phase, pressure,
                                                           mfunction);
              if (!status)
                throw std::runtime_error(
                    "Failed to initialise pressure constraint at node");
            }
          },
          tbb::simple_partitioner());
    } else
      throw std::runtime_error(
          "Nodal assignment of pressure constraint failed");

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assign nodal pressure constraints to nodes
template <unsigned Tdim>
bool mpm::Constraints<Tdim>::assign_nodal_pressure_constraints(
    const unsigned phase,
    const std::vector<std::tuple<mpm::Index, double>>& pressure_constraints) {
  bool status = false;
  try {
    for (const auto& pressure_constraint : pressure_constraints) {
      // Node id
      mpm::Index nid = std::get<0>(pressure_constraint);
      // Pressure
      double pressure = std::get<1>(pressure_constraint);

      // Apply constraint
      status = mesh_->node(nid)->assign_pressure_constraint(phase, pressure,
                                                            nullptr);

      if (!status)
        throw std::runtime_error("Node or pressure constraint is invalid");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

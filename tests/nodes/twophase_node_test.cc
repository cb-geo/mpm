#include <cmath>
#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

#include "function_base.h"
#include "geometry.h"
#include "node.h"

// Check node class for 1D case
TEST_CASE("Twophase Node is checked for 1D case", "[node][1D][2P]") {
  const unsigned Dim = 1;
  const unsigned Dof = 1;
  const unsigned Nphases = 2;
  Eigen::Matrix<double, 1, 1> coords;
  coords.setZero();

  // Check for id = 0
  SECTION("Node id is zero") {
    mpm::Index id = 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);
    REQUIRE(node->id() == 0);
  }

  // Check for id is a positive value
  SECTION("Node id is positive") {
    mpm::Index id = std::numeric_limits<mpm::Index>::max();
    std::shared_ptr<mpm::NodeBase<Dim>> node =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);
    REQUIRE(node->id() == std::numeric_limits<mpm::Index>::max());
  }

  // Check for degrees of freedom
  SECTION("Check degrees of freedom") {
    mpm::Index id = 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);
    REQUIRE(node->dof() == 1);
  }

  // Check status
  SECTION("Check status") {
    mpm::Index id = 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);
    REQUIRE(node->status() == false);
    node->assign_status(true);
    REQUIRE(node->status() == true);
  }

  SECTION("Boundary ghost id") {
    mpm::Index id = 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);
    node->ghost_id(5);
    REQUIRE(node->ghost_id() == 5);
  }

  // Check MPI Rank
  SECTION("Check MPI Rank") {
    mpm::Index id = 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);
    REQUIRE(node->id() == 0);

    // Assign MPI ranks
    node->mpi_rank(0);
    node->mpi_rank(0);
    node->mpi_rank(1);

    std::set<unsigned> ranks = node->mpi_ranks();
    REQUIRE(ranks.size() == 2);
    std::vector<unsigned> mpi_ranks = {0, 1};
    unsigned i = 0;
    for (auto it = ranks.begin(); it != ranks.end(); ++it, ++i)
      REQUIRE(*it == mpi_ranks.at(i));
  }

  // Test coordinates function
  SECTION("coordinates function is checked") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;

    std::shared_ptr<mpm::NodeBase<Dim>> node =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);

    // Check for coordinates being zero
    auto coordinates = node->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));
    REQUIRE(coordinates.size() == Dim);

    // Check for negative value of coordinates
    for (unsigned i = 0; i < coordinates.size(); ++i)
      coords(i) = -1. * std::numeric_limits<double>::max();
    node->assign_coordinates(coords);
    coordinates = node->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    REQUIRE(coordinates.size() == Dim);

    // Check for positive value of coordinates
    for (unsigned i = 0; i < coordinates.size(); ++i)
      coords(i) = std::numeric_limits<double>::max();
    node->assign_coordinates(coords);
    coordinates = node->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    REQUIRE(coordinates.size() == Dim);
  }

  SECTION("Check nodal properties") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    std::shared_ptr<mpm::NodeBase<Dim>> node =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);

    // Check mass
    REQUIRE(node->mass(mpm::NodePhase::nSolid) ==
            Approx(0.0).epsilon(Tolerance));
    REQUIRE(node->mass(mpm::NodePhase::nLiquid) ==
            Approx(0.0).epsilon(Tolerance));
    double solid_mass = 100.5;
    double liquid_mass = 200.5;
    // Update mass to 100.5 and 200.5
    REQUIRE_NOTHROW(
        node->update_mass(true, mpm::NodePhase::nSolid, solid_mass));
    REQUIRE_NOTHROW(
        node->update_mass(true, mpm::NodePhase::nLiquid, liquid_mass));
    REQUIRE(node->mass(mpm::NodePhase::nSolid) ==
            Approx(100.5).epsilon(Tolerance));
    REQUIRE(node->mass(mpm::NodePhase::nLiquid) ==
            Approx(200.5).epsilon(Tolerance));
    // Update mass to 201 and 401
    REQUIRE_NOTHROW(
        node->update_mass(true, mpm::NodePhase::nSolid, solid_mass));
    REQUIRE_NOTHROW(
        node->update_mass(true, mpm::NodePhase::nLiquid, liquid_mass));
    REQUIRE(node->mass(mpm::NodePhase::nSolid) ==
            Approx(201.0).epsilon(Tolerance));
    REQUIRE(node->mass(mpm::NodePhase::nLiquid) ==
            Approx(401.0).epsilon(Tolerance));
    // Assign mass to 100 and 200
    solid_mass = 100.;
    liquid_mass = 200.;
    REQUIRE_NOTHROW(
        node->update_mass(false, mpm::NodePhase::nSolid, solid_mass));
    REQUIRE_NOTHROW(
        node->update_mass(false, mpm::NodePhase::nLiquid, liquid_mass));
    REQUIRE(node->mass(mpm::NodePhase::nSolid) ==
            Approx(100.0).epsilon(Tolerance));
    REQUIRE(node->mass(mpm::NodePhase::nLiquid) ==
            Approx(200.0).epsilon(Tolerance));

    SECTION("Check nodal pressure") {
      // Check pressure
      REQUIRE(node->pressure(mpm::NodePhase::nSolid) ==
              Approx(0.0).epsilon(Tolerance));
      REQUIRE(node->pressure(mpm::NodePhase::nLiquid) ==
              Approx(0.0).epsilon(Tolerance));
      double pressure = 1000.7;
      double pore_pressure = 2000.7;
      // Update pressure to 1000.7 and 2000.7
      REQUIRE_NOTHROW(node->update_mass_pressure(mpm::NodePhase::nSolid,
                                                 solid_mass * pressure));
      REQUIRE_NOTHROW(node->update_mass_pressure(mpm::NodePhase::nLiquid,
                                                 liquid_mass * pore_pressure));
      REQUIRE(node->pressure(mpm::NodePhase::nSolid) ==
              Approx(1000.7).epsilon(Tolerance));
      REQUIRE(node->pressure(mpm::NodePhase::nLiquid) ==
              Approx(2000.7).epsilon(Tolerance));
      // Update pressure to 2001.4 and 4001.4
      REQUIRE_NOTHROW(node->update_mass_pressure(mpm::NodePhase::nSolid,
                                                 solid_mass * pressure));
      REQUIRE_NOTHROW(node->update_mass_pressure(mpm::NodePhase::nLiquid,
                                                 liquid_mass * pore_pressure));
      REQUIRE(node->pressure(mpm::NodePhase::nSolid) ==
              Approx(2001.4).epsilon(Tolerance));
      REQUIRE(node->pressure(mpm::NodePhase::nLiquid) ==
              Approx(4001.4).epsilon(Tolerance));
      // Assign pressure to 1000 and 2000
      pressure = 1000.;
      pore_pressure = 2000.;
      node->assign_pressure(mpm::NodePhase::nSolid, pressure);
      node->assign_pressure(mpm::NodePhase::nLiquid, pore_pressure);
      REQUIRE(node->pressure(mpm::NodePhase::nSolid) ==
              Approx(1000.0).epsilon(Tolerance));
      REQUIRE(node->pressure(mpm::NodePhase::nLiquid) ==
              Approx(2000.0).epsilon(Tolerance));
      // Assign mass to 0
      solid_mass = 0.;
      liquid_mass = 0.;
      REQUIRE_NOTHROW(
          node->update_mass(false, mpm::NodePhase::nSolid, solid_mass));
      REQUIRE_NOTHROW(
          node->update_mass(false, mpm::NodePhase::nLiquid, liquid_mass));
      REQUIRE(node->mass(mpm::NodePhase::nSolid) ==
              Approx(0.0).epsilon(Tolerance));
      REQUIRE(node->mass(mpm::NodePhase::nLiquid) ==
              Approx(0.0).epsilon(Tolerance));
      // Try to update pressure to 2000, should throw and keep to 1000.
      node->assign_pressure(mpm::NodePhase::nSolid, pressure);
      node->assign_pressure(mpm::NodePhase::nLiquid, pore_pressure);
      REQUIRE(node->pressure(mpm::NodePhase::nSolid) ==
              Approx(1000.0).epsilon(Tolerance));
      REQUIRE(node->pressure(mpm::NodePhase::nLiquid) ==
              Approx(2000.0).epsilon(Tolerance));
    }

    SECTION("Check external force") {
      // Create a force vector
      Eigen::Matrix<double, Dim, 1> force;
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 10.;

      // Check current external force is zero
      for (unsigned i = 0; i < force.size(); ++i) {
        REQUIRE(node->external_force(mpm::NodePhase::nMixture)(i) ==
                Approx(0.).epsilon(Tolerance));
        REQUIRE(node->external_force(mpm::NodePhase::nLiquid)(i) ==
                Approx(0.).epsilon(Tolerance));
      }

      // Update force to 10.0
      REQUIRE_NOTHROW(
          node->update_external_force(true, mpm::NodePhase::nMixture, force));
      REQUIRE_NOTHROW(node->update_external_force(true, mpm::NodePhase::nLiquid,
                                                  0.5 * force));
      for (unsigned i = 0; i < force.size(); ++i) {
        REQUIRE(node->external_force(mpm::NodePhase::nMixture)(i) ==
                Approx(10.).epsilon(Tolerance));
        REQUIRE(node->external_force(mpm::NodePhase::nLiquid)(i) ==
                Approx(5.).epsilon(Tolerance));
      }

      // Update force to 20.0
      REQUIRE_NOTHROW(
          node->update_external_force(true, mpm::NodePhase::nMixture, force));
      REQUIRE_NOTHROW(node->update_external_force(true, mpm::NodePhase::nLiquid,
                                                  0.5 * force));
      for (unsigned i = 0; i < force.size(); ++i) {
        REQUIRE(node->external_force(mpm::NodePhase::nMixture)(i) ==
                Approx(20.).epsilon(Tolerance));
        REQUIRE(node->external_force(mpm::NodePhase::nLiquid)(i) ==
                Approx(10.).epsilon(Tolerance));
      }

      // Assign force as 10.0
      REQUIRE_NOTHROW(
          node->update_external_force(false, mpm::NodePhase::nMixture, force));
      REQUIRE_NOTHROW(node->update_external_force(
          false, mpm::NodePhase::nLiquid, 0.5 * force));
      for (unsigned i = 0; i < force.size(); ++i) {
        REQUIRE(node->external_force(mpm::NodePhase::nMixture)(i) ==
                Approx(10.).epsilon(Tolerance));
        REQUIRE(node->external_force(mpm::NodePhase::nLiquid)(i) ==
                Approx(5.).epsilon(Tolerance));
      }

      SECTION("Check concentrated force") {
        // Set external force to zero
        force.setZero();
        REQUIRE_NOTHROW(node->update_external_force(
            false, mpm::NodePhase::nMixture, force));

        // concentrated force
        std::shared_ptr<mpm::FunctionBase> ffunction = nullptr;
        double concentrated_force = 65.32;
        const unsigned Direction = 0;
        // Check external force
        for (unsigned i = 0; i < Dim; ++i)
          REQUIRE(node->external_force(mpm::NodePhase::nMixture)(i) ==
                  Approx(0.).epsilon(Tolerance));

        REQUIRE(node->assign_concentrated_force(mpm::NodePhase::nMixture,
                                                Direction, concentrated_force,
                                                ffunction) == true);

        double current_time = 0.0;
        node->apply_concentrated_force(mpm::NodePhase::nMixture, current_time);

        for (unsigned i = 0; i < Dim; ++i) {
          if (i == Direction)
            REQUIRE(node->external_force(mpm::NodePhase::nMixture)(i) ==
                    Approx(concentrated_force).epsilon(Tolerance));
          else
            REQUIRE(node->external_force(mpm::NodePhase::nMixture)(i) ==
                    Approx(0.).epsilon(Tolerance));
        }

        // Check for incorrect direction / phase
        const unsigned wrong_dir = 4;
        REQUIRE(node->assign_concentrated_force(mpm::NodePhase::nMixture,
                                                wrong_dir, concentrated_force,
                                                ffunction) == false);

        // Check again to ensure value hasn't been updated
        for (unsigned i = 0; i < Dim; ++i) {
          if (i == Direction)
            REQUIRE(node->external_force(mpm::NodePhase::nMixture)(i) ==
                    Approx(concentrated_force).epsilon(Tolerance));
          else
            REQUIRE(node->external_force(mpm::NodePhase::nMixture)(i) ==
                    Approx(0.).epsilon(Tolerance));
        }
      }
    }

    SECTION("Check internal force") {
      // Create a force vector
      Eigen::Matrix<double, Dim, 1> force;
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 10.;

      // Check current internal force is zero
      for (unsigned i = 0; i < force.size(); ++i) {
        REQUIRE(node->internal_force(mpm::NodePhase::nMixture)(i) ==
                Approx(0.).epsilon(Tolerance));
        REQUIRE(node->internal_force(mpm::NodePhase::nLiquid)(i) ==
                Approx(0.).epsilon(Tolerance));
      }

      // Update force to 10.0
      REQUIRE_NOTHROW(
          node->update_internal_force(true, mpm::NodePhase::nMixture, force));
      REQUIRE_NOTHROW(node->update_internal_force(true, mpm::NodePhase::nLiquid,
                                                  0.5 * force));
      for (unsigned i = 0; i < force.size(); ++i) {
        REQUIRE(node->internal_force(mpm::NodePhase::nMixture)(i) ==
                Approx(10.).epsilon(Tolerance));
        REQUIRE(node->internal_force(mpm::NodePhase::nLiquid)(i) ==
                Approx(5.).epsilon(Tolerance));
      }

      // Update force to 20.0
      REQUIRE_NOTHROW(
          node->update_internal_force(true, mpm::NodePhase::nMixture, force));
      REQUIRE_NOTHROW(node->update_internal_force(true, mpm::NodePhase::nLiquid,
                                                  0.5 * force));
      for (unsigned i = 0; i < force.size(); ++i) {
        REQUIRE(node->internal_force(mpm::NodePhase::nMixture)(i) ==
                Approx(20.).epsilon(Tolerance));
        REQUIRE(node->internal_force(mpm::NodePhase::nLiquid)(i) ==
                Approx(10.).epsilon(Tolerance));
      }

      // Assign force as 10.0
      REQUIRE_NOTHROW(
          node->update_internal_force(false, mpm::NodePhase::nMixture, force));
      REQUIRE_NOTHROW(node->update_internal_force(
          false, mpm::NodePhase::nLiquid, 0.5 * force));
      for (unsigned i = 0; i < force.size(); ++i) {
        REQUIRE(node->internal_force(mpm::NodePhase::nMixture)(i) ==
                Approx(10.).epsilon(Tolerance));
        REQUIRE(node->internal_force(mpm::NodePhase::nLiquid)(i) ==
                Approx(5.).epsilon(Tolerance));
      }
    }

    SECTION("Check drag force coefficient") {
      // Create a force vector
      Eigen::Matrix<double, Dim, 1> drag_force_coefficient;
      for (unsigned i = 0; i < drag_force_coefficient.size(); ++i)
        drag_force_coefficient(i) = 10.;

      // Check current drag force coefficient is zero
      for (unsigned i = 0; i < drag_force_coefficient.size(); ++i)
        REQUIRE(node->drag_force_coefficient()(i) ==
                Approx(0.).epsilon(Tolerance));

      // Update drag force coefficient to 10.0
      REQUIRE_NOTHROW(
          node->update_drag_force_coefficient(true, drag_force_coefficient));

      for (unsigned i = 0; i < drag_force_coefficient.size(); ++i)
        REQUIRE(node->drag_force_coefficient()(i) ==
                Approx(10.).epsilon(Tolerance));

      // Update force to 20.0
      REQUIRE_NOTHROW(
          node->update_drag_force_coefficient(true, drag_force_coefficient));

      for (unsigned i = 0; i < drag_force_coefficient.size(); ++i)
        REQUIRE(node->drag_force_coefficient()(i) ==
                Approx(20.).epsilon(Tolerance));

      // Assign force as 10.0
      REQUIRE_NOTHROW(
          node->update_drag_force_coefficient(false, drag_force_coefficient));
      for (unsigned i = 0; i < drag_force_coefficient.size(); ++i)
        REQUIRE(node->drag_force_coefficient()(i) ==
                Approx(10.).epsilon(Tolerance));
    }

    SECTION("Check compute acceleration and velocity") {
      // Time step
      const double dt = 0.1;

      // Nodal mass
      double solid_mass = 100.;
      double liquid_mass = 100.;
      // Update mass to 100.5
      REQUIRE_NOTHROW(
          node->update_mass(false, mpm::NodePhase::nSolid, solid_mass));
      REQUIRE(node->mass(mpm::NodePhase::nSolid) ==
              Approx(solid_mass).epsilon(Tolerance));

      // Check internal force
      // Create a force vector
      Eigen::Matrix<double, Dim, 1> force;
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 10. * i;
      // Update force to 10.0
      REQUIRE_NOTHROW(
          node->update_internal_force(false, mpm::NodePhase::nMixture, force));
      REQUIRE_NOTHROW(node->update_internal_force(
          false, mpm::NodePhase::nLiquid, 0.5 * force));
      // Internal force
      for (unsigned i = 0; i < force.size(); ++i) {
        REQUIRE(node->internal_force(mpm::NodePhase::nMixture)(i) ==
                Approx(force(i)).epsilon(Tolerance));
        REQUIRE(node->internal_force(mpm::NodePhase::nLiquid)(i) ==
                Approx(0.5 * force(i)).epsilon(Tolerance));
      }

      // External force
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 5. * i;
      // Update force to 10.0
      REQUIRE_NOTHROW(
          node->update_external_force(false, mpm::NodePhase::nMixture, force));
      REQUIRE_NOTHROW(node->update_external_force(
          false, mpm::NodePhase::nLiquid, 0.5 * force));
      for (unsigned i = 0; i < force.size(); ++i) {
        REQUIRE(node->external_force(mpm::NodePhase::nMixture)(i) ==
                Approx(force(i)).epsilon(Tolerance));
        REQUIRE(node->external_force(mpm::NodePhase::nLiquid)(i) ==
                Approx(0.5 * force(i)).epsilon(Tolerance));
      }

      // Drag force
      Eigen::Matrix<double, Dim, 1> drag_force_coefficient;
      for (unsigned i = 0; i < drag_force_coefficient.size(); ++i)
        drag_force_coefficient(i) = 5. * i;
      // Update force to 10.0
      REQUIRE_NOTHROW(
          node->update_drag_force_coefficient(false, drag_force_coefficient));
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->drag_force_coefficient()(i) ==
                Approx(drag_force_coefficient(i)).epsilon(Tolerance));

      REQUIRE(node->compute_acceleration_velocity_twophase_explicit(dt) ==
              true);

      // Check acceleration
      Eigen::Matrix<double, Dim, 1> liquid_acceleration;
      liquid_acceleration << 0.;
      Eigen::Matrix<double, Dim, 1> solid_acceleration;
      solid_acceleration << 0.;

      for (unsigned i = 0; i < solid_acceleration.size(); ++i) {
        REQUIRE(node->acceleration(mpm::NodePhase::nSolid)(i) ==
                Approx(solid_acceleration(i)).epsilon(Tolerance));
        REQUIRE(node->acceleration(mpm::NodePhase::nLiquid)(i) ==
                Approx(liquid_acceleration(i)).epsilon(Tolerance));
      }

      // Check velocity
      Eigen::Matrix<double, Dim, 1> solid_velocity = solid_acceleration * dt;
      Eigen::Matrix<double, Dim, 1> liquid_velocity = liquid_acceleration * dt;
      for (unsigned i = 0; i < solid_velocity.size(); ++i) {
        REQUIRE(node->velocity(mpm::NodePhase::nSolid)(i) ==
                Approx(solid_velocity(i)).epsilon(Tolerance));
        REQUIRE(node->velocity(mpm::NodePhase::nLiquid)(i) ==
                Approx(liquid_velocity(i)).epsilon(Tolerance));
      }

      // Apply friction constraints
      REQUIRE(node->assign_friction_constraint(0, 1., 0.5) == true);
      // Apply friction constraints
      REQUIRE(node->assign_friction_constraint(-1, 1., 0.5) == false);
      REQUIRE(node->assign_friction_constraint(3, 1., 0.5) == false);

      // Test acceleration with constraints
      solid_acceleration[0] = 0.5 * solid_acceleration[0];
      liquid_acceleration[0] = 0.5 * liquid_acceleration[0];
      for (unsigned i = 0; i < solid_acceleration.size(); ++i) {
        REQUIRE(node->acceleration(mpm::NodePhase::nSolid)(i) ==
                Approx(solid_acceleration(i)).epsilon(Tolerance));
        REQUIRE(node->acceleration(mpm::NodePhase::nLiquid)(i) ==
                Approx(liquid_acceleration(i)).epsilon(Tolerance));
      }

      // Apply cundall damping when calculating acceleration
      REQUIRE(node->compute_acceleration_velocity_twophase_explicit_cundall(
                  dt, 0.05) == true);

      // Test acceleration with cundall damping
      solid_acceleration[0] = 0.;
      liquid_acceleration[0] = 0.;
      for (unsigned i = 0; i < solid_acceleration.size(); ++i) {
        REQUIRE(node->acceleration(mpm::NodePhase::nSolid)(i) ==
                Approx(solid_acceleration(i)).epsilon(Tolerance));
        REQUIRE(node->acceleration(mpm::NodePhase::nLiquid)(i) ==
                Approx(liquid_acceleration(i)).epsilon(Tolerance));
      }

      // Apply velocity constraints
      REQUIRE(node->assign_velocity_constraint(0, 10.5) == true);
      REQUIRE(node->assign_velocity_constraint(1, 10.5) == true);
      REQUIRE(node->compute_acceleration_velocity_twophase_explicit(dt) ==
              true);

      // Test velocity with constraints
      solid_velocity[0] = 10.5;
      liquid_velocity[0] = 10.5;
      for (unsigned i = 0; i < solid_velocity.size(); ++i) {
        REQUIRE(node->velocity(mpm::NodePhase::nSolid)(i) ==
                Approx(solid_velocity(i)).epsilon(Tolerance));
        REQUIRE(node->velocity(mpm::NodePhase::nLiquid)(i) ==
                Approx(liquid_velocity(i)).epsilon(Tolerance));
      }

      // Test acceleration with constraints
      solid_acceleration[0] = 0.;
      liquid_acceleration[0] = 0.;
      for (unsigned i = 0; i < solid_acceleration.size(); ++i) {
        REQUIRE(node->acceleration(mpm::NodePhase::nSolid)(i) ==
                Approx(solid_acceleration(i)).epsilon(Tolerance));
        REQUIRE(node->acceleration(mpm::NodePhase::nLiquid)(i) ==
                Approx(liquid_acceleration(i)).epsilon(Tolerance));
      }

      // Exception check when mass is zero
      // Update mass to 0.
      REQUIRE_NOTHROW(node->update_mass(false, mpm::NodePhase::nSolid, 0.));
      REQUIRE_NOTHROW(node->update_mass(false, mpm::NodePhase::nLiquid, 0.));
      REQUIRE(node->mass(mpm::NodePhase::nSolid) ==
              Approx(0.).epsilon(Tolerance));
      REQUIRE(node->mass(mpm::NodePhase::nLiquid) ==
              Approx(0.).epsilon(Tolerance));
      REQUIRE(node->compute_acceleration_velocity_twophase_explicit(dt) ==
              false);
    }

    SECTION("Check momentum and velocity") {
      // Check momentum
      Eigen::Matrix<double, Dim, 1> momentum;
      for (unsigned i = 0; i < momentum.size(); ++i) momentum(i) = 10.;

      // Check initial momentum
      for (unsigned i = 0; i < momentum.size(); ++i) {
        REQUIRE(node->momentum(mpm::NodePhase::nSolid)(i) ==
                Approx(0.).epsilon(Tolerance));
        REQUIRE(node->momentum(mpm::NodePhase::nLiquid)(i) ==
                Approx(0.).epsilon(Tolerance));
      }

      // Check update momentum to 10
      REQUIRE_NOTHROW(
          node->update_momentum(true, mpm::NodePhase::nSolid, momentum));
      REQUIRE_NOTHROW(
          node->update_momentum(true, mpm::NodePhase::nLiquid, momentum));
      for (unsigned i = 0; i < momentum.size(); ++i) {
        REQUIRE(node->momentum(mpm::NodePhase::nSolid)(i) ==
                Approx(10.).epsilon(Tolerance));
        REQUIRE(node->momentum(mpm::NodePhase::nLiquid)(i) ==
                Approx(10.).epsilon(Tolerance));
      }

      // Check update momentum to 20
      REQUIRE_NOTHROW(
          node->update_momentum(true, mpm::NodePhase::nSolid, momentum));
      REQUIRE_NOTHROW(
          node->update_momentum(true, mpm::NodePhase::nLiquid, momentum));
      for (unsigned i = 0; i < momentum.size(); ++i) {
        REQUIRE(node->momentum(mpm::NodePhase::nSolid)(i) ==
                Approx(20.).epsilon(Tolerance));
        REQUIRE(node->momentum(mpm::NodePhase::nLiquid)(i) ==
                Approx(20.).epsilon(Tolerance));
      }

      // Check assign momentum to 10
      REQUIRE_NOTHROW(
          node->update_momentum(false, mpm::NodePhase::nSolid, momentum));
      REQUIRE_NOTHROW(
          node->update_momentum(false, mpm::NodePhase::nLiquid, momentum));
      for (unsigned i = 0; i < momentum.size(); ++i) {
        REQUIRE(node->momentum(mpm::NodePhase::nSolid)(i) ==
                Approx(10.).epsilon(Tolerance));
        REQUIRE(node->momentum(mpm::NodePhase::nLiquid)(i) ==
                Approx(10.).epsilon(Tolerance));
      }

      // Check zero velocity
      for (unsigned i = 0; i < Dim; ++i) {
        REQUIRE(node->velocity(mpm::NodePhase::nSolid)(i) ==
                Approx(0.).epsilon(Tolerance));
        REQUIRE(node->velocity(mpm::NodePhase::nLiquid)(i) ==
                Approx(0.).epsilon(Tolerance));
      }

      // Check mass
      double mass = 0.;
      REQUIRE_NOTHROW(node->update_mass(false, mpm::NodePhase::nSolid, mass));
      REQUIRE(node->mass(mpm::NodePhase::nSolid) ==
              Approx(0.0).epsilon(Tolerance));
      REQUIRE_NOTHROW(node->update_mass(false, mpm::NodePhase::nLiquid, mass));
      REQUIRE(node->mass(mpm::NodePhase::nLiquid) ==
              Approx(0.0).epsilon(Tolerance));
      // Compute and check velocity this should throw zero mass
      node->compute_velocity();

      mass = 100.;
      // Update mass to 100.
      REQUIRE_NOTHROW(node->update_mass(true, mpm::NodePhase::nSolid, mass));
      REQUIRE(node->mass(mpm::NodePhase::nSolid) ==
              Approx(100.).epsilon(Tolerance));
      REQUIRE_NOTHROW(node->update_mass(true, mpm::NodePhase::nLiquid, mass));
      REQUIRE(node->mass(mpm::NodePhase::nLiquid) ==
              Approx(100.).epsilon(Tolerance));

      // Compute and check velocity
      node->compute_velocity();
      for (unsigned i = 0; i < Dim; ++i) {
        REQUIRE(node->velocity(mpm::NodePhase::nSolid)(i) ==
                Approx(0.1).epsilon(Tolerance));
        REQUIRE(node->velocity(mpm::NodePhase::nLiquid)(i) ==
                Approx(0.1).epsilon(Tolerance));
      }

      // Apply velocity constraints
      REQUIRE(node->assign_velocity_constraint(0, 10.5) == true);
      // Check out of bounds condition
      REQUIRE(node->assign_velocity_constraint(1, 10.5) == true);

      // Check velocity before constraints
      Eigen::Matrix<double, Dim, 1> velocity;
      velocity << 0.1;
      for (unsigned i = 0; i < velocity.size(); ++i) {
        REQUIRE(node->velocity(mpm::NodePhase::nSolid)(i) ==
                Approx(velocity(i)).epsilon(Tolerance));
        REQUIRE(node->velocity(mpm::NodePhase::nLiquid)(i) ==
                Approx(velocity(i)).epsilon(Tolerance));
      }

      // Apply constraints
      node->apply_velocity_constraints();

      // Check apply constraints
      velocity << 10.5;
      for (unsigned i = 0; i < velocity.size(); ++i) {
        REQUIRE(node->velocity(mpm::NodePhase::nSolid)(i) ==
                Approx(velocity(i)).epsilon(Tolerance));
        REQUIRE(node->velocity(mpm::NodePhase::nLiquid)(i) ==
                Approx(velocity(i)).epsilon(Tolerance));
      }
    }

    SECTION("Check acceleration") {
      // Check acceleration
      Eigen::Matrix<double, Dim, 1> acceleration;
      for (unsigned i = 0; i < acceleration.size(); ++i) acceleration(i) = 5.;

      for (unsigned i = 0; i < acceleration.size(); ++i) {
        REQUIRE(node->acceleration(mpm::NodePhase::nSolid)(i) ==
                Approx(0.).epsilon(Tolerance));
        REQUIRE(node->acceleration(mpm::NodePhase::nLiquid)(i) ==
                Approx(0.).epsilon(Tolerance));
      }

      REQUIRE_NOTHROW(node->update_acceleration(true, mpm::NodePhase::nSolid,
                                                acceleration));
      REQUIRE_NOTHROW(node->update_acceleration(true, mpm::NodePhase::nLiquid,
                                                0.5 * acceleration));

      for (unsigned i = 0; i < acceleration.size(); ++i) {
        REQUIRE(node->acceleration(mpm::NodePhase::nSolid)(i) ==
                Approx(5.).epsilon(Tolerance));
        REQUIRE(node->acceleration(mpm::NodePhase::nLiquid)(i) ==
                Approx(2.5).epsilon(Tolerance));
      }

      // Apply velocity constraints
      REQUIRE(node->assign_velocity_constraint(0, 10.5) == true);
      // Check out of bounds condition
      REQUIRE(node->assign_velocity_constraint(1, 12.5) == true);

      // Check acceleration before constraints
      acceleration.resize(Dim);
      acceleration << 5.;
      for (unsigned i = 0; i < acceleration.size(); ++i) {
        REQUIRE(node->acceleration(mpm::NodePhase::nSolid)(i) ==
                Approx(acceleration(i)).epsilon(Tolerance));
        REQUIRE(node->acceleration(mpm::NodePhase::nLiquid)(i) ==
                Approx(0.5 * acceleration(i)).epsilon(Tolerance));
      }

      // Apply constraints
      node->apply_velocity_constraints();

      // Check apply constraints
      acceleration << 0.0;
      for (unsigned i = 0; i < acceleration.size(); ++i) {
        REQUIRE(node->acceleration(mpm::NodePhase::nSolid)(i) ==
                Approx(acceleration(i)).epsilon(Tolerance));
        REQUIRE(node->acceleration(mpm::NodePhase::nLiquid)(i) ==
                Approx(acceleration(i)).epsilon(Tolerance));
      }
    }

    SECTION("Check node material ids") {
      // Add material to nodes
      node->append_material_id(0);
      node->append_material_id(1);
      node->append_material_id(4);
      node->append_material_id(0);
      node->append_material_id(2);

      // Check size of material_ids
      REQUIRE(node->material_ids().size() == 4);

      // Check elements of material_ids
      std::vector<unsigned> material_ids = {0, 1, 2, 4};
      auto mat_ids = node->material_ids();
      unsigned i = 0;
      for (auto mitr = mat_ids.begin(); mitr != mat_ids.end(); ++mitr, ++i)
        REQUIRE(*mitr == material_ids.at(i));
    }
  }
}

// \brief Check node class for 2D case
TEST_CASE("Twophase Node is checked for 2D case", "[node][2D][2P]") {
  const unsigned Dim = 2;
  const unsigned Dof = 2;
  const unsigned Nphases = 2;

  Eigen::Vector2d coords;
  coords.setZero();

  // Check for id = 0
  SECTION("Node id is zero") {
    mpm::Index id = 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);
    REQUIRE(node->id() == 0);
  }

  // Check for id is a positive value
  SECTION("Node id is positive") {
    mpm::Index id = std::numeric_limits<mpm::Index>::max();
    std::shared_ptr<mpm::NodeBase<Dim>> node =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);
    REQUIRE(node->id() == std::numeric_limits<mpm::Index>::max());
  }

  // Check for degrees of freedom
  SECTION("Check degrees of freedom") {
    mpm::Index id = 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);
    REQUIRE(node->dof() == 2);
  }

  // Check status
  SECTION("Check status") {
    mpm::Index id = 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);
    REQUIRE(node->status() == false);
    node->assign_status(true);
    REQUIRE(node->status() == true);
  }

  SECTION("Boundary ghost id") {
    mpm::Index id = 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);
    node->ghost_id(5);
    REQUIRE(node->ghost_id() == 5);
  }

  // Check MPI Rank
  SECTION("Check MPI Rank") {
    mpm::Index id = 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);
    REQUIRE(node->id() == 0);

    // Assign MPI ranks
    node->mpi_rank(0);
    node->mpi_rank(0);
    node->mpi_rank(1);

    std::set<unsigned> ranks = node->mpi_ranks();
    REQUIRE(ranks.size() == 2);
    std::vector<unsigned> mpi_ranks = {0, 1};
    unsigned i = 0;
    for (auto it = ranks.begin(); it != ranks.end(); ++it, ++i)
      REQUIRE(*it == mpi_ranks.at(i));
  }

  // Test coordinates function
  SECTION("coordinates function is checked") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;

    std::shared_ptr<mpm::NodeBase<Dim>> node =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);

    // Check for coordinates being zero
    auto coordinates = node->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));
    REQUIRE(coordinates.size() == Dim);

    // Check for negative value of coordinates
    for (unsigned i = 0; i < coordinates.size(); ++i)
      coords(i) = -1. * std::numeric_limits<double>::max();
    node->assign_coordinates(coords);
    coordinates = node->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    REQUIRE(coordinates.size() == Dim);

    // Check for positive value of coordinates
    for (unsigned i = 0; i < coordinates.size(); ++i)
      coords(i) = std::numeric_limits<double>::max();
    node->assign_coordinates(coords);
    coordinates = node->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    REQUIRE(coordinates.size() == Dim);
  }

  SECTION("Check nodal properties") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    std::shared_ptr<mpm::NodeBase<Dim>> node =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);

    // Check mass
    REQUIRE(node->mass(mpm::NodePhase::nSolid) ==
            Approx(0.0).epsilon(Tolerance));
    REQUIRE(node->mass(mpm::NodePhase::nLiquid) ==
            Approx(0.0).epsilon(Tolerance));
    double solid_mass = 100.5;
    double liquid_mass = 200.5;
    // Update mass to 100.5 and 200.5
    REQUIRE_NOTHROW(
        node->update_mass(true, mpm::NodePhase::nSolid, solid_mass));
    REQUIRE_NOTHROW(
        node->update_mass(true, mpm::NodePhase::nLiquid, liquid_mass));
    REQUIRE(node->mass(mpm::NodePhase::nSolid) ==
            Approx(100.5).epsilon(Tolerance));
    REQUIRE(node->mass(mpm::NodePhase::nLiquid) ==
            Approx(200.5).epsilon(Tolerance));
    // Update mass to 201 and 401
    REQUIRE_NOTHROW(
        node->update_mass(true, mpm::NodePhase::nSolid, solid_mass));
    REQUIRE_NOTHROW(
        node->update_mass(true, mpm::NodePhase::nLiquid, liquid_mass));
    REQUIRE(node->mass(mpm::NodePhase::nSolid) ==
            Approx(201.0).epsilon(Tolerance));
    REQUIRE(node->mass(mpm::NodePhase::nLiquid) ==
            Approx(401.0).epsilon(Tolerance));
    // Assign mass to 100 and 200
    solid_mass = 100.;
    liquid_mass = 200.;
    REQUIRE_NOTHROW(
        node->update_mass(false, mpm::NodePhase::nSolid, solid_mass));
    REQUIRE_NOTHROW(
        node->update_mass(false, mpm::NodePhase::nLiquid, liquid_mass));
    REQUIRE(node->mass(mpm::NodePhase::nSolid) ==
            Approx(100.0).epsilon(Tolerance));
    REQUIRE(node->mass(mpm::NodePhase::nLiquid) ==
            Approx(200.0).epsilon(Tolerance));

    SECTION("Check nodal pressure") {
      // Check pressure
      REQUIRE(node->pressure(mpm::NodePhase::nSolid) ==
              Approx(0.0).epsilon(Tolerance));
      REQUIRE(node->pressure(mpm::NodePhase::nLiquid) ==
              Approx(0.0).epsilon(Tolerance));
      double pressure = 1000.7;
      double pore_pressure = 2000.7;
      // Update pressure to 1000.7 and 2000.7
      REQUIRE_NOTHROW(node->update_mass_pressure(mpm::NodePhase::nSolid,
                                                 solid_mass * pressure));
      REQUIRE_NOTHROW(node->update_mass_pressure(mpm::NodePhase::nLiquid,
                                                 liquid_mass * pore_pressure));
      REQUIRE(node->pressure(mpm::NodePhase::nSolid) ==
              Approx(1000.7).epsilon(Tolerance));
      REQUIRE(node->pressure(mpm::NodePhase::nLiquid) ==
              Approx(2000.7).epsilon(Tolerance));
      // Update pressure to 2001.4 and 4001.4
      REQUIRE_NOTHROW(node->update_mass_pressure(mpm::NodePhase::nSolid,
                                                 solid_mass * pressure));
      REQUIRE_NOTHROW(node->update_mass_pressure(mpm::NodePhase::nLiquid,
                                                 liquid_mass * pore_pressure));
      REQUIRE(node->pressure(mpm::NodePhase::nSolid) ==
              Approx(2001.4).epsilon(Tolerance));
      REQUIRE(node->pressure(mpm::NodePhase::nLiquid) ==
              Approx(4001.4).epsilon(Tolerance));
      // Assign pressure to 1000 and 2000
      pressure = 1000.;
      pore_pressure = 2000.;
      node->assign_pressure(mpm::NodePhase::nSolid, pressure);
      node->assign_pressure(mpm::NodePhase::nLiquid, pore_pressure);
      REQUIRE(node->pressure(mpm::NodePhase::nSolid) ==
              Approx(1000.0).epsilon(Tolerance));
      REQUIRE(node->pressure(mpm::NodePhase::nLiquid) ==
              Approx(2000.0).epsilon(Tolerance));
      // Assign mass to 0
      solid_mass = 0.;
      liquid_mass = 0.;
      REQUIRE_NOTHROW(
          node->update_mass(false, mpm::NodePhase::nSolid, solid_mass));
      REQUIRE_NOTHROW(
          node->update_mass(false, mpm::NodePhase::nLiquid, liquid_mass));
      REQUIRE(node->mass(mpm::NodePhase::nSolid) ==
              Approx(0.0).epsilon(Tolerance));
      REQUIRE(node->mass(mpm::NodePhase::nLiquid) ==
              Approx(0.0).epsilon(Tolerance));
      // Try to update pressure to 2000, should throw and keep to 1000.
      node->assign_pressure(mpm::NodePhase::nSolid, pressure);
      node->assign_pressure(mpm::NodePhase::nLiquid, pore_pressure);
      REQUIRE(node->pressure(mpm::NodePhase::nSolid) ==
              Approx(1000.0).epsilon(Tolerance));
      REQUIRE(node->pressure(mpm::NodePhase::nLiquid) ==
              Approx(2000.0).epsilon(Tolerance));
    }

    SECTION("Check volume") {
      // Check volume
      REQUIRE(node->volume(mpm::NodePhase::nMixture) ==
              Approx(0.0).epsilon(Tolerance));
      double volume = 100.5;
      // Update volume to 100.5
      REQUIRE_NOTHROW(
          node->update_volume(true, mpm::NodePhase::nMixture, volume));
      REQUIRE(node->volume(mpm::NodePhase::nMixture) ==
              Approx(100.5).epsilon(Tolerance));
      // Update volume to 201
      REQUIRE_NOTHROW(
          node->update_volume(true, mpm::NodePhase::nMixture, volume));
      REQUIRE(node->volume(mpm::NodePhase::nMixture) ==
              Approx(201.0).epsilon(Tolerance));
      // Assign volume to 100
      volume = 100.;
      REQUIRE_NOTHROW(
          node->update_volume(false, mpm::NodePhase::nMixture, volume));
      REQUIRE(node->volume(mpm::NodePhase::nMixture) ==
              Approx(100.0).epsilon(Tolerance));
    }

    SECTION("Check external force") {
      // Create a force vector
      Eigen::Matrix<double, Dim, 1> force;
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 10.;

      // Check current external force is zero
      for (unsigned i = 0; i < force.size(); ++i) {
        REQUIRE(node->external_force(mpm::NodePhase::nMixture)(i) ==
                Approx(0.).epsilon(Tolerance));
        REQUIRE(node->external_force(mpm::NodePhase::nLiquid)(i) ==
                Approx(0.).epsilon(Tolerance));
      }

      // Update force to 10.0
      REQUIRE_NOTHROW(
          node->update_external_force(true, mpm::NodePhase::nMixture, force));
      REQUIRE_NOTHROW(node->update_external_force(true, mpm::NodePhase::nLiquid,
                                                  0.5 * force));
      for (unsigned i = 0; i < force.size(); ++i) {
        REQUIRE(node->external_force(mpm::NodePhase::nMixture)(i) ==
                Approx(10.).epsilon(Tolerance));
        REQUIRE(node->external_force(mpm::NodePhase::nLiquid)(i) ==
                Approx(5.).epsilon(Tolerance));
      }

      // Update force to 20.0
      REQUIRE_NOTHROW(
          node->update_external_force(true, mpm::NodePhase::nMixture, force));
      REQUIRE_NOTHROW(node->update_external_force(true, mpm::NodePhase::nLiquid,
                                                  0.5 * force));
      for (unsigned i = 0; i < force.size(); ++i) {
        REQUIRE(node->external_force(mpm::NodePhase::nMixture)(i) ==
                Approx(20.).epsilon(Tolerance));
        REQUIRE(node->external_force(mpm::NodePhase::nLiquid)(i) ==
                Approx(10.).epsilon(Tolerance));
      }

      // Assign force as 10.0
      REQUIRE_NOTHROW(
          node->update_external_force(false, mpm::NodePhase::nMixture, force));
      REQUIRE_NOTHROW(node->update_external_force(
          false, mpm::NodePhase::nLiquid, 0.5 * force));
      for (unsigned i = 0; i < force.size(); ++i) {
        REQUIRE(node->external_force(mpm::NodePhase::nMixture)(i) ==
                Approx(10.).epsilon(Tolerance));
        REQUIRE(node->external_force(mpm::NodePhase::nLiquid)(i) ==
                Approx(5.).epsilon(Tolerance));
      }

      SECTION("Check concentrated force") {
        // Set external force to zero
        force.setZero();
        REQUIRE_NOTHROW(node->update_external_force(
            false, mpm::NodePhase::nMixture, force));

        // Concentrated force
        std::shared_ptr<mpm::FunctionBase> ffunction = nullptr;
        double concentrated_force = 65.32;
        const unsigned Direction = 0;
        // Check traction
        for (unsigned i = 0; i < Dim; ++i)
          REQUIRE(node->external_force(mpm::NodePhase::nMixture)(i) ==
                  Approx(0.).epsilon(Tolerance));

        REQUIRE(node->assign_concentrated_force(mpm::NodePhase::nMixture,
                                                Direction, concentrated_force,
                                                ffunction) == true);
        double current_time = 0.0;
        node->apply_concentrated_force(mpm::NodePhase::nMixture, current_time);

        for (unsigned i = 0; i < Dim; ++i) {
          if (i == Direction)
            REQUIRE(node->external_force(mpm::NodePhase::nMixture)(i) ==
                    Approx(concentrated_force).epsilon(Tolerance));

          else
            REQUIRE(node->external_force(mpm::NodePhase::nMixture)(i) ==
                    Approx(0.).epsilon(Tolerance));
        }

        // Check for incorrect direction / phase
        const unsigned wrong_dir = 4;
        REQUIRE(node->assign_concentrated_force(mpm::NodePhase::nMixture,
                                                wrong_dir, concentrated_force,
                                                ffunction) == false);

        // Check again to ensure value hasn't been updated
        for (unsigned i = 0; i < Dim; ++i) {
          if (i == Direction)
            REQUIRE(node->external_force(mpm::NodePhase::nMixture)(i) ==
                    Approx(concentrated_force).epsilon(Tolerance));
          else
            REQUIRE(node->external_force(mpm::NodePhase::nMixture)(i) ==
                    Approx(0.).epsilon(Tolerance));
        }
      }
    }

    SECTION("Check internal force") {
      // Create a force vector
      Eigen::Matrix<double, Dim, 1> force;
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 10.;

      // Check current internal force is zero
      for (unsigned i = 0; i < force.size(); ++i) {
        REQUIRE(node->internal_force(mpm::NodePhase::nMixture)(i) ==
                Approx(0.).epsilon(Tolerance));
        REQUIRE(node->internal_force(mpm::NodePhase::nLiquid)(i) ==
                Approx(0.).epsilon(Tolerance));
      }

      // Update force to 10.0
      REQUIRE_NOTHROW(
          node->update_internal_force(true, mpm::NodePhase::nMixture, force));
      REQUIRE_NOTHROW(node->update_internal_force(true, mpm::NodePhase::nLiquid,
                                                  0.5 * force));
      for (unsigned i = 0; i < force.size(); ++i) {
        REQUIRE(node->internal_force(mpm::NodePhase::nMixture)(i) ==
                Approx(10.).epsilon(Tolerance));
        REQUIRE(node->internal_force(mpm::NodePhase::nLiquid)(i) ==
                Approx(5.).epsilon(Tolerance));
      }

      // Update force to 20.0
      REQUIRE_NOTHROW(
          node->update_internal_force(true, mpm::NodePhase::nMixture, force));
      REQUIRE_NOTHROW(node->update_internal_force(true, mpm::NodePhase::nLiquid,
                                                  0.5 * force));
      for (unsigned i = 0; i < force.size(); ++i) {
        REQUIRE(node->internal_force(mpm::NodePhase::nMixture)(i) ==
                Approx(20.).epsilon(Tolerance));
        REQUIRE(node->internal_force(mpm::NodePhase::nLiquid)(i) ==
                Approx(10.).epsilon(Tolerance));
      }

      // Assign force as 10.0
      REQUIRE_NOTHROW(
          node->update_internal_force(false, mpm::NodePhase::nMixture, force));
      REQUIRE_NOTHROW(node->update_internal_force(
          false, mpm::NodePhase::nLiquid, 0.5 * force));
      for (unsigned i = 0; i < force.size(); ++i) {
        REQUIRE(node->internal_force(mpm::NodePhase::nMixture)(i) ==
                Approx(10.).epsilon(Tolerance));
        REQUIRE(node->internal_force(mpm::NodePhase::nLiquid)(i) ==
                Approx(5.).epsilon(Tolerance));
      }
    }

    SECTION("Check drag force coefficient") {
      // Create a force vector
      Eigen::Matrix<double, Dim, 1> drag_force_coefficient;
      for (unsigned i = 0; i < drag_force_coefficient.size(); ++i)
        drag_force_coefficient(i) = 10.;

      // Check current drag force coefficient is zero
      for (unsigned i = 0; i < drag_force_coefficient.size(); ++i)
        REQUIRE(node->drag_force_coefficient()(i) ==
                Approx(0.).epsilon(Tolerance));

      // Update drag force coefficient to 10.0
      REQUIRE_NOTHROW(
          node->update_drag_force_coefficient(true, drag_force_coefficient));

      for (unsigned i = 0; i < drag_force_coefficient.size(); ++i)
        REQUIRE(node->drag_force_coefficient()(i) ==
                Approx(10.).epsilon(Tolerance));

      // Update force to 20.0
      REQUIRE_NOTHROW(
          node->update_drag_force_coefficient(true, drag_force_coefficient));

      for (unsigned i = 0; i < drag_force_coefficient.size(); ++i)
        REQUIRE(node->drag_force_coefficient()(i) ==
                Approx(20.).epsilon(Tolerance));

      // Assign force as 10.0
      REQUIRE_NOTHROW(
          node->update_drag_force_coefficient(false, drag_force_coefficient));
      for (unsigned i = 0; i < drag_force_coefficient.size(); ++i)
        REQUIRE(node->drag_force_coefficient()(i) ==
                Approx(10.).epsilon(Tolerance));
    }

    SECTION("Check compute acceleration and velocity") {
      // Time step
      const double dt = 0.1;

      // Nodal mass
      double solid_mass = 100.;
      double liquid_mass = 100.;
      // Update mass to 100.
      REQUIRE_NOTHROW(
          node->update_mass(false, mpm::NodePhase::nSolid, solid_mass));
      REQUIRE_NOTHROW(
          node->update_mass(false, mpm::NodePhase::nLiquid, liquid_mass));
      REQUIRE(node->mass(mpm::NodePhase::nSolid) ==
              Approx(solid_mass).epsilon(Tolerance));
      REQUIRE(node->mass(mpm::NodePhase::nLiquid) ==
              Approx(liquid_mass).epsilon(Tolerance));

      // Check internal force
      // Create a force vector
      Eigen::Matrix<double, Dim, 1> force;
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 10. * i;
      // Update force to 10.0
      REQUIRE_NOTHROW(
          node->update_internal_force(false, mpm::NodePhase::nMixture, force));
      REQUIRE_NOTHROW(node->update_internal_force(
          false, mpm::NodePhase::nLiquid, 0.5 * force));
      // Internal force
      for (unsigned i = 0; i < force.size(); ++i) {
        REQUIRE(node->internal_force(mpm::NodePhase::nMixture)(i) ==
                Approx(force(i)).epsilon(Tolerance));
        REQUIRE(node->internal_force(mpm::NodePhase::nLiquid)(i) ==
                Approx(0.5 * force(i)).epsilon(Tolerance));
      }

      // External force
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 5. * i;
      // Update force to 10.0
      REQUIRE_NOTHROW(
          node->update_external_force(false, mpm::NodePhase::nMixture, force));
      REQUIRE_NOTHROW(node->update_external_force(
          false, mpm::NodePhase::nLiquid, 0.5 * force));
      for (unsigned i = 0; i < force.size(); ++i) {
        REQUIRE(node->external_force(mpm::NodePhase::nMixture)(i) ==
                Approx(force(i)).epsilon(Tolerance));
        REQUIRE(node->external_force(mpm::NodePhase::nLiquid)(i) ==
                Approx(0.5 * force(i)).epsilon(Tolerance));
      }

      // Drag force
      Eigen::Matrix<double, Dim, 1> drag_force_coefficient;
      for (unsigned i = 0; i < drag_force_coefficient.size(); ++i)
        drag_force_coefficient(i) = 5. * i;
      // Update force to 10.0
      REQUIRE_NOTHROW(
          node->update_drag_force_coefficient(false, drag_force_coefficient));
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->drag_force_coefficient()(i) ==
                Approx(drag_force_coefficient(i)).epsilon(Tolerance));

      REQUIRE(node->compute_acceleration_velocity_twophase_explicit(dt) ==
              true);

      // Check acceleration
      Eigen::Matrix<double, Dim, 1> liquid_acceleration;
      liquid_acceleration << 0., 0.075;
      Eigen::Matrix<double, Dim, 1> solid_acceleration;
      solid_acceleration << 0., 0.075;

      for (unsigned i = 0; i < solid_acceleration.size(); ++i) {
        REQUIRE(node->acceleration(mpm::NodePhase::nSolid)(i) ==
                Approx(solid_acceleration(i)).epsilon(Tolerance));
        REQUIRE(node->acceleration(mpm::NodePhase::nLiquid)(i) ==
                Approx(liquid_acceleration(i)).epsilon(Tolerance));
      }

      // Check velocity
      Eigen::Matrix<double, Dim, 1> solid_velocity = solid_acceleration * dt;
      Eigen::Matrix<double, Dim, 1> liquid_velocity = liquid_acceleration * dt;
      for (unsigned i = 0; i < solid_velocity.size(); ++i) {
        REQUIRE(node->velocity(mpm::NodePhase::nSolid)(i) ==
                Approx(solid_velocity(i)).epsilon(Tolerance));
        REQUIRE(node->velocity(mpm::NodePhase::nLiquid)(i) ==
                Approx(liquid_velocity(i)).epsilon(Tolerance));
      }

      // Apply velocity constraints
      REQUIRE(node->assign_velocity_constraint(0, 10.5) == true);
      REQUIRE(node->assign_velocity_constraint(1, 0.03) == true);
      REQUIRE(node->assign_velocity_constraint(2, 20.5) == true);
      REQUIRE(node->assign_velocity_constraint(3, 1.03) == true);
      REQUIRE(node->compute_acceleration_velocity_twophase_explicit(dt) ==
              true);

      // Test velocity with constraints
      solid_velocity << 10.5, 0.03;
      liquid_velocity << 20.5, 1.03;
      for (unsigned i = 0; i < solid_velocity.size(); ++i) {
        REQUIRE(node->velocity(mpm::NodePhase::nSolid)(i) ==
                Approx(solid_velocity(i)).epsilon(Tolerance));
        REQUIRE(node->velocity(mpm::NodePhase::nLiquid)(i) ==
                Approx(liquid_velocity(i)).epsilon(Tolerance));
      }

      // Test acceleration with constraints
      solid_acceleration.setZero();
      liquid_acceleration.setZero();
      for (unsigned i = 0; i < solid_acceleration.size(); ++i) {
        REQUIRE(node->acceleration(mpm::NodePhase::nSolid)(i) ==
                Approx(solid_acceleration(i)).epsilon(Tolerance));
        REQUIRE(node->acceleration(mpm::NodePhase::nLiquid)(i) ==
                Approx(liquid_acceleration(i)).epsilon(Tolerance));
      }

      // Apply cundall damping when calculating acceleration
      REQUIRE(node->compute_acceleration_velocity_twophase_explicit_cundall(
                  dt, 0.05) == true);

      // Test acceleration with cundall damping
      solid_acceleration << 0., 0.;
      liquid_acceleration << 0., 0.;
      for (unsigned i = 0; i < solid_acceleration.size(); ++i) {
        REQUIRE(node->acceleration(mpm::NodePhase::nSolid)(i) ==
                Approx(solid_acceleration(i)).epsilon(Tolerance));
        REQUIRE(node->acceleration(mpm::NodePhase::nLiquid)(i) ==
                Approx(liquid_acceleration(i)).epsilon(Tolerance));
      }

      // Exception check when mass is zero
      REQUIRE_NOTHROW(node->update_mass(false, mpm::NodePhase::nSolid, 0.));
      REQUIRE_NOTHROW(node->update_mass(false, mpm::NodePhase::nLiquid, 0.));
      REQUIRE(node->mass(mpm::NodePhase::nSolid) ==
              Approx(0.).epsilon(Tolerance));
      REQUIRE(node->mass(mpm::NodePhase::nLiquid) ==
              Approx(0.).epsilon(Tolerance));
      REQUIRE(node->compute_acceleration_velocity_twophase_explicit(dt) ==
              false);
    }

    SECTION("Check momentum, velocity and acceleration") {
      // Time step
      const double dt = 0.1;

      // Check momentum
      Eigen::Matrix<double, Dim, 1> momentum;
      for (unsigned i = 0; i < momentum.size(); ++i) momentum(i) = 10.;

      // Check initial momentum
      for (unsigned i = 0; i < momentum.size(); ++i) {
        REQUIRE(node->momentum(mpm::NodePhase::nSolid)(i) ==
                Approx(0.).epsilon(Tolerance));
        REQUIRE(node->momentum(mpm::NodePhase::nLiquid)(i) ==
                Approx(0.).epsilon(Tolerance));
      }

      // Check update momentum to 10
      REQUIRE_NOTHROW(
          node->update_momentum(true, mpm::NodePhase::nSolid, momentum));
      REQUIRE_NOTHROW(
          node->update_momentum(true, mpm::NodePhase::nLiquid, momentum));
      for (unsigned i = 0; i < momentum.size(); ++i) {
        REQUIRE(node->momentum(mpm::NodePhase::nSolid)(i) ==
                Approx(10.).epsilon(Tolerance));
        REQUIRE(node->momentum(mpm::NodePhase::nLiquid)(i) ==
                Approx(10.).epsilon(Tolerance));
      }

      // Check update momentum to 20
      REQUIRE_NOTHROW(
          node->update_momentum(true, mpm::NodePhase::nSolid, momentum));
      REQUIRE_NOTHROW(
          node->update_momentum(true, mpm::NodePhase::nLiquid, momentum));
      for (unsigned i = 0; i < momentum.size(); ++i) {
        REQUIRE(node->momentum(mpm::NodePhase::nSolid)(i) ==
                Approx(20.).epsilon(Tolerance));
        REQUIRE(node->momentum(mpm::NodePhase::nLiquid)(i) ==
                Approx(20.).epsilon(Tolerance));
      }

      // Check assign momentum to 10
      REQUIRE_NOTHROW(
          node->update_momentum(false, mpm::NodePhase::nSolid, momentum));
      REQUIRE_NOTHROW(
          node->update_momentum(false, mpm::NodePhase::nLiquid, momentum));
      for (unsigned i = 0; i < momentum.size(); ++i) {
        REQUIRE(node->momentum(mpm::NodePhase::nSolid)(i) ==
                Approx(10.).epsilon(Tolerance));
        REQUIRE(node->momentum(mpm::NodePhase::nLiquid)(i) ==
                Approx(10.).epsilon(Tolerance));
      }

      // Check zero velocity
      for (unsigned i = 0; i < Dim; ++i) {
        REQUIRE(node->velocity(mpm::NodePhase::nSolid)(i) ==
                Approx(0.).epsilon(Tolerance));
        REQUIRE(node->velocity(mpm::NodePhase::nLiquid)(i) ==
                Approx(0.).epsilon(Tolerance));
      }

      // Check mass
      double mass = 0.;
      REQUIRE_NOTHROW(node->update_mass(false, mpm::NodePhase::nSolid, mass));
      REQUIRE(node->mass(mpm::NodePhase::nSolid) ==
              Approx(0.0).epsilon(Tolerance));
      REQUIRE_NOTHROW(node->update_mass(false, mpm::NodePhase::nLiquid, mass));
      REQUIRE(node->mass(mpm::NodePhase::nLiquid) ==
              Approx(0.0).epsilon(Tolerance));
      // Compute and check velocity this should throw zero mass
      node->compute_velocity();

      mass = 100.;
      // Update mass to 100.5
      REQUIRE_NOTHROW(node->update_mass(true, mpm::NodePhase::nSolid, mass));
      REQUIRE(node->mass(mpm::NodePhase::nSolid) ==
              Approx(100.).epsilon(Tolerance));
      REQUIRE_NOTHROW(node->update_mass(true, mpm::NodePhase::nLiquid, mass));
      REQUIRE(node->mass(mpm::NodePhase::nLiquid) ==
              Approx(100.).epsilon(Tolerance));

      // Compute and check velocity
      node->compute_velocity();
      for (unsigned i = 0; i < Dim; ++i) {
        REQUIRE(node->velocity(mpm::NodePhase::nSolid)(i) ==
                Approx(0.1).epsilon(Tolerance));
        REQUIRE(node->velocity(mpm::NodePhase::nLiquid)(i) ==
                Approx(0.1).epsilon(Tolerance));
      }

      // Check acceleration
      Eigen::Matrix<double, Dim, 1> acceleration;
      for (unsigned i = 0; i < acceleration.size(); ++i) acceleration(i) = 5.;

      for (unsigned i = 0; i < acceleration.size(); ++i) {
        REQUIRE(node->acceleration(mpm::NodePhase::nSolid)(i) ==
                Approx(0.).epsilon(Tolerance));
        REQUIRE(node->acceleration(mpm::NodePhase::nLiquid)(i) ==
                Approx(0.).epsilon(Tolerance));
      }

      REQUIRE_NOTHROW(node->update_acceleration(true, mpm::NodePhase::nSolid,
                                                acceleration));
      REQUIRE_NOTHROW(node->update_acceleration(true, mpm::NodePhase::nLiquid,
                                                acceleration));
      for (unsigned i = 0; i < acceleration.size(); ++i) {
        REQUIRE(node->acceleration(mpm::NodePhase::nSolid)(i) ==
                Approx(5.).epsilon(Tolerance));
        REQUIRE(node->acceleration(mpm::NodePhase::nLiquid)(i) ==
                Approx(5.).epsilon(Tolerance));
      }

      // Check if exception is handled
      Eigen::Matrix<double, Dim, 1> acceleration_bad;
      for (unsigned i = 0; i < acceleration_bad.size(); ++i)
        acceleration_bad(i) = 10.;

      // Check velocity before constraints
      Eigen::Matrix<double, Dim, 1> velocity;
      velocity << 0.1, 0.1;
      for (unsigned i = 0; i < velocity.size(); ++i) {
        REQUIRE(node->velocity(mpm::NodePhase::nSolid)(i) ==
                Approx(velocity(i)).epsilon(Tolerance));
        REQUIRE(node->velocity(mpm::NodePhase::nLiquid)(i) ==
                Approx(velocity(i)).epsilon(Tolerance));
      }

      // Check acceleration before constraints
      acceleration << 5., 5.;
      for (unsigned i = 0; i < acceleration.size(); ++i) {
        REQUIRE(node->acceleration(mpm::NodePhase::nSolid)(i) ==
                Approx(acceleration(i)).epsilon(Tolerance));
        REQUIRE(node->acceleration(mpm::NodePhase::nLiquid)(i) ==
                Approx(acceleration(i)).epsilon(Tolerance));
      }

      SECTION("Check Cartesian velocity constraints") {
        // Apply velocity constraints
        REQUIRE(node->assign_velocity_constraint(0, -12.5) == true);
        // Check out of bounds condition
        REQUIRE(node->assign_velocity_constraint(2, -12.5) == true);

        // Apply constraints
        node->apply_velocity_constraints();

        // Check apply constraints
        velocity << -12.5, 0.1;
        for (unsigned i = 0; i < velocity.size(); ++i) {
          REQUIRE(node->velocity(mpm::NodePhase::nSolid)(i) ==
                  Approx(velocity(i)).epsilon(Tolerance));
          REQUIRE(node->velocity(mpm::NodePhase::nLiquid)(i) ==
                  Approx(velocity(i)).epsilon(Tolerance));
        }

        acceleration << 0., 5.;
        for (unsigned i = 0; i < acceleration.size(); ++i) {
          REQUIRE(node->acceleration(mpm::NodePhase::nSolid)(i) ==
                  Approx(acceleration(i)).epsilon(Tolerance));
          REQUIRE(node->acceleration(mpm::NodePhase::nLiquid)(i) ==
                  Approx(acceleration(i)).epsilon(Tolerance));
        }
      }

      SECTION("Check general velocity constraints in 1 direction") {
        // Apply velocity constraints
        REQUIRE(node->assign_velocity_constraint(0, -12.5) == true);
        REQUIRE(node->assign_velocity_constraint(2, -12.5) == true);

        // Apply rotation matrix with Euler angles alpha = 10 deg, beta
        // = 30 deg
        Eigen::Matrix<double, Dim, 1> euler_angles;
        euler_angles << 10. * M_PI / 180, 30. * M_PI / 180;
        const auto rotation_matrix =
            mpm::geometry::rotation_matrix(euler_angles);
        node->assign_rotation_matrix(rotation_matrix);
        const auto inverse_rotation_matrix = rotation_matrix.inverse();

        // Apply inclined velocity constraints
        node->apply_velocity_constraints();

        // Check applied velocity constraints in the global coordinates
        velocity << -9.583478335521184, -8.025403099849004;
        for (unsigned i = 0; i < Dim; ++i) {
          REQUIRE(node->velocity(mpm::NodePhase::nSolid)(i) ==
                  Approx(velocity(i)).epsilon(Tolerance));
          REQUIRE(node->velocity(mpm::NodePhase::nLiquid)(i) ==
                  Approx(velocity(i)).epsilon(Tolerance));
        }

        // Check that the velocity is as specified in local coordinate
        REQUIRE((inverse_rotation_matrix *
                 node->velocity(mpm::NodePhase::nSolid))(0) ==
                Approx(-12.5).epsilon(Tolerance));
        REQUIRE((inverse_rotation_matrix *
                 node->velocity(mpm::NodePhase::nLiquid))(0) ==
                Approx(-12.5).epsilon(Tolerance));

        // Check applied constraints on acceleration in the global coordinates
        acceleration << -0.396139826697847, 0.472101061636807;
        for (unsigned i = 0; i < Dim; ++i) {
          REQUIRE(node->acceleration(mpm::NodePhase::nSolid)(i) ==
                  Approx(acceleration(i)).epsilon(Tolerance));
          REQUIRE(node->acceleration(mpm::NodePhase::nLiquid)(i) ==
                  Approx(acceleration(i)).epsilon(Tolerance));
        }

        // Check that the acceleration is 0 in local coordinate
        REQUIRE((inverse_rotation_matrix *
                 node->acceleration(mpm::NodePhase::nSolid))(0) ==
                Approx(0).epsilon(Tolerance));
        REQUIRE((inverse_rotation_matrix *
                 node->acceleration(mpm::NodePhase::nLiquid))(0) ==
                Approx(0).epsilon(Tolerance));
      }

      SECTION("Check general velocity constraints in all directions") {
        // Apply velocity constraints
        REQUIRE(node->assign_velocity_constraint(0, -12.5) == true);
        REQUIRE(node->assign_velocity_constraint(1, 7.5) == true);
        REQUIRE(node->assign_velocity_constraint(2, -12.5) == true);
        REQUIRE(node->assign_velocity_constraint(3, 7.5) == true);

        // Apply rotation matrix with Euler angles alpha = -10 deg, beta = 30
        // deg
        Eigen::Matrix<double, Dim, 1> euler_angles;
        euler_angles << -10. * M_PI / 180, 30. * M_PI / 180;
        const auto rotation_matrix =
            mpm::geometry::rotation_matrix(euler_angles);
        node->assign_rotation_matrix(rotation_matrix);
        const auto inverse_rotation_matrix = rotation_matrix.inverse();

        // Apply inclined velocity constraints
        node->apply_velocity_constraints();

        // Check applied velocity constraints in the global coordinates
        velocity << -14.311308834766370, 2.772442864323454;
        for (unsigned i = 0; i < Dim; ++i) {
          REQUIRE(node->velocity(mpm::NodePhase::nSolid)(i) ==
                  Approx(velocity(i)).epsilon(Tolerance));
          REQUIRE(node->velocity(mpm::NodePhase::nLiquid)(i) ==
                  Approx(velocity(i)).epsilon(Tolerance));
        }

        // Check that the velocity is as specified in local coordinate
        REQUIRE((inverse_rotation_matrix *
                 node->velocity(mpm::NodePhase::nSolid))(0) ==
                Approx(-12.5).epsilon(Tolerance));
        REQUIRE((inverse_rotation_matrix *
                 node->velocity(mpm::NodePhase::nSolid))(1) ==
                Approx(7.5).epsilon(Tolerance));
        REQUIRE((inverse_rotation_matrix *
                 node->velocity(mpm::NodePhase::nLiquid))(0) ==
                Approx(-12.5).epsilon(Tolerance));
        REQUIRE((inverse_rotation_matrix *
                 node->velocity(mpm::NodePhase::nLiquid))(1) ==
                Approx(7.5).epsilon(Tolerance));

        // Check applied constraints on acceleration in the global coordinates
        acceleration << 0, 0;
        for (unsigned i = 0; i < Dim; ++i) {
          REQUIRE(node->acceleration(mpm::NodePhase::nSolid)(i) ==
                  Approx(acceleration(i)).epsilon(Tolerance));
          REQUIRE(node->acceleration(mpm::NodePhase::nLiquid)(i) ==
                  Approx(acceleration(i)).epsilon(Tolerance));
        }

        // Check that the acceleration is 0 in local coordinate
        REQUIRE((inverse_rotation_matrix *
                 node->acceleration(mpm::NodePhase::nSolid))(0) ==
                Approx(0).epsilon(Tolerance));
        REQUIRE((inverse_rotation_matrix *
                 node->acceleration(mpm::NodePhase::nSolid))(1) ==
                Approx(0).epsilon(Tolerance));
        REQUIRE((inverse_rotation_matrix *
                 node->acceleration(mpm::NodePhase::nLiquid))(0) ==
                Approx(0).epsilon(Tolerance));
        REQUIRE((inverse_rotation_matrix *
                 node->acceleration(mpm::NodePhase::nLiquid))(1) ==
                Approx(0).epsilon(Tolerance));
      }

      SECTION("Check Cartesian friction constraints") {
        // Apply friction constraints
        REQUIRE(node->assign_friction_constraint(1, 1, 0.2) == true);
        // Check out of bounds condition
        REQUIRE(node->assign_friction_constraint(2, 1, 0.2) == false);

        // Apply friction constraints
        node->apply_friction_constraints(dt);

        // Check apply constraints
        acceleration << 4., 5.;
        for (unsigned i = 0; i < acceleration.size(); ++i)
          REQUIRE(node->acceleration(mpm::NodePhase::nSolid)(i) ==
                  Approx(acceleration(i)).epsilon(Tolerance));
      }

      SECTION("Check general friction constraints in 1 direction") {
        // Apply friction constraints
        REQUIRE(node->assign_friction_constraint(1, 1, 0.2) == true);

        // Apply rotation matrix with Euler angles alpha = 10 deg, beta = 30 deg
        Eigen::Matrix<double, Dim, 1> euler_angles;
        euler_angles << 10. * M_PI / 180, 30. * M_PI / 180;
        const auto rotation_matrix =
            mpm::geometry::rotation_matrix(euler_angles);
        node->assign_rotation_matrix(rotation_matrix);
        const auto inverse_rotation_matrix = rotation_matrix.inverse();

        // Apply general friction constraints
        node->apply_friction_constraints(dt);

        // Check applied constraints on acceleration in the global coordinates
        acceleration << 4.905579787672637, 4.920772034660430;
        for (unsigned i = 0; i < Dim; ++i)
          REQUIRE(node->acceleration(mpm::NodePhase::nSolid)(i) ==
                  Approx(acceleration(i)).epsilon(Tolerance));

        // Check the acceleration in local coordinate
        acceleration << 6.920903430595146, 0.616284167162194;
        for (unsigned i = 0; i < Dim; ++i)
          REQUIRE((inverse_rotation_matrix *
                   node->acceleration(mpm::NodePhase::nSolid))(i) ==
                  Approx(acceleration(i)).epsilon(Tolerance));
      }
    }

    SECTION("Check node material ids") {
      // Add material to nodes
      node->append_material_id(0);
      node->append_material_id(1);
      node->append_material_id(4);
      node->append_material_id(0);
      node->append_material_id(2);

      // Check size of material_ids
      REQUIRE(node->material_ids().size() == 4);

      // Check elements of material_ids
      std::vector<unsigned> material_ids = {0, 1, 2, 4};
      auto mat_ids = node->material_ids();
      unsigned i = 0;
      for (auto mitr = mat_ids.begin(); mitr != mat_ids.end(); ++mitr, ++i)
        REQUIRE(*mitr == material_ids.at(i));
    }
  }
}

// \brief Check node class for 3D case
TEST_CASE("Twophase Node is checked for 3D case", "[node][3D][2P]") {
  const unsigned Dim = 3;
  const unsigned Dof = 3;
  const unsigned Nphases = 2;

  Eigen::Vector3d coords;
  coords.setZero();

  // Check for id = 0
  SECTION("Node id is zero") {
    mpm::Index id = 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);
    REQUIRE(node->id() == 0);
  }

  // Check for id is a positive value
  SECTION("Node id is positive") {
    mpm::Index id = std::numeric_limits<mpm::Index>::max();
    std::shared_ptr<mpm::NodeBase<Dim>> node =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);
    REQUIRE(node->id() == std::numeric_limits<mpm::Index>::max());
  }

  // Check for degrees of freedom
  SECTION("Check degrees of freedom") {
    mpm::Index id = 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);
    REQUIRE(node->dof() == 3);
  }

  // Check status
  SECTION("Check status") {
    mpm::Index id = 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);
    REQUIRE(node->status() == false);
    node->assign_status(true);
    REQUIRE(node->status() == true);
  }

  SECTION("Boundary ghost id") {
    mpm::Index id = 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);
    node->ghost_id(5);
    REQUIRE(node->ghost_id() == 5);
  }

  // Check MPI Rank
  SECTION("Check MPI Rank") {
    mpm::Index id = 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);
    REQUIRE(node->id() == 0);

    // Assign MPI ranks
    node->mpi_rank(0);
    node->mpi_rank(0);
    node->mpi_rank(1);

    std::set<unsigned> ranks = node->mpi_ranks();
    REQUIRE(ranks.size() == 2);
    std::vector<unsigned> mpi_ranks = {0, 1};
    unsigned i = 0;
    for (auto it = ranks.begin(); it != ranks.end(); ++it, ++i)
      REQUIRE(*it == mpi_ranks.at(i));
  }

  // Test coordinates function
  SECTION("coordinates function is checked") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;

    std::shared_ptr<mpm::NodeBase<Dim>> node =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);

    // Check for coordinates being zero
    auto coordinates = node->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));
    REQUIRE(coordinates.size() == Dim);

    // Check for negative value of coordinates
    for (unsigned i = 0; i < coordinates.size(); ++i)
      coords(i) = -1. * std::numeric_limits<double>::max();
    node->assign_coordinates(coords);
    coordinates = node->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    REQUIRE(coordinates.size() == Dim);

    // Check for positive value of coordinates
    for (unsigned i = 0; i < coordinates.size(); ++i)
      coords(i) = std::numeric_limits<double>::max();
    node->assign_coordinates(coords);
    coordinates = node->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    REQUIRE(coordinates.size() == Dim);
  }

  SECTION("Check nodal properties") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    std::shared_ptr<mpm::NodeBase<Dim>> node =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);

    // Check mass
    REQUIRE(node->mass(mpm::NodePhase::nSolid) ==
            Approx(0.0).epsilon(Tolerance));
    REQUIRE(node->mass(mpm::NodePhase::nLiquid) ==
            Approx(0.0).epsilon(Tolerance));
    double solid_mass = 100.5;
    double liquid_mass = 200.5;
    // Update mass to 100.5 and 200.5
    REQUIRE_NOTHROW(
        node->update_mass(true, mpm::NodePhase::nSolid, solid_mass));
    REQUIRE_NOTHROW(
        node->update_mass(true, mpm::NodePhase::nLiquid, liquid_mass));
    REQUIRE(node->mass(mpm::NodePhase::nSolid) ==
            Approx(100.5).epsilon(Tolerance));
    REQUIRE(node->mass(mpm::NodePhase::nLiquid) ==
            Approx(200.5).epsilon(Tolerance));
    // Update mass to 201 and 401
    REQUIRE_NOTHROW(
        node->update_mass(true, mpm::NodePhase::nSolid, solid_mass));
    REQUIRE_NOTHROW(
        node->update_mass(true, mpm::NodePhase::nLiquid, liquid_mass));
    REQUIRE(node->mass(mpm::NodePhase::nSolid) ==
            Approx(201.0).epsilon(Tolerance));
    REQUIRE(node->mass(mpm::NodePhase::nLiquid) ==
            Approx(401.0).epsilon(Tolerance));
    // Assign mass to 100 and 200
    solid_mass = 100.;
    liquid_mass = 200.;
    REQUIRE_NOTHROW(
        node->update_mass(false, mpm::NodePhase::nSolid, solid_mass));
    REQUIRE_NOTHROW(
        node->update_mass(false, mpm::NodePhase::nLiquid, liquid_mass));
    REQUIRE(node->mass(mpm::NodePhase::nSolid) ==
            Approx(100.0).epsilon(Tolerance));
    REQUIRE(node->mass(mpm::NodePhase::nLiquid) ==
            Approx(200.0).epsilon(Tolerance));

    SECTION("Check nodal pressure") {
      // Check pressure
      REQUIRE(node->pressure(mpm::NodePhase::nSolid) ==
              Approx(0.0).epsilon(Tolerance));
      REQUIRE(node->pressure(mpm::NodePhase::nLiquid) ==
              Approx(0.0).epsilon(Tolerance));
      double pressure = 1000.7;
      double pore_pressure = 2000.7;
      // Update pressure to 1000.7 and 2000.7
      REQUIRE_NOTHROW(node->update_mass_pressure(mpm::NodePhase::nSolid,
                                                 solid_mass * pressure));
      REQUIRE_NOTHROW(node->update_mass_pressure(mpm::NodePhase::nLiquid,
                                                 liquid_mass * pore_pressure));
      REQUIRE(node->pressure(mpm::NodePhase::nSolid) ==
              Approx(1000.7).epsilon(Tolerance));
      REQUIRE(node->pressure(mpm::NodePhase::nLiquid) ==
              Approx(2000.7).epsilon(Tolerance));
      // Update pressure to 2001.4 and 4001.4
      REQUIRE_NOTHROW(node->update_mass_pressure(mpm::NodePhase::nSolid,
                                                 solid_mass * pressure));
      REQUIRE_NOTHROW(node->update_mass_pressure(mpm::NodePhase::nLiquid,
                                                 liquid_mass * pore_pressure));
      REQUIRE(node->pressure(mpm::NodePhase::nSolid) ==
              Approx(2001.4).epsilon(Tolerance));
      REQUIRE(node->pressure(mpm::NodePhase::nLiquid) ==
              Approx(4001.4).epsilon(Tolerance));
      // Assign pressure to 1000 and 2000
      pressure = 1000.;
      pore_pressure = 2000.;
      node->assign_pressure(mpm::NodePhase::nSolid, pressure);
      node->assign_pressure(mpm::NodePhase::nLiquid, pore_pressure);
      REQUIRE(node->pressure(mpm::NodePhase::nSolid) ==
              Approx(1000.0).epsilon(Tolerance));
      REQUIRE(node->pressure(mpm::NodePhase::nLiquid) ==
              Approx(2000.0).epsilon(Tolerance));
      // Assign mass to 0
      solid_mass = 0.;
      liquid_mass = 0.;
      REQUIRE_NOTHROW(
          node->update_mass(false, mpm::NodePhase::nSolid, solid_mass));
      REQUIRE_NOTHROW(
          node->update_mass(false, mpm::NodePhase::nLiquid, liquid_mass));
      REQUIRE(node->mass(mpm::NodePhase::nSolid) ==
              Approx(0.0).epsilon(Tolerance));
      REQUIRE(node->mass(mpm::NodePhase::nLiquid) ==
              Approx(0.0).epsilon(Tolerance));
      // Try to update pressure to 2000, should throw and keep to 1000.
      node->assign_pressure(mpm::NodePhase::nSolid, pressure);
      node->assign_pressure(mpm::NodePhase::nLiquid, pore_pressure);
      REQUIRE(node->pressure(mpm::NodePhase::nSolid) ==
              Approx(1000.0).epsilon(Tolerance));
      REQUIRE(node->pressure(mpm::NodePhase::nLiquid) ==
              Approx(2000.0).epsilon(Tolerance));
    }

    SECTION("Check external force") {
      // Create a force vector
      Eigen::Matrix<double, Dim, 1> force;
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 10.;

      // Check current external force is zero
      for (unsigned i = 0; i < force.size(); ++i) {
        REQUIRE(node->external_force(mpm::NodePhase::nMixture)(i) ==
                Approx(0.).epsilon(Tolerance));
        REQUIRE(node->external_force(mpm::NodePhase::nLiquid)(i) ==
                Approx(0.).epsilon(Tolerance));
      }

      // Update force to 10.0
      REQUIRE_NOTHROW(
          node->update_external_force(true, mpm::NodePhase::nMixture, force));
      REQUIRE_NOTHROW(node->update_external_force(true, mpm::NodePhase::nLiquid,
                                                  0.5 * force));
      for (unsigned i = 0; i < force.size(); ++i) {
        REQUIRE(node->external_force(mpm::NodePhase::nMixture)(i) ==
                Approx(10.).epsilon(Tolerance));
        REQUIRE(node->external_force(mpm::NodePhase::nLiquid)(i) ==
                Approx(5.).epsilon(Tolerance));
      }

      // Update force to 20.0
      REQUIRE_NOTHROW(
          node->update_external_force(true, mpm::NodePhase::nMixture, force));
      REQUIRE_NOTHROW(node->update_external_force(true, mpm::NodePhase::nLiquid,
                                                  0.5 * force));
      for (unsigned i = 0; i < force.size(); ++i) {
        REQUIRE(node->external_force(mpm::NodePhase::nMixture)(i) ==
                Approx(20.).epsilon(Tolerance));
        REQUIRE(node->external_force(mpm::NodePhase::nLiquid)(i) ==
                Approx(10.).epsilon(Tolerance));
      }

      // Assign force as 10.0
      REQUIRE_NOTHROW(
          node->update_external_force(false, mpm::NodePhase::nMixture, force));
      REQUIRE_NOTHROW(node->update_external_force(
          false, mpm::NodePhase::nLiquid, 0.5 * force));
      for (unsigned i = 0; i < force.size(); ++i) {
        REQUIRE(node->external_force(mpm::NodePhase::nMixture)(i) ==
                Approx(10.).epsilon(Tolerance));
        REQUIRE(node->external_force(mpm::NodePhase::nLiquid)(i) ==
                Approx(5.).epsilon(Tolerance));
      }

      SECTION("Check concentrated force") {
        // Set external force to zero
        force.setZero();
        REQUIRE_NOTHROW(node->update_external_force(
            false, mpm::NodePhase::nMixture, force));

        // Concentrated force
        std::shared_ptr<mpm::FunctionBase> ffunction = nullptr;
        double concentrated_force = 65.32;
        const unsigned Direction = 0;
        // Check traction
        for (unsigned i = 0; i < Dim; ++i)
          REQUIRE(node->external_force(mpm::NodePhase::nMixture)(i) ==
                  Approx(0.).epsilon(Tolerance));

        REQUIRE(node->assign_concentrated_force(mpm::NodePhase::nMixture,
                                                Direction, concentrated_force,
                                                ffunction) == true);
        double current_time = 0.0;
        node->apply_concentrated_force(mpm::NodePhase::nMixture, current_time);

        for (unsigned i = 0; i < Dim; ++i) {
          if (i == Direction)
            REQUIRE(node->external_force(mpm::NodePhase::nMixture)(i) ==
                    Approx(concentrated_force).epsilon(Tolerance));

          else
            REQUIRE(node->external_force(mpm::NodePhase::nMixture)(i) ==
                    Approx(0.).epsilon(Tolerance));
        }

        // Check for incorrect direction / phase
        const unsigned wrong_dir = 4;
        REQUIRE(node->assign_concentrated_force(mpm::NodePhase::nMixture,
                                                wrong_dir, concentrated_force,
                                                ffunction) == false);

        // Check again to ensure value hasn't been updated
        for (unsigned i = 0; i < Dim; ++i) {
          if (i == Direction)
            REQUIRE(node->external_force(mpm::NodePhase::nMixture)(i) ==
                    Approx(concentrated_force).epsilon(Tolerance));
          else
            REQUIRE(node->external_force(mpm::NodePhase::nMixture)(i) ==
                    Approx(0.).epsilon(Tolerance));
        }
      }
    }

    SECTION("Check internal force") {
      // Create a force vector
      Eigen::Matrix<double, Dim, 1> force;
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 10.;

      // Check current internal force is zero
      for (unsigned i = 0; i < force.size(); ++i) {
        REQUIRE(node->internal_force(mpm::NodePhase::nMixture)(i) ==
                Approx(0.).epsilon(Tolerance));
        REQUIRE(node->internal_force(mpm::NodePhase::nLiquid)(i) ==
                Approx(0.).epsilon(Tolerance));
      }

      // Update force to 10.0
      REQUIRE_NOTHROW(
          node->update_internal_force(true, mpm::NodePhase::nMixture, force));
      REQUIRE_NOTHROW(node->update_internal_force(true, mpm::NodePhase::nLiquid,
                                                  0.5 * force));
      for (unsigned i = 0; i < force.size(); ++i) {
        REQUIRE(node->internal_force(mpm::NodePhase::nMixture)(i) ==
                Approx(10.).epsilon(Tolerance));
        REQUIRE(node->internal_force(mpm::NodePhase::nLiquid)(i) ==
                Approx(5.).epsilon(Tolerance));
      }

      // Update force to 20.0
      REQUIRE_NOTHROW(
          node->update_internal_force(true, mpm::NodePhase::nMixture, force));
      REQUIRE_NOTHROW(node->update_internal_force(true, mpm::NodePhase::nLiquid,
                                                  0.5 * force));
      for (unsigned i = 0; i < force.size(); ++i) {
        REQUIRE(node->internal_force(mpm::NodePhase::nMixture)(i) ==
                Approx(20.).epsilon(Tolerance));
        REQUIRE(node->internal_force(mpm::NodePhase::nLiquid)(i) ==
                Approx(10.).epsilon(Tolerance));
      }

      // Assign force as 10.0
      REQUIRE_NOTHROW(
          node->update_internal_force(false, mpm::NodePhase::nMixture, force));
      REQUIRE_NOTHROW(node->update_internal_force(
          false, mpm::NodePhase::nLiquid, 0.5 * force));
      for (unsigned i = 0; i < force.size(); ++i) {
        REQUIRE(node->internal_force(mpm::NodePhase::nMixture)(i) ==
                Approx(10.).epsilon(Tolerance));
        REQUIRE(node->internal_force(mpm::NodePhase::nLiquid)(i) ==
                Approx(5.).epsilon(Tolerance));
      }
    }

    SECTION("Check drag force coefficient") {
      // Create a force vector
      Eigen::Matrix<double, Dim, 1> drag_force_coefficient;
      for (unsigned i = 0; i < drag_force_coefficient.size(); ++i)
        drag_force_coefficient(i) = 10.;

      // Check current drag force coefficient is zero
      for (unsigned i = 0; i < drag_force_coefficient.size(); ++i)
        REQUIRE(node->drag_force_coefficient()(i) ==
                Approx(0.).epsilon(Tolerance));

      // Update drag force coefficient to 10.0
      REQUIRE_NOTHROW(
          node->update_drag_force_coefficient(true, drag_force_coefficient));

      for (unsigned i = 0; i < drag_force_coefficient.size(); ++i)
        REQUIRE(node->drag_force_coefficient()(i) ==
                Approx(10.).epsilon(Tolerance));

      // Update force to 20.0
      REQUIRE_NOTHROW(
          node->update_drag_force_coefficient(true, drag_force_coefficient));

      for (unsigned i = 0; i < drag_force_coefficient.size(); ++i)
        REQUIRE(node->drag_force_coefficient()(i) ==
                Approx(20.).epsilon(Tolerance));

      // Assign force as 10.0
      REQUIRE_NOTHROW(
          node->update_drag_force_coefficient(false, drag_force_coefficient));
      for (unsigned i = 0; i < drag_force_coefficient.size(); ++i)
        REQUIRE(node->drag_force_coefficient()(i) ==
                Approx(10.).epsilon(Tolerance));
    }

    SECTION("Check compute acceleration and velocity") {
      // Time step
      const double dt = 0.1;

      // Nodal mass
      double solid_mass = 100.;
      double liquid_mass = 100.;
      // Update mass to 100.
      REQUIRE_NOTHROW(
          node->update_mass(false, mpm::NodePhase::nSolid, solid_mass));
      REQUIRE_NOTHROW(
          node->update_mass(false, mpm::NodePhase::nLiquid, liquid_mass));
      REQUIRE(node->mass(mpm::NodePhase::nSolid) ==
              Approx(solid_mass).epsilon(Tolerance));
      REQUIRE(node->mass(mpm::NodePhase::nLiquid) ==
              Approx(liquid_mass).epsilon(Tolerance));

      // Check internal force
      // Create a force vector
      Eigen::Matrix<double, Dim, 1> force;
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 10. * i;
      // Update force to 10.0
      REQUIRE_NOTHROW(
          node->update_internal_force(false, mpm::NodePhase::nMixture, force));
      REQUIRE_NOTHROW(node->update_internal_force(
          false, mpm::NodePhase::nLiquid, 0.5 * force));
      // Internal force
      for (unsigned i = 0; i < force.size(); ++i) {
        REQUIRE(node->internal_force(mpm::NodePhase::nMixture)(i) ==
                Approx(force(i)).epsilon(Tolerance));
        REQUIRE(node->internal_force(mpm::NodePhase::nLiquid)(i) ==
                Approx(0.5 * force(i)).epsilon(Tolerance));
      }

      // External force
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 5. * i;
      // Update force to 10.0
      REQUIRE_NOTHROW(
          node->update_external_force(false, mpm::NodePhase::nMixture, force));
      REQUIRE_NOTHROW(node->update_external_force(
          false, mpm::NodePhase::nLiquid, 0.5 * force));
      for (unsigned i = 0; i < force.size(); ++i) {
        REQUIRE(node->external_force(mpm::NodePhase::nMixture)(i) ==
                Approx(force(i)).epsilon(Tolerance));
        REQUIRE(node->external_force(mpm::NodePhase::nLiquid)(i) ==
                Approx(0.5 * force(i)).epsilon(Tolerance));
      }

      // Drag force
      Eigen::Matrix<double, Dim, 1> drag_force_coefficient;
      for (unsigned i = 0; i < drag_force_coefficient.size(); ++i)
        drag_force_coefficient(i) = 5. * i;
      // Update force to 10.0
      REQUIRE_NOTHROW(
          node->update_drag_force_coefficient(false, drag_force_coefficient));
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->drag_force_coefficient()(i) ==
                Approx(drag_force_coefficient(i)).epsilon(Tolerance));

      REQUIRE(node->compute_acceleration_velocity_twophase_explicit(dt) ==
              true);

      // Check acceleration
      Eigen::Matrix<double, Dim, 1> solid_acceleration;
      solid_acceleration << 0., 0.075, 0.15;
      Eigen::Matrix<double, Dim, 1> liquid_acceleration;
      liquid_acceleration << 0., 0.075, 0.15;

      // Check acceleration
      for (unsigned i = 0; i < solid_acceleration.size(); ++i) {
        REQUIRE(node->acceleration(mpm::NodePhase::nSolid)(i) ==
                Approx(solid_acceleration(i)).epsilon(Tolerance));
        REQUIRE(node->acceleration(mpm::NodePhase::nLiquid)(i) ==
                Approx(liquid_acceleration(i)).epsilon(Tolerance));
      }

      // Check velocity
      Eigen::Matrix<double, Dim, 1> solid_velocity = solid_acceleration * dt;
      Eigen::Matrix<double, Dim, 1> liquid_velocity = liquid_acceleration * dt;
      for (unsigned i = 0; i < solid_velocity.size(); ++i) {
        REQUIRE(node->velocity(mpm::NodePhase::nSolid)(i) ==
                Approx(solid_velocity(i)).epsilon(Tolerance));
        REQUIRE(node->velocity(mpm::NodePhase::nLiquid)(i) ==
                Approx(liquid_velocity(i)).epsilon(Tolerance));
      }

      // Apply velocity constraints
      REQUIRE(node->assign_velocity_constraint(0, 10.5) == true);
      REQUIRE(node->assign_velocity_constraint(1, 0.03) == true);
      REQUIRE(node->assign_velocity_constraint(2, 5.13) == true);
      REQUIRE(node->assign_velocity_constraint(3, 20.5) == true);
      REQUIRE(node->assign_velocity_constraint(4, 1.03) == true);
      REQUIRE(node->assign_velocity_constraint(5, 7.13) == true);
      REQUIRE(node->compute_acceleration_velocity_twophase_explicit(dt) ==
              true);

      // Test velocity with constraints
      solid_velocity << 10.5, 0.03, 5.13;
      liquid_velocity << 20.5, 1.03, 7.13;
      for (unsigned i = 0; i < solid_velocity.size(); ++i) {
        REQUIRE(node->velocity(mpm::NodePhase::nSolid)(i) ==
                Approx(solid_velocity(i)).epsilon(Tolerance));
        REQUIRE(node->velocity(mpm::NodePhase::nLiquid)(i) ==
                Approx(liquid_velocity(i)).epsilon(Tolerance));
      }

      // Test acceleration with constraints
      solid_acceleration.setZero();
      liquid_acceleration.setZero();
      for (unsigned i = 0; i < solid_acceleration.size(); ++i) {
        REQUIRE(node->acceleration(mpm::NodePhase::nSolid)(i) ==
                Approx(solid_acceleration(i)).epsilon(Tolerance));
        REQUIRE(node->acceleration(mpm::NodePhase::nLiquid)(i) ==
                Approx(liquid_acceleration(i)).epsilon(Tolerance));
      }

      // Apply cundall damping when calculating acceleration
      REQUIRE(node->compute_acceleration_velocity_twophase_explicit_cundall(
                  dt, 0.05) == true);

      // Test acceleration with cundall damping
      solid_acceleration << 0., 0., 0.;
      liquid_acceleration << 0., 0., 0.;
      for (unsigned i = 0; i < solid_acceleration.size(); ++i) {
        REQUIRE(node->acceleration(mpm::NodePhase::nSolid)(i) ==
                Approx(solid_acceleration(i)).epsilon(Tolerance));
        REQUIRE(node->acceleration(mpm::NodePhase::nLiquid)(i) ==
                Approx(liquid_acceleration(i)).epsilon(Tolerance));
      }

      // Apply velocity constraints
      REQUIRE(node->assign_velocity_constraint(0, 10.5) == true);
      REQUIRE(node->assign_velocity_constraint(1, 20.5) == true);
      REQUIRE(node->compute_acceleration_velocity_twophase_explicit(dt) ==
              true);

      // Exception check when mass is zero
      REQUIRE_NOTHROW(node->update_mass(false, mpm::NodePhase::nSolid, 0.));
      REQUIRE_NOTHROW(node->update_mass(false, mpm::NodePhase::nLiquid, 0.));
      REQUIRE(node->mass(mpm::NodePhase::nSolid) ==
              Approx(0.).epsilon(Tolerance));
      REQUIRE(node->mass(mpm::NodePhase::nLiquid) ==
              Approx(0.).epsilon(Tolerance));
      REQUIRE(node->compute_acceleration_velocity_twophase_explicit(dt) ==
              false);
    }

    SECTION("Check momentum, velocity and acceleration") {
      // Time step
      const double dt = 0.1;

      // Check momentum
      Eigen::Matrix<double, Dim, 1> momentum;
      for (unsigned i = 0; i < momentum.size(); ++i) momentum(i) = 10.;

      // Check initial momentum
      for (unsigned i = 0; i < momentum.size(); ++i) {
        REQUIRE(node->momentum(mpm::NodePhase::nSolid)(i) ==
                Approx(0.).epsilon(Tolerance));
        REQUIRE(node->momentum(mpm::NodePhase::nLiquid)(i) ==
                Approx(0.).epsilon(Tolerance));
      }

      // Check update momentum to 10
      REQUIRE_NOTHROW(
          node->update_momentum(true, mpm::NodePhase::nSolid, momentum));
      REQUIRE_NOTHROW(
          node->update_momentum(true, mpm::NodePhase::nLiquid, momentum));
      for (unsigned i = 0; i < momentum.size(); ++i) {
        REQUIRE(node->momentum(mpm::NodePhase::nSolid)(i) ==
                Approx(10.).epsilon(Tolerance));
        REQUIRE(node->momentum(mpm::NodePhase::nLiquid)(i) ==
                Approx(10.).epsilon(Tolerance));
      }

      // Check update momentum to 20
      REQUIRE_NOTHROW(
          node->update_momentum(true, mpm::NodePhase::nSolid, momentum));
      REQUIRE_NOTHROW(
          node->update_momentum(true, mpm::NodePhase::nLiquid, momentum));
      for (unsigned i = 0; i < momentum.size(); ++i) {
        REQUIRE(node->momentum(mpm::NodePhase::nSolid)(i) ==
                Approx(20.).epsilon(Tolerance));
        REQUIRE(node->momentum(mpm::NodePhase::nLiquid)(i) ==
                Approx(20.).epsilon(Tolerance));
      }

      // Check assign momentum to 10
      REQUIRE_NOTHROW(
          node->update_momentum(false, mpm::NodePhase::nSolid, momentum));
      REQUIRE_NOTHROW(
          node->update_momentum(false, mpm::NodePhase::nLiquid, momentum));
      for (unsigned i = 0; i < momentum.size(); ++i) {
        REQUIRE(node->momentum(mpm::NodePhase::nSolid)(i) ==
                Approx(10.).epsilon(Tolerance));
        REQUIRE(node->momentum(mpm::NodePhase::nLiquid)(i) ==
                Approx(10.).epsilon(Tolerance));
      }

      // Check zero velocity
      for (unsigned i = 0; i < Dim; ++i) {
        REQUIRE(node->velocity(mpm::NodePhase::nSolid)(i) ==
                Approx(0.).epsilon(Tolerance));
        REQUIRE(node->velocity(mpm::NodePhase::nLiquid)(i) ==
                Approx(0.).epsilon(Tolerance));
      }

      // Check mass
      double mass = 0.;
      REQUIRE_NOTHROW(node->update_mass(false, mpm::NodePhase::nSolid, mass));
      REQUIRE(node->mass(mpm::NodePhase::nSolid) ==
              Approx(0.0).epsilon(Tolerance));
      REQUIRE_NOTHROW(node->update_mass(false, mpm::NodePhase::nLiquid, mass));
      REQUIRE(node->mass(mpm::NodePhase::nLiquid) ==
              Approx(0.0).epsilon(Tolerance));
      // Compute and check velocity this should throw zero mass
      node->compute_velocity();

      mass = 100.;
      // Update mass to 100.5
      REQUIRE_NOTHROW(node->update_mass(true, mpm::NodePhase::nSolid, mass));
      REQUIRE(node->mass(mpm::NodePhase::nSolid) ==
              Approx(100.).epsilon(Tolerance));
      REQUIRE_NOTHROW(node->update_mass(true, mpm::NodePhase::nLiquid, mass));
      REQUIRE(node->mass(mpm::NodePhase::nLiquid) ==
              Approx(100.).epsilon(Tolerance));

      // Check zero velocity
      for (unsigned i = 0; i < Dim; ++i) {
        REQUIRE(node->velocity(mpm::NodePhase::nSolid)(i) ==
                Approx(0.).epsilon(Tolerance));
        REQUIRE(node->velocity(mpm::NodePhase::nLiquid)(i) ==
                Approx(0.).epsilon(Tolerance));
      }

      // Compute and check velocity
      node->compute_velocity();
      for (unsigned i = 0; i < Dim; ++i) {
        REQUIRE(node->velocity(mpm::NodePhase::nSolid)(i) ==
                Approx(0.1).epsilon(Tolerance));
        REQUIRE(node->velocity(mpm::NodePhase::nLiquid)(i) ==
                Approx(0.1).epsilon(Tolerance));
      }

      // Check acceleration
      Eigen::Matrix<double, Dim, 1> acceleration;
      for (unsigned i = 0; i < acceleration.size(); ++i) acceleration(i) = 5.;

      for (unsigned i = 0; i < acceleration.size(); ++i) {
        REQUIRE(node->acceleration(mpm::NodePhase::nSolid)(i) ==
                Approx(0.).epsilon(Tolerance));
        REQUIRE(node->acceleration(mpm::NodePhase::nLiquid)(i) ==
                Approx(0.).epsilon(Tolerance));
      }

      REQUIRE_NOTHROW(node->update_acceleration(true, mpm::NodePhase::nSolid,
                                                acceleration));
      REQUIRE_NOTHROW(node->update_acceleration(true, mpm::NodePhase::nLiquid,
                                                acceleration));
      for (unsigned i = 0; i < acceleration.size(); ++i) {
        REQUIRE(node->acceleration(mpm::NodePhase::nSolid)(i) ==
                Approx(5.).epsilon(Tolerance));
        REQUIRE(node->acceleration(mpm::NodePhase::nLiquid)(i) ==
                Approx(5.).epsilon(Tolerance));
      }

      // Check velocity before constraints
      Eigen::Matrix<double, Dim, 1> velocity;
      velocity << 0.1, 0.1, 0.1;
      for (unsigned i = 0; i < velocity.size(); ++i) {
        REQUIRE(node->velocity(mpm::NodePhase::nSolid)(i) ==
                Approx(velocity(i)).epsilon(Tolerance));
        REQUIRE(node->velocity(mpm::NodePhase::nLiquid)(i) ==
                Approx(velocity(i)).epsilon(Tolerance));
      }

      // Check acceleration before constraints
      acceleration << 5., 5.;
      for (unsigned i = 0; i < acceleration.size(); ++i) {
        REQUIRE(node->acceleration(mpm::NodePhase::nSolid)(i) ==
                Approx(acceleration(i)).epsilon(Tolerance));
        REQUIRE(node->acceleration(mpm::NodePhase::nLiquid)(i) ==
                Approx(acceleration(i)).epsilon(Tolerance));
      }
      SECTION("Check Cartesian velocity constraints") {
        // Apply velocity constraints
        REQUIRE(node->assign_velocity_constraint(0, -12.5) == true);
        // Check out of bounds condition
        REQUIRE(node->assign_velocity_constraint(3, -12.5) == true);

        // Apply constraints
        node->apply_velocity_constraints();

        // Check apply constraints
        velocity << -12.5, 0.1;
        for (unsigned i = 0; i < velocity.size(); ++i) {
          REQUIRE(node->velocity(mpm::NodePhase::nSolid)(i) ==
                  Approx(velocity(i)).epsilon(Tolerance));
          REQUIRE(node->velocity(mpm::NodePhase::nLiquid)(i) ==
                  Approx(velocity(i)).epsilon(Tolerance));
        }

        acceleration << 0., 5.;
        for (unsigned i = 0; i < acceleration.size(); ++i) {
          REQUIRE(node->acceleration(mpm::NodePhase::nSolid)(i) ==
                  Approx(acceleration(i)).epsilon(Tolerance));
          REQUIRE(node->acceleration(mpm::NodePhase::nLiquid)(i) ==
                  Approx(acceleration(i)).epsilon(Tolerance));
        }
      }
      SECTION("Check general velocity constraints in 2 directions") {
        // Apply velocity constraints
        REQUIRE(node->assign_velocity_constraint(0, -12.5) == true);
        REQUIRE(node->assign_velocity_constraint(3, -12.5) == true);

        // Apply rotation matrix with Euler angles alpha = 10 deg, beta
        // = 30 deg
        Eigen::Matrix<double, Dim, 1> euler_angles;
        euler_angles << 10. * M_PI / 180, 30. * M_PI / 180;
        const auto rotation_matrix =
            mpm::geometry::rotation_matrix(euler_angles);
        node->assign_rotation_matrix(rotation_matrix);
        const auto inverse_rotation_matrix = rotation_matrix.inverse();

        // Apply inclined velocity constraints
        node->apply_velocity_constraints();

        // Check apply constraints
        velocity << -9.5834783355, -8.0254030998, 0.1;
        for (unsigned i = 0; i < Dim; ++i) {
          REQUIRE(node->velocity(mpm::NodePhase::nSolid)(i) ==
                  Approx(velocity(i)).epsilon(Tolerance));
          REQUIRE(node->velocity(mpm::NodePhase::nLiquid)(i) ==
                  Approx(velocity(i)).epsilon(Tolerance));
        }

        // Check that the velocity is as specified in local coordinate
        REQUIRE((inverse_rotation_matrix *
                 node->velocity(mpm::NodePhase::nSolid))(0) ==
                Approx(-12.5).epsilon(Tolerance));
        REQUIRE((inverse_rotation_matrix *
                 node->velocity(mpm::NodePhase::nLiquid))(0) ==
                Approx(-12.5).epsilon(Tolerance));

        // Check applied constraints on acceleration in the global coordinates
        acceleration << -0.3961398267, 0.4721010616, 5.0;
        for (unsigned i = 0; i < Dim; ++i) {
          REQUIRE(node->acceleration(mpm::NodePhase::nSolid)(i) ==
                  Approx(acceleration(i)).epsilon(Tolerance));
          REQUIRE(node->acceleration(mpm::NodePhase::nLiquid)(i) ==
                  Approx(acceleration(i)).epsilon(Tolerance));
        }

        // Check that the acceleration is 0 in local coordinate
        REQUIRE((inverse_rotation_matrix *
                 node->acceleration(mpm::NodePhase::nSolid))(0) ==
                Approx(0).epsilon(Tolerance));
        REQUIRE((inverse_rotation_matrix *
                 node->acceleration(mpm::NodePhase::nLiquid))(0) ==
                Approx(0).epsilon(Tolerance));
      }

      SECTION("Check general velocity constraints in all directions") {
        // Apply velocity constraints
        REQUIRE(node->assign_velocity_constraint(0, 10.5) == true);
        REQUIRE(node->assign_velocity_constraint(1, -12.5) == true);
        REQUIRE(node->assign_velocity_constraint(2, 7.5) == true);
        REQUIRE(node->assign_velocity_constraint(3, 10.5) == true);
        REQUIRE(node->assign_velocity_constraint(4, -12.5) == true);
        REQUIRE(node->assign_velocity_constraint(5, 7.5) == true);

        // Apply rotation matrix with Euler angles alpha = -10 deg, beta = 20,
        // deg and gamma = -30 deg
        Eigen::Matrix<double, Dim, 1> euler_angles;
        euler_angles << -10. * M_PI / 180, 20. * M_PI / 180, -30. * M_PI / 180;
        const auto rotation_matrix =
            mpm::geometry::rotation_matrix(euler_angles);
        node->assign_rotation_matrix(rotation_matrix);
        const auto inverse_rotation_matrix = rotation_matrix.inverse();

        // Apply constraints
        node->apply_velocity_constraints();

        // Check apply constraints
        velocity << 13.351984588153375, -5.717804716697730, 10.572663655835457;
        for (unsigned i = 0; i < Dim; ++i) {
          REQUIRE(node->velocity(mpm::NodePhase::nSolid)(i) ==
                  Approx(velocity(i)).epsilon(Tolerance));
          REQUIRE(node->velocity(mpm::NodePhase::nLiquid)(i) ==
                  Approx(velocity(i)).epsilon(Tolerance));
        }

        // Check that the velocity is as specified in local coordinate
        REQUIRE((inverse_rotation_matrix *
                 node->velocity(mpm::NodePhase::nSolid))(0) ==
                Approx(10.5).epsilon(Tolerance));
        REQUIRE((inverse_rotation_matrix *
                 node->velocity(mpm::NodePhase::nSolid))(1) ==
                Approx(-12.5).epsilon(Tolerance));
        REQUIRE((inverse_rotation_matrix *
                 node->velocity(mpm::NodePhase::nSolid))(2) ==
                Approx(7.5).epsilon(Tolerance));
        REQUIRE((inverse_rotation_matrix *
                 node->velocity(mpm::NodePhase::nLiquid))(0) ==
                Approx(10.5).epsilon(Tolerance));
        REQUIRE((inverse_rotation_matrix *
                 node->velocity(mpm::NodePhase::nLiquid))(1) ==
                Approx(-12.5).epsilon(Tolerance));
        REQUIRE((inverse_rotation_matrix *
                 node->velocity(mpm::NodePhase::nLiquid))(2) ==
                Approx(7.5).epsilon(Tolerance));

        // Check apply constraints
        acceleration << 0, 0, 0;
        for (unsigned i = 0; i < Dim; ++i) {
          REQUIRE(node->acceleration(mpm::NodePhase::nSolid)(i) ==
                  Approx(acceleration(i)).epsilon(Tolerance));
          REQUIRE(node->acceleration(mpm::NodePhase::nLiquid)(i) ==
                  Approx(acceleration(i)).epsilon(Tolerance));
        }

        // Check that the acceleration is 0 in local coordinate
        REQUIRE((inverse_rotation_matrix *
                 node->acceleration(mpm::NodePhase::nSolid))(0) ==
                Approx(0).epsilon(Tolerance));
        REQUIRE((inverse_rotation_matrix *
                 node->acceleration(mpm::NodePhase::nSolid))(1) ==
                Approx(0).epsilon(Tolerance));
        REQUIRE((inverse_rotation_matrix *
                 node->acceleration(mpm::NodePhase::nLiquid))(0) ==
                Approx(0).epsilon(Tolerance));
        REQUIRE((inverse_rotation_matrix *
                 node->acceleration(mpm::NodePhase::nLiquid))(1) ==
                Approx(0).epsilon(Tolerance));
      }

      SECTION("Check Cartesian friction constraints") {
        // Apply friction constraints
        REQUIRE(node->assign_friction_constraint(2, 2, 0.3) == true);
        // Check out of bounds condition
        REQUIRE(node->assign_friction_constraint(4, 1, 0.2) == false);

        // Apply constraints
        node->apply_friction_constraints(dt);

        // Check apply constraints
        acceleration << 3.939339828220179, 3.939339828220179, 5.;
        for (unsigned i = 0; i < acceleration.size(); ++i)
          REQUIRE(node->acceleration(mpm::NodePhase::nSolid)(i) ==
                  Approx(acceleration(i)).epsilon(Tolerance));
      }

      SECTION("Check general friction constraints in 1 direction") {
        // Apply friction constraints
        REQUIRE(node->assign_friction_constraint(2, 2, 0.3) == true);

        // Apply rotation matrix with Euler angles alpha = 10 deg, beta = 20 deg
        // and gamma = 30 deg
        Eigen::Matrix<double, Dim, 1> euler_angles;
        euler_angles << 10. * M_PI / 180, 20. * M_PI / 180, 30. * M_PI / 180;
        const auto rotation_matrix =
            mpm::geometry::rotation_matrix(euler_angles);
        node->assign_rotation_matrix(rotation_matrix);
        const auto inverse_rotation_matrix = rotation_matrix.inverse();

        // Apply inclined velocity constraints
        node->apply_friction_constraints(dt);

        // Check applied constraints on acceleration in the global coordinates
        acceleration << 4.602895052828914, 4.492575657560740, 4.751301246937935;
        for (unsigned i = 0; i < Dim; ++i)
          REQUIRE(node->acceleration(mpm::NodePhase::nSolid)(i) ==
                  Approx(acceleration(i)).epsilon(Tolerance));

        // Check the acceleration in local coordinate
        acceleration << 6.878925666702865, 3.365244416454818, 2.302228080558999;
        for (unsigned i = 0; i < Dim; ++i)
          REQUIRE((inverse_rotation_matrix *
                   node->acceleration(mpm::NodePhase::nSolid))(i) ==
                  Approx(acceleration(i)).epsilon(Tolerance));
      }
    }

    SECTION("Check node material ids") {
      // Add material to nodes
      node->append_material_id(0);
      node->append_material_id(1);
      node->append_material_id(4);
      node->append_material_id(0);
      node->append_material_id(2);

      // Check size of material_ids
      REQUIRE(node->material_ids().size() == 4);

      // Check elements of material_ids
      std::vector<unsigned> material_ids = {0, 1, 2, 4};
      auto mat_ids = node->material_ids();
      unsigned i = 0;
      for (auto mitr = mat_ids.begin(); mitr != mat_ids.end(); ++mitr, ++i)
        REQUIRE(*mitr == material_ids.at(i));
    }
  }
}

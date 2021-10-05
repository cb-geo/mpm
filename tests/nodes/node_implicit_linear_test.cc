#include <cmath>
#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

#include "function_base.h"
#include "geometry.h"
#include "node.h"

// Check node class for 1D case
TEST_CASE("Implicit Linear Node is checked for 1D case",
          "[node][1D][ImplicitLinear]") {
  const unsigned Dim = 1;
  const unsigned Dof = 1;
  const unsigned Nphases = 1;
  const unsigned Nphase = 0;
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
    REQUIRE(node->mass(Nphase) == Approx(0.0).epsilon(Tolerance));
    double mass = 100.5;
    // Update mass to 100.5
    REQUIRE_NOTHROW(node->update_mass(true, Nphase, mass));
    REQUIRE(node->mass(Nphase) == Approx(100.5).epsilon(Tolerance));
    // Update mass to 201
    REQUIRE_NOTHROW(node->update_mass(true, Nphase, mass));
    REQUIRE(node->mass(Nphase) == Approx(201.0).epsilon(Tolerance));
    // Assign mass to 100
    mass = 100.;
    REQUIRE_NOTHROW(node->update_mass(false, Nphase, mass));
    REQUIRE(node->mass(Nphase) == Approx(100.0).epsilon(Tolerance));

    SECTION("Check nodal pressure") {
      // Check pressure
      REQUIRE(node->pressure(Nphase) == Approx(0.0).epsilon(Tolerance));
      double pressure = 1000.7;
      // Update pressure to 1000.7
      REQUIRE_NOTHROW(node->update_mass_pressure(Nphase, mass * pressure));
      REQUIRE(node->pressure(Nphase) == Approx(1000.7).epsilon(Tolerance));
      // Update pressure to 2001.4
      REQUIRE_NOTHROW(node->update_mass_pressure(Nphase, mass * pressure));
      REQUIRE(node->pressure(Nphase) == Approx(2001.4).epsilon(Tolerance));
      // Assign pressure to 1000
      pressure = 1000.;
      node->assign_pressure(Nphase, pressure);
      REQUIRE(node->pressure(Nphase) == Approx(1000.0).epsilon(Tolerance));
      // Assign mass to 0
      mass = 0.;
      REQUIRE_NOTHROW(node->update_mass(false, Nphase, mass));
      REQUIRE(node->mass(Nphase) == Approx(0.0).epsilon(Tolerance));
      // Try to update pressure to 2000, should throw and keep to 1000.
      node->assign_pressure(Nphase, pressure);
      REQUIRE(node->pressure(Nphase) == Approx(1000.0).epsilon(Tolerance));
      // Check pressure constraints
      SECTION("Check nodal pressure constraints") {
        // Check assign pressure constraint
        REQUIRE(node->assign_pressure_constraint(mpm::NodePhase::NSolid, 8000,
                                                 nullptr) == true);
        // Check apply pressure constraint
        REQUIRE_NOTHROW(
            node->apply_pressure_constraint(mpm::NodePhase::NSolid));
        // Check pressure
        REQUIRE(node->pressure(mpm::NodePhase::NSolid) ==
                Approx(8000).epsilon(Tolerance));
      }
    }

    SECTION("Check external force") {
      // Create a force vector
      Eigen::Matrix<double, Dim, 1> force;
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 10.;

      // Check current external force is zero
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(0.).epsilon(Tolerance));

      // Update force to 10.0
      REQUIRE_NOTHROW(node->update_external_force(true, Nphase, force));
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));

      // Update force to 20.0
      REQUIRE_NOTHROW(node->update_external_force(true, Nphase, force));
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(20.).epsilon(Tolerance));

      // Assign force as 10.0
      REQUIRE_NOTHROW(node->update_external_force(false, Nphase, force));
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));

      SECTION("Check concentrated force") {
        // Set external force to zero
        force.setZero();
        REQUIRE_NOTHROW(node->update_external_force(false, Nphase, force));

        // concentrated force
        std::shared_ptr<mpm::FunctionBase> ffunction = nullptr;
        double concentrated_force = 65.32;
        const unsigned Direction = 0;
        // Check external force
        for (unsigned i = 0; i < Dim; ++i)
          REQUIRE(node->external_force(Nphase)(i) ==
                  Approx(0.).epsilon(Tolerance));

        REQUIRE(node->assign_concentrated_force(
                    Nphase, Direction, concentrated_force, ffunction) == true);

        double current_time = 0.0;
        node->apply_concentrated_force(Nphase, current_time);

        for (unsigned i = 0; i < Dim; ++i) {
          if (i == Direction)
            REQUIRE(node->external_force(Nphase)(i) ==
                    Approx(concentrated_force).epsilon(Tolerance));
          else
            REQUIRE(node->external_force(Nphase)(i) ==
                    Approx(0.).epsilon(Tolerance));
        }

        // Check for incorrect direction / phase
        const unsigned wrong_dir = 4;
        REQUIRE(node->assign_concentrated_force(
                    Nphase, wrong_dir, concentrated_force, ffunction) == false);

        // Check again to ensure value hasn't been updated
        for (unsigned i = 0; i < Dim; ++i) {
          if (i == Direction)
            REQUIRE(node->external_force(Nphase)(i) ==
                    Approx(concentrated_force).epsilon(Tolerance));
          else
            REQUIRE(node->external_force(Nphase)(i) ==
                    Approx(0.).epsilon(Tolerance));
        }
      }
    }

    SECTION("Check internal force") {
      // Create a force vector
      Eigen::Matrix<double, Dim, 1> force;
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 10.;

      // Check current internal force is zero
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(0.).epsilon(Tolerance));

      // Update force to 10.0
      REQUIRE_NOTHROW(node->update_internal_force(true, Nphase, force));
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));

      // Update force to 20.0
      REQUIRE_NOTHROW(node->update_internal_force(true, Nphase, force));
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(20.).epsilon(Tolerance));

      // Assign force as 10.0
      REQUIRE_NOTHROW(node->update_internal_force(false, Nphase, force));
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));
    }

    SECTION("Check update of velocity and acceleration using Newmark scheme") {
      // Time step
      const double dt = 0.1;
      // Parameters of Newmark scheme
      const double newmark_beta = 0.25;
      const double newmark_gamma = 0.50;

      // Initialization
      node->initialise_implicit();

      // Nodal displacement
      // clang-format off
      Eigen::VectorXd nodal_displacement;
      nodal_displacement.resize(Dim);
      nodal_displacement << 0.1;
      // clang-format on
      // Assign nodal active_id
      node->assign_active_id(0);
      // Assign nodal displacement
      node->update_displacement_increment(nodal_displacement, Nphase, 1);
      // Check nodal displacement
      REQUIRE(node->displacement(Nphase)(0) ==
              Approx(nodal_displacement(0)).epsilon(Tolerance));

      // Nodal velocity
      Eigen::VectorXd velocity;
      velocity.resize(Dim);
      for (unsigned i = 0; i < velocity.size(); ++i) velocity(i) = i;
      for (unsigned i = 0; i < Dim; ++i)
        node->assign_velocity_constraint(i, velocity(i));
      node->apply_velocity_constraints();

      // Update nodal velocity and acceleration using Newmark scheme
      node->update_velocity_acceleration_newmark(Nphase, newmark_beta,
                                                 newmark_gamma, dt);

      // Values of nodal velocity
      Eigen::Matrix<double, 1, 1> nodal_velocity;
      // clang-format off
      nodal_velocity << 2.;
      // clang-format on
      // Check nodal velocity
      for (unsigned i = 0; i < nodal_velocity.cols(); ++i)
        REQUIRE(node->velocity(Nphase)(i) ==
                Approx(nodal_velocity(0, i)).epsilon(Tolerance));

      // Values of nodal acceleration
      Eigen::Matrix<double, 1, 1> nodal_acceleration;
      // clang-format off
      nodal_acceleration << 40.;
      // clang-format on
      // Check nodal acceleration
      for (unsigned i = 0; i < nodal_acceleration.cols(); ++i)
        REQUIRE(node->acceleration(Nphase)(i) ==
                Approx(nodal_acceleration(0, i)).epsilon(Tolerance));

      // Exception check when mass is zero
      mass = 0.;
      // Update mass to 0.
      REQUIRE_NOTHROW(node->update_mass(false, Nphase, mass));
      REQUIRE(node->mass(Nphase) == Approx(mass).epsilon(Tolerance));
      REQUIRE(node->compute_acceleration_velocity(Nphase, dt) == false);
    }

    SECTION("Check momentum, inertia, velocity and acceleration") {
      // Initialization
      node->initialise_implicit();

      // Check momentum
      Eigen::Matrix<double, Dim, 1> momentum;
      for (unsigned i = 0; i < momentum.size(); ++i) momentum(i) = 10.;

      // Check initial momentum
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(0.).epsilon(Tolerance));

      // Check update momentum to 10
      REQUIRE_NOTHROW(node->update_momentum(true, Nphase, momentum));
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(10.).epsilon(Tolerance));

      // Check update momentum to 20
      REQUIRE_NOTHROW(node->update_momentum(true, Nphase, momentum));
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(20.).epsilon(Tolerance));

      // Check assign momentum to 10
      REQUIRE_NOTHROW(node->update_momentum(false, Nphase, momentum));
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(10.).epsilon(Tolerance));

      // Check inertia
      Eigen::Matrix<double, Dim, 1> inertia;
      for (unsigned i = 0; i < inertia.size(); ++i) inertia(i) = 10.;

      // Check initial inertia
      for (unsigned i = 0; i < inertia.size(); ++i)
        REQUIRE(node->inertia(Nphase)(i) == Approx(0.).epsilon(Tolerance));

      // Check update inertia to 10
      REQUIRE_NOTHROW(node->update_inertia(true, Nphase, inertia));
      for (unsigned i = 0; i < inertia.size(); ++i)
        REQUIRE(node->inertia(Nphase)(i) == Approx(10.).epsilon(Tolerance));

      // Check update inertia to 20
      REQUIRE_NOTHROW(node->update_inertia(true, Nphase, inertia));
      for (unsigned i = 0; i < inertia.size(); ++i)
        REQUIRE(node->inertia(Nphase)(i) == Approx(20.).epsilon(Tolerance));

      // Check assign inertia to 10
      REQUIRE_NOTHROW(node->update_inertia(false, Nphase, inertia));
      for (unsigned i = 0; i < inertia.size(); ++i)
        REQUIRE(node->inertia(Nphase)(i) == Approx(10.).epsilon(Tolerance));

      // Check zero velocity
      for (unsigned i = 0; i < Dim; ++i)
        REQUIRE(node->velocity(Nphase)(i) == Approx(0.).epsilon(Tolerance));

      // Check zero acceleration
      for (unsigned i = 0; i < Dim; ++i)
        REQUIRE(node->acceleration(Nphase)(i) == Approx(0.).epsilon(Tolerance));

      // Check mass
      double mass = 0.;
      REQUIRE_NOTHROW(node->update_mass(false, Nphase, mass));
      REQUIRE(node->mass(Nphase) == Approx(0.0).epsilon(Tolerance));
      // Compute and check velocity and acceleration this should throw zero mass
      node->compute_velocity_acceleration();

      mass = 100.;
      // Update mass to 100.5
      REQUIRE_NOTHROW(node->update_mass(true, Nphase, mass));
      REQUIRE(node->mass(Nphase) == Approx(100.).epsilon(Tolerance));

      // Compute and check velocity and acceleration
      node->compute_velocity_acceleration();
      for (unsigned i = 0; i < Dim; ++i)
        REQUIRE(node->velocity(Nphase)(i) == Approx(0.1).epsilon(Tolerance));
      for (unsigned i = 0; i < Dim; ++i)
        REQUIRE(node->acceleration(Nphase)(i) ==
                Approx(0.1).epsilon(Tolerance));
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
TEST_CASE("Implicit Linear Node is checked for 2D case",
          "[node][2D][ImplicitLinear]") {
  const unsigned Dim = 2;
  const unsigned Dof = 2;
  const unsigned Nphases = 1;
  const unsigned Nphase = 0;
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
    REQUIRE(node->mass(Nphase) == Approx(0.0).epsilon(Tolerance));
    double mass = 100.5;
    // Update mass to 100.5
    REQUIRE_NOTHROW(node->update_mass(true, Nphase, mass));
    REQUIRE(node->mass(Nphase) == Approx(100.5).epsilon(Tolerance));
    // Update mass to 201
    REQUIRE_NOTHROW(node->update_mass(true, Nphase, mass));
    REQUIRE(node->mass(Nphase) == Approx(201.0).epsilon(Tolerance));
    // Assign mass to 100
    mass = 100.;
    REQUIRE_NOTHROW(node->update_mass(false, Nphase, mass));
    REQUIRE(node->mass(Nphase) == Approx(100.0).epsilon(Tolerance));

    SECTION("Check nodal pressure") {
      // Check pressure
      REQUIRE(node->pressure(Nphase) == Approx(0.0).epsilon(Tolerance));
      double pressure = 1000.7;
      // Update pressure to 1000.7
      REQUIRE_NOTHROW(node->update_mass_pressure(Nphase, mass * pressure));
      REQUIRE(node->pressure(Nphase) == Approx(1000.7).epsilon(Tolerance));
      // Update pressure to 2001.4
      REQUIRE_NOTHROW(node->update_mass_pressure(Nphase, mass * pressure));
      REQUIRE(node->pressure(Nphase) == Approx(2001.4).epsilon(Tolerance));
      // Assign pressure to 1000
      pressure = 1000.;
      node->assign_pressure(Nphase, pressure);
      REQUIRE(node->pressure(Nphase) == Approx(1000.0).epsilon(Tolerance));
      // Assign mass to 0
      mass = 0.;
      REQUIRE_NOTHROW(node->update_mass(false, Nphase, mass));
      REQUIRE(node->mass(Nphase) == Approx(0.0).epsilon(Tolerance));
      // Try to update pressure to 2000, should throw and keep to 1000.
      node->assign_pressure(Nphase, pressure);
      REQUIRE(node->pressure(Nphase) == Approx(1000.0).epsilon(Tolerance));
      // Check pressure constraints
      SECTION("Check nodal pressure constraints") {
        // Check assign pressure constraint
        REQUIRE(node->assign_pressure_constraint(mpm::NodePhase::NSolid, 8000,
                                                 nullptr) == true);
        // Check apply pressure constraint
        REQUIRE_NOTHROW(
            node->apply_pressure_constraint(mpm::NodePhase::NSolid));
        // Check pressure
        REQUIRE(node->pressure(mpm::NodePhase::NSolid) ==
                Approx(8000).epsilon(Tolerance));
      }
    }

    SECTION("Check volume") {
      // Check volume
      REQUIRE(node->volume(Nphase) == Approx(0.0).epsilon(Tolerance));
      double volume = 100.5;
      // Update volume to 100.5
      REQUIRE_NOTHROW(node->update_volume(true, Nphase, volume));
      REQUIRE(node->volume(Nphase) == Approx(100.5).epsilon(Tolerance));
      // Update volume to 201
      REQUIRE_NOTHROW(node->update_volume(true, Nphase, volume));
      REQUIRE(node->volume(Nphase) == Approx(201.0).epsilon(Tolerance));
      // Assign volume to 100
      volume = 100.;
      REQUIRE_NOTHROW(node->update_volume(false, Nphase, volume));
      REQUIRE(node->volume(Nphase) == Approx(100.0).epsilon(Tolerance));
    }

    SECTION("Check external force") {
      // Create a force vector
      Eigen::Matrix<double, Dim, 1> force;
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 10.;

      // Check current external force is zero
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(0.).epsilon(Tolerance));

      // Update force to 10.0
      REQUIRE_NOTHROW(node->update_external_force(true, Nphase, force));
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));

      // Update force to 20.0
      REQUIRE_NOTHROW(node->update_external_force(true, Nphase, force));
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(20.).epsilon(Tolerance));

      // Assign force as 10.0
      REQUIRE_NOTHROW(node->update_external_force(false, Nphase, force));
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));

      // Check if exception is handled
      Eigen::Matrix<double, Dim, 1> force_bad;
      for (unsigned i = 0; i < force_bad.size(); ++i) force_bad(i) = 10.;

      // Exception handling invalid force dimension
      // TODO Assert:    REQUIRE_NOTHROW(node->update_external_force(true, 1,
      // force_bad));

      // Exception handling invalid force dimension
      // TODO Assert: REQUIRE_NOTHROW(node->update_external_force(false, 1,
      // force_bad));

      SECTION("Check concentrated force") {
        // Set external force to zero
        force.setZero();
        REQUIRE_NOTHROW(node->update_external_force(false, Nphase, force));

        // Concentrated force
        std::shared_ptr<mpm::FunctionBase> ffunction = nullptr;
        double concentrated_force = 65.32;
        const unsigned Direction = 0;
        // Check traction
        for (unsigned i = 0; i < Dim; ++i)
          REQUIRE(node->external_force(Nphase)(i) ==
                  Approx(0.).epsilon(Tolerance));

        REQUIRE(node->assign_concentrated_force(
                    Nphase, Direction, concentrated_force, ffunction) == true);
        double current_time = 0.0;
        node->apply_concentrated_force(Nphase, current_time);

        for (unsigned i = 0; i < Dim; ++i) {
          if (i == Direction)
            REQUIRE(node->external_force(Nphase)(i) ==
                    Approx(concentrated_force).epsilon(Tolerance));

          else
            REQUIRE(node->external_force(Nphase)(i) ==
                    Approx(0.).epsilon(Tolerance));
        }

        // Check for incorrect direction / phase
        const unsigned wrong_dir = 4;
        REQUIRE(node->assign_concentrated_force(
                    Nphase, wrong_dir, concentrated_force, ffunction) == false);

        // Check again to ensure value hasn't been updated
        for (unsigned i = 0; i < Dim; ++i) {
          if (i == Direction)
            REQUIRE(node->external_force(Nphase)(i) ==
                    Approx(concentrated_force).epsilon(Tolerance));
          else
            REQUIRE(node->external_force(Nphase)(i) ==
                    Approx(0.).epsilon(Tolerance));
        }
      }
    }

    SECTION("Check internal force") {
      // Create a force vector
      Eigen::Matrix<double, Dim, 1> force;
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 10.;

      // Check current internal force is zero
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(0.).epsilon(Tolerance));

      // Update force to 10.0
      REQUIRE_NOTHROW(node->update_internal_force(true, Nphase, force));
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));

      // Update force to 20.0
      REQUIRE_NOTHROW(node->update_internal_force(true, Nphase, force));
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(20.).epsilon(Tolerance));

      // Assign force as 10.0
      REQUIRE_NOTHROW(node->update_internal_force(false, Nphase, force));
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));
    }

    SECTION("Check update of velocity and acceleration using Newmark scheme") {
      // Time step
      const double dt = 0.1;
      // Parameters of Newmark scheme
      const double newmark_beta = 0.25;
      const double newmark_gamma = 0.50;

      // Initialization
      node->initialise_implicit();

      // Nodal displacement
      // clang-format off
      Eigen::VectorXd nodal_displacement;
      nodal_displacement.resize(Dim);
      nodal_displacement << 0.,
                            0.1;
      // clang-format on
      // Assign nodal active_id
      node->assign_active_id(0);
      // Assign nodal displacement
      node->update_displacement_increment(nodal_displacement, Nphase, 1);
      // Check nodal displacement
      for (unsigned i = 0; i < nodal_displacement.size(); ++i)
        REQUIRE(node->displacement(Nphase)(i) ==
                Approx(nodal_displacement(i)).epsilon(Tolerance));

      // Nodal velocity
      Eigen::VectorXd velocity;
      velocity.resize(Dim);
      for (unsigned i = 0; i < velocity.size(); ++i) velocity(i) = i;
      for (unsigned i = 0; i < Dim; ++i)
        node->assign_velocity_constraint(i, velocity(i));
      node->apply_velocity_constraints();

      // Update nodal velocity and acceleration using Newmark scheme
      node->update_velocity_acceleration_newmark(Nphase, newmark_beta,
                                                 newmark_gamma, dt);

      // Values of nodal velocity
      Eigen::Matrix<double, 2, 1> nodal_velocity;
      // clang-format off
      nodal_velocity << 0.,
                        1.;
      // clang-format on
      // Check nodal velocity
      for (unsigned i = 0; i < nodal_velocity.cols(); ++i)
        REQUIRE(node->velocity(Nphase)(i) ==
                Approx(nodal_velocity(0, i)).epsilon(Tolerance));

      // Values of nodal acceleration
      Eigen::Matrix<double, 2, 1> nodal_acceleration;
      // clang-format off
      nodal_acceleration << 0.,
                            0.;
      // clang-format on
      // Check nodal acceleration
      for (unsigned i = 0; i < nodal_acceleration.cols(); ++i)
        REQUIRE(node->acceleration(Nphase)(i) ==
                Approx(nodal_acceleration(0, i)).epsilon(Tolerance));

      // Exception check when mass is zero
      mass = 0.;
      // Update mass to 0.
      REQUIRE_NOTHROW(node->update_mass(false, Nphase, mass));
      REQUIRE(node->mass(Nphase) == Approx(mass).epsilon(Tolerance));
      REQUIRE(node->compute_acceleration_velocity(Nphase, dt) == false);
    }

    SECTION("Check momentum, inertia, velocity and acceleration") {
      // Time step
      const double dt = 0.1;
      // Initialization
      node->initialise_implicit();

      // Check momentum
      Eigen::Matrix<double, Dim, 1> momentum;
      for (unsigned i = 0; i < momentum.size(); ++i) momentum(i) = 10.;

      // Check initial momentum
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(0.).epsilon(Tolerance));

      // Check update momentum to 10
      REQUIRE_NOTHROW(node->update_momentum(true, Nphase, momentum));
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(10.).epsilon(Tolerance));

      // Check update momentum to 20
      REQUIRE_NOTHROW(node->update_momentum(true, Nphase, momentum));
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(20.).epsilon(Tolerance));

      // Check assign momentum to 10
      REQUIRE_NOTHROW(node->update_momentum(false, Nphase, momentum));
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(10.).epsilon(Tolerance));

      // Check inertia
      Eigen::Matrix<double, Dim, 1> inertia;
      for (unsigned i = 0; i < inertia.size(); ++i) inertia(i) = 10.;

      // Check initial inertia
      for (unsigned i = 0; i < inertia.size(); ++i)
        REQUIRE(node->inertia(Nphase)(i) == Approx(0.).epsilon(Tolerance));

      // Check update inertia to 10
      REQUIRE_NOTHROW(node->update_inertia(true, Nphase, inertia));
      for (unsigned i = 0; i < inertia.size(); ++i)
        REQUIRE(node->inertia(Nphase)(i) == Approx(10.).epsilon(Tolerance));

      // Check update inertia to 20
      REQUIRE_NOTHROW(node->update_inertia(true, Nphase, inertia));
      for (unsigned i = 0; i < inertia.size(); ++i)
        REQUIRE(node->inertia(Nphase)(i) == Approx(20.).epsilon(Tolerance));

      // Check assign inertia to 10
      REQUIRE_NOTHROW(node->update_inertia(false, Nphase, inertia));
      for (unsigned i = 0; i < inertia.size(); ++i)
        REQUIRE(node->inertia(Nphase)(i) == Approx(10.).epsilon(Tolerance));

      // Check zero velocity
      for (unsigned i = 0; i < Dim; ++i)
        REQUIRE(node->velocity(Nphase)(i) == Approx(0.).epsilon(Tolerance));

      // Check zero acceleration
      for (unsigned i = 0; i < Dim; ++i)
        REQUIRE(node->acceleration(Nphase)(i) == Approx(0.).epsilon(Tolerance));

      // Check mass
      double mass = 0.;
      REQUIRE_NOTHROW(node->update_mass(false, Nphase, mass));
      REQUIRE(node->mass(Nphase) == Approx(0.0).epsilon(Tolerance));
      // Compute and check velocity this should throw zero mass
      node->compute_velocity();

      mass = 100.;
      // Update mass to 100.5
      REQUIRE_NOTHROW(node->update_mass(true, Nphase, mass));
      REQUIRE(node->mass(Nphase) == Approx(100.).epsilon(Tolerance));

      // Compute and check velocity and acceleration
      node->compute_velocity_acceleration();
      for (unsigned i = 0; i < Dim; ++i)
        REQUIRE(node->velocity(Nphase)(i) == Approx(0.1).epsilon(Tolerance));
      for (unsigned i = 0; i < Dim; ++i)
        REQUIRE(node->acceleration(Nphase)(i) ==
                Approx(0.1).epsilon(Tolerance));

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
}

// \brief Check node class for 3D case
TEST_CASE("Implicit Linear Node is checked for 3D case",
          "[node][3D][ImplicitLinear]") {
  const unsigned Dim = 3;
  const unsigned Dof = 3;
  const unsigned Nphases = 1;
  const unsigned Nphase = 0;

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
    REQUIRE(node->mass(Nphase) == Approx(0.0).epsilon(Tolerance));
    double mass = 100.5;
    // Update mass to 100.5
    REQUIRE_NOTHROW(node->update_mass(true, Nphase, mass));
    REQUIRE(node->mass(Nphase) == Approx(100.5).epsilon(Tolerance));
    // Update mass to 201
    REQUIRE_NOTHROW(node->update_mass(true, Nphase, mass));
    REQUIRE(node->mass(Nphase) == Approx(201.0).epsilon(Tolerance));
    // Assign mass to 100
    mass = 100.;
    REQUIRE_NOTHROW(node->update_mass(false, Nphase, mass));
    REQUIRE(node->mass(Nphase) == Approx(100.0).epsilon(Tolerance));

    SECTION("Check nodal pressure") {
      // Check pressure
      REQUIRE(node->pressure(Nphase) == Approx(0.0).epsilon(Tolerance));
      double pressure = 1000.7;
      // Update pressure to 1000.7
      REQUIRE_NOTHROW(node->update_mass_pressure(Nphase, mass * pressure));
      REQUIRE(node->pressure(Nphase) == Approx(1000.7).epsilon(Tolerance));
      // Update pressure to 2001.4
      REQUIRE_NOTHROW(node->update_mass_pressure(Nphase, mass * pressure));
      REQUIRE(node->pressure(Nphase) == Approx(2001.4).epsilon(Tolerance));
      // Assign pressure to 1000
      pressure = 1000.;
      node->assign_pressure(Nphase, pressure);
      REQUIRE(node->pressure(Nphase) == Approx(1000.0).epsilon(Tolerance));
      // Assign mass to 0
      mass = 0.;
      REQUIRE_NOTHROW(node->update_mass(false, Nphase, mass));
      REQUIRE(node->mass(Nphase) == Approx(0.0).epsilon(Tolerance));
      // Try to update pressure to 2000, should throw and keep to 1000.
      node->assign_pressure(Nphase, pressure);
      REQUIRE(node->pressure(Nphase) == Approx(1000.0).epsilon(Tolerance));
      // Check pressure constraints
      SECTION("Check nodal pressure constraints") {
        // Check assign pressure constraint
        REQUIRE(node->assign_pressure_constraint(mpm::NodePhase::NSolid, 8000,
                                                 nullptr) == true);
        // Check apply pressure constraint
        REQUIRE_NOTHROW(
            node->apply_pressure_constraint(mpm::NodePhase::NSolid));
        // Check pressure
        REQUIRE(node->pressure(mpm::NodePhase::NSolid) ==
                Approx(8000).epsilon(Tolerance));
      }
    }

    SECTION("Check external force") {
      // Create a force vector
      Eigen::Matrix<double, Dim, 1> force;
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 10.;

      // Check current external force
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(0.).epsilon(Tolerance));

      // Update force to 10.0
      REQUIRE_NOTHROW(node->update_external_force(true, Nphase, force));
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));

      // Update force to 2.0
      REQUIRE_NOTHROW(node->update_external_force(true, Nphase, force));
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(20.).epsilon(Tolerance));

      // Assign force as 1.0
      REQUIRE_NOTHROW(node->update_external_force(false, Nphase, force));
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));

      SECTION("Check concentrated force") {
        // Set external force to zero
        force.setZero();
        REQUIRE_NOTHROW(node->update_external_force(false, Nphase, force));

        // Concentrated froce
        std::shared_ptr<mpm::FunctionBase> ffunction = nullptr;
        double concentrated_force = 65.32;
        const unsigned Direction = 0;
        // Check external force
        for (unsigned i = 0; i < Dim; ++i)
          REQUIRE(node->external_force(Nphase)(i) ==
                  Approx(0.).epsilon(Tolerance));

        REQUIRE(node->assign_concentrated_force(
                    Nphase, Direction, concentrated_force, ffunction) == true);
        double current_time = 0.0;
        node->apply_concentrated_force(Nphase, current_time);

        for (unsigned i = 0; i < Dim; ++i) {
          if (i == Direction)
            REQUIRE(node->external_force(Nphase)(i) ==
                    Approx(concentrated_force).epsilon(Tolerance));
          else
            REQUIRE(node->external_force(Nphase)(i) ==
                    Approx(0.).epsilon(Tolerance));
        }

        // Check for incorrect direction / phase
        const unsigned wrong_dir = 4;
        REQUIRE(node->assign_concentrated_force(
                    Nphase, wrong_dir, concentrated_force, ffunction) == false);

        // Check again to ensure value hasn't been updated
        for (unsigned i = 0; i < Dim; ++i) {
          if (i == Direction)
            REQUIRE(node->external_force(Nphase)(i) ==
                    Approx(concentrated_force).epsilon(Tolerance));
          else
            REQUIRE(node->external_force(Nphase)(i) ==
                    Approx(0.).epsilon(Tolerance));
        }
      }
    }

    SECTION("Check internal force") {
      // Create a force vector
      Eigen::Matrix<double, Dim, 1> force;
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 10.;

      // Check current internal force is zero
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(0.).epsilon(Tolerance));

      // Update force to 10.0
      REQUIRE_NOTHROW(node->update_internal_force(true, Nphase, force));
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));

      // Update force to 20.0
      REQUIRE_NOTHROW(node->update_internal_force(true, Nphase, force));
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(20.).epsilon(Tolerance));

      // Assign force as 10.0
      REQUIRE_NOTHROW(node->update_internal_force(false, Nphase, force));
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));
    }

    SECTION("Check update of velocity and acceleration using Newmark scheme") {
      // Time step
      const double dt = 0.1;
      // Parameters of Newmark scheme
      const double newmark_beta = 0.25;
      const double newmark_gamma = 0.50;

      // Initialization
      node->initialise_implicit();

      // Nodal displacement
      // clang-format off
      Eigen::VectorXd nodal_displacement;
      nodal_displacement.resize(Dim);
      nodal_displacement << 0.,
                            0.1,
                            0.2;
      // clang-format on
      // Assign nodal active_id
      node->assign_active_id(0);
      // Assign nodal displacement
      node->update_displacement_increment(nodal_displacement, Nphase, 1);
      // Check nodal displacement
      for (unsigned i = 0; i < nodal_displacement.size(); ++i)
        REQUIRE(node->displacement(Nphase)(i) ==
                Approx(nodal_displacement(i)).epsilon(Tolerance));

      // Nodal velocity
      Eigen::VectorXd velocity;
      velocity.resize(Dim);
      for (unsigned i = 0; i < velocity.size(); ++i) velocity(i) = i;
      for (unsigned i = 0; i < Dim; ++i)
        node->assign_velocity_constraint(i, velocity(i));
      node->apply_velocity_constraints();

      // Update nodal velocity and acceleration using Newmark scheme
      node->update_velocity_acceleration_newmark(Nphase, newmark_beta,
                                                 newmark_gamma, dt);

      // Values of nodal velocity
      Eigen::Matrix<double, 3, 1> nodal_velocity;
      // clang-format off
      nodal_velocity << 0.,
                        1.,
                        2.;
      // clang-format on
      // Check nodal velocity
      for (unsigned i = 0; i < nodal_velocity.cols(); ++i)
        REQUIRE(node->velocity(Nphase)(i) ==
                Approx(nodal_velocity(0, i)).epsilon(Tolerance));

      // Values of nodal acceleration
      Eigen::Matrix<double, 3, 1> nodal_acceleration;
      // clang-format off
      nodal_acceleration << 0.,
                            0.,
                            0.;
      // clang-format on
      // Check nodal acceleration
      for (unsigned i = 0; i < nodal_acceleration.cols(); ++i)
        REQUIRE(node->acceleration(Nphase)(i) ==
                Approx(nodal_acceleration(0, i)).epsilon(Tolerance));

      // Exception check when mass is zero
      mass = 0.;
      // Update mass to 0.
      REQUIRE_NOTHROW(node->update_mass(false, Nphase, mass));
      REQUIRE(node->mass(Nphase) == Approx(mass).epsilon(Tolerance));
      REQUIRE(node->compute_acceleration_velocity(Nphase, dt) == false);
    }

    SECTION("Check momentum, inertia, velocity and acceleration") {
      // Time step
      const double dt = 0.1;
      // Initialization
      node->initialise_implicit();

      // Check momentum
      Eigen::Matrix<double, Dim, 1> momentum;
      for (unsigned i = 0; i < momentum.size(); ++i) momentum(i) = 10.;

      // Check initial momentum
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(0.).epsilon(Tolerance));

      // Check update momentum to 10
      REQUIRE_NOTHROW(node->update_momentum(true, Nphase, momentum));
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(10.).epsilon(Tolerance));

      // Check update momentum to 20
      REQUIRE_NOTHROW(node->update_momentum(true, Nphase, momentum));
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(20.).epsilon(Tolerance));

      // Check assign momentum to 10
      REQUIRE_NOTHROW(node->update_momentum(false, Nphase, momentum));
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(10.).epsilon(Tolerance));

      // Check inertia
      Eigen::Matrix<double, Dim, 1> inertia;
      for (unsigned i = 0; i < inertia.size(); ++i) inertia(i) = 10.;

      // Check initial inertia
      for (unsigned i = 0; i < inertia.size(); ++i)
        REQUIRE(node->inertia(Nphase)(i) == Approx(0.).epsilon(Tolerance));

      // Check update inertia to 10
      REQUIRE_NOTHROW(node->update_inertia(true, Nphase, inertia));
      for (unsigned i = 0; i < inertia.size(); ++i)
        REQUIRE(node->inertia(Nphase)(i) == Approx(10.).epsilon(Tolerance));

      // Check update inertia to 20
      REQUIRE_NOTHROW(node->update_inertia(true, Nphase, inertia));
      for (unsigned i = 0; i < inertia.size(); ++i)
        REQUIRE(node->inertia(Nphase)(i) == Approx(20.).epsilon(Tolerance));

      // Check assign inertia to 10
      REQUIRE_NOTHROW(node->update_inertia(false, Nphase, inertia));
      for (unsigned i = 0; i < inertia.size(); ++i)
        REQUIRE(node->inertia(Nphase)(i) == Approx(10.).epsilon(Tolerance));

      // Check mass
      double mass = 0.;
      REQUIRE_NOTHROW(node->update_mass(false, Nphase, mass));
      REQUIRE(node->mass(Nphase) == Approx(0.0).epsilon(Tolerance));
      // Compute and check velocity this should throw zero mass
      node->compute_velocity();

      mass = 100.;
      // Update mass to 100.5
      REQUIRE_NOTHROW(node->update_mass(true, Nphase, mass));
      REQUIRE(node->mass(Nphase) == Approx(100.).epsilon(Tolerance));

      // Check zero velocity
      for (unsigned i = 0; i < Dim; ++i)
        REQUIRE(node->velocity(Nphase)(i) == Approx(0.).epsilon(Tolerance));
      // Check zero acceleration
      for (unsigned i = 0; i < Dim; ++i)
        REQUIRE(node->acceleration(Nphase)(i) == Approx(0.).epsilon(Tolerance));

      // Compute and check velocity and acceleration
      node->compute_velocity_acceleration();
      for (unsigned i = 0; i < Dim; ++i)
        REQUIRE(node->velocity(Nphase)(i) == Approx(0.1).epsilon(Tolerance));
      for (unsigned i = 0; i < Dim; ++i)
        REQUIRE(node->acceleration(Nphase)(i) ==
                Approx(0.1).epsilon(Tolerance));
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

#include <cmath>
#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

#include "geometry.h"
#include "node.h"

// Check node class for 1D case
TEST_CASE("Node is checked for 1D case", "[node][1D]") {
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
    node->update_mass(true, Nphase, mass);
    REQUIRE(node->mass(Nphase) == Approx(100.5).epsilon(Tolerance));
    // Update mass to 201
    node->update_mass(true, Nphase, mass);
    REQUIRE(node->mass(Nphase) == Approx(201.0).epsilon(Tolerance));
    // Assign mass to 100
    mass = 100.;
    node->update_mass(false, Nphase, mass);
    REQUIRE(node->mass(Nphase) == Approx(100.0).epsilon(Tolerance));

    SECTION("Check external force") {
      // Create a force vector
      Eigen::Matrix<double, Dim, 1> force;
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 10.;

      // Check current external force is zero
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(0.).epsilon(Tolerance));

      // Update force to 10.0
      REQUIRE(node->update_external_force(true, Nphase, force) == true);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));

      // Update force to 20.0
      REQUIRE(node->update_external_force(true, Nphase, force) == true);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(20.).epsilon(Tolerance));

      // Assign force as 10.0
      REQUIRE(node->update_external_force(false, Nphase, force) == true);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));

      // Check if exception is handled
      unsigned bad_phase = 4;
      // Exception handling invalid force dimension
      REQUIRE(node->update_external_force(true, bad_phase, force) == false);
      // Exception handling invalid force dimension
      REQUIRE(node->update_external_force(false, bad_phase, force) == false);

      SECTION("Check traction") {
        // External force
        force.setZero();
        REQUIRE(node->update_external_force(false, Nphase, force) == true);

        // Traction
        double traction = 65.32;
        const unsigned Direction = 0;
        // Check traction
        for (unsigned i = 0; i < Dim; ++i)
          REQUIRE(node->external_force(Nphase)(i) ==
                  Approx(0.).epsilon(Tolerance));

        REQUIRE(node->assign_traction_force(Nphase, Direction, traction) ==
                true);

        for (unsigned i = 0; i < Dim; ++i) {
          if (i == Direction)
            REQUIRE(node->external_force(Nphase)(i) ==
                    Approx(traction).epsilon(Tolerance));
          else
            REQUIRE(node->external_force(Nphase)(i) ==
                    Approx(0.).epsilon(Tolerance));
        }

        // Check for incorrect direction / phase
        const unsigned wrong_dir = 4;
        REQUIRE(node->assign_traction_force(Nphase, wrong_dir, traction) ==
                false);

        // Check again to ensure value hasn't been updated
        for (unsigned i = 0; i < Dim; ++i) {
          if (i == Direction)
            REQUIRE(node->external_force(Nphase)(i) ==
                    Approx(traction).epsilon(Tolerance));
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
      REQUIRE(node->update_internal_force(true, Nphase, force) == true);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));

      // Update force to 20.0
      REQUIRE(node->update_internal_force(true, Nphase, force) == true);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(20.).epsilon(Tolerance));

      // Assign force as 10.0
      REQUIRE(node->update_internal_force(false, Nphase, force) == true);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));

      // Check if exception is handled
      unsigned bad_phase = 4;
      // Exception handling invalid force dimension
      REQUIRE(node->update_internal_force(true, bad_phase, force) == false);
      // Exception handling invalid force dimension
      REQUIRE(node->update_internal_force(false, bad_phase, force) == false);
    }

    SECTION("Check compute acceleration and velocity") {
      // Time step
      const double dt = 0.1;

      // Nodal mass
      double mass = 100.;
      // Update mass to 100.5
      node->update_mass(false, Nphase, mass);
      REQUIRE(node->mass(Nphase) == Approx(mass).epsilon(Tolerance));

      // Check internal force
      // Create a force vector
      Eigen::Matrix<double, Dim, 1> force;
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 10. * i;
      // Update force to 10.0
      node->update_internal_force(false, Nphase, force);
      // Internal force
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(force(i)).epsilon(Tolerance));

      // External force
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 5. * i;
      // Update force to 10.0
      node->update_external_force(false, Nphase, force);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(force(i)).epsilon(Tolerance));

      REQUIRE(node->compute_acceleration_velocity(Nphase, dt) == true);

      for (unsigned i = 0; i < force.size(); ++i) force(i) = 15. * i;

      // Check acceleration
      Eigen::Matrix<double, Dim, 1> acceleration = force / mass;
      for (unsigned i = 0; i < acceleration.size(); ++i)
        REQUIRE(node->acceleration(Nphase)(i) ==
                Approx(acceleration(i)).epsilon(Tolerance));

      // Check velocity
      Eigen::Matrix<double, Dim, 1> velocity = force / mass * dt;
      for (unsigned i = 0; i < velocity.size(); ++i)
        REQUIRE(node->velocity(Nphase)(i) ==
                Approx(velocity(i)).epsilon(Tolerance));

      // Apply velocity constraints
      REQUIRE(node->assign_velocity_constraint(0, 10.5) == true);
      REQUIRE(node->compute_acceleration_velocity(Nphase, dt) == true);

      // Test velocity with constraints
      velocity[0] = 10.5;
      for (unsigned i = 0; i < velocity.size(); ++i)
        REQUIRE(node->velocity(Nphase)(i) ==
                Approx(velocity(i)).epsilon(Tolerance));

      // Test acceleration with constraints
      acceleration[0] = 0.;
      for (unsigned i = 0; i < acceleration.size(); ++i)
        REQUIRE(node->acceleration(Nphase)(i) ==
                Approx(acceleration(i)).epsilon(Tolerance));

      // Exception check when mass is zero
      mass = 0.;
      // Update mass to 0.
      node->update_mass(false, Nphase, mass);
      REQUIRE(node->mass(Nphase) == Approx(mass).epsilon(Tolerance));
      REQUIRE(node->compute_acceleration_velocity(Nphase, dt) == false);
    }

    SECTION("Check momentum and velocity") {
      // Check momentum
      Eigen::Matrix<double, Dim, 1> momentum;
      for (unsigned i = 0; i < momentum.size(); ++i) momentum(i) = 10.;

      // Check initial momentum
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(0.).epsilon(Tolerance));

      // Check update momentum to 10
      REQUIRE(node->update_momentum(true, Nphase, momentum) == true);
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(10.).epsilon(Tolerance));

      // Check update momentum to 20
      REQUIRE(node->update_momentum(true, Nphase, momentum) == true);
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(20.).epsilon(Tolerance));

      // Check assign momentum to 10
      REQUIRE(node->update_momentum(false, Nphase, momentum) == true);
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(10.).epsilon(Tolerance));

      // Check zero velocity
      for (unsigned i = 0; i < Dim; ++i)
        REQUIRE(node->velocity(Nphase)(i) == Approx(0.).epsilon(Tolerance));

      // Check mass
      double mass = 0.;
      node->update_mass(false, Nphase, mass);
      REQUIRE(node->mass(Nphase) == Approx(0.0).epsilon(Tolerance));
      // Compute and check velocity this should throw zero mass
      node->compute_velocity();

      mass = 100.;
      // Update mass to 100.5
      node->update_mass(true, Nphase, mass);
      REQUIRE(node->mass(Nphase) == Approx(100.).epsilon(Tolerance));

      // Compute and check velocity
      node->compute_velocity();
      for (unsigned i = 0; i < Dim; ++i)
        REQUIRE(node->velocity(Nphase)(i) == Approx(0.1).epsilon(Tolerance));

      // Check if exception is handled
      unsigned bad_phase = 1;
      // Exception handling invalid momentum dimension
      REQUIRE(node->update_momentum(true, bad_phase, momentum) == false);
      // Exception handling invalid momentum dimension
      REQUIRE(node->update_momentum(false, bad_phase, momentum) == false);

      // Apply velocity constraints
      REQUIRE(node->assign_velocity_constraint(0, 10.5) == true);
      // Check out of bounds condition
      REQUIRE(node->assign_velocity_constraint(1, 0) == false);

      // Check velocity before constraints
      Eigen::Matrix<double, Dim, 1> velocity;
      velocity << 0.1;
      for (unsigned i = 0; i < velocity.size(); ++i)
        REQUIRE(node->velocity(Nphase)(i) ==
                Approx(velocity(i)).epsilon(Tolerance));

      // Apply constraints
      node->apply_velocity_constraints();

      // Check apply constraints
      velocity << 10.5;
      for (unsigned i = 0; i < velocity.size(); ++i)
        REQUIRE(node->velocity(Nphase)(i) ==
                Approx(velocity(i)).epsilon(Tolerance));
    }

    SECTION("Check acceleration") {
      // Check acceleration
      Eigen::Matrix<double, Dim, 1> acceleration;
      for (unsigned i = 0; i < acceleration.size(); ++i) acceleration(i) = 5.;

      for (unsigned i = 0; i < acceleration.size(); ++i)
        REQUIRE(node->acceleration(Nphase)(i) == Approx(0.).epsilon(Tolerance));

      REQUIRE(node->update_acceleration(true, Nphase, acceleration) == true);
      for (unsigned i = 0; i < acceleration.size(); ++i)
        REQUIRE(node->acceleration(Nphase)(i) == Approx(5.).epsilon(Tolerance));

      // Check if exception is handled
      unsigned bad_phase = 1;
      // Exception handling invalid acceleration dimension
      REQUIRE(node->update_acceleration(true, bad_phase, acceleration) ==
              false);

      // Exception handling invalid acceleration dimension
      REQUIRE(node->update_acceleration(false, bad_phase, acceleration) ==
              false);

      // Apply velocity constraints
      REQUIRE(node->assign_velocity_constraint(0, 10.5) == true);
      // Check out of bounds condition
      REQUIRE(node->assign_velocity_constraint(1, 12.5) == false);

      // Check acceleration before constraints
      acceleration.resize(Dim);
      acceleration << 5.;
      for (unsigned i = 0; i < acceleration.size(); ++i)
        REQUIRE(node->acceleration(Nphase)(i) ==
                Approx(acceleration(i)).epsilon(Tolerance));

      // Apply constraints
      node->apply_velocity_constraints();

      // Check apply constraints
      acceleration << 0.0;
      for (unsigned i = 0; i < acceleration.size(); ++i)
        REQUIRE(node->acceleration(Nphase)(i) ==
                Approx(acceleration(i)).epsilon(Tolerance));
    }
  }
}

// \brief Check node class for 2D case
TEST_CASE("Node is checked for 2D case", "[node][2D]") {
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
    node->update_mass(true, Nphase, mass);
    REQUIRE(node->mass(Nphase) == Approx(100.5).epsilon(Tolerance));
    // Update mass to 201
    node->update_mass(true, Nphase, mass);
    REQUIRE(node->mass(Nphase) == Approx(201.0).epsilon(Tolerance));
    // Assign mass to 100
    mass = 100.;
    node->update_mass(false, Nphase, mass);
    REQUIRE(node->mass(Nphase) == Approx(100.0).epsilon(Tolerance));

    SECTION("Check volume") {
      // Check volume
      REQUIRE(node->volume(Nphase) == Approx(0.0).epsilon(Tolerance));
      double volume = 100.5;
      // Update volume to 100.5
      node->update_volume(true, Nphase, volume);
      REQUIRE(node->volume(Nphase) == Approx(100.5).epsilon(Tolerance));
      // Update volume to 201
      node->update_volume(true, Nphase, volume);
      REQUIRE(node->volume(Nphase) == Approx(201.0).epsilon(Tolerance));
      // Assign volume to 100
      volume = 100.;
      node->update_volume(false, Nphase, volume);
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
      REQUIRE(node->update_external_force(true, Nphase, force) == true);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));

      // Update force to 20.0
      REQUIRE(node->update_external_force(true, Nphase, force) == true);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(20.).epsilon(Tolerance));

      // Assign force as 10.0
      REQUIRE(node->update_external_force(false, Nphase, force) == true);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));

      // Check if exception is handled
      Eigen::Matrix<double, Dim, 1> force_bad;
      for (unsigned i = 0; i < force_bad.size(); ++i) force_bad(i) = 10.;

      // Exception handling invalid force dimension
      REQUIRE(node->update_external_force(true, 1, force_bad) == false);
      // Exception handling invalid force dimension
      REQUIRE(node->update_external_force(false, 1, force_bad) == false);

      SECTION("Check traction") {
        // External force
        force.setZero();
        REQUIRE(node->update_external_force(false, Nphase, force) == true);

        // Traction
        double traction = 65.32;
        const unsigned Direction = 0;
        // Check traction
        for (unsigned i = 0; i < Dim; ++i)
          REQUIRE(node->external_force(Nphase)(i) ==
                  Approx(0.).epsilon(Tolerance));

        REQUIRE(node->assign_traction_force(Nphase, Direction, traction) ==
                true);

        for (unsigned i = 0; i < Dim; ++i) {
          if (i == Direction)
            REQUIRE(node->external_force(Nphase)(i) ==
                    Approx(traction).epsilon(Tolerance));
          else
            REQUIRE(node->external_force(Nphase)(i) ==
                    Approx(0.).epsilon(Tolerance));
        }

        // Check for incorrect direction / phase
        const unsigned wrong_dir = 4;
        REQUIRE(node->assign_traction_force(Nphase, wrong_dir, traction) ==
                false);

        // Check again to ensure value hasn't been updated
        for (unsigned i = 0; i < Dim; ++i) {
          if (i == Direction)
            REQUIRE(node->external_force(Nphase)(i) ==
                    Approx(traction).epsilon(Tolerance));
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
      REQUIRE(node->update_internal_force(true, Nphase, force) == true);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));

      // Update force to 20.0
      REQUIRE(node->update_internal_force(true, Nphase, force) == true);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(20.).epsilon(Tolerance));

      // Assign force as 10.0
      REQUIRE(node->update_internal_force(false, Nphase, force) == true);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));

      // Check if exception is handled
      unsigned bad_phase = 1;
      // Exception handling invalid force dimension
      REQUIRE(node->update_internal_force(true, bad_phase, force) == false);
      // Exception handling invalid force dimension
      REQUIRE(node->update_internal_force(false, bad_phase, force) == false);
    }

    SECTION("Check compute acceleration and velocity") {
      // Time step
      const double dt = 0.1;

      // Nodal mass
      double mass = 100.;
      // Update mass to 100.5
      node->update_mass(false, Nphase, mass);
      REQUIRE(node->mass(Nphase) == Approx(mass).epsilon(Tolerance));

      // Check internal force
      // Create a force vector
      Eigen::Matrix<double, Dim, 1> force;
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 10. * i;
      // Update force to 10.0
      node->update_internal_force(false, Nphase, force);
      // Internal force
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(force(i)).epsilon(Tolerance));

      // External force
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 5. * i;
      // Update force to 10.0
      node->update_external_force(false, Nphase, force);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(force(i)).epsilon(Tolerance));

      REQUIRE(node->compute_acceleration_velocity(Nphase, dt) == true);

      for (unsigned i = 0; i < force.size(); ++i) force(i) = 15. * i;

      // Check acceleration
      Eigen::Matrix<double, Dim, 1> acceleration = force / mass;
      for (unsigned i = 0; i < acceleration.size(); ++i)
        REQUIRE(node->acceleration(Nphase)(i) ==
                Approx(acceleration(i)).epsilon(Tolerance));

      // Check velocity
      Eigen::Matrix<double, Dim, 1> velocity = force / mass * dt;
      for (unsigned i = 0; i < velocity.size(); ++i)
        REQUIRE(node->velocity(Nphase)(i) ==
                Approx(velocity(i)).epsilon(Tolerance));

      // Apply velocity constraints
      REQUIRE(node->assign_velocity_constraint(0, 10.5) == true);
      REQUIRE(node->compute_acceleration_velocity(Nphase, dt) == true);

      // Test velocity with constraints
      velocity << 10.5, 0.03;
      for (unsigned i = 0; i < velocity.size(); ++i)
        REQUIRE(node->velocity(Nphase)(i) ==
                Approx(velocity(i)).epsilon(Tolerance));

      // Test acceleration with constraints
      acceleration[0] = 0.;
      for (unsigned i = 0; i < acceleration.size(); ++i)
        REQUIRE(node->acceleration(Nphase)(i) ==
                Approx(acceleration(i)).epsilon(Tolerance));

      // Exception check when mass is zero
      mass = 0.;
      // Update mass to 0.
      node->update_mass(false, Nphase, mass);
      REQUIRE(node->mass(Nphase) == Approx(mass).epsilon(Tolerance));
      REQUIRE(node->compute_acceleration_velocity(Nphase, dt) == false);
    }

    SECTION("Check momentum and velocity") {
      // Check momentum
      Eigen::Matrix<double, Dim, 1> momentum;
      for (unsigned i = 0; i < momentum.size(); ++i) momentum(i) = 10.;

      // Check initial momentum
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(0.).epsilon(Tolerance));

      // Check update momentum to 10
      REQUIRE(node->update_momentum(true, Nphase, momentum) == true);
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(10.).epsilon(Tolerance));

      // Check update momentum to 20
      REQUIRE(node->update_momentum(true, Nphase, momentum) == true);
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(20.).epsilon(Tolerance));

      // Check assign momentum to 10
      REQUIRE(node->update_momentum(false, Nphase, momentum) == true);
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(10.).epsilon(Tolerance));

      // Check zero velocity
      for (unsigned i = 0; i < Dim; ++i)
        REQUIRE(node->velocity(Nphase)(i) == Approx(0.).epsilon(Tolerance));

      // Check mass
      double mass = 0.;
      node->update_mass(false, Nphase, mass);
      REQUIRE(node->mass(Nphase) == Approx(0.0).epsilon(Tolerance));
      // Compute and check velocity this should throw zero mass
      node->compute_velocity();

      mass = 100.;
      // Update mass to 100.5
      node->update_mass(true, Nphase, mass);
      REQUIRE(node->mass(Nphase) == Approx(100.).epsilon(Tolerance));

      // Compute and check velocity
      node->compute_velocity();
      for (unsigned i = 0; i < Dim; ++i)
        REQUIRE(node->velocity(Nphase)(i) == Approx(0.1).epsilon(Tolerance));

      // Check if exception is handled
      unsigned bad_phase = 1;
      // Exception handling invalid momentum dimension
      REQUIRE(node->update_momentum(true, bad_phase, momentum) == false);
      // Exception handling invalid momentum dimension
      REQUIRE(node->update_momentum(false, bad_phase, momentum) == false);

      // Apply velocity constraints
      REQUIRE(node->assign_velocity_constraint(0, -12.5) == true);
      // Check out of bounds condition
      REQUIRE(node->assign_velocity_constraint(2, 0.) == false);

      // Check velocity before constraints
      Eigen::Matrix<double, Dim, 1> velocity;
      velocity << 0.1, 0.1;
      for (unsigned i = 0; i < velocity.size(); ++i)
        REQUIRE(node->velocity(Nphase)(i) ==
                Approx(velocity(i)).epsilon(Tolerance));

      // Apply constraints
      node->apply_velocity_constraints();

      // Check apply constraints
      velocity << -12.5, 0.1;
      for (unsigned i = 0; i < velocity.size(); ++i)
        REQUIRE(node->velocity(Nphase)(i) ==
                Approx(velocity(i)).epsilon(Tolerance));
    }

    SECTION("Check acceleration") {
      // Check acceleration
      Eigen::Matrix<double, Dim, 1> acceleration;
      for (unsigned i = 0; i < acceleration.size(); ++i) acceleration(i) = 5.;

      for (unsigned i = 0; i < acceleration.size(); ++i)
        REQUIRE(node->acceleration(Nphase)(i) == Approx(0.).epsilon(Tolerance));

      REQUIRE(node->update_acceleration(true, Nphase, acceleration) == true);
      for (unsigned i = 0; i < acceleration.size(); ++i)
        REQUIRE(node->acceleration(Nphase)(i) == Approx(5.).epsilon(Tolerance));

      // Check if exception is handled
      Eigen::Matrix<double, Dim, 1> acceleration_bad;
      for (unsigned i = 0; i < acceleration_bad.size(); ++i)
        acceleration_bad(i) = 10.;

      unsigned bad_phase = 3;
      // Exception handling invalid acceleration dimension
      REQUIRE(node->update_acceleration(true, bad_phase, acceleration_bad) ==
              false);

      // Apply velocity constraints
      REQUIRE(node->assign_velocity_constraint(0, -12.5) == true);
      // Check out of bounds condition
      REQUIRE(node->assign_velocity_constraint(2, 0.) == false);

      // Check acceleration before constraints
      acceleration << 5., 5.;
      for (unsigned i = 0; i < acceleration.size(); ++i)
        REQUIRE(node->acceleration(Nphase)(i) ==
                Approx(acceleration(i)).epsilon(Tolerance));

      // Apply constraints
      node->apply_velocity_constraints();

      // Check apply constraints
      acceleration << 0., 5.;
      for (unsigned i = 0; i < acceleration.size(); ++i)
        REQUIRE(node->acceleration(Nphase)(i) ==
                Approx(acceleration(i)).epsilon(Tolerance));
    }

    SECTION("Check general velocity constraints") {
      // Check momentum
      Eigen::Matrix<double, Dim, 1> momentum;
      for (unsigned i = 0; i < momentum.size(); ++i) momentum(i) = 10.;

      // Check initial momentum
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(0.).epsilon(Tolerance));

      // Check assign momentum to 10
      REQUIRE(node->update_momentum(false, Nphase, momentum) == true);
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(10.).epsilon(Tolerance));

      // Check zero velocity
      for (unsigned i = 0; i < Dim; ++i)
        REQUIRE(node->velocity(Nphase)(i) == Approx(0.).epsilon(Tolerance));

      // Check mass
      mass = 100.;
      node->update_mass(false, Nphase, mass);
      REQUIRE(node->mass(Nphase) == Approx(100.).epsilon(Tolerance));

      // Compute and check velocity
      node->compute_velocity();
      for (unsigned i = 0; i < Dim; ++i)
        REQUIRE(node->velocity(Nphase)(i) == Approx(0.1).epsilon(Tolerance));

      // Check acceleration
      Eigen::Matrix<double, Dim, 1> acceleration;
      for (unsigned i = 0; i < acceleration.size(); ++i) acceleration(i) = 5.;
      REQUIRE(node->update_acceleration(false, Nphase, acceleration) == true);
      for (unsigned i = 0; i < acceleration.size(); ++i)
        REQUIRE(node->acceleration(Nphase)(i) == Approx(5.).epsilon(Tolerance));

      // Check velocity before constraints
      Eigen::Matrix<double, Dim, 1> velocity;
      velocity << 0.1, 0.1;
      for (unsigned i = 0; i < velocity.size(); ++i)
        REQUIRE(node->velocity(Nphase)(i) ==
                Approx(velocity(i)).epsilon(Tolerance));

      // Check acceleration before constraints
      acceleration << 5., 5.;
      for (unsigned i = 0; i < acceleration.size(); ++i)
        REQUIRE(node->acceleration(Nphase)(i) ==
                Approx(acceleration(i)).epsilon(Tolerance));

      // Apply velocity constraints
      REQUIRE(node->assign_velocity_constraint(0, -12.5) == true);
      // Check out of bounds condition
      REQUIRE(node->assign_velocity_constraint(2, 0.) == false);

      // Apply rotation matrix with Euler angles alpha = 10 deg, beta = 30 deg
      Eigen::Matrix<double, Dim, 1> euler_angles;
      euler_angles << 10. * M_PI / 180, 30. * M_PI / 180;
      const auto rotation_matrix = mpm::geometry::rotation_matrix(euler_angles);
      node->assign_rotation_matrix(rotation_matrix);
      const auto inverse_rotation_matrix = rotation_matrix.inverse();

      // Apply constraints
      node->apply_velocity_constraints();

      // Check apply constraints
      velocity << -9.583478335521184, -8.025403099849004;
      for (unsigned i = 0; i < Dim; ++i)
        REQUIRE(node->velocity(Nphase)(i) ==
                Approx(velocity(i)).epsilon(Tolerance));

      // Check that the velocity in 0 direction is -12.5 in local coordinate
      REQUIRE((inverse_rotation_matrix * node->velocity(Nphase))(0) ==
              Approx(-12.5).epsilon(Tolerance));

      // Check apply constraints
      acceleration << -0.396139826697847, 0.472101061636807;
      for (unsigned i = 0; i < Dim; ++i)
        REQUIRE(node->acceleration(Nphase)(i) ==
                Approx(acceleration(i)).epsilon(Tolerance));

      // Check that the acceleration in 0 direction is 0 in local coordinate
      REQUIRE((inverse_rotation_matrix * node->acceleration(Nphase))(0) ==
              Approx(0).epsilon(Tolerance));
    }
  }
}

// \brief Check node class for 3D case
TEST_CASE("Node is checked for 3D case", "[node][3D]") {
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
    node->update_mass(true, Nphase, mass);
    REQUIRE(node->mass(Nphase) == Approx(100.5).epsilon(Tolerance));
    // Update mass to 201
    node->update_mass(true, Nphase, mass);
    REQUIRE(node->mass(Nphase) == Approx(201.0).epsilon(Tolerance));
    // Assign mass to 100
    mass = 100.;
    node->update_mass(false, Nphase, mass);
    REQUIRE(node->mass(Nphase) == Approx(100.0).epsilon(Tolerance));

    SECTION("Check external force") {
      // Create a force vector
      Eigen::Matrix<double, Dim, 1> force;
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 10.;

      // Check current external force
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(0.).epsilon(Tolerance));

      // Update force to 10.0
      node->update_external_force(true, Nphase, force);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));

      // Update force to 2.0
      node->update_external_force(true, Nphase, force);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(20.).epsilon(Tolerance));

      // Assign force as 1.0
      node->update_external_force(false, Nphase, force);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));

      SECTION("Check traction") {
        // External force
        force.setZero();
        REQUIRE(node->update_external_force(false, Nphase, force) == true);

        // Traction
        double traction = 65.32;
        const unsigned Direction = 0;
        // Check traction
        for (unsigned i = 0; i < Dim; ++i)
          REQUIRE(node->external_force(Nphase)(i) ==
                  Approx(0.).epsilon(Tolerance));

        REQUIRE(node->assign_traction_force(Nphase, Direction, traction) ==
                true);

        for (unsigned i = 0; i < Dim; ++i) {
          if (i == Direction)
            REQUIRE(node->external_force(Nphase)(i) ==
                    Approx(traction).epsilon(Tolerance));
          else
            REQUIRE(node->external_force(Nphase)(i) ==
                    Approx(0.).epsilon(Tolerance));
        }

        // Check for incorrect direction / phase
        const unsigned wrong_dir = 4;
        REQUIRE(node->assign_traction_force(Nphase, wrong_dir, traction) ==
                false);

        // Check again to ensure value hasn't been updated
        for (unsigned i = 0; i < Dim; ++i) {
          if (i == Direction)
            REQUIRE(node->external_force(Nphase)(i) ==
                    Approx(traction).epsilon(Tolerance));
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
      node->update_internal_force(true, Nphase, force);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));

      // Update force to 20.0
      node->update_internal_force(true, Nphase, force);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(20.).epsilon(Tolerance));

      // Assign force as 10.0
      node->update_internal_force(false, Nphase, force);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));
    }

    SECTION("Check compute acceleration and velocity") {
      // Time step
      const double dt = 0.1;

      // Nodal mass
      double mass = 100.;
      // Update mass to 100.5
      node->update_mass(false, Nphase, mass);
      REQUIRE(node->mass(Nphase) == Approx(mass).epsilon(Tolerance));

      // Check internal force
      // Create a force vector
      Eigen::Matrix<double, Dim, 1> force;
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 10. * i;
      // Update force to 10.0
      node->update_internal_force(false, Nphase, force);
      // Internal force
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(force(i)).epsilon(Tolerance));

      // External force
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 5. * i;
      // Update force to 10.0
      node->update_external_force(false, Nphase, force);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(force(i)).epsilon(Tolerance));

      REQUIRE(node->compute_acceleration_velocity(Nphase, dt) == true);

      for (unsigned i = 0; i < force.size(); ++i) force(i) = 15. * i;

      // Check acceleration
      Eigen::Matrix<double, Dim, 1> acceleration = force / mass;
      for (unsigned i = 0; i < acceleration.size(); ++i)
        REQUIRE(node->acceleration(Nphase)(i) ==
                Approx(acceleration(i)).epsilon(Tolerance));

      // Check velocity
      Eigen::Matrix<double, Dim, 1> velocity = force / mass * dt;
      for (unsigned i = 0; i < velocity.size(); ++i)
        REQUIRE(node->velocity(Nphase)(i) ==
                Approx(velocity(i)).epsilon(Tolerance));

      // Apply velocity constraints
      REQUIRE(node->assign_velocity_constraint(0, 10.5) == true);
      REQUIRE(node->compute_acceleration_velocity(Nphase, dt) == true);

      // Test velocity with constraints
      velocity << 10.5, 0.03, 0.06;
      for (unsigned i = 0; i < velocity.size(); ++i)
        REQUIRE(node->velocity(Nphase)(i) ==
                Approx(velocity(i)).epsilon(Tolerance));

      // Test acceleration with constraints
      acceleration[0] = 0.;
      for (unsigned i = 0; i < acceleration.size(); ++i)
        REQUIRE(node->acceleration(Nphase)(i) ==
                Approx(acceleration(i)).epsilon(Tolerance));

      // Exception check when mass is zero
      mass = 0.;
      // Update mass to 0.
      node->update_mass(false, Nphase, mass);
      REQUIRE(node->mass(Nphase) == Approx(mass).epsilon(Tolerance));
      REQUIRE(node->compute_acceleration_velocity(Nphase, dt) == false);
    }

    SECTION("Check momentum and velocity") {
      // Check momentum
      Eigen::Matrix<double, Dim, 1> momentum;
      for (unsigned i = 0; i < momentum.size(); ++i) momentum(i) = 10.;

      // Check initial momentum
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(0.).epsilon(Tolerance));

      // Check update momentum to 10
      node->update_momentum(true, Nphase, momentum);
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(10.).epsilon(Tolerance));

      // Check update momentum to 20
      node->update_momentum(true, Nphase, momentum);
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(20.).epsilon(Tolerance));

      // Check assign momentum to 10
      node->update_momentum(false, Nphase, momentum);
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(10.).epsilon(Tolerance));

      // Check mass
      double mass = 0.;
      node->update_mass(false, Nphase, mass);
      REQUIRE(node->mass(Nphase) == Approx(0.0).epsilon(Tolerance));
      // Compute and check velocity this should throw zero mass
      node->compute_velocity();

      mass = 100.;
      // Update mass to 100.5
      node->update_mass(true, Nphase, mass);
      REQUIRE(node->mass(Nphase) == Approx(100.).epsilon(Tolerance));

      // Check zero velocity
      for (unsigned i = 0; i < Dim; ++i)
        REQUIRE(node->velocity(Nphase)(i) == Approx(0.).epsilon(Tolerance));

      // Compute and check velocity
      node->compute_velocity();
      for (unsigned i = 0; i < Dim; ++i)
        REQUIRE(node->velocity(Nphase)(i) == Approx(0.1).epsilon(Tolerance));

      // Apply velocity constraints
      REQUIRE(node->assign_velocity_constraint(0, 10.5) == true);
      REQUIRE(node->assign_velocity_constraint(1, -12.5) == true);
      // Check out of bounds condition
      REQUIRE(node->assign_velocity_constraint(4, 0.) == false);

      // Check velocity before constraints
      Eigen::Matrix<double, Dim, 1> velocity;
      velocity << 0.1, 0.1, 0.1;
      for (unsigned i = 0; i < velocity.size(); ++i)
        REQUIRE(node->velocity(Nphase)(i) ==
                Approx(velocity(i)).epsilon(Tolerance));

      // Apply constraints
      node->apply_velocity_constraints();

      // Check apply constraints
      velocity << 10.5, -12.5, 0.1;
      for (unsigned i = 0; i < velocity.size(); ++i)
        REQUIRE(node->velocity(Nphase)(i) ==
                Approx(velocity(i)).epsilon(Tolerance));

      // Apply rotation matrix with Euler angles alpha = 10 deg, beta = 20 deg
      // and gamma = 30 deg
      Eigen::Matrix<double, Dim, 1> euler_angles;
      euler_angles << 10. * M_PI / 180, 20. * M_PI / 180, 30. * M_PI / 180;
      const auto rotation_matrix = mpm::geometry::rotation_matrix(euler_angles);
      node->assign_rotation_matrix(rotation_matrix);

      // Apply constraints
      node->apply_velocity_constraints();

      // Check apply constraints
      velocity << 15.774625541122909, -8.021674875339704, 2.117469796447889;
      for (unsigned i = 0; i < Dim; ++i)
        REQUIRE(node->velocity(Nphase)(i) ==
                Approx(velocity(i)).epsilon(Tolerance));
    }

    SECTION("Check acceleration") {
      // Check acceleration
      Eigen::Matrix<double, Dim, 1> acceleration;
      for (unsigned i = 0; i < acceleration.size(); ++i) acceleration(i) = 5.;

      for (unsigned i = 0; i < acceleration.size(); ++i)
        REQUIRE(node->acceleration(Nphase)(i) == Approx(0.).epsilon(Tolerance));

      REQUIRE(node->update_acceleration(true, Nphase, acceleration) == true);
      for (unsigned i = 0; i < acceleration.size(); ++i)
        REQUIRE(node->acceleration(Nphase)(i) == Approx(5.).epsilon(Tolerance));

      // Check if exception is handled
      unsigned bad_phase = 1;
      // Exception handling invalid acceleration dimension
      REQUIRE(node->update_acceleration(true, bad_phase, acceleration) ==
              false);
      // Exception handling invalid acceleration dimension
      REQUIRE(node->update_acceleration(false, bad_phase, acceleration) ==
              false);

      // Apply velocity constraints
      REQUIRE(node->assign_velocity_constraint(0, 10.5) == true);
      REQUIRE(node->assign_velocity_constraint(1, -12.5) == true);
      // Check out of bounds condition
      REQUIRE(node->assign_velocity_constraint(4, 0.) == false);

      // Check acceleration before constraints
      acceleration.resize(Dim);
      acceleration << 5., 5., 5.;
      for (unsigned i = 0; i < acceleration.size(); ++i)
        REQUIRE(node->acceleration(Nphase)(i) ==
                Approx(acceleration(i)).epsilon(Tolerance));

      // Apply constraints
      node->apply_velocity_constraints();

      // Check apply constraints
      acceleration << 0.0, 0.0, 5.;
      for (unsigned i = 0; i < acceleration.size(); ++i)
        REQUIRE(node->acceleration(Nphase)(i) ==
                Approx(acceleration(i)).epsilon(Tolerance));

      // Apply rotation matrix with Euler angles alpha = 10 deg, beta = 20 deg
      // and gamma = 30 deg
      Eigen::Matrix<double, Dim, 1> euler_angles;
      euler_angles << 10. * M_PI / 180, 20. * M_PI / 180, 30. * M_PI / 180;
      const auto rotation_matrix = mpm::geometry::rotation_matrix(euler_angles);
      node->assign_rotation_matrix(rotation_matrix);

      // Apply constraints
      node->apply_velocity_constraints();

      // Check apply constraints
      acceleration << 0.375959332951088, -2.132171329881108, 3.75;
      for (unsigned i = 0; i < Dim; ++i)
        REQUIRE(node->acceleration(Nphase)(i) ==
                Approx(acceleration(i)).epsilon(Tolerance));
    }
  }
}

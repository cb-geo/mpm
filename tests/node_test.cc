#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

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
      Eigen::VectorXd force;
      force.resize(Dim);
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 10.;

      // Check current external force is zero
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(0.).epsilon(Tolerance));

      // Update force to 10.0
      bool status = node->update_external_force(true, Nphase, force);
      REQUIRE(status == true);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));

      // Update force to 20.0
      status = node->update_external_force(true, Nphase, force);
      REQUIRE(status == true);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(20.).epsilon(Tolerance));

      // Assign force as 10.0
      status = node->update_external_force(false, Nphase, force);
      REQUIRE(status == true);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));

      // Check if exception is handled
      force.resize(Dim * 2);
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 10.;

      // Exception handling invalid force dimension
      status = node->update_external_force(true, Nphase, force);
      REQUIRE(status == false);
      // Exception handling invalid force dimension
      status = node->update_external_force(false, Nphase, force);
      REQUIRE(status == false);
    }

    SECTION("Check internal force") {
      // Create a force vector
      Eigen::VectorXd force;
      force.resize(Dim);
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 10.;

      // Check current internal force is zero
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(0.).epsilon(Tolerance));

      // Update force to 10.0
      bool status = node->update_internal_force(true, Nphase, force);
      REQUIRE(status == true);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));

      // Update force to 20.0
      status = node->update_internal_force(true, Nphase, force);
      REQUIRE(status == true);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(20.).epsilon(Tolerance));

      // Assign force as 10.0
      status = node->update_internal_force(false, Nphase, force);
      REQUIRE(status == true);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));

      // Check if exception is handled
      force.resize(Dim * 2);
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 10.;

      // Exception handling invalid force dimension
      status = node->update_internal_force(true, Nphase, force);
      REQUIRE(status == false);
      // Exception handling invalid force dimension
      status = node->update_internal_force(false, Nphase, force);
      REQUIRE(status == false);
    }

    SECTION("Check momentum and velocity") {
      // Check momentum
      Eigen::VectorXd momentum;
      momentum.resize(Dim);
      for (unsigned i = 0; i < momentum.size(); ++i) momentum(i) = 10.;

      // Check initial momentum
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(0.).epsilon(Tolerance));

      // Check update momentum to 10
      bool status = node->update_momentum(true, Nphase, momentum);
      REQUIRE(status == true);
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(10.).epsilon(Tolerance));

      // Check update momentum to 20
      status = node->update_momentum(true, Nphase, momentum);
      REQUIRE(status == true);
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(20.).epsilon(Tolerance));

      // Check assign momentum to 10
      status = node->update_momentum(false, Nphase, momentum);
      REQUIRE(status == true);
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
      momentum.resize(Dim * 2);
      for (unsigned i = 0; i < momentum.size(); ++i) momentum(i) = 10.;

      // Exception handling invalid momentum dimension
      status = node->update_momentum(true, Nphase, momentum);
      REQUIRE(status == false);
      // Exception handling invalid momentum dimension
      status = node->update_momentum(false, Nphase, momentum);
      REQUIRE(status == false);

      // Apply velocity constraints
      std::map<unsigned, double> vel_constraints;
      vel_constraints[0] = 10.5;

      REQUIRE(node->assign_velocity_constraints(vel_constraints) == true);
      // Check out of bounds condition
      vel_constraints[1] = 0.;
      REQUIRE(node->assign_velocity_constraints(vel_constraints) == false);

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
      Eigen::VectorXd acceleration;
      acceleration.resize(Dim);
      for (unsigned i = 0; i < acceleration.size(); ++i) acceleration(i) = 5.;

      for (unsigned i = 0; i < acceleration.size(); ++i)
        REQUIRE(node->acceleration(Nphase)(i) == Approx(0.).epsilon(Tolerance));

      bool status = node->update_acceleration(true, Nphase, acceleration);
      REQUIRE(status == true);
      for (unsigned i = 0; i < acceleration.size(); ++i)
        REQUIRE(node->acceleration(Nphase)(i) == Approx(5.).epsilon(Tolerance));

      // Check if exception is handled
      acceleration.resize(Dim * 2);
      for (unsigned i = 0; i < acceleration.size(); ++i) acceleration(i) = 10.;

      // Exception handling invalid acceleration dimension
      status = node->update_acceleration(true, Nphase, acceleration);
      REQUIRE(status == false);

      // Exception handling invalid acceleration dimension
      status = node->update_acceleration(false, Nphase, acceleration);
      REQUIRE(status == false);

      // Apply velocity constraints
      std::map<unsigned, double> vel_constraints;
      vel_constraints[0] = 10.5;

      REQUIRE(node->assign_velocity_constraints(vel_constraints) == true);
      // Check out of bounds condition
      vel_constraints[1] = 0.;
      REQUIRE(node->assign_velocity_constraints(vel_constraints) == false);

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
      Eigen::VectorXd force;
      force.resize(Dof);
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 10.;

      // Check current external force is zero
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(0.).epsilon(Tolerance));

      // Update force to 10.0
      bool status = node->update_external_force(true, Nphase, force);
      REQUIRE(status == true);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));

      // Update force to 20.0
      status = node->update_external_force(true, Nphase, force);
      REQUIRE(status == true);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(20.).epsilon(Tolerance));

      // Assign force as 10.0
      status = node->update_external_force(false, Nphase, force);
      REQUIRE(status == true);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->external_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));

      // Check if exception is handled
      force.resize(Dim * 2);
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 10.;

      // Exception handling invalid force dimension
      status = node->update_external_force(true, Nphase, force);
      REQUIRE(status == false);
      // Exception handling invalid force dimension
      status = node->update_external_force(false, Nphase, force);
      REQUIRE(status == false);
    }

    SECTION("Check internal force") {
      // Create a force vector
      Eigen::VectorXd force;
      force.resize(Dof);
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 10.;

      // Check current internal force is zero
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(0.).epsilon(Tolerance));

      // Update force to 10.0
      bool status = node->update_internal_force(true, Nphase, force);
      REQUIRE(status == true);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));

      // Update force to 20.0
      status = node->update_internal_force(true, Nphase, force);
      REQUIRE(status == true);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(20.).epsilon(Tolerance));

      // Assign force as 10.0
      status = node->update_internal_force(false, Nphase, force);
      REQUIRE(status == true);
      for (unsigned i = 0; i < force.size(); ++i)
        REQUIRE(node->internal_force(Nphase)(i) ==
                Approx(10.).epsilon(Tolerance));

      // Check if exception is handled
      force.resize(Dim * 2);
      for (unsigned i = 0; i < force.size(); ++i) force(i) = 10.;

      // Exception handling invalid force dimension
      status = node->update_internal_force(true, Nphase, force);
      REQUIRE(status == false);
      // Exception handling invalid force dimension
      status = node->update_internal_force(false, Nphase, force);
      REQUIRE(status == false);
    }

    SECTION("Check momentum and velocity") {
      // Check momentum
      Eigen::VectorXd momentum;
      momentum.resize(Dim);
      for (unsigned i = 0; i < momentum.size(); ++i) momentum(i) = 10.;

      // Check initial momentum
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(0.).epsilon(Tolerance));

      // Check update momentum to 10
      bool status = node->update_momentum(true, Nphase, momentum);
      REQUIRE(status == true);
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(10.).epsilon(Tolerance));

      // Check update momentum to 20
      node->update_momentum(true, Nphase, momentum);
      REQUIRE(status == true);
      for (unsigned i = 0; i < momentum.size(); ++i)
        REQUIRE(node->momentum(Nphase)(i) == Approx(20.).epsilon(Tolerance));

      // Check assign momentum to 10
      node->update_momentum(false, Nphase, momentum);
      REQUIRE(status == true);
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
      momentum.resize(Dim * 2);
      for (unsigned i = 0; i < momentum.size(); ++i) momentum(i) = 10.;

      // Exception handling invalid momentum dimension
      status = node->update_momentum(true, Nphase, momentum);
      REQUIRE(status == false);
      // Exception handling invalid momentum dimension
      status = node->update_momentum(false, Nphase, momentum);
      REQUIRE(status == false);

      // Apply velocity constraints
      std::map<unsigned, double> vel_constraints;
      vel_constraints[0] = -12.5;

      REQUIRE(node->assign_velocity_constraints(vel_constraints) == true);
      // Check out of bounds condition
      vel_constraints[2] = 0.;
      REQUIRE(node->assign_velocity_constraints(vel_constraints) == false);

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
      Eigen::VectorXd acceleration;
      acceleration.resize(Dim);
      for (unsigned i = 0; i < acceleration.size(); ++i) acceleration(i) = 5.;

      for (unsigned i = 0; i < acceleration.size(); ++i)
        REQUIRE(node->acceleration(Nphase)(i) == Approx(0.).epsilon(Tolerance));

      bool status = node->update_acceleration(true, Nphase, acceleration);
      REQUIRE(status == true);
      for (unsigned i = 0; i < acceleration.size(); ++i)
        REQUIRE(node->acceleration(Nphase)(i) == Approx(5.).epsilon(Tolerance));

      // Check if exception is handled
      acceleration.resize(Dim * 2);
      for (unsigned i = 0; i < acceleration.size(); ++i) acceleration(i) = 10.;

      // Exception handling invalid acceleration dimension
      status = node->update_acceleration(true, Nphase, acceleration);
      REQUIRE(status == false);

      // Apply velocity constraints
      std::map<unsigned, double> vel_constraints;
      vel_constraints[0] = -12.5;

      REQUIRE(node->assign_velocity_constraints(vel_constraints) == true);
      // Check out of bounds condition
      vel_constraints[2] = 0.;
      REQUIRE(node->assign_velocity_constraints(vel_constraints) == false);

      // Check acceleration before constraints
      acceleration.resize(Dim);
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
      Eigen::VectorXd force;
      force.resize(Dof);
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
    }

    SECTION("Check internal force") {
      // Create a force vector
      Eigen::VectorXd force;
      force.resize(Dof);
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
      for (unsigned i = 0; i < velocity.size(); ++i) {
        std::cout << "\nvelocity : " << node->velocity(Nphase)(i) << "\t";
        // REQUIRE(node->velocity(Nphase)(i) ==
        //      Approx(velocity(i)).epsilon(Tolerance));
      }

      // Exception check when mass is zero
      mass = 0.;
      // Update mass to 0.
      node->update_mass(false, Nphase, mass);
      REQUIRE(node->mass(Nphase) == Approx(mass).epsilon(Tolerance));
      REQUIRE(node->compute_acceleration_velocity(Nphase, dt) == false);
    }

    SECTION("Check momentum and velocity") {
      // Check momentum
      Eigen::VectorXd momentum;
      momentum.resize(Dim);
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
      std::map<unsigned, double> vel_constraints;
      vel_constraints[0] = 10.5;
      vel_constraints[1] = -12.5;

      REQUIRE(node->assign_velocity_constraints(vel_constraints) == true);
      // Check out of bounds condition
      vel_constraints[4] = 0.;
      REQUIRE(node->assign_velocity_constraints(vel_constraints) == false);

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
    }

    SECTION("Check acceleration") {
      // Check acceleration
      Eigen::VectorXd acceleration;
      acceleration.resize(Dim);
      for (unsigned i = 0; i < acceleration.size(); ++i) acceleration(i) = 5.;

      for (unsigned i = 0; i < acceleration.size(); ++i)
        REQUIRE(node->acceleration(Nphase)(i) == Approx(0.).epsilon(Tolerance));

      bool status = node->update_acceleration(true, Nphase, acceleration);
      REQUIRE(status == true);
      for (unsigned i = 0; i < acceleration.size(); ++i)
        REQUIRE(node->acceleration(Nphase)(i) == Approx(5.).epsilon(Tolerance));

      // Check if exception is handled
      acceleration.resize(Dim * 2);
      for (unsigned i = 0; i < acceleration.size(); ++i) acceleration(i) = 10.;

      // Exception handling invalid acceleration dimension
      status = node->update_acceleration(true, Nphase, acceleration);
      REQUIRE(status == false);

      // Exception handling invalid acceleration dimension
      status = node->update_acceleration(false, Nphase, acceleration);
      REQUIRE(status == false);

      // Apply velocity constraints
      std::map<unsigned, double> vel_constraints;
      vel_constraints[0] = 10.5;
      vel_constraints[1] = -12.5;

      REQUIRE(node->assign_velocity_constraints(vel_constraints) == true);
      // Check out of bounds condition
      vel_constraints[4] = 0.;
      REQUIRE(node->assign_velocity_constraints(vel_constraints) == false);

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
    }
  }
}

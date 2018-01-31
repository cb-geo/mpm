#include <limits>

#include "Eigen/Dense"
#include "catch.hpp"

#include "node.h"

//! \brief Check node class for 1D case
TEST_CASE("Node is checked for 1D case", "[node][1D]") {
  const unsigned Dim = 1;
  const unsigned Dof = 1;
  const unsigned Nphases = 1;
  const unsigned Nphase = 0;
  Eigen::Matrix<double, 1, 1> coords;
  coords.setZero();

  //! Check for id = 0
  SECTION("Node id is zero") {
    mpm::Index id = 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node = std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);
    REQUIRE(node->id() == 0);
  }

  SECTION("Node id is positive") {
    //! Check for id is a positive value
    mpm::Index id = std::numeric_limits<mpm::Index>::max();
    std::shared_ptr<mpm::NodeBase<Dim>> node = std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);
    REQUIRE(node->id() == std::numeric_limits<mpm::Index>::max());
  }

  SECTION("Check degrees of freedom") {
    //! Check for degrees of freedom
    mpm::Index id = 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node = std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);
    REQUIRE(node->dof() == 1);
  }

  //! Test coordinates function
  SECTION("coordinates function is checked") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;

    std::shared_ptr<mpm::NodeBase<Dim>> node = std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);

    //! Check for coordinates being zero
    auto coordinates = node->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));
    REQUIRE(coordinates.size() == Dim);

    //! Check for negative value of coordinates
    for (unsigned i = 0; i < coordinates.size(); ++i)
      coords(i) = -1. * std::numeric_limits<double>::max();
    node->coordinates(coords);
    coordinates = node->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    REQUIRE(coordinates.size() == Dim);

    //! Check for positive value of coordinates
    for (unsigned i = 0; i < coordinates.size(); ++i)
      coords(i) = std::numeric_limits<double>::max();
    node->coordinates(coords);
    coordinates = node->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    REQUIRE(coordinates.size() == Dim);
  }

  SECTION("Check nodal properties") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    std::shared_ptr<mpm::NodeBase<Dim>> node = std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);

    // Check mass
    REQUIRE(node->mass() == std::numeric_limits<double>::max());
    double mass = 100.5;
    node->assign_mass(mass);
    REQUIRE(node->mass() == Approx(100.5).epsilon(Tolerance));

    // Check force
    Eigen::VectorXd force;
    force.resize(Dof);
    for (unsigned i = 0; i < force.size(); ++i) force(i) = 1.;

    for (unsigned i = 0; i < force.size(); ++i)
      REQUIRE(node->force(Nphase)(i) == Approx(0.).epsilon(Tolerance));

    node->assign_force(Nphase, force);
    for (unsigned i = 0; i < force.size(); ++i)
      REQUIRE(node->force(Nphase)(i) == Approx(1.).epsilon(Tolerance));

    // Check velocity
    Eigen::VectorXd velocity;
    velocity.resize(Dof);
    for (unsigned i = 0; i < velocity.size(); ++i) velocity(i) = 1.;

    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(node->velocity(Nphase)(i) == Approx(0.).epsilon(Tolerance));

    node->assign_velocity(Nphase, velocity);
    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(node->velocity(Nphase)(i) == Approx(1.).epsilon(Tolerance));

    // Check momentum
    Eigen::VectorXd momentum;
    momentum.resize(Dof);
    for (unsigned i = 0; i < momentum.size(); ++i) momentum(i) = 1.;

    for (unsigned i = 0; i < momentum.size(); ++i)
      REQUIRE(node->momentum(Nphase)(i) == Approx(0.).epsilon(Tolerance));

    node->assign_momentum(Nphase, momentum);
    for (unsigned i = 0; i < momentum.size(); ++i)
      REQUIRE(node->momentum(Nphase)(i) == Approx(1.).epsilon(Tolerance));

    /*
    // Check acceleration
    Eigen::VectorXd acceleration;
    acceleration.resize(Dof);
    for (unsigned i = 0; i < acceleration.size(); ++i) acceleration(i) = 1.;

    for (unsigned i = 0; i < acceleration.size(); ++i)
      REQUIRE(node->acceleration(Nphase)(i) == Approx(0.).epsilon(Tolerance));

    node->assign_acceleration(Nphase, acceleration);
    for (unsigned i = 0; i < acceleration.size(); ++i)
      REQUIRE(node->acceleration(Nphase)(i) == Approx(1.).epsilon(Tolerance));
    */
  }
}

//! \brief Check node class for 2D case
TEST_CASE("Node is checked for 2D case", "[node][2D]") {
  const unsigned Dim = 2;
  const unsigned Dof = 2;
  const unsigned Nphases = 1;
  const unsigned Nphase = 0;  
  Eigen::Vector2d coords;
  coords.setZero();

  //! Check for id = 0
  SECTION("Node id is zero") {
    mpm::Index id = 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node = std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);
    REQUIRE(node->id() == 0);
  }

  SECTION("Node id is positive") {
    //! Check for id is a positive value
    mpm::Index id = std::numeric_limits<mpm::Index>::max();
    std::shared_ptr<mpm::NodeBase<Dim>> node = std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);
    REQUIRE(node->id() == std::numeric_limits<mpm::Index>::max());
  }

  SECTION("Check degrees of freedom") {
    //! Check for degrees of freedom
    mpm::Index id = 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node = std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);
    REQUIRE(node->dof() == 2);
  }

  //! Test coordinates function
  SECTION("coordinates function is checked") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;

    std::shared_ptr<mpm::NodeBase<Dim>> node = std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);

    //! Check for coordinates being zero
    auto coordinates = node->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));
    REQUIRE(coordinates.size() == Dim);

    //! Check for negative value of coordinates
    for (unsigned i = 0; i < coordinates.size(); ++i)
      coords(i) = -1. * std::numeric_limits<double>::max();
    node->coordinates(coords);
    coordinates = node->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    REQUIRE(coordinates.size() == Dim);

    //! Check for positive value of coordinates
    for (unsigned i = 0; i < coordinates.size(); ++i)
      coords(i) = std::numeric_limits<double>::max();
    node->coordinates(coords);
    coordinates = node->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    REQUIRE(coordinates.size() == Dim);
  }

  SECTION("Check nodal properties") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    std::shared_ptr<mpm::NodeBase<Dim>> node = std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);

    // Check mass
    REQUIRE(node->mass() == std::numeric_limits<double>::max());
    double mass = 100.5;
    node->assign_mass(mass);
    REQUIRE(node->mass() == Approx(100.5).epsilon(Tolerance));

    // Check force
    Eigen::VectorXd force;
    force.resize(Dof);
    for (unsigned i = 0; i < force.size(); ++i) force(i) = 1.;

    for (unsigned i = 0; i < force.size(); ++i)
      REQUIRE(node->force(Nphase)(i) == Approx(0.).epsilon(Tolerance));

    node->assign_force(Nphase, force);
    for (unsigned i = 0; i < force.size(); ++i)
      REQUIRE(node->force(Nphase)(i) == Approx(1.).epsilon(Tolerance));

    // Check velocity
    Eigen::VectorXd velocity;
    velocity.resize(Dof);
    for (unsigned i = 0; i < velocity.size(); ++i) velocity(i) = 1.;

    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(node->velocity(Nphase)(i) == Approx(0.).epsilon(Tolerance));

    node->assign_velocity(Nphase, velocity);
    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(node->velocity(Nphase)(i) == Approx(1.).epsilon(Tolerance));

    // Check momentum
    Eigen::VectorXd momentum;
    momentum.resize(Dof);
    for (unsigned i = 0; i < momentum.size(); ++i) momentum(i) = 1.;

    for (unsigned i = 0; i < momentum.size(); ++i)
      REQUIRE(node->momentum(Nphase)(i) == Approx(0.).epsilon(Tolerance));

    node->assign_momentum(Nphase, momentum);
    for (unsigned i = 0; i < momentum.size(); ++i)
      REQUIRE(node->momentum(Nphase)(i) == Approx(1.).epsilon(Tolerance));

    // Check acceleration
    Eigen::VectorXd acceleration;
    acceleration.resize(Dof);
    for (unsigned i = 0; i < acceleration.size(); ++i) acceleration(i) = 1.;

    for (unsigned i = 0; i < acceleration.size(); ++i)
      REQUIRE(node->acceleration(Nphase)(i) == Approx(0.).epsilon(Tolerance));

    node->assign_acceleration(Nphase, acceleration);
    for (unsigned i = 0; i < acceleration.size(); ++i)
      REQUIRE(node->acceleration(Nphase)(i) == Approx(1.).epsilon(Tolerance));
  }
}

//! \brief Check node class for 3D case
TEST_CASE("Node is checked for 3D case", "[node][3D]") {
  const unsigned Dim = 3;
  const unsigned Dof = 3;
  const unsigned Nphases = 1;
  const unsigned Nphase = 0;

  Eigen::Vector3d coords;
  coords.setZero();

  //! Check for id = 0
  SECTION("Node id is zero") {
    mpm::Index id = 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node = std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);
    REQUIRE(node->id() == 0);
  }

  SECTION("Node id is positive") {
    //! Check for id is a positive value
    mpm::Index id = std::numeric_limits<mpm::Index>::max();
    std::shared_ptr<mpm::NodeBase<Dim>> node = std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);
    REQUIRE(node->id() == std::numeric_limits<mpm::Index>::max());
  }

  SECTION("Check degrees of freedom") {
    //! Check for degrees of freedom
    mpm::Index id = 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node = std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);
    REQUIRE(node->dof() == 3);
  }

  //! Test coordinates function
  SECTION("coordinates function is checked") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;

    std::shared_ptr<mpm::NodeBase<Dim>> node = std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);

    //! Check for coordinates being zero
    auto coordinates = node->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));
    REQUIRE(coordinates.size() == Dim);

    //! Check for negative value of coordinates
    for (unsigned i = 0; i < coordinates.size(); ++i)
      coords(i) = -1. * std::numeric_limits<double>::max();
    node->coordinates(coords);
    coordinates = node->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    REQUIRE(coordinates.size() == Dim);

    //! Check for positive value of coordinates
    for (unsigned i = 0; i < coordinates.size(); ++i)
      coords(i) = std::numeric_limits<double>::max();
    node->coordinates(coords);
    coordinates = node->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    REQUIRE(coordinates.size() == Dim);
  }

  SECTION("Check nodal properties") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    std::shared_ptr<mpm::NodeBase<Dim>> node = std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);

    // Check mass
    REQUIRE(node->mass() == std::numeric_limits<double>::max());
    double mass = 100.5;
    node->assign_mass(mass);
    REQUIRE(node->mass() == Approx(100.5).epsilon(Tolerance));

    // Check force
    Eigen::VectorXd force;
    force.resize(Dof);
    for (unsigned i = 0; i < force.size(); ++i) force(i) = 1.;

    for (unsigned i = 0; i < force.size(); ++i)
      REQUIRE(node->force(Nphase)(i) == Approx(0.).epsilon(Tolerance));

    node->assign_force(Nphase, force);
    for (unsigned i = 0; i < force.size(); ++i)
      REQUIRE(node->force(Nphase)(i) == Approx(1.).epsilon(Tolerance));

    // Check velocity
    Eigen::VectorXd velocity;
    velocity.resize(Dof);
    for (unsigned i = 0; i < velocity.size(); ++i) velocity(i) = 1.;

    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(node->velocity(Nphase)(i) == Approx(0.).epsilon(Tolerance));

    node->assign_velocity(Nphase, velocity);
    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(node->velocity(Nphase)(i) == Approx(1.).epsilon(Tolerance));

    // Check momentum
    Eigen::VectorXd momentum;
    momentum.resize(Dof);
    for (unsigned i = 0; i < momentum.size(); ++i) momentum(i) = 1.;

    for (unsigned i = 0; i < momentum.size(); ++i)
      REQUIRE(node->momentum(Nphase)(i) == Approx(0.).epsilon(Tolerance));

    node->assign_momentum(Nphase, momentum);
    for (unsigned i = 0; i < momentum.size(); ++i)
      REQUIRE(node->momentum(Nphase)(i) == Approx(1.).epsilon(Tolerance));

    // Check acceleration
    Eigen::VectorXd acceleration;
    acceleration.resize(Dof);
    for (unsigned i = 0; i < acceleration.size(); ++i) acceleration(i) = 1.;

    for (unsigned i = 0; i < acceleration.size(); ++i)
      REQUIRE(node->acceleration(Nphase)(i) == Approx(0.).epsilon(Tolerance));

    node->assign_acceleration(Nphase, acceleration);
    for (unsigned i = 0; i < acceleration.size(); ++i)
      REQUIRE(node->acceleration(Nphase)(i) == Approx(1.).epsilon(Tolerance));
  }
}

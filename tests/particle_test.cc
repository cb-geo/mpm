#include <iostream>
#include <limits>

#include "serialize.h"

#include "catch.hpp"

#include "cell.h"
#include "node.h"
#include "particle.h"

//! \brief Check particle class for 1D case
TEST_CASE("Particle is checked for 1D case", "[particle][1D]") {
  // Dimension
  const unsigned Dim = 1;
  // Dimension
  const unsigned Dof = 1;
  // Phases
  const unsigned Nphases = 1;
  // Phase
  const unsigned Phase = 0;

  // Coordinates
  Eigen::Matrix<double, 1, 1> coords;
  coords.setZero();

  //! Check for id = 0
  SECTION("Particle id is zero") {
    mpm::Index id = 0;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);
    REQUIRE(particle->id() == 0);
    REQUIRE(particle->status() == true);
  }

  SECTION("Particle id is positive") {
    //! Check for id is a positive value
    mpm::Index id = std::numeric_limits<mpm::Index>::max();
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);
    REQUIRE(particle->id() == std::numeric_limits<mpm::Index>::max());
    REQUIRE(particle->status() == true);
  }

  //! Construct with id, coordinates and status
  SECTION("Particle with id, coordinates, and status") {
    mpm::Index id = 0;
    bool status = true;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords, status);
    REQUIRE(particle->id() == 0);
    REQUIRE(particle->status() == true);
    particle->assign_status(false);
    REQUIRE(particle->status() == false);
  }

  //! Test coordinates function
  SECTION("coordinates function is checked") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;

    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);

    // Check for coordinates being zero
    auto coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));
    REQUIRE(coordinates.size() == Dim);

    // Check for negative value of coordinates
    for (unsigned i = 0; i < coordinates.size(); ++i)
      coords(i) = -1. * std::numeric_limits<double>::max();
    particle->coordinates(coords);
    coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    REQUIRE(coordinates.size() == Dim);

    // Check for positive value of coordinates
    for (unsigned i = 0; i < coordinates.size(); ++i)
      coords(i) = std::numeric_limits<double>::max();
    particle->coordinates(coords);
    coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    REQUIRE(coordinates.size() == Dim);
  }

  SECTION("Check particle properties") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);

    // Check mass
    REQUIRE(particle->mass(Phase) == Approx(0.0).epsilon(Tolerance));
    double mass = 100.5;
    particle->assign_mass(Phase, mass);
    REQUIRE(particle->mass(Phase) == Approx(100.5).epsilon(Tolerance));

    // Check stress
    Eigen::VectorXd stress;
    stress.resize(Dim);
    for (unsigned i = 0; i < stress.size(); ++i) stress(i) = 1.;

    for (unsigned i = 0; i < stress.size(); ++i)
      REQUIRE(particle->stress(Phase)(i) == Approx(0.).epsilon(Tolerance));

    particle->assign_stress(Phase, stress);
    for (unsigned i = 0; i < stress.size(); ++i)
      REQUIRE(particle->stress(Phase)(i) == Approx(1.).epsilon(Tolerance));

    // Check velocity
    Eigen::VectorXd velocity;
    velocity.resize(Dim);
    for (unsigned i = 0; i < velocity.size(); ++i) velocity(i) = 1.;

    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(particle->velocity(Phase)(i) == Approx(0.).epsilon(Tolerance));

    particle->assign_velocity(Phase, velocity);
    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(particle->velocity(Phase)(i) == Approx(1.).epsilon(Tolerance));

    // Check momentum
    Eigen::VectorXd momentum;
    momentum.resize(Dim);
    for (unsigned i = 0; i < momentum.size(); ++i) momentum(i) = 1.;

    for (unsigned i = 0; i < momentum.size(); ++i)
      REQUIRE(particle->momentum(Phase)(i) == Approx(0.).epsilon(Tolerance));

    particle->assign_momentum(Phase, momentum);
    for (unsigned i = 0; i < momentum.size(); ++i)
      REQUIRE(particle->momentum(Phase)(i) == Approx(1.).epsilon(Tolerance));

    // Check acceleration
    Eigen::VectorXd acceleration;
    acceleration.resize(Dim);
    for (unsigned i = 0; i < acceleration.size(); ++i) acceleration(i) = 1.;

    for (unsigned i = 0; i < acceleration.size(); ++i)
      REQUIRE(particle->acceleration(Phase)(i) ==
              Approx(0.).epsilon(Tolerance));

    particle->assign_acceleration(Phase, acceleration);
    for (unsigned i = 0; i < acceleration.size(); ++i)
      REQUIRE(particle->acceleration(Phase)(i) ==
              Approx(1.).epsilon(Tolerance));
  }

  //! Test serialize function
  SECTION("Serialisation is checked") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;

    // Check for negative value of coordinates
    for (unsigned i = 0; i < coords.size(); ++i) coords(i) = i;

    // Create a string stream
    std::stringstream ss;
    // save data to archive
    {
      auto particle = std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);
      boost::archive::text_oarchive oa(ss);
      oa << *particle;
    }
    // load data from archive
    {
      mpm::Index id = 1;
      // Coordinates
      Eigen::Matrix<double, 1, 1> coordinates;
      coordinates.setZero();

      auto particle =
          std::make_shared<mpm::Particle<Dim, Nphases>>(id, coordinates);
      REQUIRE(particle->id() == 1);

      // Load from archive
      boost::archive::text_iarchive(ss) >> *particle;
      REQUIRE(particle->id() == 0);
      coordinates = particle->coordinates();
      for (unsigned i = 0; i < coordinates.size(); ++i)
        REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));
    }
  }
}

//! \brief Check particle class for 2D case
TEST_CASE("Particle is checked for 2D case", "[particle][2D]") {
  // Dimension
  const unsigned Dim = 2;
  // Degree of freedom
  const unsigned Dof = 2;
  // Number of phases
  const unsigned Nphases = 1;
  // Phase
  const unsigned Phase = 0;  
  // Number of nodes per cell
  const unsigned Nnodes = 4;
  // Tolerance
  const double Tolerance = 1.E-7;
  // Coordinates
  Eigen::Vector2d coords;
  coords.setZero();

  //! Check for id = 0
  SECTION("Particle id is zero") {
    mpm::Index id = 0;
    auto particle = std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);
    REQUIRE(particle->id() == 0);
  }

  SECTION("Particle id is positive") {
    //! Check for id is a positive value
    mpm::Index id = std::numeric_limits<mpm::Index>::max();
    auto particle = std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);
    REQUIRE(particle->id() == std::numeric_limits<mpm::Index>::max());
  }

  //! Test coordinates function
  SECTION("coordinates function is checked") {
    mpm::Index id = 0;
    auto particle = std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);

    //! Check for coordinates being zero
    auto coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));
    REQUIRE(coordinates.size() == Dim);

    //! Check for negative value of coordinates
    for (unsigned i = 0; i < coordinates.size(); ++i)
      coords(i) = -1. * std::numeric_limits<double>::max();
    particle->coordinates(coords);
    coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    REQUIRE(coordinates.size() == Dim);

    //! Check for positive value of coordinates
    for (unsigned i = 0; i < coordinates.size(); ++i)
      coords(i) = std::numeric_limits<double>::max();
    particle->coordinates(coords);
    coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    REQUIRE(coordinates.size() == Dim);
  }

  SECTION("Add a pointer to a cell to particle") {
    // Add particle
    mpm::Index id = 0;
    auto particle = std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);
    // Create cell
    auto cell = std::make_shared<mpm::Cell<Dim>>(0, Nnodes);
    // Add nodes
    std::shared_ptr<mpm::NodeBase<Dim>> node0 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

    coords << 0, 1;
    std::shared_ptr<mpm::NodeBase<Dim>> node1 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);

    coords << 1, 1;
    std::shared_ptr<mpm::NodeBase<Dim>> node2 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);

    coords << 1, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node3 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);
    cell->add_node(0, node0);
    cell->add_node(1, node1);
    cell->add_node(2, node2);
    cell->add_node(3, node3);
    REQUIRE(cell->nnodes() == 4);

    // Add cell to particle
    REQUIRE(cell->status() == false);
    particle->assign_cell(cell);
    REQUIRE(cell->status() == true);
  }

  SECTION("Check particle properties") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);

    // Check mass
    REQUIRE(particle->mass(Phase) == Approx(0.0).epsilon(Tolerance));
    double mass = 100.5;
    particle->assign_mass(Phase, mass);
    REQUIRE(particle->mass(Phase) == Approx(100.5).epsilon(Tolerance));

    // Check stress
    Eigen::VectorXd stress;
    stress.resize(Dim);
    for (unsigned i = 0; i < stress.size(); ++i) stress(i) = 1.;

    for (unsigned i = 0; i < stress.size(); ++i)
      REQUIRE(particle->stress(Phase)(i) == Approx(0.).epsilon(Tolerance));

    particle->assign_stress(Phase, stress);
    for (unsigned i = 0; i < stress.size(); ++i)
      REQUIRE(particle->stress(Phase)(i) == Approx(1.).epsilon(Tolerance));

    // Check velocity
    Eigen::VectorXd velocity;
    velocity.resize(Dim);
    for (unsigned i = 0; i < velocity.size(); ++i) velocity(i) = 1.;

    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(particle->velocity(Phase)(i) == Approx(0.).epsilon(Tolerance));

    particle->assign_velocity(Phase, velocity);
    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(particle->velocity(Phase)(i) == Approx(1.).epsilon(Tolerance));

    // Check momentum
    Eigen::VectorXd momentum;
    momentum.resize(Dim);
    for (unsigned i = 0; i < momentum.size(); ++i) momentum(i) = 1.;

    for (unsigned i = 0; i < momentum.size(); ++i)
      REQUIRE(particle->momentum(Phase)(i) == Approx(0.).epsilon(Tolerance));

    particle->assign_momentum(Phase, momentum);
    for (unsigned i = 0; i < momentum.size(); ++i)
      REQUIRE(particle->momentum(Phase)(i) == Approx(1.).epsilon(Tolerance));

    // Check acceleration
    Eigen::VectorXd acceleration;
    acceleration.resize(Dim);
    for (unsigned i = 0; i < acceleration.size(); ++i) acceleration(i) = 1.;

    for (unsigned i = 0; i < acceleration.size(); ++i)
      REQUIRE(particle->acceleration(Phase)(i) ==
              Approx(0.).epsilon(Tolerance));

    particle->assign_acceleration(Phase, acceleration);
    for (unsigned i = 0; i < acceleration.size(); ++i)
      REQUIRE(particle->acceleration(Phase)(i) ==
              Approx(1.).epsilon(Tolerance));
  }

  //! Test serialize function
  SECTION("Serialisation is checked") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;

    // Check for negative value of coordinates
    for (unsigned i = 0; i < coords.size(); ++i) coords(i) = i;

    // Create a string stream
    std::stringstream ss;
    // save data to archive
    {
      auto particle = std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);
      boost::archive::text_oarchive oa(ss);
      oa << *particle;
    }
    // load data from archive
    {
      mpm::Index id = 1;
      // Coordinates
      Eigen::Vector2d coordinates;
      coordinates.setZero();

      auto particle =
          std::make_shared<mpm::Particle<Dim, Nphases>>(id, coordinates);
      REQUIRE(particle->id() == 1);

      // Load from archive
      boost::archive::text_iarchive(ss) >> *particle;
      REQUIRE(particle->id() == 0);
      coordinates = particle->coordinates();
      for (unsigned i = 0; i < coordinates.size(); ++i)
        REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));
    }
  }
}

//! \brief Check particle class for 3D case
TEST_CASE("Particle is checked for 3D case", "[particle][3D]") {
  // Dimension
  const unsigned Dim = 3;
  // Dimension
  const unsigned Dof = 6;
  // Nnumber of phases
  const unsigned Nphases = 1;
  // Phase
  const unsigned Phase = 0;
  // Number of nodes per cell
  const unsigned Nnodes = 8;
  // Tolerance
  const double Tolerance = 1.E-7;

  // Coordinates
  Eigen::Vector3d coords;
  coords.setZero();

  //! Check for id = 0
  SECTION("Particle id is zero") {
    mpm::Index id = 0;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);
    REQUIRE(particle->id() == 0);
    REQUIRE(particle->status() == true);
  }

  SECTION("Particle id is positive") {
    //! Check for id is a positive value
    mpm::Index id = std::numeric_limits<mpm::Index>::max();
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);
    REQUIRE(particle->id() == std::numeric_limits<mpm::Index>::max());
    REQUIRE(particle->status() == true);
  }

  //! Construct with id, coordinates and status
  SECTION("Particle with id, coordinates, and status") {
    mpm::Index id = 0;
    bool status = true;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords, status);
    REQUIRE(particle->id() == 0);
    REQUIRE(particle->status() == true);
    particle->assign_status(false);
    REQUIRE(particle->status() == false);
  }

  //! Test coordinates function
  SECTION("coordinates function is checked") {
    mpm::Index id = 0;
    // Create particle
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);

    //! Check for coordinates being zero
    auto coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));
    REQUIRE(coordinates.size() == Dim);

    //! Check for negative value of coordinates
    for (unsigned i = 0; i < coordinates.size(); ++i)
      coords(i) = -1. * std::numeric_limits<double>::max();
    particle->coordinates(coords);
    coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    REQUIRE(coordinates.size() == Dim);

    //! Check for positive value of coordinates
    for (unsigned i = 0; i < coordinates.size(); ++i)
      coords(i) = std::numeric_limits<double>::max();
    particle->coordinates(coords);
    coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    REQUIRE(coordinates.size() == Dim);
  }

  //! Test assign cell pointer to particle
  SECTION("Add a pointer to a cell to particle") {
    // Add particle
    mpm::Index id = 0;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);
    // Create cell
    auto cell = std::make_shared<mpm::Cell<Dim>>(0, Nnodes);
    // Add nodes
    coords << 0, 0, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node0 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

    coords << 1, 0, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node1 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);

    coords << 0, 1, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node2 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);

    coords << 1, 1, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node3 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);

    coords << 0, 0, 1;
    std::shared_ptr<mpm::NodeBase<Dim>> node4 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(4, coords);

    coords << 1, 0, 1;
    std::shared_ptr<mpm::NodeBase<Dim>> node5 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(5, coords);

    coords << 0, 1, 1;
    std::shared_ptr<mpm::NodeBase<Dim>> node6 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(6, coords);

    coords << 1, 1, 1;
    std::shared_ptr<mpm::NodeBase<Dim>> node7 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(7, coords);

    cell->add_node(0, node0);
    cell->add_node(1, node1);
    cell->add_node(2, node2);
    cell->add_node(3, node3);
    cell->add_node(4, node4);
    cell->add_node(5, node5);
    cell->add_node(6, node6);
    cell->add_node(7, node7);
    REQUIRE(cell->nnodes() == 8);

    // Add cell to particle
    REQUIRE(cell->status() == false);
    particle->assign_cell(cell);
    REQUIRE(cell->status() == true);
  }

  SECTION("Check particle properties") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);

    // Check mass
    REQUIRE(particle->mass(Phase) == Approx(0.0).epsilon(Tolerance));
    double mass = 100.5;
    particle->assign_mass(Phase, mass);
    REQUIRE(particle->mass(Phase) == Approx(100.5).epsilon(Tolerance));

    // Check stress
    Eigen::VectorXd stress;
    stress.resize(Dim);
    for (unsigned i = 0; i < stress.size(); ++i) stress(i) = 1.;

    for (unsigned i = 0; i < stress.size(); ++i)
      REQUIRE(particle->stress(Phase)(i) == Approx(0.).epsilon(Tolerance));

    particle->assign_stress(Phase, stress);
    for (unsigned i = 0; i < stress.size(); ++i)
      REQUIRE(particle->stress(Phase)(i) == Approx(1.).epsilon(Tolerance));

    // Check velocity
    Eigen::VectorXd velocity;
    velocity.resize(Dim);
    for (unsigned i = 0; i < velocity.size(); ++i) velocity(i) = 1.;

    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(particle->velocity(Phase)(i) == Approx(0.).epsilon(Tolerance));

    particle->assign_velocity(Phase, velocity);
    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(particle->velocity(Phase)(i) == Approx(1.).epsilon(Tolerance));

    // Check momentum
    Eigen::VectorXd momentum;
    momentum.resize(Dim);
    for (unsigned i = 0; i < momentum.size(); ++i) momentum(i) = 1.;

    for (unsigned i = 0; i < momentum.size(); ++i)
      REQUIRE(particle->momentum(Phase)(i) == Approx(0.).epsilon(Tolerance));

    particle->assign_momentum(Phase, momentum);
    for (unsigned i = 0; i < momentum.size(); ++i)
      REQUIRE(particle->momentum(Phase)(i) == Approx(1.).epsilon(Tolerance));

    // Check acceleration
    Eigen::VectorXd acceleration;
    acceleration.resize(Dim);
    for (unsigned i = 0; i < acceleration.size(); ++i) acceleration(i) = 1.;

    for (unsigned i = 0; i < acceleration.size(); ++i)
      REQUIRE(particle->acceleration(Phase)(i) ==
              Approx(0.).epsilon(Tolerance));

    particle->assign_acceleration(Phase, acceleration);
    for (unsigned i = 0; i < acceleration.size(); ++i)
      REQUIRE(particle->acceleration(Phase)(i) ==
              Approx(1.).epsilon(Tolerance));
  }

  //! Test serialize function
  SECTION("Serialisation is checked") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;

    // Check for negative value of coordinates
    for (unsigned i = 0; i < coords.size(); ++i) coords(i) = i;

    // Create a string stream
    std::stringstream ss;
    // save data to archive
    {
      auto particle = std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);
      boost::archive::text_oarchive oa(ss);
      oa << *particle;
    }
    // load data from archive
    {
      mpm::Index id = 1;
      // Coordinates
      Eigen::Vector3d coordinates;
      coordinates.setZero();

      auto particle =
          std::make_shared<mpm::Particle<Dim, Nphases>>(id, coordinates);
      REQUIRE(particle->id() == 1);

      // Load from archive
      boost::archive::text_iarchive(ss) >> *particle;
      REQUIRE(particle->id() == 0);
      coordinates = particle->coordinates();
      for (unsigned i = 0; i < coordinates.size(); ++i)
        REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));
    }
  }
}

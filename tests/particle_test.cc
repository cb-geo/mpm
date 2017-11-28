#include <limits>

#include "Eigen/Dense"
#include "catch.hpp"

#include "cell.h"
#include "particle.h"

//! \brief Check particle class for 1D case
TEST_CASE("Particle is checked for 1D case", "[particle][1D]") {
  // Dimension
  const unsigned Dim = 1;
  // Coordinates
  Eigen::Matrix<double, 1, 1> coords;
  coords.setZero();

  //! Check for id = 0
  SECTION("Particle id is zero") {
    mpm::Index id = 0;
    auto particle = std::make_shared<mpm::Particle<Dim>>(id, coords);
    REQUIRE(particle->id() == 0);
  }

  SECTION("Particle id is positive") {
    //! Check for id is a positive value
    mpm::Index id = std::numeric_limits<mpm::Index>::max();
    auto particle = std::make_shared<mpm::Particle<Dim>>(id, coords);
    REQUIRE(particle->id() == std::numeric_limits<mpm::Index>::max());
  }

  //! Test coordinates function
  SECTION("coordinates function is checked") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;

    auto particle = std::make_shared<mpm::Particle<Dim>>(id, coords);

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
}

//! \brief Check particle class for 2D case
TEST_CASE("Particle is checked for 2D case", "[particle][2D]") {
  // Dimension
  const unsigned Dim = 2;
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
    auto particle = std::make_shared<mpm::Particle<Dim>>(id, coords);
    REQUIRE(particle->id() == 0);
  }

  SECTION("Particle id is positive") {
    //! Check for id is a positive value
    mpm::Index id = std::numeric_limits<mpm::Index>::max();
    auto particle = std::make_shared<mpm::Particle<Dim>>(id, coords);
    REQUIRE(particle->id() == std::numeric_limits<mpm::Index>::max());
  }

  //! Test coordinates function
  SECTION("coordinates function is checked") {
    mpm::Index id = 0;
    auto particle = std::make_shared<mpm::Particle<Dim>>(id, coords);

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
    auto particle = std::make_shared<mpm::Particle<Dim>>(id, coords);
    // Create cell
    auto cell = std::make_shared<mpm::Cell<Dim>>(0, Nnodes);
    // Add nodes
    auto node0 = std::make_shared<mpm::Node<Dim>>(0, coords);

    coords << 0, 1;
    auto node1 = std::make_shared<mpm::Node<Dim>>(1, coords);

    coords << 1, 1;
    auto node2 = std::make_shared<mpm::Node<Dim>>(2, coords);
    
    coords << 1, 0;
    auto node3 = std::make_shared<mpm::Node<Dim>>(3, coords);
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
}

//! \brief Check particle class for 3D case
TEST_CASE("Particle is checked for 3D case", "[particle][3D]") {
  const unsigned Dim = 3;
  Eigen::Vector3d coords;
  coords.setZero();

  //! Check for id = 0
  SECTION("Particle id is zero") {
    mpm::Index id = 0;
    auto particle = std::make_shared<mpm::Particle<Dim>>(id, coords);
    REQUIRE(particle->id() == 0);
  }

  SECTION("Particle id is positive") {
    //! Check for id is a positive value
    mpm::Index id = std::numeric_limits<mpm::Index>::max();
    auto particle = std::make_shared<mpm::Particle<Dim>>(id, coords);
    REQUIRE(particle->id() == std::numeric_limits<mpm::Index>::max());
  }

  //! Test coordinates function
  SECTION("coordinates function is checked") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;

    auto particle = std::make_shared<mpm::Particle<Dim>>(id, coords);

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
}

#include <limits>

#include "node.h"

#include "catch.hpp"
#include "Eigen/Dense"

//! \brief Check node class for 1D case
TEST_CASE("Node is checked for 1D case", "[node][1D]") {
  const unsigned Dim = 1;
  Eigen::Matrix<double, 1, 1> coords;
  coords.setZero();

  //! Check for id = 0
  SECTION("Node id is zero") {
    mpm::Index id = 0;
    auto node = std::make_shared<mpm::Node<Dim>>(id, coords);
    REQUIRE(node->id() == 0);
  }

  SECTION("Node id is negative") {
    //! Check for negative node id
    mpm::Index id = std::numeric_limits<mpm::Index>::min();
    auto node = std::make_shared<mpm::Node<Dim>>(id, coords);
    REQUIRE(node->id() == std::numeric_limits<mpm::Index>::min());
  }

  SECTION("Node id is positive") {
    //! Check for id is a positive value
    mpm::Index id = std::numeric_limits<mpm::Index>::max();
    auto node = std::make_shared<mpm::Node<Dim>>(id, coords);
    REQUIRE(node->id() == std::numeric_limits<mpm::Index>::max());
  }

  //! Test coordinates function
  SECTION("coordinates function is checked") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;

    auto node = std::make_shared<mpm::Node<Dim>>(id, coords);

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
}

//! \brief Check node class for 2D case
TEST_CASE("Node is checked for 2D case", "[node][2D]") {
  const unsigned Dim = 2;
  Eigen::Vector2d coords;
  coords.setZero();

  //! Check for id = 0
  SECTION("Node id is zero") {
    mpm::Index id = 0;
    auto node = std::make_shared<mpm::Node<Dim>>(id, coords);
    REQUIRE(node->id() == 0);
  }

  SECTION("Node id is negative") {
    //! Check for negative node id
    mpm::Index id = std::numeric_limits<mpm::Index>::min();
    auto node = std::make_shared<mpm::Node<Dim>>(id, coords);
    REQUIRE(node->id() == std::numeric_limits<mpm::Index>::min());
  }

  SECTION("Node id is positive") {
    //! Check for id is a positive value
    mpm::Index id = std::numeric_limits<mpm::Index>::max();
    auto node = std::make_shared<mpm::Node<Dim>>(id, coords);
    REQUIRE(node->id() == std::numeric_limits<mpm::Index>::max());
  }

  //! Test coordinates function
  SECTION("coordinates function is checked") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;

    auto node = std::make_shared<mpm::Node<Dim>>(id, coords);

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
}

//! \brief Check node class for 3D case
TEST_CASE("Node is checked for 3D case", "[node][3D]") {
  const unsigned Dim = 3;
  Eigen::Vector3d  coords;
  coords.setZero();

  //! Check for id = 0
  SECTION("Node id is zero") {
    mpm::Index id = 0;
    auto node = std::make_shared<mpm::Node<Dim>>(id, coords);
    REQUIRE(node->id() == 0);
  }

  SECTION("Node id is negative") {
    //! Check for negative node id
    mpm::Index id = std::numeric_limits<mpm::Index>::min();
    auto node = std::make_shared<mpm::Node<Dim>>(id, coords);
    REQUIRE(node->id() == std::numeric_limits<mpm::Index>::min());
  }

  SECTION("Node id is positive") {
    //! Check for id is a positive value
    mpm::Index id = std::numeric_limits<mpm::Index>::max();
    auto node = std::make_shared<mpm::Node<Dim>>(id, coords);
    REQUIRE(node->id() == std::numeric_limits<mpm::Index>::max());
  }

  //! Test coordinates function
  SECTION("coordinates function is checked") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;

    auto node = std::make_shared<mpm::Node<Dim>>(id, coords);

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
}

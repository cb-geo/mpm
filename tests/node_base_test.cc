#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

#include "node_base.h"

//! \brief Check nodebasebase class for 1D case
TEST_CASE("NodeBase is checked for 1D case", "[nodebase][1D]") {
  const unsigned Dim = 1;
  Eigen::Matrix<double, 1, 1> coords;
  coords.setZero();

  //! Check for id = 0
  SECTION("NodeBase id is zero") {
    mpm::Index id = 0;
    auto nodebase = std::make_shared<mpm::NodeBase<Dim>>(id, coords);
    REQUIRE(nodebase->id() == 0);
  }

  SECTION("NodeBase id is positive") {
    //! Check for id is a positive value
    mpm::Index id = std::numeric_limits<mpm::Index>::max();
    auto nodebase = std::make_shared<mpm::NodeBase<Dim>>(id, coords);
    REQUIRE(nodebase->id() == std::numeric_limits<mpm::Index>::max());
  }

  //! Test coordinates function
  SECTION("coordinates function is checked") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;

    auto nodebase = std::make_shared<mpm::NodeBase<Dim>>(id, coords);

    //! Check for coordinates being zero
    auto coordinates = nodebase->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));
    REQUIRE(coordinates.size() == Dim);

    //! Check for negative value of coordinates
    for (unsigned i = 0; i < coordinates.size(); ++i)
      coords(i) = -1. * std::numeric_limits<double>::max();
    nodebase->coordinates(coords);
    coordinates = nodebase->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    REQUIRE(coordinates.size() == Dim);

    //! Check for positive value of coordinates
    for (unsigned i = 0; i < coordinates.size(); ++i)
      coords(i) = std::numeric_limits<double>::max();
    nodebase->coordinates(coords);
    coordinates = nodebase->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    REQUIRE(coordinates.size() == Dim);
  }
}

//! \brief Check nodebase class for 2D case
TEST_CASE("NodeBase is checked for 2D case", "[nodebase][2D]") {
  const unsigned Dim = 2;
  Eigen::Vector2d coords;
  coords.setZero();

  //! Check for id = 0
  SECTION("NodeBase id is zero") {
    mpm::Index id = 0;
    auto nodebase = std::make_shared<mpm::NodeBase<Dim>>(id, coords);
    REQUIRE(nodebase->id() == 0);
  }

  SECTION("NodeBase id is positive") {
    //! Check for id is a positive value
    mpm::Index id = std::numeric_limits<mpm::Index>::max();
    auto nodebase = std::make_shared<mpm::NodeBase<Dim>>(id, coords);
    REQUIRE(nodebase->id() == std::numeric_limits<mpm::Index>::max());
  }

  //! Test coordinates function
  SECTION("coordinates function is checked") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;

    auto nodebase = std::make_shared<mpm::NodeBase<Dim>>(id, coords);

    //! Check for coordinates being zero
    auto coordinates = nodebase->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));
    REQUIRE(coordinates.size() == Dim);

    //! Check for negative value of coordinates
    for (unsigned i = 0; i < coordinates.size(); ++i)
      coords(i) = -1. * std::numeric_limits<double>::max();
    nodebase->coordinates(coords);
    coordinates = nodebase->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    REQUIRE(coordinates.size() == Dim);

    //! Check for positive value of coordinates
    for (unsigned i = 0; i < coordinates.size(); ++i)
      coords(i) = std::numeric_limits<double>::max();
    nodebase->coordinates(coords);
    coordinates = nodebase->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    REQUIRE(coordinates.size() == Dim);
  }
}

//! \brief Check nodebase class for 3D case
TEST_CASE("NodeBase is checked for 3D case", "[nodebase][3D]") {
  const unsigned Dim = 3;
  Eigen::Vector3d coords;
  coords.setZero();

  //! Check for id = 0
  SECTION("NodeBase id is zero") {
    mpm::Index id = 0;
    auto nodebase = std::make_shared<mpm::NodeBase<Dim>>(id, coords);
    REQUIRE(nodebase->id() == 0);
  }

  SECTION("NodeBase id is positive") {
    //! Check for id is a positive value
    mpm::Index id = std::numeric_limits<mpm::Index>::max();
    auto nodebase = std::make_shared<mpm::NodeBase<Dim>>(id, coords);
    REQUIRE(nodebase->id() == std::numeric_limits<mpm::Index>::max());
  }

  //! Test coordinates function
  SECTION("coordinates function is checked") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;

    auto nodebase = std::make_shared<mpm::NodeBase<Dim>>(id, coords);

    //! Check for coordinates being zero
    auto coordinates = nodebase->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));
    REQUIRE(coordinates.size() == Dim);

    //! Check for negative value of coordinates
    for (unsigned i = 0; i < coordinates.size(); ++i)
      coords(i) = -1. * std::numeric_limits<double>::max();
    nodebase->coordinates(coords);
    coordinates = nodebase->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    REQUIRE(coordinates.size() == Dim);

    //! Check for positive value of coordinates
    for (unsigned i = 0; i < coordinates.size(); ++i)
      coords(i) = std::numeric_limits<double>::max();
    nodebase->coordinates(coords);
    coordinates = nodebase->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    REQUIRE(coordinates.size() == Dim);
  }
}

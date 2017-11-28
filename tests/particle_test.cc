#include <limits>

#include "Eigen/Dense"
#include "catch.hpp"

#include "particle.h"

//! \brief Check particle class for 1D case
TEST_CASE("Particle is checked for 1D case", "[particle][1D]") {
  const unsigned Dim = 1;
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
  const unsigned Dim = 2;
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

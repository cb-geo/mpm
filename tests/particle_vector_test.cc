#include <functional>
#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

#include "particle.h"
#include "vector.h"

//! \brief Check particle vector class for 2D case
TEST_CASE("Particle vector is checked for 2D case", "[particlevector][2D]") {
  // Dimension
  const unsigned Dim = 2;
  // Tolerance
  const double Tolerance = 1.E-7;
  // Particle 1
  mpm::Index id1 = 0;
  Eigen::Vector2d coords;
  coords.setZero();
  std::shared_ptr<mpm::Particle<Dim>> particle1 =
      std::make_shared<mpm::Particle<Dim>>(id1, coords);

  // Particle 2
  mpm::Index id2 = 1;
  std::shared_ptr<mpm::Particle<Dim>> particle2 =
      std::make_shared<mpm::Particle<Dim>>(id2, coords);

  // Particle vector
  auto particlevector = std::make_shared<mpm::Vector<mpm::Particle<Dim>>>();

  // Check add particle
  SECTION("Check add particle functionality") {
    // Add particle 1 and check
    REQUIRE(particlevector->add(particle1, false) == true);
    // Add particle 2 and check
    REQUIRE(particlevector->add(particle2, false) == true);
    // Try and particle 2 again and check
    REQUIRE(particlevector->add(particle2) == false);

    // Check size of particle hanlder
    REQUIRE(particlevector->size() == 2);

    // Remove particle 2 and check
    REQUIRE(particlevector->remove(particle2) == true);
    // Try and remove particle 2 again and check
    REQUIRE(particlevector->remove(particle2) == false);
    // Check size of particle hanlder
    REQUIRE(particlevector->size() == 1);

    // Clear particle vector
    particlevector->clear();
    // Check size of particle hanlder
    REQUIRE(particlevector->size() == 0);
  }

  // Check iterator
  SECTION("Check particle range iterator") {
    // Add particle 1
    particlevector->add(particle1);
    // Add particle 2
    particlevector->add(particle2);
    // Check size of particle hanlder
    std::size_t counter = 0;
    for (auto itr = particlevector->cbegin(); itr != particlevector->cend();
         ++itr) {
      auto coords = (*itr)->coordinates();
      // Check if coordinates for each particle is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(0.).epsilon(Tolerance));
      ++counter;
    }

    // Iterate over particles and check if the number of particles is good
    REQUIRE(counter == 2);
  }

  // Check for_each
  SECTION("Check particle for_each") {
    // Add particle 1
    particlevector->add(particle1);
    // Add particle 2
    particlevector->add(particle2);
    // Check size of particle hanlder
    REQUIRE(particlevector->size() == 2);

    // Check coordinates before updating
    for (auto itr = particlevector->cbegin(); itr != particlevector->cend();
         ++itr) {
      auto coords = (*itr)->coordinates();
      // Check if coordinates for each particle is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(0.).epsilon(Tolerance));
    }

    // Set coordinates to unity
    coords << 1., 1.;

    // Iterate through particle vector to update coordinaates
    particlevector->for_each(  // function structure
        std::bind(&mpm::Particle<Dim>::assign_coordinates,
                  std::placeholders::_1, coords));

    // Check if update has gone through
    for (auto itr = particlevector->cbegin(); itr != particlevector->cend();
         ++itr) {
      auto coords = (*itr)->coordinates();
      // Check if coordinates for each particle is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(1.).epsilon(Tolerance));
    }
  }
}

//! \brief Check particle vector class for 3D case
TEST_CASE("Particle vector is checked for 3D case", "[particlevector][3D]") {
  // Dimension
  const unsigned Dim = 3;
  // Phases
  const unsigned Nphases = 1;
  // Tolerance
  const double Tolerance = 1.E-7;

  // Particle 1
  mpm::Index id1 = 0;
  Eigen::Vector3d coords;
  coords.setZero();
  std::shared_ptr<mpm::Particle<Dim>> particle1 =
      std::make_shared<mpm::Particle<Dim>>(id1, coords);

  // Particle 2
  mpm::Index id2 = 1;
  std::shared_ptr<mpm::Particle<Dim>> particle2 =
      std::make_shared<mpm::Particle<Dim>>(id2, coords);

  // Particle vector
  auto particlevector = std::make_shared<mpm::Vector<mpm::Particle<Dim>>>();

  // Check add particle
  SECTION("Check add particle functionality") {
    // Add particle 1 and check
    REQUIRE(particlevector->add(particle1, false) == true);
    // Add particle 2 and check
    REQUIRE(particlevector->add(particle2, false) == true);
    // Try and particle 2 again and check
    REQUIRE(particlevector->add(particle2, true) == false);

    // Check size of particle hanlder
    REQUIRE(particlevector->size() == 2);

    // Remove particle 2 and check
    REQUIRE(particlevector->remove(particle2) == true);
    // Try and remove particle 2 again and check
    REQUIRE(particlevector->remove(particle2) == false);
    // Check size of particle hanlder
    REQUIRE(particlevector->size() == 1);

    // Clear particle vector
    particlevector->clear();
    // Check size of particle hanlder
    REQUIRE(particlevector->size() == 0);
  }

  // Check iterator
  SECTION("Check particle range iterator") {
    // Add particle 1
    particlevector->add(particle1);
    // Add particle 2
    particlevector->add(particle2);
    // Check size of particle hanlder
    std::size_t counter = 0;
    for (auto itr = particlevector->cbegin(); itr != particlevector->cend();
         ++itr) {
      auto coords = (*itr)->coordinates();
      // Check if coordinates for each particle is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(0.).epsilon(Tolerance));

      ++counter;
    }
    // Iterate over particles and check if the number of particles is good
    REQUIRE(counter == 2);
  }

  // Check for_each
  SECTION("Check particle for_each") {
    // Add particle 1
    particlevector->add(particle1);
    // Add particle 2
    particlevector->add(particle2);
    // Check size of particle hanlder
    REQUIRE(particlevector->size() == 2);

    for (auto itr = particlevector->cbegin(); itr != particlevector->cend();
         ++itr) {
      auto coords = (*itr)->coordinates();
      // Check if coordinates for each particle is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(0.).epsilon(Tolerance));
    }

    // Set coordinates to unity
    coords << 1., 1., 1.;

    // Iterate through particle vector to update coordinaates
    particlevector->for_each(std::bind(&mpm::Particle<Dim>::assign_coordinates,
                                       std::placeholders::_1, coords));

    // Check if update has gone through
    for (auto itr = particlevector->cbegin(); itr != particlevector->cend();
         ++itr) {
      auto coords = (*itr)->coordinates();
      // Check if coordinates for each particle is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(1.).epsilon(Tolerance));
    }
  }
}

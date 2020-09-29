#include <functional>
#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

#include "node.h"
#include "vector.h"

//! \brief Check node vector class for 2D case
TEST_CASE("Node vector is checked for 2D case", "[nodevector][2D]") {
  // Dimension
  const unsigned Dim = 2;
  // Degrees of freedom
  const unsigned Dof = 2;
  // Number of phases
  const unsigned Nphases = 1;
  // Tolerance
  const double Tolerance = 1.E-7;

  // Node 1
  mpm::Index id1 = 0;
  Eigen::Vector2d coords;
  coords.setZero();
  std::shared_ptr<mpm::NodeBase<Dim>> node1 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id1, coords);

  // Node 2
  mpm::Index id2 = 1;
  std::shared_ptr<mpm::NodeBase<Dim>> node2 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id2, coords);

  // Node vector
  auto nodevector = std::make_shared<mpm::Vector<mpm::NodeBase<Dim>>>();

  // Check add node
  SECTION("Check add node functionality") {
    // Add node 1 and check
    REQUIRE(nodevector->add(node1, true) == true);
    // Add node 2 and check
    REQUIRE(nodevector->add(node2, false) == true);
    // Try and node 2 again and check
    REQUIRE(nodevector->add(node2) == false);

    // Check size of node hanlder
    REQUIRE(nodevector->size() == 2);

    // Remove node 2 and check
    REQUIRE(nodevector->remove(node2) == true);
    // Try and remove node 2 again and check
    REQUIRE(nodevector->remove(node2) == false);
    // Check size of node hanlder
    REQUIRE(nodevector->size() == 1);

    // Clear node vector
    nodevector->clear();
    // Check size of node hanlder
    REQUIRE(nodevector->size() == 0);
  }

  // Check iterator
  SECTION("Check node range iterator") {
    // Add node 1
    nodevector->add(node1);
    // Add node 2
    nodevector->add(node2);
    // Check size of node hanlder
    std::size_t counter = 0;
    for (auto itr = nodevector->cbegin(); itr != nodevector->cend(); ++itr) {
      auto coords = (*itr)->coordinates();
      // Check if coordinates for each node is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(0.).epsilon(Tolerance));
      ++counter;
    }

    // Iterate over nodes and check if the number of nodes is good
    REQUIRE(counter == 2);
  }

  // Check for_each
  SECTION("Check node for_each") {
    // Add node 1
    nodevector->add(node1);
    // Add node 2
    nodevector->add(node2);
    // Check size of node hanlder
    REQUIRE(nodevector->size() == 2);

    // Check coordinates before updating
    for (auto itr = nodevector->cbegin(); itr != nodevector->cend(); ++itr) {
      auto coords = (*itr)->coordinates();
      // Check if coordinates for each node is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(0.).epsilon(Tolerance));
    }

    // Set coordinates to unity
    coords << 1., 1.;

    // Iterate through node vector to update coordinaates
    nodevector->for_each(  // function structure
        std::bind(&mpm::NodeBase<Dim>::assign_coordinates,
                  std::placeholders::_1, coords));

    // Check if update has gone through
    for (auto itr = nodevector->cbegin(); itr != nodevector->cend(); ++itr) {
      auto coords = (*itr)->coordinates();
      // Check if coordinates for each node is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(1.).epsilon(Tolerance));
    }
  }
}

//! \brief Check node vector class for 3D case
TEST_CASE("Node vector is checked for 3D case", "[nodevector][3D]") {
  // Dimension
  const unsigned Dim = 3;
  // Degrees of freedom
  const unsigned Dof = 6;
  // Number of phases
  const unsigned Nphases = 1;
  // Tolerance
  const double Tolerance = 1.E-7;

  // Node 1
  mpm::Index id1 = 0;
  Eigen::Vector3d coords;
  coords.setZero();
  std::shared_ptr<mpm::NodeBase<Dim>> node1 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id1, coords);

  // Node 2
  mpm::Index id2 = 1;
  std::shared_ptr<mpm::NodeBase<Dim>> node2 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id2, coords);

  // Node vector
  auto nodevector = std::make_shared<mpm::Vector<mpm::NodeBase<Dim>>>();

  // Check add node
  SECTION("Check add node functionality") {
    // Add node 1 and check
    REQUIRE(nodevector->add(node1) == true);
    // Add node 2 and check
    REQUIRE(nodevector->add(node2) == true);
    // Try and node 2 again and check
    REQUIRE(nodevector->add(node2) == false);

    // Check size of node hanlder
    REQUIRE(nodevector->size() == 2);

    // Remove node 2 and check
    REQUIRE(nodevector->remove(node2) == true);
    // Try and remove node 2 again and check
    REQUIRE(nodevector->remove(node2) == false);
    // Check size of node hanlder
    REQUIRE(nodevector->size() == 1);

    // Clear node vector
    nodevector->clear();
    // Check size of node hanlder
    REQUIRE(nodevector->size() == 0);
  }

  // Check iterator
  SECTION("Check node range iterator") {
    // Add node 1
    nodevector->add(node1);
    // Add node 2
    nodevector->add(node2);
    // Check size of node hanlder
    std::size_t counter = 0;
    for (auto itr = nodevector->cbegin(); itr != nodevector->cend(); ++itr) {
      auto coords = (*itr)->coordinates();
      // Check if coordinates for each node is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(0.).epsilon(Tolerance));

      ++counter;
    }
    // Iterate over nodes and check if the number of nodes is good
    REQUIRE(counter == 2);
  }

  // Check for_each
  SECTION("Check node for_each") {
    // Add node 1
    nodevector->add(node1);
    // Add node 2
    nodevector->add(node2);
    // Check size of node hanlder
    REQUIRE(nodevector->size() == 2);

    for (auto itr = nodevector->cbegin(); itr != nodevector->cend(); ++itr) {
      auto coords = (*itr)->coordinates();
      // Check if coordinates for each node is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(0.).epsilon(Tolerance));
    }

    // Set coordinates to unity
    coords << 1., 1., 1.;

    // Iterate through node vector to update coordinaates
    nodevector->for_each(
        // function structure
        std::bind(&mpm::NodeBase<Dim>::assign_coordinates,
                  std::placeholders::_1, coords));

    // Check if update has gone through
    for (auto itr = nodevector->cbegin(); itr != nodevector->cend(); ++itr) {
      auto coords = (*itr)->coordinates();
      // Check if coordinates for each node is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(1.).epsilon(Tolerance));
    }
  }
}

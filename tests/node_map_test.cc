// No-debase test
#include <functional>
#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

#include "map.h"
#include "node.h"

//! \brief Check node map class for 2D case
TEST_CASE("Node map is checked for 2D case", "[nodemap][2D]") {
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

  // Node 3
  mpm::Index id3 = 2;
  std::shared_ptr<mpm::NodeBase<Dim>> node3 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id3, coords);

  // Node map
  auto nodemap = std::make_shared<mpm::Map<mpm::NodeBase<Dim>>>();

  // Check insert and remove node
  SECTION("Check insert and remove node functionality") {
    // Insert node 1 and check
    REQUIRE(nodemap->insert(node1) == true);
    // Insert node 2 and check
    REQUIRE(nodemap->insert(node2->id(), node2) == true);
    // Insert node 3 and check
    REQUIRE(nodemap->insert(node3->id(), node3) == true);
    // Check size of node hanlder
    REQUIRE(nodemap->size() == 3);
    // Remove node 3 and check
    REQUIRE(nodemap->remove(node3->id()) == true);
    // Try to remove node 3 again
    REQUIRE(nodemap->remove(node3->id()) == false);
    // Check size of node hanlder
    REQUIRE(nodemap->size() == 2);
  }

  SECTION("Check operator []") {
    // Insert node 1
    nodemap->insert(id1, node1);
    // Insert node 2
    nodemap->insert(id2, node2);

    // Check operator []
    REQUIRE((*nodemap)[0]->id() == id1);
    REQUIRE((*nodemap)[1]->id() == id2);
  }

  // Check iterator
  SECTION("Check node range iterator") {
    // Insert node 1
    nodemap->insert(node1);
    // Insert node 2
    nodemap->insert(node2);

    // Check find iterator
    REQUIRE(nodemap->find(id1) != nodemap->end());
    REQUIRE(nodemap->find(id2) != nodemap->end());
    REQUIRE(nodemap->find(501) == nodemap->end());

    // Check size of node hanlder
    std::size_t counter = 0;
    for (auto itr = nodemap->begin(); itr != nodemap->end(); ++itr) {
      auto coords = ((*itr).second)->coordinates();
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
    // Insert node 1
    nodemap->insert(node1);
    // Insert node 2
    nodemap->insert(node2);
    // Check size of node hanlder
    REQUIRE(nodemap->size() == 2);

    // Check coordinates before updating
    for (auto itr = nodemap->begin(); itr != nodemap->end(); ++itr) {
      auto coords = ((*itr).second)->coordinates();
      // Check if coordinates for each node is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(0.).epsilon(Tolerance));
    }

    // Set coordinates to unity
    coords << 1., 1.;

    // Iterate through node map to update coordinaates
    nodemap->for_each(  // function structure
        std::bind(&mpm::NodeBase<Dim>::assign_coordinates,
                  std::placeholders::_1, coords));

    // Check if update has gone through
    for (auto itr = nodemap->begin(); itr != nodemap->end(); ++itr) {
      auto coords = ((*itr).second)->coordinates();
      // Check if coordinates for each node is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(1.).epsilon(Tolerance));
    }
  }
}

//! \brief Check node map class for 3D case
TEST_CASE("Node map is checked for 3D case", "[nodemap][3D]") {
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

  // Node 3
  mpm::Index id3 = 2;
  std::shared_ptr<mpm::NodeBase<Dim>> node3 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id3, coords);

  // Node map
  auto nodemap = std::make_shared<mpm::Map<mpm::NodeBase<Dim>>>();

  // Check insert node
  SECTION("Check insert node functionality") {
    // Insert node 1 and check
    REQUIRE(nodemap->insert(node1) == true);
    // Insert node 2 and check
    REQUIRE(nodemap->insert(node2->id(), node2) == true);
    // Insert node 3 and check
    REQUIRE(nodemap->insert(node3->id(), node3) == true);
    // Check size of node hanlder
    REQUIRE(nodemap->size() == 3);
    // Remove node 3 and check
    REQUIRE(nodemap->remove(node3->id()) == true);
    // Try to remove node 3 again
    REQUIRE(nodemap->remove(node3->id()) == false);
    // Check size of node hanlder
    REQUIRE(nodemap->size() == 2);
  }

  SECTION("Check operator []") {
    // Insert node 1
    nodemap->insert(id1, node1);
    // Insert node 2
    nodemap->insert(id2, node2);

    // Check operator []
    REQUIRE((*nodemap)[0]->id() == id1);
    REQUIRE((*nodemap)[1]->id() == id2);
  }

  // Check iterator
  SECTION("Check node range iterator") {
    // Insert node 1
    nodemap->insert(node1);
    // Insert node 2
    nodemap->insert(node2);

    // Check find iterator
    REQUIRE(nodemap->find(id1) != nodemap->end());
    REQUIRE(nodemap->find(id2) != nodemap->end());
    REQUIRE(nodemap->find(501) == nodemap->end());

    // Check size of node hanlder
    std::size_t counter = 0;
    for (auto itr = nodemap->begin(); itr != nodemap->end(); ++itr) {
      auto coords = ((*itr).second)->coordinates();
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
    // Insert node 1
    nodemap->insert(node1);
    // Insert node 2
    nodemap->insert(node2);
    // Check size of node hanlder
    REQUIRE(nodemap->size() == 2);

    for (auto itr = nodemap->begin(); itr != nodemap->end(); ++itr) {
      auto coords = ((*itr).second)->coordinates();
      // Check if coordinates for each node is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(0.).epsilon(Tolerance));
    }

    // Set coordinates to unity
    coords << 1., 1., 1.;

    // Iterate through node map to update coordinaates
    nodemap->for_each(
        // function structure
        std::bind(&mpm::NodeBase<Dim>::assign_coordinates,
                  std::placeholders::_1, coords));

    // Check if update has gone through
    for (auto itr = nodemap->begin(); itr != nodemap->end(); ++itr) {
      auto coords = ((*itr).second)->coordinates();
      // Check if coordinates for each node is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(1.).epsilon(Tolerance));
    }
  }
}

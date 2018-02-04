#include <functional>
#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

#include "container.h"
#include "node.h"

//! \brief Check node container class for 2D case
TEST_CASE("Node container is checked for 2D case", "[nodecontainer][2D]") {
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

  // Node container
  auto nodecontainer = std::make_shared<mpm::Container<mpm::NodeBase<Dim>>>();

  // Check add node
  SECTION("Check add node functionality") {
    // Add node 1 and check status
    bool status1 = nodecontainer->add(node1);
    REQUIRE(status1 == true);
    // Add node 2 and check status
    bool status2 = nodecontainer->add(node2);
    REQUIRE(status2 == true);
    // Try and node 2 again and check status
    bool status3 = nodecontainer->add(node2);
    REQUIRE(status3 == false);

    // Check size of node hanlder
    REQUIRE(nodecontainer->size() == 2);

    // Remove node 2 and check status
    bool remove_status = nodecontainer->remove(node2);
    REQUIRE(remove_status == true);
    // Try and remove node 2 again and check status
    bool remove_status1 = nodecontainer->remove(node2);
    REQUIRE(remove_status1 == false);
    // Check size of node hanlder
    REQUIRE(nodecontainer->size() == 1);

    // Clear node container
    nodecontainer->clear();
    // Check size of node hanlder
    REQUIRE(nodecontainer->size() == 0);
  }

  // Check iterator
  SECTION("Check node range iterator") {
    // Add node 1
    nodecontainer->add(node1);
    // Add node 2
    nodecontainer->add(node2);
    // Check size of node hanlder
    std::size_t counter = 0;
    for (auto itr = nodecontainer->cbegin(); itr != nodecontainer->cend();
         ++itr) {
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
    nodecontainer->add(node1);
    // Add node 2
    nodecontainer->add(node2);
    // Check size of node hanlder
    REQUIRE(nodecontainer->size() == 2);

    // Check coordinates before updating
    for (auto itr = nodecontainer->cbegin(); itr != nodecontainer->cend();
         ++itr) {
      auto coords = (*itr)->coordinates();
      // Check if coordinates for each node is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(0.).epsilon(Tolerance));
    }

    // Set coordinates to unity
    coords << 1., 1.;

    // Iterate through node container to update coordinaates
    nodecontainer->for_each(  // function structure
        std::bind(static_cast<void (mpm::NodeBase<Dim>::*)(
                      const Eigen::Matrix<double, Dim, 1>&)>(
                      // function
                      &mpm::NodeBase<Dim>::coordinates),
                  // arguments
                  std::placeholders::_1, coords));

    // Check if update has gone through
    for (auto itr = nodecontainer->cbegin(); itr != nodecontainer->cend();
         ++itr) {
      auto coords = (*itr)->coordinates();
      // Check if coordinates for each node is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(1.).epsilon(Tolerance));
    }
  }
}

//! \brief Check node container class for 3D case
TEST_CASE("Node container is checked for 3D case", "[nodecontainer][3D]") {
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

  // Node container
  auto nodecontainer = std::make_shared<mpm::Container<mpm::NodeBase<Dim>>>();

  // Check add node
  SECTION("Check add node functionality") {
    // Add node 1 and check status
    bool status1 = nodecontainer->add(node1);
    REQUIRE(status1 == true);
    // Add node 2 and check status
    bool status2 = nodecontainer->add(node2);
    REQUIRE(status2 == true);
    // Try and node 2 again and check status
    bool status3 = nodecontainer->add(node2);
    REQUIRE(status3 == false);

    // Check size of node hanlder
    REQUIRE(nodecontainer->size() == 2);

    // Remove node 2 and check status
    bool remove_status = nodecontainer->remove(node2);
    REQUIRE(remove_status == true);
    // Try and remove node 2 again and check status
    bool remove_status1 = nodecontainer->remove(node2);
    REQUIRE(remove_status1 == false);
    // Check size of node hanlder
    REQUIRE(nodecontainer->size() == 1);

    // Clear node container
    nodecontainer->clear();
    // Check size of node hanlder
    REQUIRE(nodecontainer->size() == 0);
  }

  // Check iterator
  SECTION("Check node range iterator") {
    // Add node 1
    nodecontainer->add(node1);
    // Add node 2
    nodecontainer->add(node2);
    // Check size of node hanlder
    std::size_t counter = 0;
    for (auto itr = nodecontainer->cbegin(); itr != nodecontainer->cend();
         ++itr) {
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
    nodecontainer->add(node1);
    // Add node 2
    nodecontainer->add(node2);
    // Check size of node hanlder
    REQUIRE(nodecontainer->size() == 2);

    for (auto itr = nodecontainer->cbegin(); itr != nodecontainer->cend();
         ++itr) {
      auto coords = (*itr)->coordinates();
      // Check if coordinates for each node is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(0.).epsilon(Tolerance));
    }

    // Set coordinates to unity
    coords << 1., 1., 1.;

    // Iterate through node container to update coordinaates
    nodecontainer->for_each(
        // function structure
        std::bind(static_cast<void (mpm::NodeBase<Dim>::*)(
                      const Eigen::Matrix<double, Dim, 1>&)>(
                      // function
                      &mpm::NodeBase<Dim>::coordinates),
                  // arguments
                  std::placeholders::_1, coords));

    // Check if update has gone through
    for (auto itr = nodecontainer->cbegin(); itr != nodecontainer->cend();
         ++itr) {
      auto coords = (*itr)->coordinates();
      // Check if coordinates for each node is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(1.).epsilon(Tolerance));
    }
  }
}

// No-debase test
#include <functional>
#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

#include "container.h"
#include "node_base.h"

//! \brief Check nodebase container class for 2D case
TEST_CASE("NodeBase container is checked for 2D case",
          "[nodebasecontainer][2D]") {
  // Dimension
  const unsigned Dim = 2;
  // Tolerance
  const double Tolerance = 1.E-7;

  // NodeBase 1
  mpm::Index id1 = 0;
  Eigen::Vector2d coords;
  coords.setZero();
  auto nodebase1 = std::make_shared<mpm::NodeBase<Dim>>(id1, coords);

  // NodeBase 2
  mpm::Index id2 = 1;
  auto nodebase2 = std::make_shared<mpm::NodeBase<Dim>>(id2, coords);

  // NodeBase container
  auto nodebasecontainer =
      std::make_shared<mpm::Container<mpm::NodeBase<Dim>>>();

  // Check add nodebase
  SECTION("Check add nodebase functionality") {
    // Add nodebase 1 and check status
    bool status1 = nodebasecontainer->add(nodebase1);
    REQUIRE(status1 == true);
    // Add nodebase 2 and check status
    bool status2 = nodebasecontainer->add(nodebase2);
    REQUIRE(status2 == true);
    // Try and nodebase 2 again and check status
    bool status3 = nodebasecontainer->add(nodebase2);
    REQUIRE(status3 == false);

    // Check size of nodebase hanlder
    REQUIRE(nodebasecontainer->size() == 2);

    // Remove nodebase 2 and check status
    bool remove_status = nodebasecontainer->remove(nodebase2);
    REQUIRE(remove_status == true);
    // Try and remove nodebase 2 again and check status
    bool remove_status1 = nodebasecontainer->remove(nodebase2);
    REQUIRE(remove_status1 == false);
    // Check size of nodebase hanlder
    REQUIRE(nodebasecontainer->size() == 1);

    // Clear nodebase container
    nodebasecontainer->clear();
    // Check size of nodebase hanlder
    REQUIRE(nodebasecontainer->size() == 0);
  }

  // Check iterator
  SECTION("Check nodebase range iterator") {
    // Add nodebase 1
    nodebasecontainer->add(nodebase1);
    // Add nodebase 2
    nodebasecontainer->add(nodebase2);
    // Check size of nodebase hanlder
    std::size_t counter = 0;
    for (auto itr = nodebasecontainer->begin(); itr != nodebasecontainer->end();
         ++itr) {
      auto coords = (*itr)->coordinates();
      // Check if coordinates for each nodebase is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(0.).epsilon(Tolerance));
      ++counter;
    }

    // Iterate over nodebases and check if the number of nodebases is good
    REQUIRE(counter == 2);
  }

  // Check for_each
  SECTION("Check nodebase for_each") {
    // Add nodebase 1
    nodebasecontainer->add(nodebase1);
    // Add nodebase 2
    nodebasecontainer->add(nodebase2);
    // Check size of nodebase hanlder
    REQUIRE(nodebasecontainer->size() == 2);

    // Check coordinates before updating
    for (auto itr = nodebasecontainer->begin(); itr != nodebasecontainer->end();
         ++itr) {
      auto coords = (*itr)->coordinates();
      // Check if coordinates for each nodebase is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(0.).epsilon(Tolerance));
    }

    // Set coordinates to unity
    coords << 1., 1.;

    // Iterate through nodebase container to update coordinaates
    nodebasecontainer->for_each(  // function structure
        std::bind(static_cast<void (mpm::NodeBase<Dim>::*)(
                      const Eigen::Matrix<double, Dim, 1>&)>(
                      // function
                      &mpm::NodeBase<Dim>::coordinates),
                  // arguments
                  std::placeholders::_1, coords));

    // Check if update has gone through
    for (auto itr = nodebasecontainer->begin(); itr != nodebasecontainer->end();
         ++itr) {
      auto coords = (*itr)->coordinates();
      // Check if coordinates for each nodebase is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(1.).epsilon(Tolerance));
    }
  }
}

//! \brief Check nodebase container class for 3D case
TEST_CASE("NodeBase container is checked for 3D case",
          "[nodebasecontainer][3D]") {
  // Dimension
  const unsigned Dim = 3;
  // Tolerance
  const double Tolerance = 1.E-7;

  // NodeBase 1
  mpm::Index id1 = 0;
  Eigen::Vector3d coords;
  coords.setZero();
  auto nodebase1 = std::make_shared<mpm::NodeBase<Dim>>(id1, coords);

  // NodeBase 2
  mpm::Index id2 = 1;
  auto nodebase2 = std::make_shared<mpm::NodeBase<Dim>>(id2, coords);

  // NodeBase container
  auto nodebasecontainer =
      std::make_shared<mpm::Container<mpm::NodeBase<Dim>>>();

  // Check add nodebase
  SECTION("Check add nodebase functionality") {
    // Add nodebase 1 and check status
    bool status1 = nodebasecontainer->add(nodebase1);
    REQUIRE(status1 == true);
    // Add nodebase 2 and check status
    bool status2 = nodebasecontainer->add(nodebase2);
    REQUIRE(status2 == true);
    // Try and nodebase 2 again and check status
    bool status3 = nodebasecontainer->add(nodebase2);
    REQUIRE(status3 == false);

    // Check size of nodebase hanlder
    REQUIRE(nodebasecontainer->size() == 2);

    // Remove nodebase 2 and check status
    bool remove_status = nodebasecontainer->remove(nodebase2);
    REQUIRE(remove_status == true);
    // Try and remove nodebase 2 again and check status
    bool remove_status1 = nodebasecontainer->remove(nodebase2);
    REQUIRE(remove_status1 == false);
    // Check size of nodebase hanlder
    REQUIRE(nodebasecontainer->size() == 1);

    // Clear nodebase container
    nodebasecontainer->clear();
    // Check size of nodebase hanlder
    REQUIRE(nodebasecontainer->size() == 0);
  }

  // Check iterator
  SECTION("Check nodebase range iterator") {
    // Add nodebase 1
    nodebasecontainer->add(nodebase1);
    // Add nodebase 2
    nodebasecontainer->add(nodebase2);
    // Check size of nodebase hanlder
    std::size_t counter = 0;
    for (auto itr = nodebasecontainer->begin(); itr != nodebasecontainer->end();
         ++itr) {
      auto coords = (*itr)->coordinates();
      // Check if coordinates for each nodebase is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(0.).epsilon(Tolerance));

      ++counter;
    }
    // Iterate over nodebases and check if the number of nodebases is good
    REQUIRE(counter == 2);
  }

  // Check for_each
  SECTION("Check nodebase for_each") {
    // Add nodebase 1
    nodebasecontainer->add(nodebase1);
    // Add nodebase 2
    nodebasecontainer->add(nodebase2);
    // Check size of nodebase hanlder
    REQUIRE(nodebasecontainer->size() == 2);

    for (auto itr = nodebasecontainer->begin(); itr != nodebasecontainer->end();
         ++itr) {
      auto coords = (*itr)->coordinates();
      // Check if coordinates for each nodebase is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(0.).epsilon(Tolerance));
    }

    // Set coordinates to unity
    coords << 1., 1., 1.;

    // Iterate through nodebase container to update coordinaates
    nodebasecontainer->for_each(
        // function structure
        std::bind(static_cast<void (mpm::NodeBase<Dim>::*)(
                      const Eigen::Matrix<double, Dim, 1>&)>(
                      // function
                      &mpm::NodeBase<Dim>::coordinates),
                  // arguments
                  std::placeholders::_1, coords));

    // Check if update has gone through
    for (auto itr = nodebasecontainer->begin(); itr != nodebasecontainer->end();
         ++itr) {
      auto coords = (*itr)->coordinates();
      // Check if coordinates for each nodebase is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(1.).epsilon(Tolerance));
    }
  }
}

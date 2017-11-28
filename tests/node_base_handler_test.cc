// No-debase test
#include <functional>
#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

#include "handler.h"
#include "node_base.h"

//! \brief Check nodebase handler class for 2D case
TEST_CASE("NodeBase handler is checked for 2D case", "[nodebasehandler][2D]") {
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

  // NodeBase handler
  auto nodebasehandler = std::make_shared<mpm::Handler<mpm::NodeBase<Dim>>>();

  // Check insert nodebase
  SECTION("Check insert nodebase functionality") {
    // Insert nodebase 1 and check status
    bool status1 = nodebasehandler->insert(nodebase1);
    REQUIRE(status1 == true);
    // Insert nodebase 2 and check status
    bool status2 = nodebasehandler->insert(nodebase2->id(), nodebase2);
    REQUIRE(status2 == true);
    // Check size of nodebase hanlder
    REQUIRE(nodebasehandler->size() == 2);
  }

  // Check iterator
  SECTION("Check nodebase range iterator") {
    // Insert nodebase 1
    nodebasehandler->insert(nodebase1);
    // Insert nodebase 2
    nodebasehandler->insert(nodebase2);
    // Check size of nodebase hanlder
    std::size_t counter = 0;
    for (auto itr = nodebasehandler->begin(); itr != nodebasehandler->end();
         ++itr) {
      auto coords = ((*itr).second)->coordinates();
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
    // Insert nodebase 1
    nodebasehandler->insert(nodebase1);
    // Insert nodebase 2
    nodebasehandler->insert(nodebase2);
    // Check size of nodebase hanlder
    REQUIRE(nodebasehandler->size() == 2);

    // Check coordinates before updating
    for (auto itr = nodebasehandler->begin(); itr != nodebasehandler->end();
         ++itr) {
      auto coords = ((*itr).second)->coordinates();
      // Check if coordinates for each nodebase is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(0.).epsilon(Tolerance));
    }

    // Set coordinates to unity
    coords << 1., 1.;

    // Iterate through nodebase handler to update coordinaates
    nodebasehandler->for_each(  // function structure
        std::bind(static_cast<void (mpm::NodeBase<Dim>::*)(
                      const Eigen::Matrix<double, Dim, 1>&)>(
                      // function
                      &mpm::NodeBase<Dim>::coordinates),
                  // arguments
                  std::placeholders::_1, coords));

    // Check if update has gone through
    for (auto itr = nodebasehandler->begin(); itr != nodebasehandler->end();
         ++itr) {
      auto coords = ((*itr).second)->coordinates();
      // Check if coordinates for each nodebase is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(1.).epsilon(Tolerance));
    }
  }

}

//! \brief Check nodebase handler class for 3D case
TEST_CASE("NodeBase handler is checked for 3D case", "[nodebasehandler][3D]") {
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

  // NodeBase handler
  auto nodebasehandler = std::make_shared<mpm::Handler<mpm::NodeBase<Dim>>>();

  // Check insert nodebase
  SECTION("Check insert nodebase functionality") {
    // Insert nodebase 1 and check status
    bool status1 = nodebasehandler->insert(nodebase1);
    REQUIRE(status1 == true);
    // Insert nodebase 2 and check status
    bool status2 = nodebasehandler->insert(nodebase2->id(), nodebase2);
    REQUIRE(status2 == true);
    // Check size of nodebase hanlder
    REQUIRE(nodebasehandler->size() == 2);
  }

  // Check iterator
  SECTION("Check nodebase range iterator") {
    // Insert nodebase 1
    nodebasehandler->insert(nodebase1);
    // Insert nodebase 2
    nodebasehandler->insert(nodebase2);
    // Check size of nodebase hanlder
    std::size_t counter = 0;
    for (auto itr = nodebasehandler->begin(); itr != nodebasehandler->end();
         ++itr) {
      auto coords = ((*itr).second)->coordinates();
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
    // Insert nodebase 1
    nodebasehandler->insert(nodebase1);
    // Insert nodebase 2
    nodebasehandler->insert(nodebase2);
    // Check size of nodebase hanlder
    REQUIRE(nodebasehandler->size() == 2);
    
    for (auto itr = nodebasehandler->begin(); itr != nodebasehandler->end();
         ++itr) {
      auto coords = ((*itr).second)->coordinates();
      // Check if coordinates for each nodebase is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(0.).epsilon(Tolerance));
    }

    // Set coordinates to unity
    coords << 1., 1., 1.;

    // Iterate through nodebase handler to update coordinaates
    nodebasehandler->for_each(
        // function structure
        std::bind(static_cast<void (mpm::NodeBase<Dim>::*)(
                      const Eigen::Matrix<double, Dim, 1>&)>(
                      // function
                      &mpm::NodeBase<Dim>::coordinates),
                  // arguments
                  std::placeholders::_1, coords));

    // Check if update has gone through
    for (auto itr = nodebasehandler->begin(); itr != nodebasehandler->end();
         ++itr) {
      auto coords = ((*itr).second)->coordinates();
      // Check if coordinates for each nodebase is zero
      for (unsigned i = 0; i < coords.size(); ++i)
        REQUIRE(coords[i] == Approx(1.).epsilon(Tolerance));
    }
  }
}

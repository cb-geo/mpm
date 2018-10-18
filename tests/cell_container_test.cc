#include <functional>
#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

#include "cell.h"
#include "container.h"
#include "element.h"
#include "hexahedron_element.h"
#include "node.h"
#include "quadrilateral_element.h"

//! \brief Check cell container class for 2D case
TEST_CASE("Cell container is checked for 2D case", "[cellcontainer][2D]") {
  // Dimension
  const unsigned Dim = 2;
  // Degrees of freedom
  const unsigned Dof = 2;
  // Number of nodes per cell
  const unsigned Nnodes = 4;
  // Number of phases
  const unsigned Nphases = 1;

  // Tolerance
  const double Tolerance = 1.E-7;

  // Element
  std::shared_ptr<mpm::Element<Dim>> element =
      std::make_shared<mpm::QuadrilateralElement<Dim, 4>>();

  // Cell 1
  mpm::Index id1 = 0;
  auto cell1 = std::make_shared<mpm::Cell<Dim>>(id1, Nnodes, element);

  // Cell 2
  mpm::Index id2 = 1;
  auto cell2 = std::make_shared<mpm::Cell<Dim>>(id2, Nnodes, element);

  // Cell container
  auto cellcontainer = std::make_shared<mpm::Container<mpm::Cell<Dim>>>();

  // Check add cell
  SECTION("Check add cell functionality") {
    // Add cell 1 and check
    REQUIRE(cellcontainer->add(cell1) == true);
    // Add cell 2 and check
    REQUIRE(cellcontainer->add(cell2) == true);
    // Try and cell 2 again and check
    REQUIRE(cellcontainer->add(cell2) == false);

    // Check size of cell hanlder
    REQUIRE(cellcontainer->size() == 2);

    // Remove cell 2 and check
    REQUIRE(cellcontainer->remove(cell2) == true);
    // Try and remove cell 2 again and check
    REQUIRE(cellcontainer->remove(cell2) == false);
    // Check size of cell hanlder
    REQUIRE(cellcontainer->size() == 1);

    // Clear cell container
    cellcontainer->clear();
    // Check size of cell hanlder
    REQUIRE(cellcontainer->size() == 0);
  }

  // Check iterator
  SECTION("Check cell range iterator") {
    // Add cell 1
    cellcontainer->add(cell1);
    // Add cell 2
    cellcontainer->add(cell2);
    // Check size of cell hanlder
    std::size_t counter = 0;
    for (auto itr = cellcontainer->cbegin(); itr != cellcontainer->cend();
         ++itr) {
      REQUIRE((*itr)->nnodes() == 0);
      ++counter;
    }

    // Iterate over cells and check if the number of cells is good
    REQUIRE(counter == 2);
  }

  // Check for_each
  SECTION("Check cell for_each") {
    // Add cell 1
    cellcontainer->add(cell1);
    // Add cell 2
    cellcontainer->add(cell2);
    // Check size of cell hanlder
    REQUIRE(cellcontainer->size() == 2);

    for (auto itr = cellcontainer->cbegin(); itr != cellcontainer->cend();
         ++itr) {
      REQUIRE((*itr)->nfunctions() == 4);
    }

    Eigen::Vector2d coords;
    coords.setZero();

    std::shared_ptr<mpm::NodeBase<Dim>> node0 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

    coords << 2., 0.;
    std::shared_ptr<mpm::NodeBase<Dim>> node1 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);

    coords << 2., 2.;
    std::shared_ptr<mpm::NodeBase<Dim>> node2 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);

    coords << 0., 2.;
    std::shared_ptr<mpm::NodeBase<Dim>> node3 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);

    coords << 4., 0.;
    std::shared_ptr<mpm::NodeBase<Dim>> node4 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(4, coords);

    coords << 4., 2.;
    std::shared_ptr<mpm::NodeBase<Dim>> node5 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(5, coords);

    REQUIRE(cell1->add_node(0, node0) == true);
    REQUIRE(cell1->add_node(1, node1) == true);
    REQUIRE(cell1->add_node(2, node2) == true);
    REQUIRE(cell1->add_node(3, node3) == true);
    REQUIRE(cell1->nnodes() == 4);

    REQUIRE(cell2->add_node(0, node1) == true);
    REQUIRE(cell2->add_node(1, node4) == true);
    REQUIRE(cell2->add_node(2, node5) == true);
    REQUIRE(cell2->add_node(3, node2) == true);
    REQUIRE(cell2->nnodes() == 4);

    // Check if cell has been initialised
    for (auto itr = cellcontainer->cbegin(); itr != cellcontainer->cend();
         ++itr)
      REQUIRE((*itr)->is_initialised() == false);

    // Iterate through cell container to initialise
    cellcontainer->for_each(
        std::bind(&mpm::Cell<Dim>::initialise, std::placeholders::_1));

    // Check if update has gone through
    for (auto itr = cellcontainer->cbegin(); itr != cellcontainer->cend();
         ++itr)
      REQUIRE((*itr)->is_initialised() == true);
  }
}

//! \brief Check cell container class for 3D case
TEST_CASE("Cell container is checked for 3D case", "[cellcontainer][3D]") {
  // Dimension
  const unsigned Dim = 3;
  // Degrees of freedom
  const unsigned Dof = 6;
  // Number of phases
  const unsigned Nphases = 1;
  // Number of nodes per cell
  const unsigned Nnodes = 8;

  // Tolerance
  const double Tolerance = 1.E-7;

  // Element
  std::shared_ptr<mpm::Element<Dim>> element =
      std::make_shared<mpm::HexahedronElement<Dim, 8>>();

  // Cell 1
  mpm::Index id1 = 0;
  auto cell1 = std::make_shared<mpm::Cell<Dim>>(id1, Nnodes, element);

  // Cell 2
  mpm::Index id2 = 1;
  auto cell2 = std::make_shared<mpm::Cell<Dim>>(id2, Nnodes, element);

  // Cell container
  auto cellcontainer = std::make_shared<mpm::Container<mpm::Cell<Dim>>>();

  // Check add cell
  SECTION("Check add cell functionality") {
    // Add cell 1 and check
    REQUIRE(cellcontainer->add(cell1) == true);
    // Add cell 2 and check
    REQUIRE(cellcontainer->add(cell2) == true);
    // Try and cell 2 again and check
    REQUIRE(cellcontainer->add(cell2) == false);

    // Check size of cell hanlder
    REQUIRE(cellcontainer->size() == 2);

    // Remove cell 2 and check
    REQUIRE(cellcontainer->remove(cell2) == true);
    // Try and remove cell 2 again and check
    REQUIRE(cellcontainer->remove(cell2) == false);
    // Check size of cell hanlder
    REQUIRE(cellcontainer->size() == 1);

    // Clear cell container
    cellcontainer->clear();
    // Check size of cell hanlder
    REQUIRE(cellcontainer->size() == 0);
  }

  // Check iterator
  SECTION("Check cell range iterator") {
    // Add cell 1
    cellcontainer->add(cell1);
    // Add cell 2
    cellcontainer->add(cell2);
    // Check size of cell hanlder
    std::size_t counter = 0;
    for (auto itr = cellcontainer->cbegin(); itr != cellcontainer->cend();
         ++itr) {
      REQUIRE((*itr)->nnodes() == 0);
      ++counter;
    }

    // Iterate over cells and check if the number of cells is good
    REQUIRE(counter == 2);
  }

  // Check for_each
  SECTION("Check cell for_each") {
    // Add cell 1
    cellcontainer->add(cell1);
    // Add cell 2
    cellcontainer->add(cell2);
    // Check size of cell hanlder
    REQUIRE(cellcontainer->size() == 2);

    for (auto itr = cellcontainer->cbegin(); itr != cellcontainer->cend();
         ++itr) {
      REQUIRE((*itr)->nfunctions() == 8);
    }

    // Coordinates
    Eigen::Vector3d coords;

    coords << 0, 0, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node0 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

    coords << 2, 0, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node1 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);

    coords << 2, 2, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node2 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);

    coords << 0, 2, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node3 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);

    coords << 0, 0, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node4 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(4, coords);

    coords << 2, 0, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node5 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(5, coords);

    coords << 2, 2, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node6 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(6, coords);

    coords << 0, 2, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node7 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(7, coords);

    coords << 4, 0, 0.;
    std::shared_ptr<mpm::NodeBase<Dim>> node8 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(8, coords);

    coords << 4, 2., 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node9 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(9, coords);

    coords << 4., 0., 2.;
    std::shared_ptr<mpm::NodeBase<Dim>> node10 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(10, coords);

    coords << 4., 2., 0.;
    std::shared_ptr<mpm::NodeBase<Dim>> node11 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(11, coords);

    // Cell 1
    REQUIRE(cell1->add_node(0, node0) == true);
    REQUIRE(cell1->add_node(1, node1) == true);
    REQUIRE(cell1->add_node(2, node2) == true);
    REQUIRE(cell1->add_node(3, node3) == true);
    REQUIRE(cell1->add_node(4, node4) == true);
    REQUIRE(cell1->add_node(5, node5) == true);
    REQUIRE(cell1->add_node(6, node6) == true);
    REQUIRE(cell1->add_node(7, node7) == true);
    REQUIRE(cell1->nnodes() == 8);

    // Cell 2
    REQUIRE(cell2->add_node(0, node1) == true);
    REQUIRE(cell2->add_node(1, node8) == true);
    REQUIRE(cell2->add_node(2, node9) == true);
    REQUIRE(cell2->add_node(3, node2) == true);
    REQUIRE(cell2->add_node(4, node5) == true);
    REQUIRE(cell2->add_node(5, node10) == true);
    REQUIRE(cell2->add_node(6, node11) == true);
    REQUIRE(cell2->add_node(7, node6) == true);
    REQUIRE(cell2->nnodes() == 8);

    // Check if cell has been initialised
    for (auto itr = cellcontainer->cbegin(); itr != cellcontainer->cend();
         ++itr)
      REQUIRE((*itr)->is_initialised() == false);

    // Iterate through cell container to initialise
    cellcontainer->for_each(
        std::bind(&mpm::Cell<Dim>::initialise, std::placeholders::_1));

    // Check if update has gone through
    for (auto itr = cellcontainer->cbegin(); itr != cellcontainer->cend();
         ++itr)
      REQUIRE((*itr)->is_initialised() == true);
  }
}

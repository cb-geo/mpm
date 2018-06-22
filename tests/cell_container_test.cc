#include <functional>
#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

#include "cell.h"
#include "container.h"
#include "hex_shapefn.h"
#include "quad_shapefn.h"
#include "shapefn.h"

//! \brief Check cell container class for 2D case
TEST_CASE("Cell container is checked for 2D case", "[cellcontainer][2D]") {
  // Dimension
  const unsigned Dim = 2;
  // Number of nodes per cell
  const unsigned Nnodes = 4;

  // Tolerance
  const double Tolerance = 1.E-7;

  // Shape function
  std::shared_ptr<mpm::ShapeFn<Dim>> shapefn =
      std::make_shared<mpm::QuadrilateralShapeFn<Dim, 4>>();

  // Cell 1
  mpm::Index id1 = 0;
  auto cell1 = std::make_shared<mpm::Cell<Dim>>(id1, Nnodes, shapefn);

  // Cell 2
  mpm::Index id2 = 1;
  auto cell2 = std::make_shared<mpm::Cell<Dim>>(id2, Nnodes, shapefn);

  // Cell container
  auto cellcontainer = std::make_shared<mpm::Container<mpm::Cell<Dim>>>();

  // Check add cell
  SECTION("Check add cell functionality") {
    // Add cell 1 and check status
    bool status1 = cellcontainer->add(cell1);
    REQUIRE(status1 == true);
    // Add cell 2 and check status
    bool status2 = cellcontainer->add(cell2);
    REQUIRE(status2 == true);
    // Try and cell 2 again and check status
    bool status3 = cellcontainer->add(cell2);
    REQUIRE(status3 == false);

    // Check size of cell hanlder
    REQUIRE(cellcontainer->size() == 2);

    // Remove cell 2 and check status
    bool remove_status = cellcontainer->remove(cell2);
    REQUIRE(remove_status == true);
    // Try and remove cell 2 again and check status
    bool remove_status1 = cellcontainer->remove(cell2);
    REQUIRE(remove_status1 == false);
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

    // Iterate through cell container to update coordinaates
    // TODO: Remove shapefn call
    cellcontainer->for_each(
        std::bind(&mpm::Cell<Dim>::shapefn, std::placeholders::_1, shapefn));

    // Check if update has gone through
    for (auto itr = cellcontainer->cbegin(); itr != cellcontainer->cend();
         ++itr) {
      REQUIRE((*itr)->nfunctions() == 4);
    }
  }
}

//! \brief Check cell container class for 3D case
TEST_CASE("Cell container is checked for 3D case", "[cellcontainer][3D]") {
  // Dimension
  const unsigned Dim = 3;  // Dimension
  // Number of nodes per cell
  const unsigned Nnodes = 8;

  // Tolerance
  const double Tolerance = 1.E-7;

  // Shape function
  std::shared_ptr<mpm::ShapeFn<Dim>> shapefn =
      std::make_shared<mpm::HexahedronShapeFn<Dim, 8>>();

  // Cell 1
  mpm::Index id1 = 0;
  auto cell1 = std::make_shared<mpm::Cell<Dim>>(id1, Nnodes, shapefn);

  // Cell 2
  mpm::Index id2 = 1;
  auto cell2 = std::make_shared<mpm::Cell<Dim>>(id2, Nnodes, shapefn);

  // Cell container
  auto cellcontainer = std::make_shared<mpm::Container<mpm::Cell<Dim>>>();

  // Check add cell
  SECTION("Check add cell functionality") {
    // Add cell 1 and check status
    bool status1 = cellcontainer->add(cell1);
    REQUIRE(status1 == true);
    // Add cell 2 and check status
    bool status2 = cellcontainer->add(cell2);
    REQUIRE(status2 == true);
    // Try and cell 2 again and check status
    bool status3 = cellcontainer->add(cell2);
    REQUIRE(status3 == false);

    // Check size of cell hanlder
    REQUIRE(cellcontainer->size() == 2);

    // Remove cell 2 and check status
    bool remove_status = cellcontainer->remove(cell2);
    REQUIRE(remove_status == true);
    // Try and remove cell 2 again and check status
    bool remove_status1 = cellcontainer->remove(cell2);
    REQUIRE(remove_status1 == false);
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

    // Iterate through cell container to update coordinaates
    // TODO: Remove functions
    cellcontainer->for_each(
        std::bind(&mpm::Cell<Dim>::shapefn, std::placeholders::_1, shapefn));

    // Check if update has gone through
    for (auto itr = cellcontainer->cbegin(); itr != cellcontainer->cend();
         ++itr) {
      REQUIRE((*itr)->nfunctions() == 8);
    }
  }
}

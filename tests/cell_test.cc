#include <limits>

#include "cell.h"

#include "catch.hpp"

//! \brief Check cell class for 2D case
TEST_CASE("Cell is checked for 2D case", "[cell][2D]") {
  // Dimension
  const unsigned Dim = 2;
  // Number of nodes per cell
  const unsigned Nnodes = 4;

  //! Check for id = 0
  SECTION("Cell id is zero") {
    mpm::Index id = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes);
    REQUIRE(cell->id() == 0);
  }

  SECTION("Cell id is positive") {
    //! Check for id is a positive value
    mpm::Index id = std::numeric_limits<mpm::Index>::max();
    auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes);
    REQUIRE(cell->id() == std::numeric_limits<mpm::Index>::max());
  }

}

//! \brief Check cell class for 3D case
TEST_CASE("Cell is checked for 3D case", "[cell][3D]") {
  // Dimension
  const unsigned Dim = 3;
  // Number of nodes per cell
  const unsigned Nnodes = 4;

  //! Check for id = 0
  SECTION("Cell id is zero") {
    mpm::Index id = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes);
    REQUIRE(cell->id() == 0);
  }

  SECTION("Cell id is positive") {
    //! Check for id is a positive value
    mpm::Index id = std::numeric_limits<mpm::Index>::max();
    auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes);
    REQUIRE(cell->id() == std::numeric_limits<mpm::Index>::max());
  }

}

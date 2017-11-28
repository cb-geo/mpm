#include <limits>
#include <memory>

#include "catch.hpp"
#include "Eigen/Dense"

#include "cell.h"
#include "hex_shapefn.h"
#include "node_base.h"
#include "quad_shapefn.h"
#include "shapefn.h"

//! \brief Check cell class for 2D case
TEST_CASE("Cell is checked for 2D case", "[cell][2D]") {
  // Dimension
  const unsigned Dim = 2;
  // Number of nodebases per cell
  const unsigned Nnodes = 4;

  Eigen::Vector2d coords;
  coords.setZero();
  
  auto nodebase0 = std::make_shared<mpm::NodeBase<Dim>>(0, coords);

  coords << 0, 1;
  auto nodebase1 = std::make_shared<mpm::NodeBase<Dim>>(1, coords);
  
  coords << 1, 1;
  auto nodebase2 = std::make_shared<mpm::NodeBase<Dim>>(2, coords);

  coords << 1, 0;
  auto nodebase3 = std::make_shared<mpm::NodeBase<Dim>>(3, coords);  

  //! Check Cell IDs
  SECTION("Check cell ids") {
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

  SECTION("Add nodebases") {
    mpm::Index id = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes);
    cell->add_node(0, nodebase0);
    cell->add_node(1, nodebase1);
    cell->add_node(2, nodebase2);
    cell->add_node(3, nodebase3);
    REQUIRE(cell->nnodes() == 4);
  }

  SECTION("Add neighbours") {
    auto cell = std::make_shared<mpm::Cell<Dim>>(0, Nnodes);
    auto neighbourcell = std::make_shared<mpm::Cell<Dim>>(1, Nnodes);
    REQUIRE(cell->nneighbours() == 0);
    cell->add_neighbour(0, neighbourcell);
    REQUIRE(cell->nneighbours() == 1);
  }
  
  SECTION("Check shape functions") {
    mpm::Index id = 0;
    auto shapefn = std::make_shared<mpm::QuadrilateralShapeFn<Dim>>(Nnodes);
    auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes,shapefn);
    REQUIRE(cell->nfunctions() == 4);
    // Check 8-nodebased function
    shapefn = std::make_shared<mpm::QuadrilateralShapeFn<Dim>>(8);
    cell->shapefn(shapefn);
    REQUIRE(cell->nfunctions() == 8);
    // Check 9-nodebased function
    shapefn = std::make_shared<mpm::QuadrilateralShapeFn<Dim>>(9);
    cell->shapefn(shapefn);
    REQUIRE(cell->nfunctions() == 9);
  }
}

//! \brief Check cell class for 3D case
TEST_CASE("Cell is checked for 3D case", "[cell][3D]") {
  // Dimension
  const unsigned Dim = 3;
  // Number of nodebases per cell
  const unsigned Nnodes = 8;

  // Coordinaates
  Eigen::Vector3d coords;

  coords << 0, 0, 0;
  auto nodebase0 = std::make_shared<mpm::NodeBase<Dim>>(0, coords);

  coords << 1, 0, 0;
  auto nodebase1 = std::make_shared<mpm::NodeBase<Dim>>(1, coords);
  
  coords << 0, 1, 0;
  auto nodebase2 = std::make_shared<mpm::NodeBase<Dim>>(2, coords);

  coords << 1, 1, 0;
  auto nodebase3 = std::make_shared<mpm::NodeBase<Dim>>(3, coords);  

  coords << 0, 0, 1;
  auto nodebase4 = std::make_shared<mpm::NodeBase<Dim>>(4, coords);

  coords << 1, 0, 1;
  auto nodebase5 = std::make_shared<mpm::NodeBase<Dim>>(5, coords);
  
  coords << 0, 1, 1;
  auto nodebase6 = std::make_shared<mpm::NodeBase<Dim>>(6, coords);

  coords << 1, 1, 1;
  auto nodebase7 = std::make_shared<mpm::NodeBase<Dim>>(7, coords);  

  //! Check Cell IDs
  SECTION("Check cell ids") {
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

  // Check nodebase additions
  SECTION("Add nodebases") {
    mpm::Index id = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes);
    cell->add_node(0, nodebase0);
    cell->add_node(1, nodebase1);
    cell->add_node(2, nodebase2);
    cell->add_node(3, nodebase3);
    cell->add_node(4, nodebase4);
    cell->add_node(5, nodebase5);
    cell->add_node(6, nodebase6);
    cell->add_node(7, nodebase7);
    REQUIRE(cell->nnodes() == 8);
  }

  SECTION("Add neighbours") {
    auto cell = std::make_shared<mpm::Cell<Dim>>(0, Nnodes);
    auto neighbourcell = std::make_shared<mpm::Cell<Dim>>(1, Nnodes);
    REQUIRE(cell->nneighbours() == 0);
    cell->add_neighbour(0, neighbourcell);
    REQUIRE(cell->nneighbours() == 1);
  }

  SECTION("Check shape functions") {
    mpm::Index id = 0;
    auto shapefn = std::make_shared<mpm::HexahedronShapeFn<Dim>>(Nnodes);
    auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes,shapefn);
    REQUIRE(cell->nfunctions() == 8);
    // Check 20-nodebased function
    shapefn = std::make_shared<mpm::HexahedronShapeFn<Dim>>(20);
    cell->shapefn(shapefn);
    REQUIRE(cell->nfunctions() == 20);
  }
}

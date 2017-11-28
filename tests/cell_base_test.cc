#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

#include "cell_base.h"
#include "hex_shapefn.h"
#include "node_base.h"
#include "quad_shapefn.h"
#include "shapefn.h"

//! \brief Check cellbase class for 2D case
TEST_CASE("CellBase is checked for 2D case", "[cellbase][2D]") {
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

  //! Check CellBase IDs
  SECTION("Check cellbase ids") {
    //! Check for id = 0
    SECTION("CellBase id is zero") {
      mpm::Index id = 0;
      auto cellbase = std::make_shared<mpm::CellBase<Dim>>(id, Nnodes);
      REQUIRE(cellbase->id() == 0);
    }

    SECTION("CellBase id is positive") {
      //! Check for id is a positive value
      mpm::Index id = std::numeric_limits<mpm::Index>::max();
      auto cellbase = std::make_shared<mpm::CellBase<Dim>>(id, Nnodes);
      REQUIRE(cellbase->id() == std::numeric_limits<mpm::Index>::max());
    }
  }

  SECTION("Add nodebases") {
    mpm::Index id = 0;
    auto cellbase = std::make_shared<mpm::CellBase<Dim>>(id, Nnodes);
    cellbase->add_node(0, nodebase0);
    cellbase->add_node(1, nodebase1);
    cellbase->add_node(2, nodebase2);
    cellbase->add_node(3, nodebase3);
    REQUIRE(cellbase->nnodes() == 4);
  }

  SECTION("Add neighbours") {
    auto cellbase = std::make_shared<mpm::CellBase<Dim>>(0, Nnodes);
    auto neighbourcellbase = std::make_shared<mpm::CellBase<Dim>>(1, Nnodes);
    REQUIRE(cellbase->nneighbours() == 0);
    cellbase->add_neighbour(0, neighbourcellbase);
    REQUIRE(cellbase->nneighbours() == 1);
  }

  SECTION("Check shape functions") {
    mpm::Index id = 0;
    auto shapefn = std::make_shared<mpm::QuadrilateralShapeFn<Dim>>(Nnodes);
    auto cellbase = std::make_shared<mpm::CellBase<Dim>>(id, Nnodes, shapefn);
    REQUIRE(cellbase->nfunctions() == 4);
    // Check 8-nodebased function
    shapefn = std::make_shared<mpm::QuadrilateralShapeFn<Dim>>(8);
    cellbase->shapefn(shapefn);
    REQUIRE(cellbase->nfunctions() == 8);
    // Check 9-nodebased function
    shapefn = std::make_shared<mpm::QuadrilateralShapeFn<Dim>>(9);
    cellbase->shapefn(shapefn);
    REQUIRE(cellbase->nfunctions() == 9);
  }
}

//! \brief Check cellbase class for 3D case
TEST_CASE("CellBase is checked for 3D case", "[cellbase][3D]") {
  // Dimension
  const unsigned Dim = 3;
  // Number of nodebases per cellbase
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

  //! Check CellBase IDs
  SECTION("Check cellbase ids") {
    //! Check for id = 0
    SECTION("CellBase id is zero") {
      mpm::Index id = 0;
      auto cellbase = std::make_shared<mpm::CellBase<Dim>>(id, Nnodes);
      REQUIRE(cellbase->id() == 0);
    }

    SECTION("CellBase id is positive") {
      //! Check for id is a positive value
      mpm::Index id = std::numeric_limits<mpm::Index>::max();
      auto cellbase = std::make_shared<mpm::CellBase<Dim>>(id, Nnodes);
      REQUIRE(cellbase->id() == std::numeric_limits<mpm::Index>::max());
    }
  }

  // Check nodebase additions
  SECTION("Add nodebases") {
    mpm::Index id = 0;
    auto cellbase = std::make_shared<mpm::CellBase<Dim>>(id, Nnodes);
    cellbase->add_node(0, nodebase0);
    cellbase->add_node(1, nodebase1);
    cellbase->add_node(2, nodebase2);
    cellbase->add_node(3, nodebase3);
    cellbase->add_node(4, nodebase4);
    cellbase->add_node(5, nodebase5);
    cellbase->add_node(6, nodebase6);
    cellbase->add_node(7, nodebase7);
    REQUIRE(cellbase->nnodes() == 8);
  }

  SECTION("Add neighbours") {
    auto cellbase = std::make_shared<mpm::CellBase<Dim>>(0, Nnodes);
    auto neighbourcellbase = std::make_shared<mpm::CellBase<Dim>>(1, Nnodes);
    REQUIRE(cellbase->nneighbours() == 0);
    cellbase->add_neighbour(0, neighbourcellbase);
    REQUIRE(cellbase->nneighbours() == 1);
  }

  SECTION("Check shape functions") {
    mpm::Index id = 0;
    auto shapefn = std::make_shared<mpm::HexahedronShapeFn<Dim>>(Nnodes);
    auto cellbase = std::make_shared<mpm::CellBase<Dim>>(id, Nnodes, shapefn);
    REQUIRE(cellbase->nfunctions() == 8);
    // Check 20-nodebased function
    shapefn = std::make_shared<mpm::HexahedronShapeFn<Dim>>(20);
    cellbase->shapefn(shapefn);
    REQUIRE(cellbase->nfunctions() == 20);
  }
}

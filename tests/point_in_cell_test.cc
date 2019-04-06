#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

#include "cell.h"
#include "element.h"
#include "factory.h"
#include "hexahedron_element.h"
#include "hexahedron_quadrature.h"
#include "node.h"
#include "quadrilateral_element.h"
#include "quadrilateral_quadrature.h"

TEST_CASE("Point in cell 2D", "[PointInCell][2D]") {
  // Dimension
  const unsigned Dim = 2;
  // Degrees of freedom
  const unsigned Dof = 2;
  // Number of phases
  const unsigned Nphases = 1;
  // Number of nodes per cell
  const unsigned Nnodes = 4;
  // Tolerance
  const double Tolerance = 1.E-7;

  SECTION("Transform real to unit cell analytical solution") {
    // Number of nodes in cell
    const unsigned Nnodes = 4;

    // Coordinates
    Eigen::Vector2d coords;

    coords << 0.656514162228664, 0.448587131356584;
    std::shared_ptr<mpm::NodeBase<Dim>> node0 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

    coords << 0.609997617675458, 0.448995487014756;
    std::shared_ptr<mpm::NodeBase<Dim>> node1 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);

    coords << 0.612187210083002, 0.414580484205138;
    std::shared_ptr<mpm::NodeBase<Dim>> node2 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);

    coords << 0.651629357356265, 0.391627886274249;
    std::shared_ptr<mpm::NodeBase<Dim>> node3 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);

    // 4-noded quadrilateral shape functions
    std::shared_ptr<mpm::Element<Dim>> element =
        Factory<mpm::Element<Dim>>::instance()->create("ED2Q4");

    mpm::Index id = 0;
    auto cell = std::make_unique<mpm::Cell<Dim>>(id, Nnodes, element);

    REQUIRE(cell->add_node(0, node0) == true);
    REQUIRE(cell->add_node(1, node1) == true);
    REQUIRE(cell->add_node(2, node2) == true);
    REQUIRE(cell->add_node(3, node3) == true);
    REQUIRE(cell->nnodes() == 4);

    // Initialise cell
    REQUIRE(cell->initialise() == true);

    // Coordinates of a point in real cell
    Eigen::Vector2d point;
    point << 0.632582, 0.425948;
    // Test if point is in cell
    REQUIRE(cell->is_point_in_cell(point) == true);

    point << 0.632585, 0.42595;
    // Test if point is in cell
    REQUIRE(cell->is_point_in_cell(point) == true);
  }

  SECTION("Transform real to unit cell Newton-Raphson") {
    // Number of nodes in cell
    const unsigned Nnodes = 4;

    // Coordinates
    Eigen::Vector2d coords;

    coords << 0.049340385470457, 0.546167667109886;
    std::shared_ptr<mpm::NodeBase<Dim>> node0 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

    coords << 0.049008570276153, 0.497592363325129;
    std::shared_ptr<mpm::NodeBase<Dim>> node1 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);

    coords << 0.097545161257934, 0.490392640151913;
    std::shared_ptr<mpm::NodeBase<Dim>> node2 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);

    coords << 0.098928337700656, 0.541016130614386;
    std::shared_ptr<mpm::NodeBase<Dim>> node3 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);

    // 4-noded quadrilateral shape functions
    std::shared_ptr<mpm::Element<Dim>> element =
        Factory<mpm::Element<Dim>>::instance()->create("ED2Q4");

    mpm::Index id = 0;
    auto cell = std::make_unique<mpm::Cell<Dim>>(id, Nnodes, element);

    REQUIRE(cell->add_node(0, node0) == true);
    REQUIRE(cell->add_node(1, node1) == true);
    REQUIRE(cell->add_node(2, node2) == true);
    REQUIRE(cell->add_node(3, node3) == true);
    REQUIRE(cell->nnodes() == 4);

    // Initialise cell
    REQUIRE(cell->initialise() == true);

    // Coordinates of a point in real cell
    Eigen::Vector2d point;
    point << 0.0597025, 0.534722;
    // Test if point is in cell
    REQUIRE(cell->is_point_in_cell(point) == true);
  }

  SECTION("Transform real to unit cell analytical solution") {
    // Number of nodes in cell
    const unsigned Nnodes = 4;

    // Coordinates
    Eigen::Vector2d coords;

    coords << 0.899835351184034, 2.41201458730691;
    std::shared_ptr<mpm::NodeBase<Dim>> node0 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

    coords << 0.900000000003745, 2.5;
    std::shared_ptr<mpm::NodeBase<Dim>> node1 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);

    coords << 0.800000000003329, 2.5;
    std::shared_ptr<mpm::NodeBase<Dim>> node2 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);

    coords << 0.804253384820319, 2.41407134328757;
    std::shared_ptr<mpm::NodeBase<Dim>> node3 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);

    // 4-noded quadrilateral shape functions
    std::shared_ptr<mpm::Element<Dim>> element =
        Factory<mpm::Element<Dim>>::instance()->create("ED2Q4");

    mpm::Index id = 0;
    auto cell = std::make_unique<mpm::Cell<Dim>>(id, Nnodes, element);

    REQUIRE(cell->add_node(0, node0) == true);
    REQUIRE(cell->add_node(1, node1) == true);
    REQUIRE(cell->add_node(2, node2) == true);
    REQUIRE(cell->add_node(3, node3) == true);
    REQUIRE(cell->nnodes() == 4);

    // Initialise cell
    REQUIRE(cell->initialise() == true);

    // Coordinates of a point in real cell
    Eigen::Vector2d point;
    point << 0.879474, 2.43095;
    // Test if point is in cell
    REQUIRE(cell->is_point_in_cell(point) == true);

    point << 0.87903, 2.4815;
    // Test if point is in cell
    REQUIRE(cell->is_point_in_cell(point) == true);

    point << 0.821834, 2.48175;
    // Test if point is in cell
    REQUIRE(cell->is_point_in_cell(point) == true);

    point << 0.823751, 2.43189;
    // Test if point is in cell
    REQUIRE(cell->is_point_in_cell(point) == true);
  }
  // Check if a point is in a skewed isoparametric element
  SECTION(
      "Point in a skewed isoparametric element with Affine transformation and "
      "Newton Raphson") {

    // Shape function
    // 4-noded quadrilateral shape functions
    std::shared_ptr<mpm::Element<Dim>> element =
        Factory<mpm::Element<Dim>>::instance()->create("ED2Q4");

    Eigen::Vector2d coords;
    coords.setZero();

    // Check point in cell
    Eigen::Vector2d point;
    point << 0.37763268079517553, 0.10469406123558715;

    // Local coordinates
    Eigen::Vector2d xi;

    // Coordinates
    coords << 0.375, 0.1013;
    std::shared_ptr<mpm::NodeBase<Dim>> node0 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

    coords << 0.38, 0.098247997;
    std::shared_ptr<mpm::NodeBase<Dim>> node1 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);

    coords << 0.38, 0.10325;
    std::shared_ptr<mpm::NodeBase<Dim>> node2 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);

    coords << 0.375, 0.1063;
    std::shared_ptr<mpm::NodeBase<Dim>> node3 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);

    //! Check Cell IDs
    mpm::Index id = 0;
    auto cell1 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

    // Check if cell is initialised, before addition of nodes
    REQUIRE(cell1->is_initialised() == false);

    cell1->add_node(0, node0);
    cell1->add_node(1, node1);
    cell1->add_node(2, node2);
    cell1->add_node(3, node3);
    REQUIRE(cell1->nnodes() == 4);

    REQUIRE(cell1->nfunctions() == 4);

    // Check using unit cell with affine transformation / Newton-Raphson
    REQUIRE(cell1->is_point_in_cell(point) == false);

    // Cell 2
    mpm::Index id2 = 1;
    auto cell2 = std::make_shared<mpm::Cell<Dim>>(id2, Nnodes, element);

    // Check if cell is initialised, before addition of nodes
    REQUIRE(cell2->is_initialised() == false);

    // Element 1
    coords << 0.375, 0.1063;
    std::shared_ptr<mpm::NodeBase<Dim>> node10 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(10, coords);

    coords << 0.38, 0.10325;
    std::shared_ptr<mpm::NodeBase<Dim>> node11 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(11, coords);

    coords << 0.38, 0.10825;
    std::shared_ptr<mpm::NodeBase<Dim>> node12 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(12, coords);

    coords << 0.375, 0.1113;
    std::shared_ptr<mpm::NodeBase<Dim>> node13 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(13, coords);

    cell2->add_node(0, node10);
    cell2->add_node(1, node11);
    cell2->add_node(2, node12);
    cell2->add_node(3, node13);
    REQUIRE(cell2->nnodes() == 4);
    REQUIRE(cell2->nfunctions() == 4);

    // Point in cell captures with tolerance
    REQUIRE(cell2->is_point_in_cell(point) == true);
  }
}

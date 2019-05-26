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
    auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

    REQUIRE(cell->add_node(0, node0) == true);
    REQUIRE(cell->add_node(1, node1) == true);
    REQUIRE(cell->add_node(2, node2) == true);
    REQUIRE(cell->add_node(3, node3) == true);
    REQUIRE(cell->nnodes() == 4);

    // Initialise cell
    REQUIRE(cell->initialise() == true);

    Eigen::Vector2d xi;
    // Coordinates of a point in real cell
    Eigen::Vector2d point;
    point << 0.632582, 0.425948;
    // Test if point is in cell
    REQUIRE(cell->is_point_in_cell(point, &xi) == true);

    point << 0.632585, 0.42595;
    // Test if point is in cell
    REQUIRE(cell->is_point_in_cell(point, &xi) == true);
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
    auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

    REQUIRE(cell->add_node(0, node0) == true);
    REQUIRE(cell->add_node(1, node1) == true);
    REQUIRE(cell->add_node(2, node2) == true);
    REQUIRE(cell->add_node(3, node3) == true);
    REQUIRE(cell->nnodes() == 4);

    // Initialise cell
    REQUIRE(cell->initialise() == true);

    Eigen::Vector2d xi;
    // Coordinates of a point in real cell
    Eigen::Vector2d point;
    point << 0.0597025, 0.534722;
    // Test if point is in cell
    REQUIRE(cell->is_point_in_cell(point, &xi) == true);
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
    auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

    REQUIRE(cell->add_node(0, node0) == true);
    REQUIRE(cell->add_node(1, node1) == true);
    REQUIRE(cell->add_node(2, node2) == true);
    REQUIRE(cell->add_node(3, node3) == true);
    REQUIRE(cell->nnodes() == 4);

    // Initialise cell
    REQUIRE(cell->initialise() == true);

    Eigen::Vector2d xi;
    // Coordinates of a point in real cell
    Eigen::Vector2d point;
    point << 0.879474, 2.43095;
    // Test if point is in cell
    REQUIRE(cell->is_point_in_cell(point, &xi) == true);

    point << 0.87903, 2.4815;
    // Test if point is in cell
    REQUIRE(cell->is_point_in_cell(point, &xi) == true);

    point << 0.821834, 2.48175;
    // Test if point is in cell
    REQUIRE(cell->is_point_in_cell(point, &xi) == true);

    point << 0.823751, 2.43189;
    // Test if point is in cell
    REQUIRE(cell->is_point_in_cell(point, &xi) == true);
  }

  SECTION("Transform real to unit cell analytical solution") {
    // Number of nodes in cell
    const unsigned Nnodes = 4;

    // Coordinates
    Eigen::Vector2d coords;

    coords << -2.0, -2.0;
    std::shared_ptr<mpm::NodeBase<Dim>> node0 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

    coords << 2.0, -2.0;
    std::shared_ptr<mpm::NodeBase<Dim>> node1 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);

    coords << 2.0, 2.0;
    std::shared_ptr<mpm::NodeBase<Dim>> node2 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);

    coords << -2.0, 2.0;
    std::shared_ptr<mpm::NodeBase<Dim>> node3 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);

    // 4-noded quadrilateral shape functions
    std::shared_ptr<mpm::Element<Dim>> element =
        Factory<mpm::Element<Dim>>::instance()->create("ED2Q4");

    mpm::Index id = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

    REQUIRE(cell->add_node(0, node0) == true);
    REQUIRE(cell->add_node(1, node1) == true);
    REQUIRE(cell->add_node(2, node2) == true);
    REQUIRE(cell->add_node(3, node3) == true);
    REQUIRE(cell->nnodes() == 4);

    // Initialise cell
    REQUIRE(cell->initialise() == true);

    Eigen::Vector2d xi;

    // Coordinates of a point in real cell
    Eigen::Vector2d point;
    point << 0.5, 0.5;

    // Test if point is in cell
    REQUIRE(cell->is_point_in_cell(point, &xi) == true);

    // Coordinates of the point in an unit cell
    Eigen::Matrix<double, 2, 1> point_unit_cell;
    point_unit_cell << 0.25, 0.25;

    // Use Newton-raphson iteration to find local coordinates
    auto local_point = cell->transform_real_to_unit_cell(point);
    for (unsigned i = 0; i < local_point.size(); ++i)
      REQUIRE(local_point[i] == Approx(point_unit_cell[i]).epsilon(Tolerance));

    // Use analytical solution
    local_point = cell->local_coordinates_point_2d(point);
    for (unsigned i = 0; i < local_point.size(); ++i)
      REQUIRE(local_point[i] == Approx(point_unit_cell[i]).epsilon(Tolerance));

    // Coordinates of a point in real cell
    point << 0., 0.;
    // Test if point is in cell
    REQUIRE(cell->is_point_in_cell(point, &xi) == true);
    // Coordinates of the point in an unit cell
    point_unit_cell << 0., 0.;
    // Use analytical solution
    local_point = cell->local_coordinates_point_2d(point);
    for (unsigned i = 0; i < local_point.size(); ++i)
      REQUIRE(local_point[i] == Approx(point_unit_cell[i]).epsilon(Tolerance));

    // Coordinates of a point in real cell
    point << -1.5, -1.5;
    // Test if point is in cell
    REQUIRE(cell->is_point_in_cell(point, &xi) == true);
    // Coordinates of the point in an unit cell
    point_unit_cell << -0.75, -0.75;
    // Use analytical solution
    local_point = cell->local_coordinates_point_2d(point);
    for (unsigned i = 0; i < local_point.size(); ++i)
      REQUIRE(local_point[i] == Approx(point_unit_cell[i]).epsilon(Tolerance));
  }

  SECTION("Transform real to unit cell analytical solution") {
    // Number of nodes in cell
    const unsigned Nnodes = 4;

    // Coordinates
    Eigen::Vector2d coords;

    coords << 2.0, 1.0;
    std::shared_ptr<mpm::NodeBase<Dim>> node0 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

    coords << 4.0, 2.0;
    std::shared_ptr<mpm::NodeBase<Dim>> node1 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);

    coords << 2.0, 4.0;
    std::shared_ptr<mpm::NodeBase<Dim>> node2 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);

    coords << 1.0, 3.0;
    std::shared_ptr<mpm::NodeBase<Dim>> node3 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);

    // 4-noded quadrilateral shape functions
    std::shared_ptr<mpm::Element<Dim>> element =
        Factory<mpm::Element<Dim>>::instance()->create("ED2Q4");

    mpm::Index id = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

    REQUIRE(cell->add_node(0, node0) == true);
    REQUIRE(cell->add_node(1, node1) == true);
    REQUIRE(cell->add_node(2, node2) == true);
    REQUIRE(cell->add_node(3, node3) == true);
    REQUIRE(cell->nnodes() == 4);

    // Initialise cell
    REQUIRE(cell->initialise() == true);

    Eigen::Vector2d xi;

    // Coordinates of a point in real cell
    Eigen::Vector2d point;
    point << 2.1875, 3.25;

    // Test if point is in cell
    REQUIRE(cell->is_point_in_cell(point, &xi) == true);

    // Coordinates of the point in an unit cell
    Eigen::Matrix<double, 2, 1> point_unit_cell;
    point_unit_cell << 0.5, 0.5;

    // Use Newton-raphson iteration to find local coordinates
    auto local_point = cell->transform_real_to_unit_cell(point);
    for (unsigned i = 0; i < local_point.size(); ++i)
      REQUIRE(local_point[i] == Approx(point_unit_cell[i]).epsilon(Tolerance));

    // Use analytical solution
    local_point = cell->local_coordinates_point_2d(point);
    for (unsigned i = 0; i < local_point.size(); ++i)
      REQUIRE(local_point[i] == Approx(point_unit_cell[i]).epsilon(Tolerance));

    // Coordinates of a point in real cell
    point << 3., 3.;
    // Coordinates of the point in an unit cell
    point_unit_cell << 1., 0.;
    // Test if point is in cell
    REQUIRE(cell->is_point_in_cell(point, &xi) == true);
    // Use analytical solution
    local_point = cell->local_coordinates_point_2d(point);
    for (unsigned i = 0; i < local_point.size(); ++i)
      REQUIRE(local_point[i] == Approx(point_unit_cell[i]).epsilon(Tolerance));
  }

  SECTION("Check point in unit cell analytical solution") {
    // Number of nodes in cell
    const unsigned Nnodes = 4;

    // Coordinates
    Eigen::Vector2d coords;

    coords << 0.375, 0.1063;
    std::shared_ptr<mpm::NodeBase<Dim>> node0 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

    coords << 0.38, 0.10325;
    std::shared_ptr<mpm::NodeBase<Dim>> node1 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);

    coords << 0.38, 0.10825;
    std::shared_ptr<mpm::NodeBase<Dim>> node2 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);

    coords << 0.375, 0.1113;
    std::shared_ptr<mpm::NodeBase<Dim>> node3 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);

    coords << 0.375, 0.1013;
    std::shared_ptr<mpm::NodeBase<Dim>> node4 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(4, coords);

    coords << 0.38, 0.098247997;
    std::shared_ptr<mpm::NodeBase<Dim>> node5 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(5, coords);

    // 4-noded quadrilateral shape functions
    std::shared_ptr<mpm::Element<Dim>> element =
        Factory<mpm::Element<Dim>>::instance()->create("ED2Q4");

    Eigen::Vector2d xi;
    // Coordinates of a point in real cell
    Eigen::Vector2d point;
    point << 0.377632680795176, 0.104694061235587;

    // Cell 1
    mpm::Index id = 0;
    auto cell1 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

    REQUIRE(cell1->add_node(0, node0) == true);
    REQUIRE(cell1->add_node(1, node1) == true);
    REQUIRE(cell1->add_node(2, node2) == true);
    REQUIRE(cell1->add_node(3, node3) == true);
    REQUIRE(cell1->nnodes() == 4);

    // Initialise cell
    REQUIRE(cell1->initialise() == true);
    // Test if point is in cell
    REQUIRE(cell1->is_point_in_cell(point, &xi) == false);

    // Cell 2
    mpm::Index id2 = 1;
    auto cell2 = std::make_shared<mpm::Cell<Dim>>(id2, Nnodes, element);

    REQUIRE(cell2->add_node(0, node3) == true);
    REQUIRE(cell2->add_node(1, node2) == true);
    REQUIRE(cell2->add_node(2, node5) == true);
    REQUIRE(cell2->add_node(3, node4) == true);
    REQUIRE(cell2->nnodes() == 4);

    // Initialise cell
    REQUIRE(cell2->initialise() == true);
    // Test if point is in cell
    REQUIRE(cell2->is_point_in_cell(point, &xi) == true);

    // Coordinates of the point in an unit cell
    Eigen::Matrix<double, 2, 1> point_unit_cell;
    point_unit_cell << 0.0530723, -0.000104758;

    // Use Newton-raphson iteration to find local coordinates
    auto local_point = cell2->transform_real_to_unit_cell(point);
    for (unsigned i = 0; i < local_point.size(); ++i)
      REQUIRE(local_point[i] == Approx(point_unit_cell[i]).epsilon(Tolerance));

    // Use analytical solution
    // REQUIRE(cell2->point_in_2d_cell(point) == true);
    local_point = cell2->local_coordinates_point_2d(point);
    for (unsigned i = 0; i < local_point.size(); ++i)
      REQUIRE(local_point[i] == Approx(point_unit_cell[i]).epsilon(Tolerance));
  }
}

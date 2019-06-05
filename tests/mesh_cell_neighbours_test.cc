#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

#include "cell.h"
#include "element.h"
#include "factory.h"
#include "hexahedron_element.h"
#include "hexahedron_quadrature.h"
#include "mesh.h"
#include "node.h"
#include "quadrilateral_element.h"
#include "quadrilateral_quadrature.h"

TEST_CASE("Mesh cell neighbours 2D", "[MeshCell][2D]") {
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

  SECTION("Mesh cell neighbours 2D") {
    // Number of nodes in cell
    const unsigned Nnodes = 4;

    // Coordinates
    Eigen::Vector2d coords;

    coords << 0., 0.;
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

    // 4-noded quadrilateral shape functions
    std::shared_ptr<mpm::Element<Dim>> element =
        Factory<mpm::Element<Dim>>::instance()->create("ED2Q4");

    mpm::Index id = 0;
    auto cell0 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

    REQUIRE(cell0->add_node(0, node0) == true);
    REQUIRE(cell0->add_node(1, node1) == true);
    REQUIRE(cell0->add_node(2, node2) == true);
    REQUIRE(cell0->add_node(3, node3) == true);
    REQUIRE(cell0->nnodes() == 4);

    // Initialise cell
    REQUIRE(cell0->initialise() == true);

    // Check face node ids of cell 0
    auto fnodes = cell0->sorted_face_node_ids();
    std::vector<std::vector<mpm::Index>> fnodes_check = {
        {0, 1}, {1, 2}, {2, 3}, {0, 3}};
    for (unsigned i = 0; i < fnodes_check.size(); ++i)
      for (unsigned j = 0; j < fnodes_check.at(i).size(); ++j)
        REQUIRE(fnodes.at(i).at(j) == fnodes_check.at(i).at(j));

    id = 1;
    auto cell1 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

    REQUIRE(cell1->add_node(0, node1) == true);
    REQUIRE(cell1->add_node(1, node4) == true);
    REQUIRE(cell1->add_node(2, node5) == true);
    REQUIRE(cell1->add_node(3, node2) == true);
    REQUIRE(cell1->nnodes() == 4);

    // Initialise cell
    REQUIRE(cell1->initialise() == true);

    // Check face node ids of cell 1
    fnodes = cell1->sorted_face_node_ids();
    fnodes_check = {{1, 4}, {4, 5}, {2, 5}, {1, 2}};
    for (unsigned i = 0; i < fnodes_check.size(); ++i)
      for (unsigned j = 0; j < fnodes_check.at(i).size(); ++j)
        REQUIRE(fnodes.at(i).at(j) == fnodes_check.at(i).at(j));
  }
}

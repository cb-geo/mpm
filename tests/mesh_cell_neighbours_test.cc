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
  // Number of nodes per cell
  const unsigned Nnodes = 4;
  // Tolerance
  const double Tolerance = 1.E-7;

  SECTION("Mesh cell neighbours 2D") {
    // Number of nodes in cell
    const unsigned Nnodes = 4;

    auto mesh = std::make_shared<mpm::Mesh<Dim>>(0);
    // Check mesh is active
    REQUIRE(mesh->status() == false);

    // Coordinates
    Eigen::Vector2d coords;

    coords << 0., 0.;
    std::shared_ptr<mpm::NodeBase<Dim>> node0 =
        std::make_shared<mpm::Node<Dim, Dof>>(0, coords);

    coords << 2., 0.;
    std::shared_ptr<mpm::NodeBase<Dim>> node1 =
        std::make_shared<mpm::Node<Dim, Dof>>(1, coords);

    coords << 2., 2.;
    std::shared_ptr<mpm::NodeBase<Dim>> node2 =
        std::make_shared<mpm::Node<Dim, Dof>>(2, coords);

    coords << 0., 2.;
    std::shared_ptr<mpm::NodeBase<Dim>> node3 =
        std::make_shared<mpm::Node<Dim, Dof>>(3, coords);

    coords << 4., 0.;
    std::shared_ptr<mpm::NodeBase<Dim>> node4 =
        std::make_shared<mpm::Node<Dim, Dof>>(4, coords);

    coords << 4., 2.;
    std::shared_ptr<mpm::NodeBase<Dim>> node5 =
        std::make_shared<mpm::Node<Dim, Dof>>(5, coords);

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

    // Initialise cell and to mesh
    REQUIRE(cell0->initialise() == true);
    REQUIRE(mesh->add_cell(cell0) == true);

    id = 1;
    auto cell1 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

    REQUIRE(cell1->add_node(0, node1) == true);
    REQUIRE(cell1->add_node(1, node4) == true);
    REQUIRE(cell1->add_node(2, node5) == true);
    REQUIRE(cell1->add_node(3, node2) == true);
    REQUIRE(cell1->nnodes() == 4);

    // Initialise cell and add to mesh
    REQUIRE(cell1->initialise() == true);
    REQUIRE(mesh->add_cell(cell1) == true);

    SECTION("Test face node neighbour ids") {
      // Check face node ids of cell 0
      auto fnodes = cell0->sorted_face_node_ids();
      REQUIRE(fnodes.size() == 4);

      std::vector<std::vector<mpm::Index>> fnodes_check = {
          {0, 1}, {1, 2}, {2, 3}, {0, 3}};
      for (unsigned i = 0; i < fnodes_check.size(); ++i)
        for (unsigned j = 0; j < fnodes_check.at(i).size(); ++j)
          REQUIRE(fnodes.at(i).at(j) == fnodes_check.at(i).at(j));

      // Check face node ids of cell 1
      fnodes = cell1->sorted_face_node_ids();
      REQUIRE(fnodes.size() == 4);

      // Assign neighbours
      fnodes_check = {{1, 4}, {4, 5}, {2, 5}, {1, 2}};
      for (unsigned i = 0; i < fnodes_check.size(); ++i)
        for (unsigned j = 0; j < fnodes_check.at(i).size(); ++j)
          REQUIRE(fnodes.at(i).at(j) == fnodes_check.at(i).at(j));
    }

    SECTION("Check add neighbours") {
      // Add neighbours to cell 0
      REQUIRE(cell0->nneighbours() == 0);
      REQUIRE(cell0->add_neighbour(1) == true);
      REQUIRE(cell0->add_neighbour(0) == false);
      REQUIRE(cell0->add_neighbour(0) == false);
      REQUIRE(cell0->nneighbours() == 1);
      for (auto n : cell0->neighbours()) REQUIRE(n == 1);

      // Add neighbours to cell 1
      REQUIRE(cell1->nneighbours() == 0);
      REQUIRE(cell1->add_neighbour(0) == true);
      REQUIRE(cell1->add_neighbour(0) == false);
      REQUIRE(cell1->add_neighbour(1) == false);
      REQUIRE(cell1->nneighbours() == 1);
      for (auto n : cell1->neighbours()) REQUIRE(n == 0);
    }

    SECTION("Compute cell neighbours") {
      mesh->compute_cell_neighbours();
      REQUIRE(cell0->nneighbours() == 1);
      REQUIRE(cell1->nneighbours() == 1);
      for (auto n : cell0->neighbours()) REQUIRE(n == 1);
      for (auto n : cell1->neighbours()) REQUIRE(n == 0);

      SECTION("Locate particles in mesh") {
        coords << 3., 1.5;
        std::shared_ptr<mpm::ParticleBase<Dim>> particle1 =
            std::make_shared<mpm::Particle<Dim>>(1, coords);
        // Add particle 1 and check
        REQUIRE(mesh->add_particle(particle1, false) == true);

        // Assign incorrect cell id of 0
        REQUIRE(particle1->assign_cell_id(0) == true);

        // Locate particles in a mesh
        auto particles = mesh->locate_particles_mesh();

        // Should find all particles in mesh
        REQUIRE(particles.size() == 0);
      }
    }
  }
}

TEST_CASE("Mesh cell neighbours 3D", "[MeshCell][3D]") {
  // Dimension
  const unsigned Dim = 3;
  // Degrees of freedom
  const unsigned Dof = 6;
  // Number of nodes per cell
  const unsigned Nnodes = 8;
  // Tolerance
  const double Tolerance = 1.E-9;

  // 8-noded hexahedron element
  std::shared_ptr<mpm::Element<Dim>> element =
      Factory<mpm::Element<Dim>>::instance()->create("ED3H8");

  SECTION("Mesh cell neighbours 3D") {

    auto mesh = std::make_shared<mpm::Mesh<Dim>>(0);
    // Check mesh is active
    REQUIRE(mesh->status() == false);

    // Define nodes
    Eigen::Vector3d coords;
    coords << 0, 0, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node0 =
        std::make_shared<mpm::Node<Dim, Dof>>(0, coords);
    REQUIRE(mesh->add_node(node0) == true);

    coords << 2, 0, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node1 =
        std::make_shared<mpm::Node<Dim, Dof>>(1, coords);
    REQUIRE(mesh->add_node(node1) == true);

    coords << 2, 2, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node2 =
        std::make_shared<mpm::Node<Dim, Dof>>(2, coords);
    REQUIRE(mesh->add_node(node2) == true);

    coords << 0, 2, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node3 =
        std::make_shared<mpm::Node<Dim, Dof>>(3, coords);
    REQUIRE(mesh->add_node(node3) == true);

    coords << 0, 0, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node4 =
        std::make_shared<mpm::Node<Dim, Dof>>(4, coords);
    REQUIRE(mesh->add_node(node4) == true);

    coords << 2, 0, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node5 =
        std::make_shared<mpm::Node<Dim, Dof>>(5, coords);
    REQUIRE(mesh->add_node(node5) == true);

    coords << 2, 2, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node6 =
        std::make_shared<mpm::Node<Dim, Dof>>(6, coords);
    REQUIRE(mesh->add_node(node6) == true);

    coords << 0, 2, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node7 =
        std::make_shared<mpm::Node<Dim, Dof>>(7, coords);
    REQUIRE(mesh->add_node(node7) == true);

    // Create cell0
    auto cell0 = std::make_shared<mpm::Cell<Dim>>(0, Nnodes, element);

    // Add nodes to cell
    cell0->add_node(0, node0);
    cell0->add_node(1, node1);
    cell0->add_node(2, node2);
    cell0->add_node(3, node3);
    cell0->add_node(4, node4);
    cell0->add_node(5, node5);
    cell0->add_node(6, node6);
    cell0->add_node(7, node7);

    REQUIRE(cell0->nnodes() == 8);

    REQUIRE(mesh->add_cell(cell0) == true);

    // Cell 1
    coords << 4, 0, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node8 =
        std::make_shared<mpm::Node<Dim, Dof>>(8, coords);
    REQUIRE(mesh->add_node(node8) == true);

    coords << 4, 2, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node9 =
        std::make_shared<mpm::Node<Dim, Dof>>(9, coords);
    REQUIRE(mesh->add_node(node9) == true);

    coords << 4, 0, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node10 =
        std::make_shared<mpm::Node<Dim, Dof>>(10, coords);
    REQUIRE(mesh->add_node(node10) == true);

    coords << 4, 2, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node11 =
        std::make_shared<mpm::Node<Dim, Dof>>(11, coords);
    REQUIRE(mesh->add_node(node11) == true);

    // Create cell0
    auto cell1 = std::make_shared<mpm::Cell<Dim>>(1, Nnodes, element);

    // Add nodes to cell
    cell1->add_node(0, node1);
    cell1->add_node(1, node8);
    cell1->add_node(2, node9);
    cell1->add_node(3, node2);
    cell1->add_node(4, node5);
    cell1->add_node(5, node10);
    cell1->add_node(6, node11);
    cell1->add_node(7, node6);

    REQUIRE(cell1->nnodes() == 8);

    REQUIRE(mesh->add_cell(cell1) == true);

    SECTION("Check face node ids") {

      // Check face node ids of cell 0
      auto fnodes = cell0->sorted_face_node_ids();
      std::vector<std::vector<mpm::Index>> fnodes_check = {
          {0, 1, 4, 5}, {1, 2, 5, 6}, {2, 3, 6, 7},
          {0, 3, 4, 7}, {0, 1, 2, 3}, {4, 5, 6, 7}};

      for (unsigned i = 0; i < fnodes.size(); ++i)
        for (unsigned j = 0; j < fnodes.at(i).size(); ++j)
          REQUIRE(fnodes.at(i).at(j) == fnodes_check.at(i).at(j));

      // Check face node ids of cell 0
      fnodes = cell1->sorted_face_node_ids();
      fnodes_check = {{1, 5, 8, 10}, {8, 9, 10, 11}, {2, 6, 9, 11},
                      {1, 2, 5, 6},  {1, 2, 8, 9},   {5, 6, 10, 11}};
      for (unsigned i = 0; i < fnodes.size(); ++i)
        for (unsigned j = 0; j < fnodes.at(i).size(); ++j)
          REQUIRE(fnodes.at(i).at(j) == fnodes_check.at(i).at(j));
    }

    SECTION("Check assign neighbours") {
      // Add neighbours to cell 0
      REQUIRE(cell0->nneighbours() == 0);
      REQUIRE(cell0->add_neighbour(1) == true);
      REQUIRE(cell0->add_neighbour(0) == false);
      REQUIRE(cell0->add_neighbour(0) == false);
      REQUIRE(cell0->nneighbours() == 1);
      for (auto n : cell0->neighbours()) REQUIRE(n == 1);

      // Add neighbours to cell 1
      REQUIRE(cell1->nneighbours() == 0);
      REQUIRE(cell1->add_neighbour(0) == true);
      REQUIRE(cell1->add_neighbour(0) == false);
      REQUIRE(cell1->add_neighbour(1) == false);
      REQUIRE(cell1->nneighbours() == 1);
      for (auto n : cell1->neighbours()) REQUIRE(n == 0);
    }
    // Compute cell neighbours
    SECTION("Compute cell neighbours") {
      mesh->compute_cell_neighbours();
      REQUIRE(cell0->nneighbours() == 1);
      REQUIRE(cell1->nneighbours() == 1);
      for (auto n : cell0->neighbours()) REQUIRE(n == 1);
      for (auto n : cell1->neighbours()) REQUIRE(n == 0);

      REQUIRE(cell0->initialise() == true);
      REQUIRE(cell1->initialise() == true);

      SECTION("Locate particles in mesh") {
        coords << 3., 1.5, 1.5;
        std::shared_ptr<mpm::ParticleBase<Dim>> particle1 =
            std::make_shared<mpm::Particle<Dim>>(1, coords);
        // Add particle 1 and check
        REQUIRE(mesh->add_particle(particle1, false) == true);

        // Assign incorrect cell id of 0
        REQUIRE(particle1->assign_cell_id(0) == true);

        // Locate particles in a mesh
        auto particles = mesh->locate_particles_mesh();

        // Should find all particles in mesh
        REQUIRE(particles.size() == 0);
      }
    }
  }
}

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

    auto mesh = std::make_shared<mpm::Mesh<Dim>>(0);
    // Check mesh is active
    REQUIRE(mesh->status() == false);

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

    coords << 4., 4.;
    std::shared_ptr<mpm::NodeBase<Dim>> node6 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(6, coords);

    coords << 2., 4.;
    std::shared_ptr<mpm::NodeBase<Dim>> node7 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(7, coords);

    coords << 0., 4.;
    std::shared_ptr<mpm::NodeBase<Dim>> node8 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(8, coords);

    coords << 6., 0.;
    std::shared_ptr<mpm::NodeBase<Dim>> node9 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(9, coords);

    coords << 6., 2.;
    std::shared_ptr<mpm::NodeBase<Dim>> node10 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(10, coords);

    coords << 6., 4.;
    std::shared_ptr<mpm::NodeBase<Dim>> node11 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(11, coords);

    coords << 6., 6.;
    std::shared_ptr<mpm::NodeBase<Dim>> node12 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(12, coords);

    coords << 4., 6.;
    std::shared_ptr<mpm::NodeBase<Dim>> node13 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(13, coords);

    coords << 2., 6.;
    std::shared_ptr<mpm::NodeBase<Dim>> node14 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(14, coords);

    coords << 0., 6.;
    std::shared_ptr<mpm::NodeBase<Dim>> node15 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(15, coords);

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

    id = 2;
    auto cell2 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

    REQUIRE(cell2->add_node(0, node2) == true);
    REQUIRE(cell2->add_node(1, node5) == true);
    REQUIRE(cell2->add_node(2, node6) == true);
    REQUIRE(cell2->add_node(3, node7) == true);
    REQUIRE(cell2->nnodes() == 4);

    // Initialise cell and add to mesh
    REQUIRE(cell2->initialise() == true);
    REQUIRE(mesh->add_cell(cell2) == true);

    id = 3;
    auto cell3 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

    REQUIRE(cell3->add_node(0, node3) == true);
    REQUIRE(cell3->add_node(1, node2) == true);
    REQUIRE(cell3->add_node(2, node7) == true);
    REQUIRE(cell3->add_node(3, node8) == true);
    REQUIRE(cell3->nnodes() == 4);

    // Initialise cell and add to mesh
    REQUIRE(cell3->initialise() == true);
    REQUIRE(mesh->add_cell(cell3) == true);

    id = 4;
    auto cell4 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

    REQUIRE(cell4->add_node(0, node4) == true);
    REQUIRE(cell4->add_node(1, node9) == true);
    REQUIRE(cell4->add_node(2, node10) == true);
    REQUIRE(cell4->add_node(3, node5) == true);
    REQUIRE(cell4->nnodes() == 4);

    // Initialise cell and add to mesh
    REQUIRE(cell4->initialise() == true);
    REQUIRE(mesh->add_cell(cell4) == true);

    id = 5;
    auto cell5 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

    REQUIRE(cell5->add_node(0, node5) == true);
    REQUIRE(cell5->add_node(1, node10) == true);
    REQUIRE(cell5->add_node(2, node11) == true);
    REQUIRE(cell5->add_node(3, node6) == true);
    REQUIRE(cell5->nnodes() == 4);

    // Initialise cell and add to mesh
    REQUIRE(cell5->initialise() == true);
    REQUIRE(mesh->add_cell(cell5) == true);

    id = 6;
    auto cell6 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

    REQUIRE(cell6->add_node(0, node6) == true);
    REQUIRE(cell6->add_node(1, node11) == true);
    REQUIRE(cell6->add_node(2, node12) == true);
    REQUIRE(cell6->add_node(3, node13) == true);
    REQUIRE(cell6->nnodes() == 4);

    // Initialise cell and add to mesh
    REQUIRE(cell6->initialise() == true);
    REQUIRE(mesh->add_cell(cell6) == true);

    id = 7;
    auto cell7 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

    REQUIRE(cell7->add_node(0, node7) == true);
    REQUIRE(cell7->add_node(1, node6) == true);
    REQUIRE(cell7->add_node(2, node13) == true);
    REQUIRE(cell7->add_node(3, node14) == true);
    REQUIRE(cell7->nnodes() == 4);

    // Initialise cell and add to mesh
    REQUIRE(cell7->initialise() == true);
    REQUIRE(mesh->add_cell(cell7) == true);

    id = 8;
    auto cell8 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

    REQUIRE(cell8->add_node(0, node8) == true);
    REQUIRE(cell8->add_node(1, node7) == true);
    REQUIRE(cell8->add_node(2, node14) == true);
    REQUIRE(cell8->add_node(3, node15) == true);
    REQUIRE(cell8->nnodes() == 4);

    // Initialise cell and add to mesh
    REQUIRE(cell8->initialise() == true);
    REQUIRE(mesh->add_cell(cell8) == true);

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
      REQUIRE(cell0->nneighbours() == 3);
      REQUIRE(cell1->nneighbours() == 5);
      REQUIRE(cell2->nneighbours() == 8);
      REQUIRE(cell3->nneighbours() == 5);
      REQUIRE(cell4->nneighbours() == 3);
      REQUIRE(cell5->nneighbours() == 5);
      REQUIRE(cell6->nneighbours() == 3);
      REQUIRE(cell7->nneighbours() == 5);
      REQUIRE(cell8->nneighbours() == 3);

      // Check solutions
      std::set<mpm::Index> n0 = {1, 2, 3};
      std::set<mpm::Index> n1 = {0, 2, 3, 4, 5};
      std::set<mpm::Index> n2 = {0, 1, 3, 4, 5, 6, 7, 8};
      std::set<mpm::Index> n3 = {0, 1, 2, 7, 8};
      std::set<mpm::Index> n4 = {1, 2, 5};
      std::set<mpm::Index> n5 = {1, 2, 4, 6, 7};
      std::set<mpm::Index> n6 = {2, 5, 7};
      std::set<mpm::Index> n7 = {2, 3, 5, 6, 8};
      std::set<mpm::Index> n8 = {2, 3, 7};

      REQUIRE(cell0->neighbours() == n0);
      REQUIRE(cell1->neighbours() == n1);
      REQUIRE(cell2->neighbours() == n2);
      REQUIRE(cell3->neighbours() == n3);
      REQUIRE(cell4->neighbours() == n4);
      REQUIRE(cell5->neighbours() == n5);
      REQUIRE(cell6->neighbours() == n6);
      REQUIRE(cell7->neighbours() == n7);
      REQUIRE(cell8->neighbours() == n8);

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
  // Number of phases
  const unsigned Nphases = 1;
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
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);
    REQUIRE(mesh->add_node(node0) == true);

    coords << 2, 0, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node1 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);
    REQUIRE(mesh->add_node(node1) == true);

    coords << 2, 2, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node2 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);
    REQUIRE(mesh->add_node(node2) == true);

    coords << 0, 2, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node3 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);
    REQUIRE(mesh->add_node(node3) == true);

    coords << 0, 0, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node4 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(4, coords);
    REQUIRE(mesh->add_node(node4) == true);

    coords << 2, 0, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node5 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(5, coords);
    REQUIRE(mesh->add_node(node5) == true);

    coords << 2, 2, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node6 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(6, coords);
    REQUIRE(mesh->add_node(node6) == true);

    coords << 0, 2, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node7 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(7, coords);
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
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(8, coords);
    REQUIRE(mesh->add_node(node8) == true);

    coords << 4, 2, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node9 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(9, coords);
    REQUIRE(mesh->add_node(node9) == true);

    coords << 4, 0, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node10 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(10, coords);
    REQUIRE(mesh->add_node(node10) == true);

    coords << 4, 2, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node11 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(11, coords);
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

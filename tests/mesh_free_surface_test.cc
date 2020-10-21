#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

#include "cell.h"
#include "element.h"
#include "factory.h"
#include "hexahedron_element.h"
#include "hexahedron_quadrature.h"
#include "material.h"
#include "mesh.h"
#include "node.h"
#include "quadrilateral_element.h"
#include "quadrilateral_quadrature.h"

TEST_CASE("Mesh free surface 2D", "[MeshCell][2D][free_surface]") {
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

  // Initialise material
  Json jmaterial;
  jmaterial["density"] = 1000.;
  jmaterial["bulk_modulus"] = 8333333.333333333;
  jmaterial["dynamic_viscosity"] = 8.9E-4;
  auto material =
      Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
          "Newtonian2D", std::move(0), jmaterial);

  auto mesh = std::make_shared<mpm::Mesh<Dim>>(0);
  // Check mesh is active
  REQUIRE(mesh->status() == false);

  // Coordinates
  Eigen::Vector2d coords;

  coords << 0., 0.;
  std::shared_ptr<mpm::NodeBase<Dim>> node0 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);
  REQUIRE(mesh->add_node(node0) == true);

  coords << 2., 0.;
  std::shared_ptr<mpm::NodeBase<Dim>> node1 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);
  REQUIRE(mesh->add_node(node1) == true);

  coords << 4., 0.;
  std::shared_ptr<mpm::NodeBase<Dim>> node2 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);
  REQUIRE(mesh->add_node(node2) == true);

  coords << 6., 0.;
  std::shared_ptr<mpm::NodeBase<Dim>> node3 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);
  REQUIRE(mesh->add_node(node3) == true);

  coords << 8., 0.;
  std::shared_ptr<mpm::NodeBase<Dim>> node4 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(4, coords);
  REQUIRE(mesh->add_node(node4) == true);

  coords << 10., 0.;
  std::shared_ptr<mpm::NodeBase<Dim>> node5 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(5, coords);
  REQUIRE(mesh->add_node(node5) == true);

  coords << 0., 2.;
  std::shared_ptr<mpm::NodeBase<Dim>> node6 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(6, coords);
  REQUIRE(mesh->add_node(node6) == true);

  coords << 2., 2.;
  std::shared_ptr<mpm::NodeBase<Dim>> node7 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(7, coords);
  REQUIRE(mesh->add_node(node7) == true);

  coords << 4., 2.;
  std::shared_ptr<mpm::NodeBase<Dim>> node8 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(8, coords);
  REQUIRE(mesh->add_node(node8) == true);

  coords << 6., 2.;
  std::shared_ptr<mpm::NodeBase<Dim>> node9 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(9, coords);
  REQUIRE(mesh->add_node(node9) == true);

  coords << 8., 2.;
  std::shared_ptr<mpm::NodeBase<Dim>> node10 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(10, coords);
  REQUIRE(mesh->add_node(node10) == true);

  coords << 10., 2.;
  std::shared_ptr<mpm::NodeBase<Dim>> node11 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(11, coords);
  REQUIRE(mesh->add_node(node11) == true);

  coords << 0., 4.;
  std::shared_ptr<mpm::NodeBase<Dim>> node12 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(12, coords);
  REQUIRE(mesh->add_node(node12) == true);

  coords << 2., 4.;
  std::shared_ptr<mpm::NodeBase<Dim>> node13 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(13, coords);
  REQUIRE(mesh->add_node(node13) == true);

  coords << 4., 4.;
  std::shared_ptr<mpm::NodeBase<Dim>> node14 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(14, coords);
  REQUIRE(mesh->add_node(node14) == true);

  coords << 6., 4.;
  std::shared_ptr<mpm::NodeBase<Dim>> node15 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(15, coords);
  REQUIRE(mesh->add_node(node15) == true);

  coords << 8., 4.;
  std::shared_ptr<mpm::NodeBase<Dim>> node16 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(16, coords);
  REQUIRE(mesh->add_node(node16) == true);

  coords << 10., 4.;
  std::shared_ptr<mpm::NodeBase<Dim>> node17 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(17, coords);
  REQUIRE(mesh->add_node(node17) == true);

  coords << 0., 6.;
  std::shared_ptr<mpm::NodeBase<Dim>> node18 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(18, coords);
  REQUIRE(mesh->add_node(node18) == true);

  coords << 2., 6.;
  std::shared_ptr<mpm::NodeBase<Dim>> node19 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(19, coords);
  REQUIRE(mesh->add_node(node19) == true);

  coords << 4., 6.;
  std::shared_ptr<mpm::NodeBase<Dim>> node20 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(20, coords);
  REQUIRE(mesh->add_node(node20) == true);

  coords << 6., 6.;
  std::shared_ptr<mpm::NodeBase<Dim>> node21 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(21, coords);
  REQUIRE(mesh->add_node(node21) == true);

  coords << 8., 6.;
  std::shared_ptr<mpm::NodeBase<Dim>> node22 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(22, coords);
  REQUIRE(mesh->add_node(node22) == true);

  coords << 10., 6.;
  std::shared_ptr<mpm::NodeBase<Dim>> node23 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(23, coords);
  REQUIRE(mesh->add_node(node23) == true);

  coords << 0., 8.;
  std::shared_ptr<mpm::NodeBase<Dim>> node24 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(24, coords);
  REQUIRE(mesh->add_node(node24) == true);

  coords << 2., 8.;
  std::shared_ptr<mpm::NodeBase<Dim>> node25 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(25, coords);
  REQUIRE(mesh->add_node(node25) == true);

  coords << 4., 8.;
  std::shared_ptr<mpm::NodeBase<Dim>> node26 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(26, coords);
  REQUIRE(mesh->add_node(node26) == true);

  coords << 6., 8.;
  std::shared_ptr<mpm::NodeBase<Dim>> node27 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(27, coords);
  REQUIRE(mesh->add_node(node27) == true);

  coords << 8., 8.;
  std::shared_ptr<mpm::NodeBase<Dim>> node28 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(28, coords);
  REQUIRE(mesh->add_node(node28) == true);

  coords << 10., 8.;
  std::shared_ptr<mpm::NodeBase<Dim>> node29 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(29, coords);
  REQUIRE(mesh->add_node(node29) == true);

  coords << 0., 10.;
  std::shared_ptr<mpm::NodeBase<Dim>> node30 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(30, coords);
  REQUIRE(mesh->add_node(node30) == true);

  coords << 2., 10.;
  std::shared_ptr<mpm::NodeBase<Dim>> node31 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(31, coords);
  REQUIRE(mesh->add_node(node31) == true);

  coords << 4., 10.;
  std::shared_ptr<mpm::NodeBase<Dim>> node32 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(32, coords);
  REQUIRE(mesh->add_node(node32) == true);

  coords << 6., 10.;
  std::shared_ptr<mpm::NodeBase<Dim>> node33 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(33, coords);
  REQUIRE(mesh->add_node(node33) == true);

  coords << 8., 10.;
  std::shared_ptr<mpm::NodeBase<Dim>> node34 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(34, coords);
  REQUIRE(mesh->add_node(node34) == true);

  coords << 10., 10.;
  std::shared_ptr<mpm::NodeBase<Dim>> node35 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(35, coords);
  REQUIRE(mesh->add_node(node35) == true);

  // 4-noded quadrilateral shape functions
  std::shared_ptr<mpm::Element<Dim>> element =
      Factory<mpm::Element<Dim>>::instance()->create("ED2Q4");

  mpm::Index id = 0;
  auto cell0 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

  REQUIRE(cell0->add_node(0, node0) == true);
  REQUIRE(cell0->add_node(1, node1) == true);
  REQUIRE(cell0->add_node(2, node7) == true);
  REQUIRE(cell0->add_node(3, node6) == true);
  REQUIRE(cell0->nnodes() == 4);

  // Initialise cell and to mesh
  REQUIRE(cell0->initialise() == true);
  REQUIRE(mesh->add_cell(cell0) == true);

  id = 1;
  auto cell1 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

  REQUIRE(cell1->add_node(0, node1) == true);
  REQUIRE(cell1->add_node(1, node2) == true);
  REQUIRE(cell1->add_node(2, node8) == true);
  REQUIRE(cell1->add_node(3, node7) == true);
  REQUIRE(cell1->nnodes() == 4);

  // Initialise cell and add to mesh
  REQUIRE(cell1->initialise() == true);
  REQUIRE(mesh->add_cell(cell1) == true);

  id = 2;
  auto cell2 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

  REQUIRE(cell2->add_node(0, node2) == true);
  REQUIRE(cell2->add_node(1, node3) == true);
  REQUIRE(cell2->add_node(2, node9) == true);
  REQUIRE(cell2->add_node(3, node8) == true);
  REQUIRE(cell2->nnodes() == 4);

  // Initialise cell and add to mesh
  REQUIRE(cell2->initialise() == true);
  REQUIRE(mesh->add_cell(cell2) == true);

  id = 3;
  auto cell3 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

  REQUIRE(cell3->add_node(0, node3) == true);
  REQUIRE(cell3->add_node(1, node4) == true);
  REQUIRE(cell3->add_node(2, node10) == true);
  REQUIRE(cell3->add_node(3, node9) == true);
  REQUIRE(cell3->nnodes() == 4);

  // Initialise cell and add to mesh
  REQUIRE(cell3->initialise() == true);
  REQUIRE(mesh->add_cell(cell3) == true);

  id = 4;
  auto cell4 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

  REQUIRE(cell4->add_node(0, node4) == true);
  REQUIRE(cell4->add_node(1, node5) == true);
  REQUIRE(cell4->add_node(2, node11) == true);
  REQUIRE(cell4->add_node(3, node10) == true);
  REQUIRE(cell4->nnodes() == 4);

  // Initialise cell and add to mesh
  REQUIRE(cell4->initialise() == true);
  REQUIRE(mesh->add_cell(cell4) == true);

  id = 5;
  auto cell5 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

  REQUIRE(cell5->add_node(0, node6) == true);
  REQUIRE(cell5->add_node(1, node7) == true);
  REQUIRE(cell5->add_node(2, node13) == true);
  REQUIRE(cell5->add_node(3, node12) == true);
  REQUIRE(cell5->nnodes() == 4);

  // Initialise cell and add to mesh
  REQUIRE(cell5->initialise() == true);
  REQUIRE(mesh->add_cell(cell5) == true);

  id = 6;
  auto cell6 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

  REQUIRE(cell6->add_node(0, node7) == true);
  REQUIRE(cell6->add_node(1, node8) == true);
  REQUIRE(cell6->add_node(2, node14) == true);
  REQUIRE(cell6->add_node(3, node13) == true);
  REQUIRE(cell6->nnodes() == 4);

  // Initialise cell and add to mesh
  REQUIRE(cell6->initialise() == true);
  REQUIRE(mesh->add_cell(cell6) == true);

  id = 7;
  auto cell7 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

  REQUIRE(cell7->add_node(0, node8) == true);
  REQUIRE(cell7->add_node(1, node9) == true);
  REQUIRE(cell7->add_node(2, node15) == true);
  REQUIRE(cell7->add_node(3, node14) == true);
  REQUIRE(cell7->nnodes() == 4);

  // Initialise cell and add to mesh
  REQUIRE(cell7->initialise() == true);
  REQUIRE(mesh->add_cell(cell7) == true);

  id = 8;
  auto cell8 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

  REQUIRE(cell8->add_node(0, node9) == true);
  REQUIRE(cell8->add_node(1, node10) == true);
  REQUIRE(cell8->add_node(2, node16) == true);
  REQUIRE(cell8->add_node(3, node15) == true);
  REQUIRE(cell8->nnodes() == 4);

  // Initialise cell and add to mesh
  REQUIRE(cell8->initialise() == true);
  REQUIRE(mesh->add_cell(cell8) == true);

  id = 9;
  auto cell9 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

  REQUIRE(cell9->add_node(0, node10) == true);
  REQUIRE(cell9->add_node(1, node11) == true);
  REQUIRE(cell9->add_node(2, node17) == true);
  REQUIRE(cell9->add_node(3, node16) == true);
  REQUIRE(cell9->nnodes() == 4);

  // Initialise cell and add to mesh
  REQUIRE(cell9->initialise() == true);
  REQUIRE(mesh->add_cell(cell9) == true);

  id = 10;
  auto cell10 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

  REQUIRE(cell10->add_node(0, node12) == true);
  REQUIRE(cell10->add_node(1, node13) == true);
  REQUIRE(cell10->add_node(2, node19) == true);
  REQUIRE(cell10->add_node(3, node18) == true);
  REQUIRE(cell10->nnodes() == 4);

  // Initialise cell and add to mesh
  REQUIRE(cell10->initialise() == true);
  REQUIRE(mesh->add_cell(cell10) == true);

  id = 11;
  auto cell11 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

  REQUIRE(cell11->add_node(0, node13) == true);
  REQUIRE(cell11->add_node(1, node14) == true);
  REQUIRE(cell11->add_node(2, node20) == true);
  REQUIRE(cell11->add_node(3, node19) == true);
  REQUIRE(cell11->nnodes() == 4);

  // Initialise cell and add to mesh
  REQUIRE(cell11->initialise() == true);
  REQUIRE(mesh->add_cell(cell11) == true);

  id = 12;
  auto cell12 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

  REQUIRE(cell12->add_node(0, node14) == true);
  REQUIRE(cell12->add_node(1, node15) == true);
  REQUIRE(cell12->add_node(2, node21) == true);
  REQUIRE(cell12->add_node(3, node20) == true);
  REQUIRE(cell12->nnodes() == 4);

  // Initialise cell and add to mesh
  REQUIRE(cell12->initialise() == true);
  REQUIRE(mesh->add_cell(cell12) == true);

  id = 13;
  auto cell13 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

  REQUIRE(cell13->add_node(0, node15) == true);
  REQUIRE(cell13->add_node(1, node16) == true);
  REQUIRE(cell13->add_node(2, node22) == true);
  REQUIRE(cell13->add_node(3, node21) == true);
  REQUIRE(cell13->nnodes() == 4);

  // Initialise cell and add to mesh
  REQUIRE(cell13->initialise() == true);
  REQUIRE(mesh->add_cell(cell13) == true);

  id = 14;
  auto cell14 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

  REQUIRE(cell14->add_node(0, node16) == true);
  REQUIRE(cell14->add_node(1, node17) == true);
  REQUIRE(cell14->add_node(2, node23) == true);
  REQUIRE(cell14->add_node(3, node22) == true);
  REQUIRE(cell14->nnodes() == 4);

  // Initialise cell and add to mesh
  REQUIRE(cell14->initialise() == true);
  REQUIRE(mesh->add_cell(cell14) == true);

  id = 15;
  auto cell15 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

  REQUIRE(cell15->add_node(0, node18) == true);
  REQUIRE(cell15->add_node(1, node19) == true);
  REQUIRE(cell15->add_node(2, node25) == true);
  REQUIRE(cell15->add_node(3, node24) == true);
  REQUIRE(cell15->nnodes() == 4);

  // Initialise cell and add to mesh
  REQUIRE(cell15->initialise() == true);
  REQUIRE(mesh->add_cell(cell15) == true);

  id = 16;
  auto cell16 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

  REQUIRE(cell16->add_node(0, node19) == true);
  REQUIRE(cell16->add_node(1, node20) == true);
  REQUIRE(cell16->add_node(2, node26) == true);
  REQUIRE(cell16->add_node(3, node25) == true);
  REQUIRE(cell16->nnodes() == 4);

  // Initialise cell and add to mesh
  REQUIRE(cell16->initialise() == true);
  REQUIRE(mesh->add_cell(cell16) == true);

  id = 17;
  auto cell17 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

  REQUIRE(cell17->add_node(0, node20) == true);
  REQUIRE(cell17->add_node(1, node21) == true);
  REQUIRE(cell17->add_node(2, node27) == true);
  REQUIRE(cell17->add_node(3, node26) == true);
  REQUIRE(cell17->nnodes() == 4);

  // Initialise cell and add to mesh
  REQUIRE(cell17->initialise() == true);
  REQUIRE(mesh->add_cell(cell17) == true);

  id = 18;
  auto cell18 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

  REQUIRE(cell18->add_node(0, node21) == true);
  REQUIRE(cell18->add_node(1, node22) == true);
  REQUIRE(cell18->add_node(2, node28) == true);
  REQUIRE(cell18->add_node(3, node27) == true);
  REQUIRE(cell18->nnodes() == 4);

  // Initialise cell and add to mesh
  REQUIRE(cell18->initialise() == true);
  REQUIRE(mesh->add_cell(cell18) == true);

  id = 19;
  auto cell19 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

  REQUIRE(cell19->add_node(0, node22) == true);
  REQUIRE(cell19->add_node(1, node23) == true);
  REQUIRE(cell19->add_node(2, node29) == true);
  REQUIRE(cell19->add_node(3, node28) == true);
  REQUIRE(cell19->nnodes() == 4);

  // Initialise cell and add to mesh
  REQUIRE(cell19->initialise() == true);
  REQUIRE(mesh->add_cell(cell19) == true);

  id = 20;
  auto cell20 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

  REQUIRE(cell20->add_node(0, node24) == true);
  REQUIRE(cell20->add_node(1, node25) == true);
  REQUIRE(cell20->add_node(2, node31) == true);
  REQUIRE(cell20->add_node(3, node30) == true);
  REQUIRE(cell20->nnodes() == 4);

  // Initialise cell and add to mesh
  REQUIRE(cell20->initialise() == true);
  REQUIRE(mesh->add_cell(cell20) == true);

  id = 21;
  auto cell21 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

  REQUIRE(cell21->add_node(0, node25) == true);
  REQUIRE(cell21->add_node(1, node26) == true);
  REQUIRE(cell21->add_node(2, node32) == true);
  REQUIRE(cell21->add_node(3, node31) == true);
  REQUIRE(cell21->nnodes() == 4);

  // Initialise cell and add to mesh
  REQUIRE(cell21->initialise() == true);
  REQUIRE(mesh->add_cell(cell21) == true);

  id = 22;
  auto cell22 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

  REQUIRE(cell22->add_node(0, node26) == true);
  REQUIRE(cell22->add_node(1, node27) == true);
  REQUIRE(cell22->add_node(2, node33) == true);
  REQUIRE(cell22->add_node(3, node32) == true);
  REQUIRE(cell22->nnodes() == 4);

  // Initialise cell and add to mesh
  REQUIRE(cell22->initialise() == true);
  REQUIRE(mesh->add_cell(cell22) == true);

  id = 23;
  auto cell23 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

  REQUIRE(cell23->add_node(0, node27) == true);
  REQUIRE(cell23->add_node(1, node28) == true);
  REQUIRE(cell23->add_node(2, node34) == true);
  REQUIRE(cell23->add_node(3, node33) == true);
  REQUIRE(cell23->nnodes() == 4);

  // Initialise cell and add to mesh
  REQUIRE(cell23->initialise() == true);
  REQUIRE(mesh->add_cell(cell23) == true);

  id = 24;
  auto cell24 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

  REQUIRE(cell24->add_node(0, node28) == true);
  REQUIRE(cell24->add_node(1, node29) == true);
  REQUIRE(cell24->add_node(2, node35) == true);
  REQUIRE(cell24->add_node(3, node34) == true);
  REQUIRE(cell24->nnodes() == 4);

  // Initialise cell and add to mesh
  REQUIRE(cell24->initialise() == true);
  REQUIRE(mesh->add_cell(cell24) == true);

  // Find cell neighbours
  mesh->find_cell_neighbours();

  std::shared_ptr<mpm::ParticleBase<Dim>> particle0, particle1, particle2,
      particle3, particle4, particle5, particle6, particle7, particle8,
      particle9, particle10, particle11, particle12, particle13, particle14,
      particle15;

  coords << 3.5, 3.5;
  particle0 = std::make_shared<mpm::Particle<Dim>>(0, coords);
  REQUIRE(mesh->add_particle(particle0, false) == true);

  coords << 4.5, 3.5;
  particle1 = std::make_shared<mpm::Particle<Dim>>(1, coords);
  REQUIRE(mesh->add_particle(particle1, false) == true);

  coords << 5.5, 3.5;
  particle2 = std::make_shared<mpm::Particle<Dim>>(2, coords);
  REQUIRE(mesh->add_particle(particle2, false) == true);

  coords << 6.5, 3.5;
  particle3 = std::make_shared<mpm::Particle<Dim>>(3, coords);
  REQUIRE(mesh->add_particle(particle3, false) == true);

  coords << 3.5, 4.5;
  particle4 = std::make_shared<mpm::Particle<Dim>>(4, coords);
  REQUIRE(mesh->add_particle(particle4, false) == true);

  coords << 4.5, 4.5;
  particle5 = std::make_shared<mpm::Particle<Dim>>(5, coords);
  REQUIRE(mesh->add_particle(particle5, false) == true);

  coords << 5.5, 4.5;
  particle6 = std::make_shared<mpm::Particle<Dim>>(6, coords);
  REQUIRE(mesh->add_particle(particle6, false) == true);

  coords << 6.5, 4.5;
  particle7 = std::make_shared<mpm::Particle<Dim>>(7, coords);
  REQUIRE(mesh->add_particle(particle7, false) == true);

  coords << 3.5, 5.5;
  particle8 = std::make_shared<mpm::Particle<Dim>>(8, coords);
  REQUIRE(mesh->add_particle(particle8, false) == true);

  coords << 4.5, 5.5;
  particle9 = std::make_shared<mpm::Particle<Dim>>(9, coords);
  REQUIRE(mesh->add_particle(particle9, false) == true);

  coords << 5.5, 5.5;
  particle10 = std::make_shared<mpm::Particle<Dim>>(10, coords);
  REQUIRE(mesh->add_particle(particle10, false) == true);

  coords << 6.5, 5.5;
  particle11 = std::make_shared<mpm::Particle<Dim>>(11, coords);
  REQUIRE(mesh->add_particle(particle11, false) == true);

  coords << 3.5, 6.5;
  particle12 = std::make_shared<mpm::Particle<Dim>>(12, coords);
  REQUIRE(mesh->add_particle(particle12, false) == true);

  coords << 4.5, 6.5;
  particle13 = std::make_shared<mpm::Particle<Dim>>(13, coords);
  REQUIRE(mesh->add_particle(particle13, false) == true);

  coords << 5.5, 6.5;
  particle14 = std::make_shared<mpm::Particle<Dim>>(14, coords);
  REQUIRE(mesh->add_particle(particle14, false) == true);

  coords << 6.5, 6.5;
  particle15 = std::make_shared<mpm::Particle<Dim>>(15, coords);
  REQUIRE(mesh->add_particle(particle15, false) == true);

  // Assign material to particles
  mesh->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Dim>::assign_material, std::placeholders::_1,
                material, 0));

  // Locate particles in a mesh
  auto particles = mesh->locate_particles_mesh();

  // Should find all particles in mesh
  REQUIRE(particles.size() == 0);

  // Check particles inside cells
  REQUIRE(cell0->nparticles() == 0);
  REQUIRE(cell1->nparticles() == 0);
  REQUIRE(cell2->nparticles() == 0);
  REQUIRE(cell3->nparticles() == 0);
  REQUIRE(cell4->nparticles() == 0);
  REQUIRE(cell5->nparticles() == 0);
  REQUIRE(cell6->nparticles() == 1);
  REQUIRE(cell7->nparticles() == 2);
  REQUIRE(cell8->nparticles() == 1);
  REQUIRE(cell9->nparticles() == 0);
  REQUIRE(cell10->nparticles() == 0);
  REQUIRE(cell11->nparticles() == 2);
  REQUIRE(cell12->nparticles() == 4);
  REQUIRE(cell13->nparticles() == 2);
  REQUIRE(cell14->nparticles() == 0);
  REQUIRE(cell15->nparticles() == 0);
  REQUIRE(cell16->nparticles() == 1);
  REQUIRE(cell17->nparticles() == 2);
  REQUIRE(cell18->nparticles() == 1);
  REQUIRE(cell19->nparticles() == 0);
  REQUIRE(cell20->nparticles() == 0);
  REQUIRE(cell21->nparticles() == 0);
  REQUIRE(cell22->nparticles() == 0);
  REQUIRE(cell23->nparticles() == 0);
  REQUIRE(cell24->nparticles() == 0);

  // Find particle neighbours
  mesh->find_particle_neighbours();

  // Check particle neighbours
  REQUIRE(particle0->nneighbours() == 8);
  REQUIRE(particle1->nneighbours() == 11);
  REQUIRE(particle2->nneighbours() == 11);
  REQUIRE(particle3->nneighbours() == 8);
  REQUIRE(particle4->nneighbours() == 11);
  REQUIRE(particle5->nneighbours() == 15);
  REQUIRE(particle6->nneighbours() == 15);
  REQUIRE(particle7->nneighbours() == 11);
  REQUIRE(particle8->nneighbours() == 11);
  REQUIRE(particle9->nneighbours() == 15);
  REQUIRE(particle10->nneighbours() == 15);
  REQUIRE(particle11->nneighbours() == 11);
  REQUIRE(particle12->nneighbours() == 8);
  REQUIRE(particle13->nneighbours() == 11);
  REQUIRE(particle14->nneighbours() == 11);
  REQUIRE(particle15->nneighbours() == 8);

  // Initialise particle variables
  // Assign particle volume
  mesh->iterate_over_particles(std::bind(&mpm::ParticleBase<Dim>::assign_volume,
                                         std::placeholders::_1, 1.));

  // Compute mass
  mesh->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Dim>::compute_mass, std::placeholders::_1));

  // Initialise nodes
  mesh->iterate_over_nodes(
      std::bind(&mpm::NodeBase<Dim>::initialise, std::placeholders::_1));

  mesh->iterate_over_cells(
      std::bind(&mpm::Cell<Dim>::activate_nodes, std::placeholders::_1));

  // Iterate over each particle to compute shapefn
  mesh->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Dim>::compute_shapefn, std::placeholders::_1));

  // Assign mass and momentum to nodes
  mesh->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Dim>::map_mass_momentum_to_nodes,
                std::placeholders::_1));

  SECTION("Mesh check initial condition") {

    // Check cell free surface
    REQUIRE(cell0->free_surface() == false);
    REQUIRE(cell1->free_surface() == false);
    REQUIRE(cell2->free_surface() == false);
    REQUIRE(cell3->free_surface() == false);
    REQUIRE(cell4->free_surface() == false);
    REQUIRE(cell5->free_surface() == false);
    REQUIRE(cell6->free_surface() == false);
    REQUIRE(cell7->free_surface() == false);
    REQUIRE(cell8->free_surface() == false);
    REQUIRE(cell9->free_surface() == false);
    REQUIRE(cell10->free_surface() == false);
    REQUIRE(cell11->free_surface() == false);
    REQUIRE(cell12->free_surface() == false);
    REQUIRE(cell13->free_surface() == false);
    REQUIRE(cell14->free_surface() == false);
    REQUIRE(cell15->free_surface() == false);
    REQUIRE(cell16->free_surface() == false);
    REQUIRE(cell17->free_surface() == false);
    REQUIRE(cell18->free_surface() == false);
    REQUIRE(cell19->free_surface() == false);
    REQUIRE(cell20->free_surface() == false);
    REQUIRE(cell21->free_surface() == false);
    REQUIRE(cell22->free_surface() == false);
    REQUIRE(cell23->free_surface() == false);
    REQUIRE(cell24->free_surface() == false);

    // Check cell status
    REQUIRE(cell0->status() == false);
    REQUIRE(cell1->status() == false);
    REQUIRE(cell2->status() == false);
    REQUIRE(cell3->status() == false);
    REQUIRE(cell4->status() == false);
    REQUIRE(cell5->status() == false);
    REQUIRE(cell6->status() == true);
    REQUIRE(cell7->status() == true);
    REQUIRE(cell8->status() == true);
    REQUIRE(cell9->status() == false);
    REQUIRE(cell10->status() == false);
    REQUIRE(cell11->status() == true);
    REQUIRE(cell12->status() == true);
    REQUIRE(cell13->status() == true);
    REQUIRE(cell14->status() == false);
    REQUIRE(cell15->status() == false);
    REQUIRE(cell16->status() == true);
    REQUIRE(cell17->status() == true);
    REQUIRE(cell18->status() == true);
    REQUIRE(cell19->status() == false);
    REQUIRE(cell20->status() == false);
    REQUIRE(cell21->status() == false);
    REQUIRE(cell22->status() == false);
    REQUIRE(cell23->status() == false);
    REQUIRE(cell24->status() == false);

    // Check node free surface
    REQUIRE(node0->free_surface() == false);
    REQUIRE(node1->free_surface() == false);
    REQUIRE(node2->free_surface() == false);
    REQUIRE(node3->free_surface() == false);
    REQUIRE(node4->free_surface() == false);
    REQUIRE(node5->free_surface() == false);
    REQUIRE(node6->free_surface() == false);
    REQUIRE(node7->free_surface() == false);
    REQUIRE(node8->free_surface() == false);
    REQUIRE(node9->free_surface() == false);
    REQUIRE(node10->free_surface() == false);
    REQUIRE(node11->free_surface() == false);
    REQUIRE(node12->free_surface() == false);
    REQUIRE(node13->free_surface() == false);
    REQUIRE(node14->free_surface() == false);
    REQUIRE(node15->free_surface() == false);
    REQUIRE(node16->free_surface() == false);
    REQUIRE(node17->free_surface() == false);
    REQUIRE(node18->free_surface() == false);
    REQUIRE(node19->free_surface() == false);
    REQUIRE(node20->free_surface() == false);
    REQUIRE(node21->free_surface() == false);
    REQUIRE(node22->free_surface() == false);
    REQUIRE(node23->free_surface() == false);
    REQUIRE(node24->free_surface() == false);
    REQUIRE(node25->free_surface() == false);
    REQUIRE(node26->free_surface() == false);
    REQUIRE(node27->free_surface() == false);
    REQUIRE(node28->free_surface() == false);
    REQUIRE(node29->free_surface() == false);
    REQUIRE(node30->free_surface() == false);
    REQUIRE(node31->free_surface() == false);
    REQUIRE(node32->free_surface() == false);
    REQUIRE(node33->free_surface() == false);
    REQUIRE(node34->free_surface() == false);
    REQUIRE(node35->free_surface() == false);

    // Check node status
    REQUIRE(node0->status() == false);
    REQUIRE(node1->status() == false);
    REQUIRE(node2->status() == false);
    REQUIRE(node3->status() == false);
    REQUIRE(node4->status() == false);
    REQUIRE(node5->status() == false);
    REQUIRE(node6->status() == false);
    REQUIRE(node7->status() == true);
    REQUIRE(node8->status() == true);
    REQUIRE(node9->status() == true);
    REQUIRE(node10->status() == true);
    REQUIRE(node11->status() == false);
    REQUIRE(node12->status() == false);
    REQUIRE(node13->status() == true);
    REQUIRE(node14->status() == true);
    REQUIRE(node15->status() == true);
    REQUIRE(node16->status() == true);
    REQUIRE(node17->status() == false);
    REQUIRE(node18->status() == false);
    REQUIRE(node19->status() == true);
    REQUIRE(node20->status() == true);
    REQUIRE(node21->status() == true);
    REQUIRE(node22->status() == true);
    REQUIRE(node23->status() == false);
    REQUIRE(node24->status() == false);
    REQUIRE(node25->status() == true);
    REQUIRE(node26->status() == true);
    REQUIRE(node27->status() == true);
    REQUIRE(node28->status() == true);
    REQUIRE(node29->status() == false);
    REQUIRE(node30->status() == false);
    REQUIRE(node31->status() == false);
    REQUIRE(node32->status() == false);
    REQUIRE(node33->status() == false);
    REQUIRE(node34->status() == false);
    REQUIRE(node35->status() == false);

    // Check particle free surface
    REQUIRE(particle0->free_surface() == false);
    REQUIRE(particle1->free_surface() == false);
    REQUIRE(particle2->free_surface() == false);
    REQUIRE(particle3->free_surface() == false);
    REQUIRE(particle4->free_surface() == false);
    REQUIRE(particle5->free_surface() == false);
    REQUIRE(particle6->free_surface() == false);
    REQUIRE(particle7->free_surface() == false);
    REQUIRE(particle8->free_surface() == false);
    REQUIRE(particle9->free_surface() == false);
    REQUIRE(particle10->free_surface() == false);
    REQUIRE(particle11->free_surface() == false);
    REQUIRE(particle12->free_surface() == false);
    REQUIRE(particle13->free_surface() == false);
    REQUIRE(particle14->free_surface() == false);
    REQUIRE(particle15->free_surface() == false);
  }

  SECTION("Mesh free surface 2D by density") {

    REQUIRE(mesh->compute_free_surface_by_density(0.25) == true);

    // Check cell volume fraction
    REQUIRE(cell0->volume_fraction() == Approx(0.0).epsilon(Tolerance));
    REQUIRE(cell1->volume_fraction() == Approx(0.0).epsilon(Tolerance));
    REQUIRE(cell2->volume_fraction() == Approx(0.0).epsilon(Tolerance));
    REQUIRE(cell3->volume_fraction() == Approx(0.0).epsilon(Tolerance));
    REQUIRE(cell4->volume_fraction() == Approx(0.0).epsilon(Tolerance));
    REQUIRE(cell5->volume_fraction() == Approx(0.0).epsilon(Tolerance));
    REQUIRE(cell6->volume_fraction() == Approx(0.25).epsilon(Tolerance));
    REQUIRE(cell7->volume_fraction() == Approx(0.5).epsilon(Tolerance));
    REQUIRE(cell8->volume_fraction() == Approx(0.25).epsilon(Tolerance));
    REQUIRE(cell9->volume_fraction() == Approx(0.0).epsilon(Tolerance));
    REQUIRE(cell10->volume_fraction() == Approx(0.0).epsilon(Tolerance));
    REQUIRE(cell11->volume_fraction() == Approx(0.5).epsilon(Tolerance));
    REQUIRE(cell12->volume_fraction() == Approx(1.0).epsilon(Tolerance));
    REQUIRE(cell13->volume_fraction() == Approx(0.5).epsilon(Tolerance));
    REQUIRE(cell14->volume_fraction() == Approx(0.0).epsilon(Tolerance));
    REQUIRE(cell15->volume_fraction() == Approx(0.0).epsilon(Tolerance));
    REQUIRE(cell16->volume_fraction() == Approx(0.25).epsilon(Tolerance));
    REQUIRE(cell17->volume_fraction() == Approx(0.5).epsilon(Tolerance));
    REQUIRE(cell18->volume_fraction() == Approx(0.25).epsilon(Tolerance));
    REQUIRE(cell19->volume_fraction() == Approx(0.0).epsilon(Tolerance));
    REQUIRE(cell20->volume_fraction() == Approx(0.0).epsilon(Tolerance));
    REQUIRE(cell21->volume_fraction() == Approx(0.0).epsilon(Tolerance));
    REQUIRE(cell22->volume_fraction() == Approx(0.0).epsilon(Tolerance));
    REQUIRE(cell23->volume_fraction() == Approx(0.0).epsilon(Tolerance));
    REQUIRE(cell24->volume_fraction() == Approx(0.0).epsilon(Tolerance));

    // Check nodal density for active nodes
    REQUIRE(node7->density(0) == Approx(62.5).epsilon(Tolerance));
    REQUIRE(node8->density(0) == Approx(218.75).epsilon(Tolerance));
    REQUIRE(node9->density(0) == Approx(218.75).epsilon(Tolerance));
    REQUIRE(node10->density(0) == Approx(62.5).epsilon(Tolerance));

    REQUIRE(node13->density(0) == Approx(218.75).epsilon(Tolerance));
    REQUIRE(node14->density(0) == Approx(765.625).epsilon(Tolerance));
    REQUIRE(node15->density(0) == Approx(765.625).epsilon(Tolerance));
    REQUIRE(node16->density(0) == Approx(218.75).epsilon(Tolerance));

    REQUIRE(node19->density(0) == Approx(218.75).epsilon(Tolerance));
    REQUIRE(node20->density(0) == Approx(765.625).epsilon(Tolerance));
    REQUIRE(node21->density(0) == Approx(765.625).epsilon(Tolerance));
    REQUIRE(node22->density(0) == Approx(218.75).epsilon(Tolerance));

    REQUIRE(node25->density(0) == Approx(62.5).epsilon(Tolerance));
    REQUIRE(node26->density(0) == Approx(218.75).epsilon(Tolerance));
    REQUIRE(node27->density(0) == Approx(218.75).epsilon(Tolerance));
    REQUIRE(node28->density(0) == Approx(62.5).epsilon(Tolerance));

    // Check solutions
    std::set<mpm::Index> fsc = {6, 7, 8, 11, 13, 16, 17, 18};
    std::set<mpm::Index> fsn = {7, 8, 9, 10, 13, 16, 19, 22, 25, 26, 27, 28};
    std::set<mpm::Index> fsp = {0, 1, 2, 3, 4, 7, 8, 11, 12, 13, 14, 15};

    // Check cell, node, and particle free surface
    REQUIRE(mesh->free_surface_cells() == fsc);
    REQUIRE(mesh->free_surface_nodes() == fsn);
    REQUIRE(mesh->free_surface_particles() == fsp);
  }

  SECTION("Mesh free surface 2D by geometry") {

    REQUIRE(mesh->compute_free_surface_by_geometry(0.25) == true);

    // Check solutions
    std::set<mpm::Index> fsc = {6, 7, 8, 11, 13, 16, 17, 18};
    std::set<mpm::Index> fsn = {7, 8, 9, 10, 13, 16, 19, 22, 25, 26, 27, 28};
    std::set<mpm::Index> fsp = {0, 1, 2, 3, 4, 7, 8, 11, 12, 13, 14, 15};

    // Check cell, node, and particle free surface
    REQUIRE(mesh->free_surface_cells() == fsc);
    REQUIRE(mesh->free_surface_nodes() == fsn);
    REQUIRE(mesh->free_surface_particles() == fsp);

    // Check particle normal vector
    REQUIRE(particle0->normal().norm() == Approx(1.0).epsilon(Tolerance));
    REQUIRE(particle1->normal().norm() == Approx(1.0).epsilon(Tolerance));
    REQUIRE(particle2->normal().norm() == Approx(1.0).epsilon(Tolerance));
    REQUIRE(particle3->normal().norm() == Approx(1.0).epsilon(Tolerance));
    REQUIRE(particle4->normal().norm() == Approx(1.0).epsilon(Tolerance));
    REQUIRE(particle5->normal().norm() == Approx(0.0).epsilon(Tolerance));
    REQUIRE(particle6->normal().norm() == Approx(0.0).epsilon(Tolerance));
    REQUIRE(particle7->normal().norm() == Approx(1.0).epsilon(Tolerance));
    REQUIRE(particle8->normal().norm() == Approx(1.0).epsilon(Tolerance));
    REQUIRE(particle9->normal().norm() == Approx(0.0).epsilon(Tolerance));
    REQUIRE(particle10->normal().norm() == Approx(0.0).epsilon(Tolerance));
    REQUIRE(particle11->normal().norm() == Approx(1.0).epsilon(Tolerance));
    REQUIRE(particle12->normal().norm() == Approx(1.0).epsilon(Tolerance));
    REQUIRE(particle13->normal().norm() == Approx(1.0).epsilon(Tolerance));
    REQUIRE(particle14->normal().norm() == Approx(1.0).epsilon(Tolerance));
    REQUIRE(particle15->normal().norm() == Approx(1.0).epsilon(Tolerance));
  }

  SECTION("Mesh free surface 2D") {

    // Check solutions
    std::set<mpm::Index> fsc = {6, 7, 8, 11, 13, 16, 17, 18};
    std::set<mpm::Index> fsn = {7, 8, 9, 10, 13, 16, 19, 22, 25, 26, 27, 28};
    std::set<mpm::Index> fsp = {0, 1, 2, 3, 4, 7, 8, 11, 12, 13, 14, 15};

    std::string method = "density";
    REQUIRE(mesh->compute_free_surface(method, 0.25) == true);

    // Check cell, node, and particle free surface
    REQUIRE(mesh->free_surface_cells() == fsc);
    REQUIRE(mesh->free_surface_nodes() == fsn);
    REQUIRE(mesh->free_surface_particles() == fsp);

    method = "geometry";
    REQUIRE(mesh->compute_free_surface(method, 0.25) == true);

    // Check cell, node, and particle free surface
    REQUIRE(mesh->free_surface_cells() == fsc);
    REQUIRE(mesh->free_surface_nodes() == fsn);
    REQUIRE(mesh->free_surface_particles() == fsp);

    method = "other";
    REQUIRE(mesh->compute_free_surface(method, 0.25) == true);
  }
}


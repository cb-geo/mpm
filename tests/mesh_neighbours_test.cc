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

TEST_CASE("Mesh cell neighbours 2D", "[MeshCell][2D][mpi]") {
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

    SECTION("Check add cell neighbours") {
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

    SECTION("Find cell neighbours") {
      mesh->find_cell_neighbours();
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

      SECTION("Find particle neighbours") {
        // Get number of MPI ranks
        int mpi_size, mpi_rank;
#ifdef USE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

        if (mpi_size == 4) {
          cell0->rank(0);
          cell1->rank(0);
          cell4->rank(0);
          cell2->rank(1);
          cell3->rank(1);
          cell5->rank(1);
          cell8->rank(2);
          cell7->rank(2);
          cell6->rank(3);

          std::shared_ptr<mpm::ParticleBase<Dim>> particle0, particle1,
              particle2, particle3, particle4, particle5, particle6, particle7,
              particle8, particle9, particle10, particle11, particle12,
              particle13, particle14, particle15, particle16, particle17;

          if (mpi_rank == 0) {
            coords << 1.0, 1.0;
            particle0 = std::make_shared<mpm::Particle<Dim>>(0, coords);
            REQUIRE(mesh->add_particle(particle0, false) == true);

            coords << 1.5, 1.5;
            particle1 = std::make_shared<mpm::Particle<Dim>>(1, coords);
            REQUIRE(mesh->add_particle(particle1, false) == true);

            coords << 2.1, 1.0;
            particle2 = std::make_shared<mpm::Particle<Dim>>(2, coords);
            REQUIRE(mesh->add_particle(particle2, false) == true);

            coords << 3.0, 1.5;
            particle3 = std::make_shared<mpm::Particle<Dim>>(3, coords);
            REQUIRE(mesh->add_particle(particle3, false) == true);

            coords << 4.1, 1.0;
            particle4 = std::make_shared<mpm::Particle<Dim>>(4, coords);
            REQUIRE(mesh->add_particle(particle4, false) == true);

            coords << 4.5, 1.5;
            particle5 = std::make_shared<mpm::Particle<Dim>>(5, coords);
            REQUIRE(mesh->add_particle(particle5, false) == true);
          }
          if (mpi_rank == 1) {
            coords << 1.0, 2.5;
            particle6 = std::make_shared<mpm::Particle<Dim>>(6, coords);
            REQUIRE(mesh->add_particle(particle6, false) == true);

            coords << 1.5, 2.5;
            particle7 = std::make_shared<mpm::Particle<Dim>>(7, coords);
            REQUIRE(mesh->add_particle(particle7, false) == true);

            coords << 2.1, 2.5;
            particle8 = std::make_shared<mpm::Particle<Dim>>(8, coords);
            REQUIRE(mesh->add_particle(particle8, false) == true);

            coords << 3.0, 2.5;
            particle9 = std::make_shared<mpm::Particle<Dim>>(9, coords);
            REQUIRE(mesh->add_particle(particle9, false) == true);

            coords << 4.1, 2.5;
            particle10 = std::make_shared<mpm::Particle<Dim>>(10, coords);
            REQUIRE(mesh->add_particle(particle10, false) == true);

            coords << 4.5, 2.5;
            particle11 = std::make_shared<mpm::Particle<Dim>>(11, coords);
            REQUIRE(mesh->add_particle(particle11, false) == true);
          }
          if (mpi_rank == 2) {
            coords << 1.0, 4.5;
            particle12 = std::make_shared<mpm::Particle<Dim>>(12, coords);
            REQUIRE(mesh->add_particle(particle12, false) == true);

            coords << 1.5, 4.5;
            particle13 = std::make_shared<mpm::Particle<Dim>>(13, coords);
            REQUIRE(mesh->add_particle(particle13, false) == true);

            coords << 2.1, 4.5;
            particle14 = std::make_shared<mpm::Particle<Dim>>(14, coords);
            REQUIRE(mesh->add_particle(particle14, false) == true);

            coords << 3.0, 4.5;
            particle15 = std::make_shared<mpm::Particle<Dim>>(15, coords);
            REQUIRE(mesh->add_particle(particle15, false) == true);
          }
          if (mpi_rank == 3) {
            coords << 4.1, 4.5;
            particle16 = std::make_shared<mpm::Particle<Dim>>(16, coords);
            REQUIRE(mesh->add_particle(particle16, false) == true);

            coords << 4.5, 4.5;
            particle17 = std::make_shared<mpm::Particle<Dim>>(17, coords);
            REQUIRE(mesh->add_particle(particle17, false) == true);
          }

          // Locate particles in a mesh
          auto particles = mesh->locate_particles_mesh();

          // Should find all particles in mesh
          REQUIRE(particles.size() == 0);

          // Check particles inside cells
          if (mpi_rank == 0) {
            REQUIRE(cell0->nparticles() == 2);
            REQUIRE(cell1->nparticles() == 2);
            REQUIRE(cell4->nparticles() == 2);
          }
          if (mpi_rank == 1) {
            REQUIRE(cell2->nparticles() == 2);
            REQUIRE(cell3->nparticles() == 2);
            REQUIRE(cell5->nparticles() == 2);
          }
          if (mpi_rank == 2) {
            REQUIRE(cell7->nparticles() == 2);
            REQUIRE(cell8->nparticles() == 2);
          }
          if (mpi_rank == 3) {
            REQUIRE(cell6->nparticles() == 2);
          }

          // Check solutions
          std::vector<mpm::Index> p0 = {0, 1};
          std::vector<mpm::Index> p1 = {2, 3};
          std::vector<mpm::Index> p2 = {8, 9};
          std::vector<mpm::Index> p3 = {6, 7};
          std::vector<mpm::Index> p4 = {4, 5};
          std::vector<mpm::Index> p5 = {10, 11};
          std::vector<mpm::Index> p6 = {16, 17};
          std::vector<mpm::Index> p7 = {14, 15};
          std::vector<mpm::Index> p8 = {12, 13};

          // Check whether particle id inside cells are as intended
          if (mpi_rank == 0) {
            REQUIRE(cell0->particles() == p0);
            REQUIRE(cell1->particles() == p1);
            REQUIRE(cell4->particles() == p4);
          }
          if (mpi_rank == 1) {
            REQUIRE(cell2->particles() == p2);
            REQUIRE(cell3->particles() == p3);
            REQUIRE(cell5->particles() == p5);
          }
          if (mpi_rank == 2) {
            REQUIRE(cell7->particles() == p7);
            REQUIRE(cell8->particles() == p8);
          }
          if (mpi_rank == 3) {
            REQUIRE(cell6->particles() == p6);
          }

          // Find particle neighbours
          mesh->find_particle_neighbours();

          // Check particle neighbours
          if (mpi_rank == 0) {
            REQUIRE(particle0->nneighbours() == 7);
            REQUIRE(particle1->nneighbours() == 7);
            REQUIRE(particle2->nneighbours() == 11);
            REQUIRE(particle3->nneighbours() == 11);
            REQUIRE(particle4->nneighbours() == 7);
            REQUIRE(particle5->nneighbours() == 7);
          }

          if (mpi_rank == 1) {
            REQUIRE(particle6->nneighbours() == 11);
            REQUIRE(particle7->nneighbours() == 11);
            REQUIRE(particle8->nneighbours() == 17);
            REQUIRE(particle9->nneighbours() == 17);
            REQUIRE(particle10->nneighbours() == 11);
            REQUIRE(particle11->nneighbours() == 11);
          }

          if (mpi_rank == 2) {
            REQUIRE(particle12->nneighbours() == 7);

            REQUIRE(particle13->nneighbours() == 7);
            REQUIRE(particle14->nneighbours() == 11);
            REQUIRE(particle15->nneighbours() == 11);
          }

          if (mpi_rank == 3) {
            REQUIRE(particle16->nneighbours() == 7);
            REQUIRE(particle17->nneighbours() == 7);
          }

          // Check solutions
          std::vector<mpm::Index> np0 = {1, 2, 3, 8, 9, 6, 7};
          std::vector<mpm::Index> np1 = {0, 2, 3, 8, 9, 6, 7};
          std::vector<mpm::Index> np2 = {3, 0, 1, 8, 9, 6, 7, 4, 5, 10, 11};
          std::vector<mpm::Index> np3 = {2, 0, 1, 8, 9, 6, 7, 4, 5, 10, 11};
          std::vector<mpm::Index> np4 = {5, 2, 3, 8, 9, 10, 11};
          std::vector<mpm::Index> np5 = {4, 2, 3, 8, 9, 10, 11};
          std::vector<mpm::Index> np6 = {7, 0, 1, 2, 3, 8, 9, 14, 15, 12, 13};
          std::vector<mpm::Index> np7 = {6, 0, 1, 2, 3, 8, 9, 14, 15, 12, 13};
          std::vector<mpm::Index> np8 = {9,  0,  1,  2,  3,  6,  7,  4, 5,
                                         10, 11, 16, 17, 14, 15, 12, 13};
          std::vector<mpm::Index> np9 = {8,  0,  1,  2,  3,  6,  7,  4, 5,
                                         10, 11, 16, 17, 14, 15, 12, 13};
          std::vector<mpm::Index> np10 = {11, 2, 3, 8, 9, 4, 5, 16, 17, 14, 15};
          std::vector<mpm::Index> np11 = {10, 2, 3, 8, 9, 4, 5, 16, 17, 14, 15};
          std::vector<mpm::Index> np12 = {13, 8, 9, 6, 7, 14, 15};
          std::vector<mpm::Index> np13 = {12, 8, 9, 6, 7, 14, 15};
          std::vector<mpm::Index> np14 = {15, 8,  9,  6,  7, 10,
                                          11, 16, 17, 12, 13};
          std::vector<mpm::Index> np15 = {14, 8,  9,  6,  7, 10,
                                          11, 16, 17, 12, 13};
          std::vector<mpm::Index> np16 = {17, 8, 9, 10, 11, 14, 15};
          std::vector<mpm::Index> np17 = {16, 8, 9, 10, 11, 14, 15};

          // Check whether particle id inside cells are as intended
          if (mpi_rank == 0) {
            REQUIRE(particle0->neighbours() == np0);
            REQUIRE(particle1->neighbours() == np1);
            REQUIRE(particle2->neighbours() == np2);
            REQUIRE(particle3->neighbours() == np3);
            REQUIRE(particle4->neighbours() == np4);
            REQUIRE(particle5->neighbours() == np5);
          }
          if (mpi_rank == 1) {
            REQUIRE(particle6->neighbours() == np6);
            REQUIRE(particle7->neighbours() == np7);
            REQUIRE(particle8->neighbours() == np8);
            REQUIRE(particle9->neighbours() == np9);
            REQUIRE(particle10->neighbours() == np10);
            REQUIRE(particle11->neighbours() == np11);
          }
          if (mpi_rank == 2) {
            REQUIRE(particle12->neighbours() == np12);
            REQUIRE(particle13->neighbours() == np13);
            REQUIRE(particle14->neighbours() == np14);
            REQUIRE(particle15->neighbours() == np15);
          }
          if (mpi_rank == 3) {
            REQUIRE(particle16->neighbours() == np16);
            REQUIRE(particle17->neighbours() == np17);
          }
#endif
        }
      }
    }
  }
}

TEST_CASE("Mesh cell neighbours 3D", "[MeshCell][3D][mpi]") {
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

    // Create cell1
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

    // Cell 2
    coords << 4, 4, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node12 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(12, coords);
    REQUIRE(mesh->add_node(node12) == true);

    coords << 4, 4, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node13 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(13, coords);
    REQUIRE(mesh->add_node(node13) == true);

    coords << 2, 4, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node14 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(14, coords);
    REQUIRE(mesh->add_node(node14) == true);

    coords << 2, 4, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node15 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(15, coords);
    REQUIRE(mesh->add_node(node15) == true);

    // Create cell2
    auto cell2 = std::make_shared<mpm::Cell<Dim>>(2, Nnodes, element);

    // Add nodes to cell
    cell2->add_node(0, node2);
    cell2->add_node(1, node9);
    cell2->add_node(2, node12);
    cell2->add_node(3, node14);
    cell2->add_node(4, node6);
    cell2->add_node(5, node11);
    cell2->add_node(6, node13);
    cell2->add_node(7, node15);

    REQUIRE(cell2->nnodes() == 8);

    REQUIRE(mesh->add_cell(cell2) == true);

    // Cell 3
    coords << 0, 4, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node16 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(16, coords);
    REQUIRE(mesh->add_node(node16) == true);

    coords << 0, 4, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node17 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(17, coords);
    REQUIRE(mesh->add_node(node17) == true);

    // Create cell3
    auto cell3 = std::make_shared<mpm::Cell<Dim>>(3, Nnodes, element);

    // Add nodes to cell
    cell3->add_node(0, node3);
    cell3->add_node(1, node2);
    cell3->add_node(2, node14);
    cell3->add_node(3, node16);
    cell3->add_node(4, node7);
    cell3->add_node(5, node6);
    cell3->add_node(6, node15);
    cell3->add_node(7, node17);

    REQUIRE(cell3->nnodes() == 8);

    REQUIRE(mesh->add_cell(cell3) == true);

    // Cell 4
    coords << 6, 0, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node18 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(18, coords);
    REQUIRE(mesh->add_node(node18) == true);

    coords << 6, 2, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node19 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(19, coords);
    REQUIRE(mesh->add_node(node19) == true);

    coords << 6, 0, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node20 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(20, coords);
    REQUIRE(mesh->add_node(node20) == true);

    coords << 6, 2, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node21 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(21, coords);
    REQUIRE(mesh->add_node(node21) == true);

    // Create cell4
    auto cell4 = std::make_shared<mpm::Cell<Dim>>(4, Nnodes, element);

    // Add nodes to cell
    cell4->add_node(0, node8);
    cell4->add_node(1, node18);
    cell4->add_node(2, node19);
    cell4->add_node(3, node9);
    cell4->add_node(4, node10);
    cell4->add_node(5, node20);
    cell4->add_node(6, node21);
    cell4->add_node(7, node11);

    REQUIRE(cell4->nnodes() == 8);

    REQUIRE(mesh->add_cell(cell4) == true);

    // Cell 5
    coords << 6, 4, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node22 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(22, coords);
    REQUIRE(mesh->add_node(node22) == true);

    coords << 6, 4, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node23 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(23, coords);
    REQUIRE(mesh->add_node(node23) == true);

    // Create cell5
    auto cell5 = std::make_shared<mpm::Cell<Dim>>(5, Nnodes, element);

    // Add nodes to cell
    cell5->add_node(0, node9);
    cell5->add_node(1, node19);
    cell5->add_node(2, node22);
    cell5->add_node(3, node12);
    cell5->add_node(4, node11);
    cell5->add_node(5, node21);
    cell5->add_node(6, node23);
    cell5->add_node(7, node13);

    REQUIRE(cell5->nnodes() == 8);

    REQUIRE(mesh->add_cell(cell5) == true);

    // Cell 6
    coords << 6, 6, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node24 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(24, coords);
    REQUIRE(mesh->add_node(node24) == true);

    coords << 6, 6, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node25 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(25, coords);
    REQUIRE(mesh->add_node(node25) == true);

    coords << 4, 6, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node26 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(26, coords);
    REQUIRE(mesh->add_node(node26) == true);

    coords << 4, 6, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node27 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(27, coords);
    REQUIRE(mesh->add_node(node27) == true);

    // Create cell6
    auto cell6 = std::make_shared<mpm::Cell<Dim>>(6, Nnodes, element);

    // Add nodes to cell
    cell6->add_node(0, node12);
    cell6->add_node(1, node22);
    cell6->add_node(2, node24);
    cell6->add_node(3, node26);
    cell6->add_node(4, node13);
    cell6->add_node(5, node23);
    cell6->add_node(6, node25);
    cell6->add_node(7, node27);

    REQUIRE(cell6->nnodes() == 8);

    REQUIRE(mesh->add_cell(cell6) == true);

    // Cell 7
    coords << 2, 6, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node28 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(28, coords);
    REQUIRE(mesh->add_node(node28) == true);

    coords << 2, 6, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node29 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(29, coords);
    REQUIRE(mesh->add_node(node29) == true);

    // Create cell7
    auto cell7 = std::make_shared<mpm::Cell<Dim>>(7, Nnodes, element);

    // Add nodes to cell
    cell7->add_node(0, node14);
    cell7->add_node(1, node12);
    cell7->add_node(2, node26);
    cell7->add_node(3, node28);
    cell7->add_node(4, node15);
    cell7->add_node(5, node13);
    cell7->add_node(6, node27);
    cell7->add_node(7, node29);

    REQUIRE(cell7->nnodes() == 8);

    REQUIRE(mesh->add_cell(cell7) == true);

    // Cell 8
    coords << 0, 6, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node30 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(30, coords);
    REQUIRE(mesh->add_node(node30) == true);

    coords << 0, 6, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node31 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(31, coords);
    REQUIRE(mesh->add_node(node31) == true);

    // Create cell8
    auto cell8 = std::make_shared<mpm::Cell<Dim>>(8, Nnodes, element);

    // Add nodes to cell
    cell8->add_node(0, node16);
    cell8->add_node(1, node14);
    cell8->add_node(2, node28);
    cell8->add_node(3, node30);
    cell8->add_node(4, node17);
    cell8->add_node(5, node15);
    cell8->add_node(6, node29);
    cell8->add_node(7, node31);

    REQUIRE(cell8->nnodes() == 8);

    REQUIRE(mesh->add_cell(cell8) == true);

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

    SECTION("Check assign cell neighbours") {
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
    SECTION("Find cell neighbours") {
      mesh->find_cell_neighbours();
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

      REQUIRE(cell0->initialise() == true);
      REQUIRE(cell1->initialise() == true);
      REQUIRE(cell2->initialise() == true);
      REQUIRE(cell3->initialise() == true);
      REQUIRE(cell4->initialise() == true);
      REQUIRE(cell5->initialise() == true);
      REQUIRE(cell6->initialise() == true);
      REQUIRE(cell7->initialise() == true);
      REQUIRE(cell8->initialise() == true);

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

      SECTION("Find particle neighbours") {
        // Get number of MPI ranks
        int mpi_size, mpi_rank;
#ifdef USE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

        if (mpi_size == 4) {
          // Assign cell ranks
          cell0->rank(0);
          cell1->rank(0);
          cell4->rank(0);
          cell2->rank(1);
          cell3->rank(1);
          cell5->rank(1);
          cell7->rank(2);
          cell8->rank(2);
          cell6->rank(3);

          std::shared_ptr<mpm::ParticleBase<Dim>> particle0, particle1,
              particle2, particle3, particle4, particle5, particle6, particle7,
              particle8, particle9, particle10, particle11, particle12,
              particle13, particle14, particle15, particle16, particle17;

          if (mpi_rank == 0) {
            coords << 1.0, 1.0, 1.5;
            particle0 = std::make_shared<mpm::Particle<Dim>>(0, coords);
            REQUIRE(mesh->add_particle(particle0, false) == true);

            coords << 1.5, 1.5, 1.5;
            particle1 = std::make_shared<mpm::Particle<Dim>>(1, coords);
            REQUIRE(mesh->add_particle(particle1, false) == true);

            coords << 2.1, 1.0, 1.5;
            particle2 = std::make_shared<mpm::Particle<Dim>>(2, coords);
            REQUIRE(mesh->add_particle(particle2, false) == true);

            coords << 3.0, 1.5, 1.5;
            particle3 = std::make_shared<mpm::Particle<Dim>>(3, coords);
            REQUIRE(mesh->add_particle(particle3, false) == true);

            coords << 4.1, 1.0, 1.5;
            particle4 = std::make_shared<mpm::Particle<Dim>>(4, coords);
            REQUIRE(mesh->add_particle(particle4, false) == true);

            coords << 4.5, 1.5, 1.5;
            particle5 = std::make_shared<mpm::Particle<Dim>>(5, coords);
            REQUIRE(mesh->add_particle(particle5, false) == true);
          }
          if (mpi_rank == 1) {
            coords << 1.0, 2.5, 1.5;
            particle6 = std::make_shared<mpm::Particle<Dim>>(6, coords);
            REQUIRE(mesh->add_particle(particle6, false) == true);

            coords << 1.5, 2.5, 1.5;
            particle7 = std::make_shared<mpm::Particle<Dim>>(7, coords);
            REQUIRE(mesh->add_particle(particle7, false) == true);

            coords << 2.1, 2.5, 1.5;
            particle8 = std::make_shared<mpm::Particle<Dim>>(8, coords);
            REQUIRE(mesh->add_particle(particle8, false) == true);

            coords << 3.0, 2.5, 1.5;
            particle9 = std::make_shared<mpm::Particle<Dim>>(9, coords);
            REQUIRE(mesh->add_particle(particle9, false) == true);

            coords << 4.1, 2.5, 1.5;
            particle10 = std::make_shared<mpm::Particle<Dim>>(10, coords);
            REQUIRE(mesh->add_particle(particle10, false) == true);

            coords << 4.5, 2.5, 1.5;
            particle11 = std::make_shared<mpm::Particle<Dim>>(11, coords);
            REQUIRE(mesh->add_particle(particle11, false) == true);
          }
          if (mpi_rank == 2) {
            coords << 1.0, 4.5, 1.5;
            particle12 = std::make_shared<mpm::Particle<Dim>>(12, coords);
            REQUIRE(mesh->add_particle(particle12, false) == true);

            coords << 1.5, 4.5, 1.5;
            particle13 = std::make_shared<mpm::Particle<Dim>>(13, coords);
            REQUIRE(mesh->add_particle(particle13, false) == true);

            coords << 2.1, 4.5, 1.5;
            particle14 = std::make_shared<mpm::Particle<Dim>>(14, coords);
            REQUIRE(mesh->add_particle(particle14, false) == true);

            coords << 3.0, 4.5, 1.5;
            particle15 = std::make_shared<mpm::Particle<Dim>>(15, coords);
            REQUIRE(mesh->add_particle(particle15, false) == true);
          }
          if (mpi_rank == 3) {
            coords << 4.1, 4.5, 1.5;
            particle16 = std::make_shared<mpm::Particle<Dim>>(16, coords);
            REQUIRE(mesh->add_particle(particle16, false) == true);

            coords << 4.5, 4.5, 1.5;
            particle17 = std::make_shared<mpm::Particle<Dim>>(17, coords);
            REQUIRE(mesh->add_particle(particle17, false) == true);
          }

          // Locate particles in a mesh
          auto particles = mesh->locate_particles_mesh();

          // Should find all particles in mesh
          REQUIRE(particles.size() == 0);

          // Check particles inside cells
          if (mpi_rank == 0) {
            REQUIRE(cell0->nparticles() == 2);
            REQUIRE(cell1->nparticles() == 2);
            REQUIRE(cell4->nparticles() == 2);
          }
          if (mpi_rank == 1) {
            REQUIRE(cell2->nparticles() == 2);
            REQUIRE(cell3->nparticles() == 2);
            REQUIRE(cell5->nparticles() == 2);
          }
          if (mpi_rank == 2) {
            REQUIRE(cell7->nparticles() == 2);
            REQUIRE(cell8->nparticles() == 2);
          }
          if (mpi_rank == 3) {
            REQUIRE(cell6->nparticles() == 2);
          }

          // Check solutions
          std::vector<mpm::Index> p0 = {0, 1};
          std::vector<mpm::Index> p1 = {2, 3};
          std::vector<mpm::Index> p2 = {8, 9};
          std::vector<mpm::Index> p3 = {6, 7};
          std::vector<mpm::Index> p4 = {4, 5};
          std::vector<mpm::Index> p5 = {10, 11};
          std::vector<mpm::Index> p6 = {16, 17};
          std::vector<mpm::Index> p7 = {14, 15};
          std::vector<mpm::Index> p8 = {12, 13};

          // Check whether particle id inside cells are as intended
          if (mpi_rank == 0) {
            REQUIRE(cell0->particles() == p0);
            REQUIRE(cell1->particles() == p1);
            REQUIRE(cell4->particles() == p4);
          }
          if (mpi_rank == 1) {
            REQUIRE(cell2->particles() == p2);
            REQUIRE(cell3->particles() == p3);
            REQUIRE(cell5->particles() == p5);
          }
          if (mpi_rank == 2) {
            REQUIRE(cell7->particles() == p7);
            REQUIRE(cell8->particles() == p8);
          }
          if (mpi_rank == 3) {
            REQUIRE(cell6->particles() == p6);
          }
          // Find particle neighbours
          mesh->find_particle_neighbours();

          // Check particle neighbours
          if (mpi_rank == 0) {
            REQUIRE(particle0->nneighbours() == 7);
            REQUIRE(particle1->nneighbours() == 7);
            REQUIRE(particle2->nneighbours() == 11);
            REQUIRE(particle3->nneighbours() == 11);
            REQUIRE(particle4->nneighbours() == 7);
            REQUIRE(particle5->nneighbours() == 7);
          }
          if (mpi_rank == 1) {
            REQUIRE(particle6->nneighbours() == 11);
            REQUIRE(particle7->nneighbours() == 11);
            REQUIRE(particle8->nneighbours() == 17);
            REQUIRE(particle9->nneighbours() == 17);
            REQUIRE(particle10->nneighbours() == 11);
            REQUIRE(particle11->nneighbours() == 11);
          }
          if (mpi_rank == 2) {
            REQUIRE(particle12->nneighbours() == 7);
            REQUIRE(particle13->nneighbours() == 7);
            REQUIRE(particle14->nneighbours() == 11);
            REQUIRE(particle15->nneighbours() == 11);
          }
          if (mpi_rank == 3) {
            REQUIRE(particle16->nneighbours() == 7);
            REQUIRE(particle17->nneighbours() == 7);
          }

          // Check solutions
          std::vector<mpm::Index> np0 = {1, 2, 3, 8, 9, 6, 7};
          std::vector<mpm::Index> np1 = {0, 2, 3, 8, 9, 6, 7};
          std::vector<mpm::Index> np2 = {3, 0, 1, 8, 9, 6, 7, 4, 5, 10, 11};
          std::vector<mpm::Index> np3 = {2, 0, 1, 8, 9, 6, 7, 4, 5, 10, 11};
          std::vector<mpm::Index> np4 = {5, 2, 3, 8, 9, 10, 11};
          std::vector<mpm::Index> np5 = {4, 2, 3, 8, 9, 10, 11};
          std::vector<mpm::Index> np6 = {7, 0, 1, 2, 3, 8, 9, 14, 15, 12, 13};
          std::vector<mpm::Index> np7 = {6, 0, 1, 2, 3, 8, 9, 14, 15, 12, 13};
          std::vector<mpm::Index> np8 = {9,  0,  1,  2,  3,  6,  7,  4, 5,
                                         10, 11, 16, 17, 14, 15, 12, 13};
          std::vector<mpm::Index> np9 = {8,  0,  1,  2,  3,  6,  7,  4, 5,
                                         10, 11, 16, 17, 14, 15, 12, 13};
          std::vector<mpm::Index> np10 = {11, 2, 3, 8, 9, 4, 5, 16, 17, 14, 15};
          std::vector<mpm::Index> np11 = {10, 2, 3, 8, 9, 4, 5, 16, 17, 14, 15};
          std::vector<mpm::Index> np12 = {13, 8, 9, 6, 7, 14, 15};
          std::vector<mpm::Index> np13 = {12, 8, 9, 6, 7, 14, 15};
          std::vector<mpm::Index> np14 = {15, 8,  9,  6,  7, 10,
                                          11, 16, 17, 12, 13};
          std::vector<mpm::Index> np15 = {14, 8,  9,  6,  7, 10,
                                          11, 16, 17, 12, 13};
          std::vector<mpm::Index> np16 = {17, 8, 9, 10, 11, 14, 15};
          std::vector<mpm::Index> np17 = {16, 8, 9, 10, 11, 14, 15};

          // Check whether particle id inside cells are as intended
          if (mpi_rank == 0) {
            REQUIRE(particle0->neighbours() == np0);
            REQUIRE(particle1->neighbours() == np1);
            REQUIRE(particle2->neighbours() == np2);
            REQUIRE(particle3->neighbours() == np3);
            REQUIRE(particle4->neighbours() == np4);
            REQUIRE(particle5->neighbours() == np5);
          }
          if (mpi_rank == 1) {
            REQUIRE(particle6->neighbours() == np6);
            REQUIRE(particle7->neighbours() == np7);
            REQUIRE(particle8->neighbours() == np8);
            REQUIRE(particle9->neighbours() == np9);
            REQUIRE(particle10->neighbours() == np10);
            REQUIRE(particle11->neighbours() == np11);
          }
          if (mpi_rank == 2) {
            REQUIRE(particle12->neighbours() == np12);
            REQUIRE(particle13->neighbours() == np13);
            REQUIRE(particle14->neighbours() == np14);
            REQUIRE(particle15->neighbours() == np15);
          }
          if (mpi_rank == 3) {
            REQUIRE(particle16->neighbours() == np16);
            REQUIRE(particle17->neighbours() == np17);
          }
        }
#endif
      }
    }
  }
}

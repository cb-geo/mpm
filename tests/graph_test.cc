#ifdef USE_PARMETIS

#include <functional>
#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

#include "cell.h"
#include "container.h"
#include "data_types.h"
#include "element.h"
#include "graph.h"
#include "hdf5_particle.h"
#include "hexahedron_element.h"
#include "material/material.h"
#include "mesh.h"
#include "mpi_datatypes.h"
#include "node.h"
#include "particle.h"
#include "quadrilateral_element.h"

// MPI
#ifdef USE_MPI
#include "mpi.h"
#endif

//! \brief Check graph class for 2D case
TEST_CASE("Graph is checked for 2D case", "[graph][2D]") {
  // Dimension
  const unsigned Dim = 2;
  // Deress of freedom
  const unsigned Dof = 2;
  // Number of nodes per cell
  const unsigned Nnodes = 4;
  // Number of phases
  const unsigned Nphases = 1;

  // Tolerance
  const double Tolerance = 1.E-7;

  //! Try to use the graph from the guide to test the graph for 2D

  // Element
  std::shared_ptr<mpm::Element<Dim>> element =
      std::make_shared<mpm::QuadrilateralElement<Dim, 4>>();

  // Cell 1
  mpm::Index id0 = 0;
  auto cell0 = std::make_shared<mpm::Cell<Dim>>(id0, Nnodes, element);

  // Cell 2
  mpm::Index id1 = 1;
  auto cell1 = std::make_shared<mpm::Cell<Dim>>(id1, Nnodes, element);

  // Cell 3
  mpm::Index id2 = 2;
  auto cell2 = std::make_shared<mpm::Cell<Dim>>(id2, Nnodes, element);

  // Cell 4
  mpm::Index id3 = 3;
  auto cell3 = std::make_shared<mpm::Cell<Dim>>(id3, Nnodes, element);

  // Cell 5
  mpm::Index id4 = 4;
  auto cell4 = std::make_shared<mpm::Cell<Dim>>(id4, Nnodes, element);

  // Cell 6
  mpm::Index id5 = 5;
  auto cell5 = std::make_shared<mpm::Cell<Dim>>(id5, Nnodes, element);

  // Cell 7
  mpm::Index id6 = 6;
  auto cell6 = std::make_shared<mpm::Cell<Dim>>(id6, Nnodes, element);

  // Cell 8
  mpm::Index id7 = 7;
  auto cell7 = std::make_shared<mpm::Cell<Dim>>(id7, Nnodes, element);

  // Cell 9
  mpm::Index id8 = 8;
  auto cell8 = std::make_shared<mpm::Cell<Dim>>(id8, Nnodes, element);

  // Cell 10
  mpm::Index id9 = 9;
  auto cell9 = std::make_shared<mpm::Cell<Dim>>(id9, Nnodes, element);

  // Cell 11
  mpm::Index id10 = 10;
  auto cell10 = std::make_shared<mpm::Cell<Dim>>(id10, Nnodes, element);

  // Cell 12
  mpm::Index id11 = 11;
  auto cell11 = std::make_shared<mpm::Cell<Dim>>(id11, Nnodes, element);

  // Cell 13
  mpm::Index id12 = 12;
  auto cell12 = std::make_shared<mpm::Cell<Dim>>(id12, Nnodes, element);

  // Cell 14
  mpm::Index id13 = 13;
  auto cell13 = std::make_shared<mpm::Cell<Dim>>(id13, Nnodes, element);

  // Cell 15
  mpm::Index id14 = 14;
  auto cell14 = std::make_shared<mpm::Cell<Dim>>(id14, Nnodes, element);

  //! create 15 cells
  //! add neighbours

  //! cell0
  cell0->add_neighbour(id1);
  cell0->add_neighbour(id5);

  //! cell1
  cell1->add_neighbour(id0);
  cell1->add_neighbour(id2);
  cell1->add_neighbour(id6);

  //! cell2
  cell2->add_neighbour(id1);
  cell2->add_neighbour(id3);
  cell2->add_neighbour(id7);

  //! cell3
  cell3->add_neighbour(id2);
  cell3->add_neighbour(id4);
  cell3->add_neighbour(id8);

  //! cell4
  cell4->add_neighbour(id3);
  cell4->add_neighbour(id9);

  //! cell5
  cell5->add_neighbour(id0);
  cell5->add_neighbour(id6);
  cell5->add_neighbour(id10);

  //! cell6
  cell6->add_neighbour(id1);
  cell6->add_neighbour(id5);
  cell6->add_neighbour(id7);
  cell6->add_neighbour(id11);

  // cell7
  cell7->add_neighbour(id2);
  cell7->add_neighbour(id6);
  cell7->add_neighbour(id8);
  cell7->add_neighbour(id12);

  // cell8
  cell8->add_neighbour(id3);
  cell8->add_neighbour(id7);
  cell8->add_neighbour(id9);
  cell8->add_neighbour(id13);

  // cell9
  cell9->add_neighbour(id4);
  cell9->add_neighbour(id8);
  cell9->add_neighbour(id14);

  // cell10
  cell10->add_neighbour(id5);
  cell10->add_neighbour(id11);

  // cell11
  cell11->add_neighbour(id6);
  cell11->add_neighbour(id10);
  cell11->add_neighbour(id12);

  // cell12
  cell12->add_neighbour(id7);
  cell12->add_neighbour(id11);
  cell12->add_neighbour(id13);

  // cell13
  cell13->add_neighbour(id8);
  cell13->add_neighbour(id12);
  cell13->add_neighbour(id14);

  // cell14
  cell14->add_neighbour(id9);
  cell14->add_neighbour(id13);

  // adding neighour complete
  // Cell container
  auto cellcontainer = std::make_shared<mpm::Container<mpm::Cell<Dim>>>();

  // add cells into container
  cellcontainer->add(cell0);
  cellcontainer->add(cell1);
  cellcontainer->add(cell2);
  cellcontainer->add(cell3);
  cellcontainer->add(cell4);
  cellcontainer->add(cell5);
  cellcontainer->add(cell6);
  cellcontainer->add(cell7);
  cellcontainer->add(cell8);
  cellcontainer->add(cell9);
  cellcontainer->add(cell10);
  cellcontainer->add(cell11);
  cellcontainer->add(cell12);
  cellcontainer->add(cell13);
  cellcontainer->add(cell14);

  // simluate number of tasks
  int num_threads = 3;

  // initialize graph
  mpm::Graph<Dim> graph1 = mpm::Graph<Dim>(*cellcontainer, num_threads, 0);
  mpm::Graph<Dim> graph2 = mpm::Graph<Dim>(*cellcontainer, num_threads, 1);
  mpm::Graph<Dim> graph3 = mpm::Graph<Dim>(*cellcontainer, num_threads, 2);

  // Check graph structure
  SECTION("Check graph initialize function") {
    // check element in xadj in graph1
    REQUIRE(graph1.xadj()[0] == 0);
    REQUIRE(graph1.xadj()[1] == 2);
    REQUIRE(graph1.xadj()[2] == 5);
    REQUIRE(graph1.xadj()[3] == 8);
    REQUIRE(graph1.xadj()[4] == 11);
    REQUIRE(graph1.xadj()[5] == 13);

    // check element in xadj in graph2
    REQUIRE(graph2.xadj()[0] == 0);
    REQUIRE(graph2.xadj()[1] == 3);
    REQUIRE(graph2.xadj()[2] == 7);
    REQUIRE(graph2.xadj()[3] == 11);
    REQUIRE(graph2.xadj()[4] == 15);
    REQUIRE(graph2.xadj()[5] == 18);

    // check element in xadj in graph3
    REQUIRE(graph3.xadj()[0] == 0);
    REQUIRE(graph3.xadj()[1] == 2);
    REQUIRE(graph3.xadj()[2] == 5);
    REQUIRE(graph3.xadj()[3] == 8);
    REQUIRE(graph3.xadj()[4] == 11);
    REQUIRE(graph3.xadj()[5] == 13);

    // check element in adjncy in graph1
    REQUIRE(graph1.adjncy()[0] == 1);
    REQUIRE(graph1.adjncy()[1] == 5);
    REQUIRE(graph1.adjncy()[2] == 0);
    REQUIRE(graph1.adjncy()[3] == 2);
    REQUIRE(graph1.adjncy()[4] == 6);
    REQUIRE(graph1.adjncy()[5] == 1);
    REQUIRE(graph1.adjncy()[6] == 3);
    REQUIRE(graph1.adjncy()[7] == 7);
    REQUIRE(graph1.adjncy()[8] == 2);
    REQUIRE(graph1.adjncy()[9] == 4);
    REQUIRE(graph1.adjncy()[10] == 8);
    REQUIRE(graph1.adjncy()[11] == 3);
    REQUIRE(graph1.adjncy()[12] == 9);

    // check element in adjncy in graph1
    REQUIRE(graph2.adjncy()[0] == 0);
    REQUIRE(graph2.adjncy()[1] == 6);
    REQUIRE(graph2.adjncy()[2] == 10);
    REQUIRE(graph2.adjncy()[3] == 1);
    REQUIRE(graph2.adjncy()[4] == 5);
    REQUIRE(graph2.adjncy()[5] == 7);
    REQUIRE(graph2.adjncy()[6] == 11);
    REQUIRE(graph2.adjncy()[7] == 2);
    REQUIRE(graph2.adjncy()[8] == 6);
    REQUIRE(graph2.adjncy()[9] == 8);
    REQUIRE(graph2.adjncy()[10] == 12);
    REQUIRE(graph2.adjncy()[11] == 3);
    REQUIRE(graph2.adjncy()[12] == 7);
    REQUIRE(graph2.adjncy()[13] == 9);
    REQUIRE(graph2.adjncy()[14] == 13);
    REQUIRE(graph2.adjncy()[15] == 4);
    REQUIRE(graph2.adjncy()[16] == 8);
    REQUIRE(graph2.adjncy()[17] == 14);

    // check element in adjncy in graph1
    REQUIRE(graph3.adjncy()[0] == 5);
    REQUIRE(graph3.adjncy()[1] == 11);
    REQUIRE(graph3.adjncy()[2] == 6);
    REQUIRE(graph3.adjncy()[3] == 10);
    REQUIRE(graph3.adjncy()[4] == 12);
    REQUIRE(graph3.adjncy()[5] == 7);
    REQUIRE(graph3.adjncy()[6] == 11);
    REQUIRE(graph3.adjncy()[7] == 13);
    REQUIRE(graph3.adjncy()[8] == 8);
    REQUIRE(graph3.adjncy()[9] == 12);
    REQUIRE(graph3.adjncy()[10] == 14);
    REQUIRE(graph3.adjncy()[11] == 9);
    REQUIRE(graph3.adjncy()[12] == 13);
  }
}

// MPI
#ifdef USE_MPI
//! Check domain decomposition
TEST_CASE("Graph Partitioning in 2D", "[mpi][graph][2D]") {
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

  SECTION("Transfer particles in mesh in nonrank MPI cells") {

    // Get number of MPI ranks
    int mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    int receiver = 1;

    // Assign material
    unsigned mid = 1;
    // Initialise material
    Json jmaterial;
    jmaterial["density"] = 1000.;
    jmaterial["youngs_modulus"] = 1.0E+7;
    jmaterial["poisson_ratio"] = 0.3;

    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "LinearElastic2D", std::move(mid), jmaterial);

    std::map<unsigned, std::shared_ptr<mpm::Material<Dim>>> materials;
    materials[mid] = material;

    // 4-noded quadrilateral element
    std::shared_ptr<mpm::Element<Dim>> element =
        Factory<mpm::Element<Dim>>::instance()->create("ED2Q4");

    // Cell 1
    mpm::Index id1 = 0;
    Eigen::Vector2d coords;
    coords.setZero();

    // Mesh
    auto mesh = std::make_shared<mpm::Mesh<Dim>>(0);
    // Check mesh is inactive (false)
    REQUIRE(mesh->status() == false);

    // Define nodes
    std::shared_ptr<mpm::NodeBase<Dim>> node0 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

    // Add node 0 and check
    REQUIRE(mesh->add_node(node0) == true);

    coords << 2., 0.;
    std::shared_ptr<mpm::NodeBase<Dim>> node1 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);

    // Add node 1 and check
    REQUIRE(mesh->add_node(node1) == true);

    coords << 2., 2.;
    std::shared_ptr<mpm::NodeBase<Dim>> node2 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);

    // Add node 2 and check
    REQUIRE(mesh->add_node(node2) == true);

    coords << 0., 2.;
    std::shared_ptr<mpm::NodeBase<Dim>> node3 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);

    // Add node 3 and check
    REQUIRE(mesh->add_node(node3) == true);

    // Create cell0
    coords.setZero();
    auto cell0 = std::make_shared<mpm::Cell<Dim>>(0, Nnodes, element);

    // Add nodes to cell
    cell0->add_node(0, node0);
    cell0->add_node(1, node1);
    cell0->add_node(2, node2);
    cell0->add_node(3, node3);

    // Initialize cell
    REQUIRE(cell0->initialise() == true);

    // Initialize material models
    mesh->initialise_material_models(materials);

    // Add cell 1 and check
    REQUIRE(mesh->add_cell(cell0) == true);

    // Cell1
    coords << 4., 0.;
    std::shared_ptr<mpm::NodeBase<Dim>> node4 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(4, coords);

    // Add node 4 and check
    REQUIRE(mesh->add_node(node4) == true);

    coords << 4., 2.;
    std::shared_ptr<mpm::NodeBase<Dim>> node5 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(5, coords);

    // Add node 5 and check
    REQUIRE(mesh->add_node(node5) == true);

    auto cell1 = std::make_shared<mpm::Cell<Dim>>(1, Nnodes, element);

    // Add nodes to cell
    cell1->add_node(0, node1);
    cell1->add_node(1, node4);
    cell1->add_node(2, node5);
    cell1->add_node(3, node2);

    // Initialize cell
    REQUIRE(cell1->initialise() == true);

    // Initialize material models
    mesh->initialise_material_models(materials);

    // Add cell 1 and check
    REQUIRE(mesh->add_cell(cell1) == true);

    // Cell2
    coords << 0., 4.;
    std::shared_ptr<mpm::NodeBase<Dim>> node6 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(6, coords);

    // Add node 6 and check
    REQUIRE(mesh->add_node(node6) == true);

    coords << 2., 4.;
    std::shared_ptr<mpm::NodeBase<Dim>> node7 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(7, coords);

    // Add node 7 and check
    REQUIRE(mesh->add_node(node7) == true);

    auto cell2 = std::make_shared<mpm::Cell<Dim>>(2, Nnodes, element);

    // Add nodes to cell
    cell2->add_node(0, node3);
    cell2->add_node(1, node2);
    cell2->add_node(2, node7);
    cell2->add_node(3, node6);

    // Initialize cell
    REQUIRE(cell2->initialise() == true);

    // Initialize material models
    mesh->initialise_material_models(materials);

    // Add cell 2 and check
    REQUIRE(mesh->add_cell(cell2) == true);

    // Cell3
    coords << 4., 4.;
    std::shared_ptr<mpm::NodeBase<Dim>> node8 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(8, coords);

    // Add node 4 and check
    REQUIRE(mesh->add_node(node8) == true);

    auto cell3 = std::make_shared<mpm::Cell<Dim>>(3, Nnodes, element);

    // Add nodes to cell
    cell3->add_node(0, node2);
    cell3->add_node(1, node5);
    cell3->add_node(2, node8);
    cell3->add_node(3, node7);

    // Initialize cell
    REQUIRE(cell3->initialise() == true);

    // Initialize material models
    mesh->initialise_material_models(materials);

    // Add cell 3 and check
    REQUIRE(mesh->add_cell(cell3) == true);

    // Compute cell neighbours
    mesh->compute_cell_neighbours();

    SECTION("Decompose mesh graph") {
      // Create graph
      auto graph =
          std::make_shared<mpm::Graph<Dim>>(mesh->cells(), mpi_size, mpi_rank);
      // Initialize MPI
      MPI_Comm comm;
      MPI_Comm_dup(MPI_COMM_WORLD, &comm);

      // Check number of ghost cells
      REQUIRE(mesh->nghost_cells() == 0);
      REQUIRE(mesh->nlocal_ghost_cells() == 0);

      // Create partition
      bool graphpartition = graph->create_partitions(&comm);

      // Collect the partitions
      graph->collect_partitions(mesh->ncells(), mpi_size, mpi_rank, &comm);

      // Delete all the particles which is not in local task parititon
      mesh->remove_all_nonrank_particles();

      // Identify shared nodes across MPI domains
      mesh->find_domain_shared_nodes();

      // Identify ghost boundary cells
      mesh->find_ghost_boundary_cells();

      if (mpi_rank == 4) {
        REQUIRE(mesh->nghost_cells() == 3);
        REQUIRE(mesh->nlocal_ghost_cells() == 1);
      }
    }
  }
}

//! Check domain decomposition
TEST_CASE("Graph Partitioning in 3D", "[mpi][graph][3D]") {
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

  SECTION("Transfer particles in mesh in nonrank MPI cells") {

    // Get number of MPI ranks
    int mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    int receiver = 1;

    // Assign material
    unsigned mid = 1;
    // Initialise material
    Json jmaterial;
    jmaterial["density"] = 1000.;
    jmaterial["youngs_modulus"] = 1.0E+7;
    jmaterial["poisson_ratio"] = 0.3;

    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "LinearElastic3D", std::move(mid), jmaterial);

    std::map<unsigned, std::shared_ptr<mpm::Material<Dim>>> materials;
    materials[mid] = material;

    // 8-noded hexahedron element
    std::shared_ptr<mpm::Element<Dim>> element =
        Factory<mpm::Element<Dim>>::instance()->create("ED3H8");

    // Cell 1
    mpm::Index id1 = 0;
    Eigen::Vector3d coords;
    coords.setZero();

    // Mesh
    auto mesh = std::make_shared<mpm::Mesh<Dim>>(0);
    // Check mesh is inactive (false)
    REQUIRE(mesh->status() == false);

    // Define nodes
    // Define nodes
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

    // Initialize cell
    REQUIRE(cell0->initialise() == true);

    // Initialize material models
    mesh->initialise_material_models(materials);

    coords << 0, 0, 4;
    std::shared_ptr<mpm::NodeBase<Dim>> node8 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(8, coords);
    REQUIRE(mesh->add_node(node8) == true);

    coords << 2, 0, 4;
    std::shared_ptr<mpm::NodeBase<Dim>> node9 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(9, coords);
    REQUIRE(mesh->add_node(node9) == true);

    coords << 2, 2, 4;
    std::shared_ptr<mpm::NodeBase<Dim>> node10 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(10, coords);
    REQUIRE(mesh->add_node(node10) == true);

    coords << 0, 2, 4;
    std::shared_ptr<mpm::NodeBase<Dim>> node11 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(11, coords);
    REQUIRE(mesh->add_node(node11) == true);

    auto cell1 = std::make_shared<mpm::Cell<Dim>>(1, Nnodes, element);

    // Add nodes to cell
    cell1->add_node(0, node4);
    cell1->add_node(1, node5);
    cell1->add_node(2, node6);
    cell1->add_node(3, node7);
    cell1->add_node(4, node8);
    cell1->add_node(5, node9);
    cell1->add_node(6, node10);
    cell1->add_node(7, node11);
    // Initialize cell
    REQUIRE(cell1->initialise() == true);

    coords << 4, 0, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node12 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(12, coords);
    REQUIRE(mesh->add_node(node12) == true);

    coords << 4, 2, 0;
    std::shared_ptr<mpm::NodeBase<Dim>> node13 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(13, coords);
    REQUIRE(mesh->add_node(node13) == true);

    coords << 4, 0, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node14 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(14, coords);
    REQUIRE(mesh->add_node(node14) == true);

    coords << 4, 2, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node15 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(15, coords);
    REQUIRE(mesh->add_node(node15) == true);

    coords << 4, 0, 4;
    std::shared_ptr<mpm::NodeBase<Dim>> node16 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(16, coords);
    REQUIRE(mesh->add_node(node16) == true);

    coords << 4, 2, 4;
    std::shared_ptr<mpm::NodeBase<Dim>> node17 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(17, coords);
    REQUIRE(mesh->add_node(node17) == true);

    auto cell2 = std::make_shared<mpm::Cell<Dim>>(2, Nnodes, element);

    // Add nodes to cell
    cell2->add_node(0, node1);
    cell2->add_node(1, node12);
    cell2->add_node(2, node13);
    cell2->add_node(3, node2);
    cell2->add_node(4, node5);
    cell2->add_node(5, node14);
    cell2->add_node(6, node15);
    cell2->add_node(7, node6);
    // Initialize cell
    REQUIRE(cell2->initialise() == true);

    auto cell3 = std::make_shared<mpm::Cell<Dim>>(3, Nnodes, element);
    // Add nodes to cell
    cell3->add_node(0, node5);
    cell3->add_node(1, node14);
    cell3->add_node(2, node15);
    cell3->add_node(3, node6);
    cell3->add_node(4, node9);
    cell3->add_node(5, node16);
    cell3->add_node(6, node17);
    cell3->add_node(7, node10);
    // Initialize cell
    REQUIRE(cell3->initialise() == true);

    // Initialize material models
    mesh->initialise_material_models(materials);

    // Add cells to mesh and check
    REQUIRE(mesh->add_cell(cell0) == true);
    REQUIRE(mesh->add_cell(cell1) == true);
    REQUIRE(mesh->add_cell(cell2) == true);
    REQUIRE(mesh->add_cell(cell3) == true);

    // Compute cell neighbours
    mesh->compute_cell_neighbours();

    SECTION("Decompose mesh graph") {
      // Create graph
      auto graph =
          std::make_shared<mpm::Graph<Dim>>(mesh->cells(), mpi_size, mpi_rank);
      // Initialize MPI
      MPI_Comm comm;
      MPI_Comm_dup(MPI_COMM_WORLD, &comm);

      // Check number of ghost cells
      REQUIRE(mesh->nghost_cells() == 0);
      REQUIRE(mesh->nlocal_ghost_cells() == 0);

      // Create partition
      bool graphpartition = graph->create_partitions(&comm);

      // Collect the partitions
      graph->collect_partitions(mesh->ncells(), mpi_size, mpi_rank, &comm);

      // Delete all the particles which is not in local task parititon
      mesh->remove_all_nonrank_particles();

      // Identify shared nodes across MPI domains
      mesh->find_domain_shared_nodes();

      // Identify ghost boundary cells
      mesh->find_ghost_boundary_cells();

      if (mpi_rank == 4) {
        REQUIRE(mesh->nghost_cells() == 3);
        REQUIRE(mesh->nlocal_ghost_cells() == 1);
      }
    }
  }
}
#endif
#endif

#include <functional>
#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

#include "cell.h"
#include "container.h"
#include "element.h"
#include "graph.h"
#include "hexahedron_element.h"
#include "node.h"
#include "quadrilateral_element.h"

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

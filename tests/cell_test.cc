#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

#include "cell.h"
#include "hex_shapefn.h"
#include "node.h"
#include "quad_shapefn.h"
#include "shapefn.h"

//! \brief Check cell class for 2D case
TEST_CASE("Cell is checked for 2D case", "[cell][2D]") {
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

  Eigen::Vector2d coords;
  coords.setZero();

  std::shared_ptr<mpm::NodeBase<Dim>> node0 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

  coords << 0, 1;
  std::shared_ptr<mpm::NodeBase<Dim>> node1 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);

  coords << 1, 1;
  std::shared_ptr<mpm::NodeBase<Dim>> node2 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);

  coords << 1, 0;
  std::shared_ptr<mpm::NodeBase<Dim>> node3 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);

  //! Check Cell IDs
  SECTION("Check cell ids") {
    //! Check for id = 0
    SECTION("Cell id is zero") {
      mpm::Index id = 0;
      auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes);
      REQUIRE(cell->id() == 0);
    }

    SECTION("Cell id is positive") {
      //! Check for id is a positive value
      mpm::Index id = std::numeric_limits<mpm::Index>::max();
      auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes);
      REQUIRE(cell->id() == std::numeric_limits<mpm::Index>::max());
    }
  }

  SECTION("Add nodes") {
    mpm::Index id = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes);
    cell->add_node(0, node0);
    cell->add_node(1, node1);
    cell->add_node(2, node2);
    cell->add_node(3, node3);
    REQUIRE(cell->nnodes() == 4);

    SECTION("Compute volume of a cell") {
      std::shared_ptr<mpm::ShapeFn<Dim>> shapefn =
          std::make_shared<mpm::QuadrilateralShapeFn<Dim, 4>>();
      cell->shapefn(shapefn);
      REQUIRE(cell->nfunctions() == 4);
      // Compute volume
      REQUIRE(cell->volume() ==
              Approx(std::numeric_limits<double>::max()).epsilon(Tolerance));
      cell->compute_volume();
      REQUIRE(cell->volume() == Approx(1.0).epsilon(Tolerance));

      SECTION("Check if a point is in a cell") {
        // Check point in cell
        Eigen::Vector2d point;
        point << 0.5, 0.5;
        REQUIRE(cell->point_in_cell(point) == true);

        // Check point on vertex
        point << 0., 0.;
        REQUIRE(cell->point_in_cell(point) == true);

        // Check point on edge
        point << 0.5, 0.;
        REQUIRE(cell->point_in_cell(point) == true);

        // Check point outside
        point << -2, 2.;
        REQUIRE(cell->point_in_cell(point) == false);
      }
    }
  }

  SECTION("Add neighbours") {
    auto cell = std::make_shared<mpm::Cell<Dim>>(0, Nnodes);
    auto neighbourcell = std::make_shared<mpm::Cell<Dim>>(1, Nnodes);
    REQUIRE(cell->nneighbours() == 0);
    cell->add_neighbour(0, neighbourcell);
    REQUIRE(cell->nneighbours() == 1);
  }

  SECTION("Check shape functions") {
    mpm::Index id = 0;
    // Check 4-noded function
    SECTION("Check 4-noded Quadrilateral") {
      std::shared_ptr<mpm::ShapeFn<Dim>> shapefn =
          std::make_shared<mpm::QuadrilateralShapeFn<Dim, 4>>();
      auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, shapefn);
      REQUIRE(cell->nfunctions() == 4);
    }
    // Check 8-noded function
    SECTION("Check 8-noded Quadrilateral") {
      std::shared_ptr<mpm::ShapeFn<Dim>> shapefn =
          std::make_shared<mpm::QuadrilateralShapeFn<Dim, 8>>();
      auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, shapefn);
      REQUIRE(cell->nfunctions() == 8);
    }
    // Check 9-noded function
    SECTION("Check 9-noded Quadrilateral") {
      std::shared_ptr<mpm::ShapeFn<Dim>> shapefn =
          std::make_shared<mpm::QuadrilateralShapeFn<Dim, 9>>();
      auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, shapefn);
      REQUIRE(cell->nfunctions() == 9);
    }
  }

  SECTION("Test particle addition deletion") {
    mpm::Index pid = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(0, Nnodes);
    REQUIRE(cell->status() == false);
    bool status = cell->add_particle_id(pid);
    REQUIRE(status == true);
    REQUIRE(cell->status() == true);
    cell->remove_particle_id(pid);
    REQUIRE(cell->status() == false);
  }

  SECTION("Test particle information mapping") {
    mpm::Index id = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes);
    cell->add_node(0, node0);
    cell->add_node(1, node3);
    cell->add_node(2, node2);
    cell->add_node(3, node1);
    Eigen::Vector2d xi;
    xi << 0., 0.;
    Eigen::Matrix<double,0,1> pmass;
    pmass << 4.;
    Eigen::Vector2d pvelocity;
    pvelocity << 1., 1.;
    SECTION("Check particle mass and velocity mapping") {
      cell->assign_mass_to_nodes(xi, pmass);
      cell->assign_momentum_to_nodes(xi, pmass, pvelocity);
      // cell->interpolate_velocity() == pvelocity
    }
  }
}



//! \brief Check cell class for 3D case
TEST_CASE("Cell is checked for 3D case", "[cell][3D]") {
  // Dimension
  const unsigned Dim = 3;
  // Degrees of freedom
  const unsigned Dof = 6;
  // Number of phases
  const unsigned Nphases = 1;
  // Number of nodes per cell
  const unsigned Nnodes = 8;
  // Tolerance
  const double Tolerance = 1.E-7;

  // Coordinaates
  Eigen::Vector3d coords;

  coords << 0, 0, 0;
  std::shared_ptr<mpm::NodeBase<Dim>> node0 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

  coords << 1, 0, 0;
  std::shared_ptr<mpm::NodeBase<Dim>> node1 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);

  coords << 1, 1, 0;
  std::shared_ptr<mpm::NodeBase<Dim>> node2 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);

  coords << 0, 1, 0;
  std::shared_ptr<mpm::NodeBase<Dim>> node3 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);

  coords << 0, 0, 1;
  std::shared_ptr<mpm::NodeBase<Dim>> node4 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(4, coords);

  coords << 1, 0, 1;
  std::shared_ptr<mpm::NodeBase<Dim>> node5 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(5, coords);

  coords << 1, 1, 1;
  std::shared_ptr<mpm::NodeBase<Dim>> node6 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(6, coords);

  coords << 0, 1, 1;
  std::shared_ptr<mpm::NodeBase<Dim>> node7 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(7, coords);

  //! Check Cell IDs
  SECTION("Check cell ids") {
    //! Check for id = 0
    SECTION("Cell id is zero") {
      mpm::Index id = 0;
      auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes);
      REQUIRE(cell->id() == 0);
    }

    SECTION("Cell id is positive") {
      //! Check for id is a positive value
      mpm::Index id = std::numeric_limits<mpm::Index>::max();
      auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes);
      REQUIRE(cell->id() == std::numeric_limits<mpm::Index>::max());
    }
  }

  // Check node additions
  SECTION("Add nodes") {
    mpm::Index id = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes);
    cell->add_node(0, node0);
    cell->add_node(1, node1);
    cell->add_node(2, node2);
    cell->add_node(3, node3);
    cell->add_node(4, node4);
    cell->add_node(5, node5);
    cell->add_node(6, node6);
    cell->add_node(7, node7);
    REQUIRE(cell->nnodes() == 8);

    SECTION("Compute volume of a cell") {
      std::shared_ptr<mpm::ShapeFn<Dim>> shapefn =
          std::make_shared<mpm::HexahedronShapeFn<Dim, 8>>();
      cell->shapefn(shapefn);
      REQUIRE(cell->nfunctions() == 8);
      // Compute volume
      REQUIRE(cell->volume() ==
              Approx(std::numeric_limits<double>::max()).epsilon(Tolerance));
      cell->compute_volume();
      REQUIRE(cell->volume() == Approx(1.0).epsilon(Tolerance));

      SECTION("Check if a point is in a cell") {
        // Check point in cell
        Eigen::Vector3d point;
        point << 0.5, 0.5, 0.5;
        REQUIRE(cell->point_in_cell(point) == true);

        // Check point on vertex
        point << 0., 0., 0.;
        REQUIRE(cell->point_in_cell(point) == true);

        // Check point on edge
        point << 0.5, 0., 0.;
        REQUIRE(cell->point_in_cell(point) == true);

        // Check point on surface
        point << 0.5, 0.5, 0.;
        REQUIRE(cell->point_in_cell(point) == true);

        // Check point outside
        point << 2., 2., 2.;
        REQUIRE(cell->point_in_cell(point) == false);
      }
    }
  }

  SECTION("Add neighbours") {
    auto cell = std::make_shared<mpm::Cell<Dim>>(0, Nnodes);
    auto neighbourcell = std::make_shared<mpm::Cell<Dim>>(1, Nnodes);
    REQUIRE(cell->nneighbours() == 0);
    cell->add_neighbour(0, neighbourcell);
    REQUIRE(cell->nneighbours() == 1);
  }

  SECTION("Check shape functions") {
    mpm::Index id = 0;
    // Check 8-noded function
    SECTION("Check 8-noded Hexahedron") {
      std::shared_ptr<mpm::ShapeFn<Dim>> shapefn =
          std::make_shared<mpm::HexahedronShapeFn<Dim, 8>>();
      auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, shapefn);
      REQUIRE(cell->nfunctions() == 8);
    }
    // Check 20-noded function
    SECTION("Check 20-noded Hexahedron") {
      std::shared_ptr<mpm::ShapeFn<Dim>> shapefn =
          std::make_shared<mpm::HexahedronShapeFn<Dim, 20>>();
      auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, shapefn);
      REQUIRE(cell->nfunctions() == 20);
    }
  }

  SECTION("Test particle addition deletion") {
    mpm::Index pid = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(0, Nnodes);
    REQUIRE(cell->status() == false);
    bool status = cell->add_particle_id(pid);
    REQUIRE(status == true);
    REQUIRE(cell->status() == true);
    cell->remove_particle_id(pid);
    REQUIRE(cell->status() == false);
  }
}

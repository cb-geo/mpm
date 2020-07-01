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

  // Shape function
  // 4-noded quadrilateral shape functions
  std::shared_ptr<mpm::Element<Dim>> element =
      Factory<mpm::Element<Dim>>::instance()->create("ED2Q4");

  Eigen::Vector2d coords;
  coords.setZero();

  std::shared_ptr<mpm::NodeBase<Dim>> node0 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(10, coords);

  coords << 2., 0.;
  std::shared_ptr<mpm::NodeBase<Dim>> node1 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(11, coords);

  coords << 2., 2.;
  std::shared_ptr<mpm::NodeBase<Dim>> node2 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(12, coords);

  coords << 0., 2.;
  std::shared_ptr<mpm::NodeBase<Dim>> node3 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(13, coords);

  //! Check Cell IDs
  SECTION("Check cell ids") {
    //! Check for id = 0
    SECTION("Cell id is zero") {
      mpm::Index id = 0;
      auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);
      REQUIRE(cell->id() == 0);
    }

    SECTION("Cell id is positive") {
      //! Check for id is a positive value
      mpm::Index id = std::numeric_limits<mpm::Index>::max();
      auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);
      REQUIRE(cell->id() == std::numeric_limits<mpm::Index>::max());
    }
  }

  //! Check cell rank
  SECTION("Check cell rank") {
    mpm::Index id = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element, true);
    REQUIRE(cell->rank() == 0);
    cell->rank(1);
    REQUIRE(cell->rank() == 1);
  }

  SECTION("Add nodes") {
    mpm::Index id = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element, true);
    // Test initialisation for failure
    REQUIRE(cell->initialise() == false);

    // Check if cell is initialised, before addition of nodes
    REQUIRE(cell->is_initialised() == false);

    REQUIRE(cell->add_node(0, node0) == true);
    REQUIRE(cell->add_node(1, node1) == true);
    REQUIRE(cell->add_node(2, node2) == true);
    REQUIRE(cell->add_node(3, node3) == true);
    REQUIRE(cell->nnodes() == 4);

    // Test failing add node
    coords << 1., 1.;
    std::shared_ptr<mpm::NodeBase<Dim>> node4 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(4, coords);
    // Fail add node
    REQUIRE(cell->add_node(4, node4) == false);

    // Check if cell is initialised, after addition of nodes
    REQUIRE(cell->is_initialised() == false);

    // Check cell length at initialisation
    REQUIRE(cell->mean_length() == std::numeric_limits<double>::max());
    // Check volume before initialisation
    REQUIRE(cell->volume() ==
            Approx(std::numeric_limits<double>::lowest()).epsilon(Tolerance));

    // Initialise cell
    REQUIRE(cell->initialise() == true);
    // Check if cell is initialised, after addition of nodes
    REQUIRE(cell->is_initialised() == true);

    // Check MPI rank
    SECTION("Assign and check MPI rank on nodes") {
      cell->rank(1);
      REQUIRE(cell->rank() == 1);
      cell->assign_mpi_rank_to_nodes();

      REQUIRE(node0->mpi_ranks().size() == 1);
      REQUIRE(node1->mpi_ranks().size() == 1);
      REQUIRE(node2->mpi_ranks().size() == 1);
      REQUIRE(node3->mpi_ranks().size() == 1);

      // Update MPI rank of cell
      REQUIRE(cell->previous_mpirank() == 0);
      REQUIRE(cell->rank() == 1);
      cell->rank(2);
      REQUIRE(cell->rank() == 2);
      REQUIRE(cell->previous_mpirank() == 1);
    }

    // Check cell length calculation
    SECTION("Compute mean length of cell") {
      // Length of the cell
      const double length = 2.0;

      // Check cell length calculations
      REQUIRE(cell->mean_length() == Approx(length).epsilon(Tolerance));
    }

    // Check cell coordinates
    SECTION("Check cell nodal coordinates") {
      Eigen::Matrix<double, 4, Dim> nodal_coords;
      // clang-format off
      nodal_coords << 0., 0.,
                      2., 0.,
                      2., 2.,
                      0., 2.;
      // clang-format on
      auto coords = cell->nodal_coordinates();

      // Check if the sizes are the same
      REQUIRE(coords.size() == nodal_coords.size());

      for (unsigned i = 0; i < 4; ++i)
        for (unsigned j = 0; j < Dim; ++j)
          REQUIRE(coords(i, j) ==
                  Approx(nodal_coords(i, j)).epsilon(Tolerance));
    }

    // Check side indices
    SECTION("Check side indices") {
      const auto side_indices = cell->side_node_pairs();
      REQUIRE(side_indices.size() == 4);
    }
    // Check shape functions
    SECTION("Check shape functions") {
      const auto sf_ptr = cell->element_ptr();
      REQUIRE(sf_ptr->nfunctions() == element->nfunctions());
    }

    // Check centroid calculation
    SECTION("Compute centroid of a cell") {
      REQUIRE(cell->nfunctions() == 4);

      // Compute centroid
      cell->compute_centroid();
      auto check_centroid = cell->centroid();

      // Centroid
      Eigen::Matrix<double, Dim, 1> centroid;
      centroid << 1.0, 1.0;

      // Check number of coordinates
      REQUIRE(check_centroid.size() == centroid.size());

      // Check centroid calculations
      for (unsigned i = 0; i < check_centroid.size(); ++i)
        REQUIRE(check_centroid[i] == Approx(centroid[i]).epsilon(Tolerance));
    }

    // Check cell volume
    SECTION("Compute volume of a cell") {
      REQUIRE(cell->nfunctions() == 4);

      REQUIRE(cell->volume() == Approx(4.0).epsilon(Tolerance));

      // Check if cell is initialised, after volume calculation
      REQUIRE(cell->is_initialised() == true);

      SECTION("Check negative or zero volume cell") {
        mpm::Index id = 0;
        auto cell1 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);
        REQUIRE(cell1->add_node(0, node1) == true);
        REQUIRE(cell1->add_node(1, node1) == true);
        REQUIRE(cell1->add_node(2, node1) == true);
        REQUIRE(cell1->add_node(3, node1) == true);
        REQUIRE(cell1->nnodes() == 4);
        cell1->compute_volume();
      }

      SECTION("Check if a point is in a cell") {
        Eigen::Vector2d xi;
        // Check point in cell
        Eigen::Vector2d point;
        point << 0.5, 0.5;

        REQUIRE(cell->point_in_cartesian_cell(point) == true);
        REQUIRE(cell->is_point_in_cell(point, &xi) == true);

        // Check point on vertex
        point << 0., 0.;
        REQUIRE(cell->point_in_cartesian_cell(point) == true);
        REQUIRE(cell->is_point_in_cell(point, &xi) == true);
        REQUIRE(xi(0) == -1. + std::numeric_limits<double>::epsilon());
        REQUIRE(xi(1) == -1. + std::numeric_limits<double>::epsilon());

        // Check point on edge
        point << 0.5, 0.;
        REQUIRE(cell->point_in_cartesian_cell(point) == true);
        REQUIRE(cell->is_point_in_cell(point, &xi) == true);
        REQUIRE(xi(0) == -0.5);
        REQUIRE(xi(1) == -1. + std::numeric_limits<double>::epsilon());

        // Check point outside
        point << -2, 2.;
        REQUIRE(cell->point_in_cartesian_cell(point) == false);
        REQUIRE(cell->is_point_in_cell(point, &xi) == false);
      }

      // Find local coordinates of a point in a cell
      SECTION("Find local coordinates of a point in cell") {
        // Coordinates of a point in real cell
        Eigen::Vector2d point;
        point << 0.5, 1.5;

        // Coordinates of the point in an unit cell
        Eigen::Matrix<double, 2, 1> point_unit_cell;
        point_unit_cell << -0.5, 0.5;

        // Get local coordinates of the point
        auto local_point = cell->local_coordinates_point(point);
        for (unsigned i = 0; i < local_point.size(); ++i)
          REQUIRE(local_point[i] ==
                  Approx(point_unit_cell[i]).epsilon(Tolerance));

        // Use Newton-raphson iteration
        local_point = cell->transform_real_to_unit_cell(point);
        for (unsigned i = 0; i < local_point.size(); ++i)
          REQUIRE(local_point[i] ==
                  Approx(point_unit_cell[i]).epsilon(Tolerance));
      }

      SECTION("Transform real to unit cell, xi is in negative coordinates") {
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

        // Coordinates of a point in real cell
        Eigen::Vector2d point;
        point << 2.1875, 3.25;

        Eigen::Vector2d xi;
        // Test if point is in cell
        REQUIRE(cell->is_point_in_cell(point, &xi) == true);

        // Coordinates of the point in an unit cell
        Eigen::Matrix<double, 2, 1> point_unit_cell;
        point_unit_cell << 0.5, 0.5;

        // Use Newton-raphson iteration to find local coordinates
        auto local_point = cell->transform_real_to_unit_cell(point);
        for (unsigned i = 0; i < local_point.size(); ++i)
          REQUIRE(local_point[i] ==
                  Approx(point_unit_cell[i]).epsilon(Tolerance));

        // Coordinates of a point in real cell
        Eigen::Vector2d point1;
        point1 << 2., 1.;

        // Test if point is in cell
        REQUIRE(cell->is_point_in_cell(point1, &xi) == true);

        // Coordinates of the point in an unit cell
        Eigen::Matrix<double, 2, 1> point_unit_cell1;
        point_unit_cell1 << -1., -1.;

        // Use Newton-raphson iteration to find local coordinates
        auto local_point1 = cell->transform_real_to_unit_cell(point1);
        for (unsigned i = 0; i < local_point1.size(); ++i)
          REQUIRE(local_point1[i] ==
                  Approx(point_unit_cell1[i]).epsilon(Tolerance));
      }

      SECTION("Transform real to unit cell, xi is (-0.25, -0.5)") {
        // Number of nodes in cell
        const unsigned Nnodes = 4;

        // Coordinates
        Eigen::Vector2d coords;

        coords << 0., 0.5;
        std::shared_ptr<mpm::NodeBase<Dim>> node0 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

        coords << 2., 0.5;
        std::shared_ptr<mpm::NodeBase<Dim>> node1 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);

        coords << 3., 1.5;
        std::shared_ptr<mpm::NodeBase<Dim>> node2 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);

        coords << 1., 1.5;
        std::shared_ptr<mpm::NodeBase<Dim>> node3 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);

        // 4-noded quadrilateral shape functions
        std::shared_ptr<mpm::Element<Dim>> element =
            Factory<mpm::Element<Dim>>::instance()->create("ED2Q4");

        mpm::Index id = 0;
        auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

        cell->add_node(0, node0);
        cell->add_node(1, node1);
        cell->add_node(2, node2);
        cell->add_node(3, node3);
        REQUIRE(cell->nnodes() == 4);

        // Initialise cell
        REQUIRE(cell->initialise() == true);

        // Coordinates of a point in real cell
        Eigen::Vector2d point;
        point << 1.0, 0.75;

        // Coordinates of the point in an unit cell
        Eigen::Matrix<double, 2, 1> point_unit_cell;
        point_unit_cell << -0.25, -0.5;

        // Use Newton-raphson iteration to find local coordinates
        auto local_point = cell->transform_real_to_unit_cell(point);
        for (unsigned i = 0; i < local_point.size(); ++i)
          REQUIRE(local_point[i] ==
                  Approx(point_unit_cell[i]).epsilon(Tolerance));
      }

      SECTION("Generate points in cell") {
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

        // Assuming a quadrature of 1x1
        auto points = cell->generate_points();
        REQUIRE(points.size() == 1);

        auto quad = std::make_shared<mpm::QuadrilateralQuadrature<Dim, 1>>();
        auto quadrature = quad->quadratures();
        // Check if the output coordinates match local quadratures
        for (unsigned i = 0; i < points.size(); ++i) {
          auto local_point = cell->transform_real_to_unit_cell(points.at(i));
          for (unsigned j = 0; j < Dim; ++j)
            REQUIRE(local_point[j] ==
                    Approx(quadrature(j, i)).epsilon(Tolerance));
        }

        // Assign quadrature 2x2
        cell->assign_quadrature(2);

        points = cell->generate_points();
        REQUIRE(points.size() == 4);

        auto quad4 = std::make_shared<mpm::QuadrilateralQuadrature<Dim, 4>>();
        quadrature = quad4->quadratures();
        // Check if the output coordinates match local quadratures
        for (unsigned i = 0; i < points.size(); ++i) {
          auto local_point = cell->transform_real_to_unit_cell(points.at(i));
          for (unsigned j = 0; j < Dim; ++j)
            REQUIRE(local_point[j] ==
                    Approx(quadrature(j, i)).epsilon(Tolerance));
        }

        // Assign quadrature 3x3
        cell->assign_quadrature(3);

        points = cell->generate_points();
        REQUIRE(points.size() == 9);

        auto quad9 = std::make_shared<mpm::QuadrilateralQuadrature<Dim, 9>>();
        quadrature = quad9->quadratures();
        // Check if the output coordinates match local quadratures
        for (unsigned i = 0; i < points.size(); ++i) {
          auto local_point = cell->transform_real_to_unit_cell(points.at(i));
          for (unsigned j = 0; j < Dim; ++j)
            REQUIRE(local_point[j] ==
                    Approx(quadrature(j, i)).epsilon(Tolerance));
        }

        // Assign quadrature 4x4
        cell->assign_quadrature(4);

        points = cell->generate_points();
        REQUIRE(points.size() == 16);

        auto quad16 = std::make_shared<mpm::QuadrilateralQuadrature<Dim, 16>>();
        quadrature = quad16->quadratures();
        // Check if the output coordinates match local quadratures
        for (unsigned i = 0; i < points.size(); ++i) {
          auto local_point = cell->transform_real_to_unit_cell(points.at(i));
          for (unsigned j = 0; j < Dim; ++j)
            REQUIRE(local_point[j] ==
                    Approx(quadrature(j, i)).epsilon(Tolerance));
        }
      }
    }
  }

  SECTION("Check shape functions") {
    mpm::Index id = 0;
    // Check 4-noded function
    SECTION("Check 4-noded Quadrilateral") {
      // 4-noded quadrilateral shape functions
      std::shared_ptr<mpm::Element<Dim>> element =
          Factory<mpm::Element<Dim>>::instance()->create("ED2Q4");
      auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);
      REQUIRE(cell->nfunctions() == 4);

      const unsigned nnodes = 8;
      auto cell1 = std::make_shared<mpm::Cell<Dim>>(id, nnodes, element);
    }
    // Check 8-noded function
    SECTION("Check 8-noded Quadrilateral") {
      // 8-noded quadrilateral shape functions
      const unsigned nnodes = 8;
      std::shared_ptr<mpm::Element<Dim>> element =
          Factory<mpm::Element<Dim>>::instance()->create("ED2Q8");
      auto cell = std::make_shared<mpm::Cell<Dim>>(id, nnodes, element);
      REQUIRE(cell->nfunctions() == 8);
    }
    // Check 9-noded function
    SECTION("Check 9-noded Quadrilateral") {
      // 9-noded quadrilateral shape functions
      const unsigned nnodes = 9;
      std::shared_ptr<mpm::Element<Dim>> element =
          Factory<mpm::Element<Dim>>::instance()->create("ED2Q9");
      auto cell = std::make_shared<mpm::Cell<Dim>>(id, nnodes, element);
      REQUIRE(cell->nfunctions() == 9);
    }
  }

  SECTION("Test particle addition deletion") {
    mpm::Index pid = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(0, Nnodes, element);
    REQUIRE(cell->status() == false);
    REQUIRE(cell->nparticles() == 0);
    REQUIRE(cell->add_particle_id(pid) == true);
    REQUIRE(cell->status() == true);
    REQUIRE(cell->nparticles() == 1);
    REQUIRE(cell->particles().size() == 1);
    REQUIRE(cell->particles().at(0) == pid);
    cell->remove_particle_id(pid);
    REQUIRE(cell->status() == false);
    REQUIRE(cell->nparticles() == 0);
    REQUIRE(cell->particles().size() == 0);
  }

  SECTION("Test node status") {
    mpm::Index pid = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(0, Nnodes, element);
    cell->add_node(0, node0);
    cell->add_node(1, node1);
    cell->add_node(2, node2);
    cell->add_node(3, node3);

    std::vector<std::shared_ptr<mpm::NodeBase<Dim>>> nodes;
    nodes.emplace_back(node0);
    nodes.emplace_back(node1);
    nodes.emplace_back(node2);
    nodes.emplace_back(node3);

    // When number of particles are zero
    REQUIRE(cell->nparticles() == 0);
    REQUIRE(cell->status() == false);
    cell->activate_nodes();
    for (const auto& node : nodes) REQUIRE(node->status() == false);

    // Add a particle
    REQUIRE(cell->add_particle_id(pid) == true);
    REQUIRE(cell->status() == true);
    REQUIRE(cell->nparticles() == 1);
    cell->activate_nodes();
    for (const auto& node : nodes) REQUIRE(node->status() == true);

    // Remove a particle
    cell->remove_particle_id(pid);
    REQUIRE(cell->status() == false);
    REQUIRE(cell->nparticles() == 0);
    // initialise nodes
    for (const auto& node : nodes) node->initialise();
    cell->activate_nodes();
    for (const auto& node : nodes) REQUIRE(node->status() == false);
  }

  SECTION("Test nglobal particles") {
    mpm::Index pid = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(0, Nnodes, element);
    REQUIRE(cell->nglobal_particles() == 0);
    cell->nglobal_particles(5);
    REQUIRE(cell->nglobal_particles() == 5);
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

  // Coordinates
  Eigen::Vector3d coords;

  coords << 0, 0, 0;
  std::shared_ptr<mpm::NodeBase<Dim>> node0 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

  coords << 2, 0, 0;
  std::shared_ptr<mpm::NodeBase<Dim>> node1 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);

  coords << 2, 2, 0;
  std::shared_ptr<mpm::NodeBase<Dim>> node2 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);

  coords << 0, 2, 0;
  std::shared_ptr<mpm::NodeBase<Dim>> node3 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);

  coords << 0, 0, 2;
  std::shared_ptr<mpm::NodeBase<Dim>> node4 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(4, coords);

  coords << 2, 0, 2;
  std::shared_ptr<mpm::NodeBase<Dim>> node5 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(5, coords);

  coords << 2, 2, 2;
  std::shared_ptr<mpm::NodeBase<Dim>> node6 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(6, coords);

  coords << 0, 2, 2;
  std::shared_ptr<mpm::NodeBase<Dim>> node7 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(7, coords);

  // 8-noded hexahedron shape functions
  std::shared_ptr<mpm::Element<Dim>> element =
      Factory<mpm::Element<Dim>>::instance()->create("ED3H8");

  //! Check Cell IDs
  SECTION("Check cell ids") {
    //! Check for id = 0
    SECTION("Cell id is zero") {
      mpm::Index id = 0;
      auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);
      REQUIRE(cell->id() == 0);
    }

    SECTION("Cell id is positive") {
      //! Check for id is a positive value
      mpm::Index id = std::numeric_limits<mpm::Index>::max();
      auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);
      REQUIRE(cell->id() == std::numeric_limits<mpm::Index>::max());
    }
  }

  //! Check cell rank
  SECTION("Check cell rank") {
    mpm::Index id = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element, true);
    REQUIRE(cell->rank() == 0);
    cell->rank(1);
    REQUIRE(cell->rank() == 1);
  }

  // Check node additions
  SECTION("Add nodes") {
    mpm::Index id = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element, true);

    // Test initialisation for failure
    REQUIRE(cell->initialise() == false);

    // Check if cell is initialised, before addition of nodes
    REQUIRE(cell->is_initialised() == false);

    REQUIRE(cell->add_node(0, node0) == true);
    REQUIRE(cell->add_node(1, node1) == true);
    REQUIRE(cell->add_node(2, node2) == true);
    REQUIRE(cell->add_node(3, node3) == true);
    REQUIRE(cell->add_node(4, node4) == true);
    REQUIRE(cell->add_node(5, node5) == true);
    REQUIRE(cell->add_node(6, node6) == true);
    REQUIRE(cell->add_node(7, node7) == true);
    REQUIRE(cell->nnodes() == 8);

    // Test failing add node
    coords << 1., 1., 1.;
    std::shared_ptr<mpm::NodeBase<Dim>> node100 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(100, coords);
    // Fail add node
    REQUIRE(cell->add_node(8, node100) == false);

    // Check if cell is initialised, after addition of nodes
    REQUIRE(cell->is_initialised() == false);

    // Check cell length at initialisation
    REQUIRE(cell->mean_length() == std::numeric_limits<double>::max());
    // Check volume before initialisation
    REQUIRE(cell->volume() ==
            Approx(std::numeric_limits<double>::lowest()).epsilon(Tolerance));

    // Initialise cell
    REQUIRE(cell->initialise() == true);
    // Check if cell is initialised, after addition of nodes
    REQUIRE(cell->is_initialised() == true);

    // Check MPI rank
    SECTION("Assign and check MPI rank on nodes") {
      cell->rank(1);
      cell->assign_mpi_rank_to_nodes();

      REQUIRE(node0->mpi_ranks().size() == 1);
      REQUIRE(node1->mpi_ranks().size() == 1);
      REQUIRE(node2->mpi_ranks().size() == 1);
      REQUIRE(node3->mpi_ranks().size() == 1);

      // Update MPI rank of cell
      REQUIRE(cell->previous_mpirank() == 0);
      REQUIRE(cell->rank() == 1);
      cell->rank(2);
      REQUIRE(cell->rank() == 2);
      REQUIRE(cell->previous_mpirank() == 1);
    }

    // Check cell length calculation
    SECTION("Compute mean length of cell") {
      // Length of the cell
      const double length = 2.0;

      // Check cell length calculations
      REQUIRE(cell->mean_length() == Approx(length).epsilon(Tolerance));
    }

    // Check cell coordinates
    SECTION("Check cell nodal coordinates") {
      Eigen::Matrix<double, 8, Dim> nodal_coords;
      // clang-format off
      nodal_coords << 0, 0, 0,
                      2, 0, 0,
                      2, 2, 0,
                      0, 2, 0,
                      0, 0, 2,
                      2, 0, 2,
                      2, 2, 2,
                      0, 2, 2;
      // clang-format on
      auto coords = cell->nodal_coordinates();

      // Check if the sizes are the same
      REQUIRE(coords.size() == nodal_coords.size());

      for (unsigned i = 0; i < 8; ++i)
        for (unsigned j = 0; j < Dim; ++j)
          REQUIRE(coords(i, j) ==
                  Approx(nodal_coords(i, j)).epsilon(Tolerance));
    }

    // Check centroid calculation
    SECTION("Compute centroid of a cell") {
      REQUIRE(cell->nfunctions() == 8);

      // Compute centroid
      auto check_centroid = cell->centroid();

      // Centroid
      Eigen::Matrix<double, Dim, 1> centroid;
      centroid << 1.0, 1.0, 1.0;

      // Check number of coordinates
      REQUIRE(check_centroid.size() == centroid.size());

      // Check centroid calculations
      for (unsigned i = 0; i < check_centroid.size(); ++i)
        REQUIRE(check_centroid[i] == Approx(centroid[i]).epsilon(Tolerance));
    }

    // Check shape functions
    SECTION("Check shape functions") {
      const auto sf_ptr = cell->element_ptr();
      REQUIRE(sf_ptr->nfunctions() == element->nfunctions());
    }

    // Check cell volume calculation
    SECTION("Compute volume of a cell") {
      REQUIRE(cell->nfunctions() == 8);

      REQUIRE(cell->volume() == Approx(8.0).epsilon(Tolerance));

      SECTION("Check negative or zero volume cell") {
        mpm::Index id = 0;
        auto cell1 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);
        REQUIRE(cell1->add_node(0, node1) == true);
        REQUIRE(cell1->add_node(1, node1) == true);
        REQUIRE(cell1->add_node(2, node1) == true);
        REQUIRE(cell1->add_node(3, node1) == true);
        REQUIRE(cell1->add_node(4, node1) == true);
        REQUIRE(cell1->add_node(5, node1) == true);
        REQUIRE(cell1->add_node(6, node1) == true);
        REQUIRE(cell1->add_node(7, node1) == true);
        REQUIRE(cell1->nnodes() == 8);
        cell1->compute_volume();
      }

      SECTION("Check if a point is in a cell") {
        Eigen::Vector3d xi;

        // Check point in cell
        Eigen::Vector3d point;
        point << 0.5, 0.5, 0.5;

        REQUIRE(cell->point_in_cartesian_cell(point) == true);
        REQUIRE(cell->is_point_in_cell(point, &xi) == true);

        // Check point on vertex
        point << 0., 0., 0.;
        REQUIRE(cell->point_in_cartesian_cell(point) == true);
        REQUIRE(cell->is_point_in_cell(point, &xi) == true);
        REQUIRE(xi(0) == -1. + std::numeric_limits<double>::epsilon());
        REQUIRE(xi(1) == -1. + std::numeric_limits<double>::epsilon());
        REQUIRE(xi(2) == -1. + std::numeric_limits<double>::epsilon());

        // Check point on edge
        point << 0.5, 0., 0.;
        REQUIRE(cell->point_in_cartesian_cell(point) == true);
        REQUIRE(cell->is_point_in_cell(point, &xi) == true);
        REQUIRE(xi(0) == -0.5);
        REQUIRE(xi(1) == -1. + std::numeric_limits<double>::epsilon());
        REQUIRE(xi(2) == -1. + std::numeric_limits<double>::epsilon());

        // Check point on surface
        point << 0.5, 0.5, 0.;
        REQUIRE(cell->point_in_cartesian_cell(point) == true);
        REQUIRE(cell->is_point_in_cell(point, &xi) == true);
        REQUIRE(xi(0) == -0.5);
        REQUIRE(xi(1) == -0.5);
        REQUIRE(xi(2) == -1. + std::numeric_limits<double>::epsilon());

        // Check point outside
        point << 2.5, 2.5, 2.5;
        REQUIRE(cell->point_in_cartesian_cell(point) == false);
        REQUIRE(cell->is_point_in_cell(point, &xi) == false);
      }

      // Find local coordinates of a point in a cell
      SECTION("Find local coordinates of a point in cell") {
        // Coordinates of a point in real cell
        Eigen::Vector3d point;
        point << 0.5, 1.0, 1.5;

        // Coordinates of the point in an unit cell
        Eigen::Matrix<double, 3, 1> point_unit_cell;
        point_unit_cell << -0.5, 0., 0.5;

        // Get local coordinates of the point
        auto local_point = cell->local_coordinates_point(point);
        for (unsigned i = 0; i < local_point.size(); ++i)
          REQUIRE(local_point[i] ==
                  Approx(point_unit_cell[i]).epsilon(Tolerance));

        // Use Newton-raphson iteration
        local_point = cell->transform_real_to_unit_cell(point);
        for (unsigned i = 0; i < local_point.size(); ++i)
          REQUIRE(local_point[i] ==
                  Approx(point_unit_cell[i]).epsilon(Tolerance));
      }

      SECTION("Transform real to unit cell, xi is in positive coordinates") {
        // Number of nodes in cell
        const unsigned Nnodes = 8;

        // Coordinates
        Eigen::Vector3d coords;

        coords << 2, 1, -1;
        std::shared_ptr<mpm::NodeBase<Dim>> node0 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

        coords << 4, 2, -1.;
        std::shared_ptr<mpm::NodeBase<Dim>> node1 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);

        coords << 2, 4, -1.;
        std::shared_ptr<mpm::NodeBase<Dim>> node2 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);

        coords << 1, 3, -1.;
        std::shared_ptr<mpm::NodeBase<Dim>> node3 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);

        coords << 2, 1, 1.;
        std::shared_ptr<mpm::NodeBase<Dim>> node4 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(4, coords);

        coords << 4, 2, 1.;
        std::shared_ptr<mpm::NodeBase<Dim>> node5 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(5, coords);

        coords << 2, 4, 1.;
        std::shared_ptr<mpm::NodeBase<Dim>> node6 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(6, coords);

        coords << 1, 3, 1.;
        std::shared_ptr<mpm::NodeBase<Dim>> node7 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(7, coords);

        // 8-noded hexahedron shape functions
        std::shared_ptr<mpm::Element<Dim>> element =
            Factory<mpm::Element<Dim>>::instance()->create("ED3H8");

        mpm::Index id = 0;
        auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

        cell->add_node(0, node0);
        cell->add_node(1, node1);
        cell->add_node(2, node2);
        cell->add_node(3, node3);
        cell->add_node(4, node4);
        cell->add_node(5, node5);
        cell->add_node(6, node6);
        cell->add_node(7, node7);
        REQUIRE(cell->nnodes() == 8);

        // Initialise cell
        REQUIRE(cell->initialise() == true);

        // Coordinates of a point in real cell
        Eigen::Vector3d point;
        point << 2.1875, 3.25, 0.;

        // Coordinates of the point in an unit cell
        Eigen::Matrix<double, 3, 1> point_unit_cell;
        point_unit_cell << 0.5, 0.5, 0.;

        // Use Newton-raphson iteration to find local coordinates
        auto local_point = cell->transform_real_to_unit_cell(point);
        for (unsigned i = 0; i < local_point.size(); ++i)
          REQUIRE(local_point[i] ==
                  Approx(point_unit_cell[i]).epsilon(Tolerance));
      }

      SECTION("Transform real to unit cell, locate points in cells") {
        // Number of nodes in cell
        const unsigned Nnodes = 8;

        // Coordinates
        Eigen::Vector3d coords;

        coords << 0, 0., 0;
        std::shared_ptr<mpm::NodeBase<Dim>> node0 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

        coords << 2, 0., 0;
        std::shared_ptr<mpm::NodeBase<Dim>> node1 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);

        coords << 2, .5, 0;
        std::shared_ptr<mpm::NodeBase<Dim>> node2 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);

        coords << 0, .5, 0;
        std::shared_ptr<mpm::NodeBase<Dim>> node3 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);

        coords << 1., 0., 2.;
        std::shared_ptr<mpm::NodeBase<Dim>> node4 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(4, coords);

        coords << 3., 0., 2.;
        std::shared_ptr<mpm::NodeBase<Dim>> node5 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(5, coords);

        coords << 3., .5, 2.;
        std::shared_ptr<mpm::NodeBase<Dim>> node6 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(6, coords);

        coords << 1., .5, 2.;
        std::shared_ptr<mpm::NodeBase<Dim>> node7 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(7, coords);

        // Cell 2
        coords << 2, 1.5, 0;
        std::shared_ptr<mpm::NodeBase<Dim>> node8 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(8, coords);

        coords << 0, 1.5, 0;
        std::shared_ptr<mpm::NodeBase<Dim>> node9 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(9, coords);

        coords << 3., 1.5, 2.;
        std::shared_ptr<mpm::NodeBase<Dim>> node10 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(10, coords);

        coords << 1., 1.5, 2.;
        std::shared_ptr<mpm::NodeBase<Dim>> node11 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(11, coords);

        // 8-noded hexahedron shape functions
        std::shared_ptr<mpm::Element<Dim>> element =
            Factory<mpm::Element<Dim>>::instance()->create("ED3H8");

        // Cell 1
        mpm::Index id1 = 0;
        auto cell1 = std::make_shared<mpm::Cell<Dim>>(id1, Nnodes, element);

        cell1->add_node(0, node0);
        cell1->add_node(1, node1);
        cell1->add_node(2, node2);
        cell1->add_node(3, node3);
        cell1->add_node(4, node4);
        cell1->add_node(5, node5);
        cell1->add_node(6, node6);
        cell1->add_node(7, node7);
        REQUIRE(cell1->nnodes() == 8);

        // Initialise cell
        REQUIRE(cell1->initialise() == true);

        // Cell 2
        mpm::Index id2 = 0;
        auto cell2 = std::make_shared<mpm::Cell<Dim>>(id2, Nnodes, element);

        cell2->add_node(0, node3);
        cell2->add_node(1, node2);
        cell2->add_node(2, node8);
        cell2->add_node(3, node9);
        cell2->add_node(4, node7);
        cell2->add_node(5, node6);
        cell2->add_node(6, node10);
        cell2->add_node(7, node11);
        REQUIRE(cell2->nnodes() == 8);

        // Initialise cell
        REQUIRE(cell2->initialise() == true);

        Eigen::Vector3d xi;
        // Check point 1
        // Coordinates of a point in real cell
        Eigen::Vector3d point1;
        point1 << 2.25, 0.375, 1.5;

        // Coordinates of the point in an unit cell
        Eigen::Matrix<double, 3, 1> point_unit_cell1;
        point_unit_cell1 << 0.5, 0.5, 0.5;

        // Check if point is in cell
        REQUIRE(cell1->is_point_in_cell(point1, &xi) == true);
        REQUIRE(cell2->is_point_in_cell(point1, &xi) == false);

        // Use Newton-raphson iteration to find local coordinates
        auto local_point1 = cell1->transform_real_to_unit_cell(point1);
        for (unsigned i = 0; i < local_point1.size(); ++i)
          REQUIRE(local_point1[i] ==
                  Approx(point_unit_cell1[i]).epsilon(Tolerance));

        // Check point 2
        // Coordinates of a point in real cell
        Eigen::Vector3d point2;
        point2 << 0.75, 0.75, 0.5;

        // Coordinates of the point in an unit cell
        Eigen::Matrix<double, 3, 1> point_unit_cell2;
        point_unit_cell2 << -0.5, -0.5, -0.5;

        // Check if point is in cell
        REQUIRE(cell1->is_point_in_cell(point2, &xi) == false);
        REQUIRE(cell2->is_point_in_cell(point2, &xi) == true);

        // Use Newton-raphson iteration to find local coordinates
        auto local_point2 = cell2->transform_real_to_unit_cell(point2);
        for (unsigned i = 0; i < local_point2.size(); ++i)
          REQUIRE(local_point2[i] ==
                  Approx(point_unit_cell2[i]).epsilon(Tolerance));

        // Check point 3
        // Coordinates of a point in real cell
        Eigen::Vector3d point3;
        point3 << 0., 0., 0.;

        // Coordinates of the point in an unit cell
        Eigen::Matrix<double, 3, 1> point_unit_cell3;
        point_unit_cell3 << -1., -1., -1.;

        // Check if point is in cell
        REQUIRE(cell1->is_point_in_cell(point3, &xi) == true);
        REQUIRE(cell2->is_point_in_cell(point3, &xi) == false);

        // Use Newton-raphson iteration to find local coordinates
        auto local_point3 = cell1->transform_real_to_unit_cell(point3);
        for (unsigned i = 0; i < local_point3.size(); ++i)
          REQUIRE(local_point3[i] ==
                  Approx(point_unit_cell3[i]).epsilon(Tolerance));

        // Check point 4
        // Coordinates of a point in real cell
        Eigen::Vector3d point4;
        point4 << 3., 1.5, 2.;

        // Coordinates of the point in an unit cell
        Eigen::Matrix<double, 3, 1> point_unit_cell4;
        point_unit_cell4 << 1., 1., 1.;

        // Check if point is in cell
        REQUIRE(cell1->is_point_in_cell(point4, &xi) == false);
        REQUIRE(cell2->is_point_in_cell(point4, &xi) == true);

        // Use Newton-raphson iteration to find local coordinates
        auto local_point4 = cell2->transform_real_to_unit_cell(point4);
        for (unsigned i = 0; i < local_point4.size(); ++i)
          REQUIRE(local_point4[i] ==
                  Approx(point_unit_cell4[i]).epsilon(Tolerance));
      }

      SECTION("Generate points in cell") {
        // Number of nodes in cell
        const unsigned Nnodes = 8;

        // Coordinates
        Eigen::Vector3d coords;

        coords << 2, 1, -1;
        std::shared_ptr<mpm::NodeBase<Dim>> node0 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

        coords << 4, 2, -1.;
        std::shared_ptr<mpm::NodeBase<Dim>> node1 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);

        coords << 2, 4, -1.;
        std::shared_ptr<mpm::NodeBase<Dim>> node2 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);

        coords << 1, 3, -1.;
        std::shared_ptr<mpm::NodeBase<Dim>> node3 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);

        coords << 2, 1, 1.;
        std::shared_ptr<mpm::NodeBase<Dim>> node4 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(4, coords);

        coords << 4, 2, 1.;
        std::shared_ptr<mpm::NodeBase<Dim>> node5 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(5, coords);

        coords << 2, 4, 1.;
        std::shared_ptr<mpm::NodeBase<Dim>> node6 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(6, coords);

        coords << 1, 3, 1.;
        std::shared_ptr<mpm::NodeBase<Dim>> node7 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(7, coords);

        // 8-noded hexahedron shape functions
        std::shared_ptr<mpm::Element<Dim>> element =
            Factory<mpm::Element<Dim>>::instance()->create("ED3H8");

        mpm::Index id = 0;
        auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

        cell->add_node(0, node0);
        cell->add_node(1, node1);
        cell->add_node(2, node2);
        cell->add_node(3, node3);
        cell->add_node(4, node4);
        cell->add_node(5, node5);
        cell->add_node(6, node6);
        cell->add_node(7, node7);
        REQUIRE(cell->nnodes() == 8);

        // Initialise cell
        REQUIRE(cell->initialise() == true);

        // Assuming a quadrature of 1x1x1
        auto points = cell->generate_points();
        REQUIRE(points.size() == 1);

        auto quad = std::make_shared<mpm::HexahedronQuadrature<Dim, 1>>();
        auto quadrature = quad->quadratures();
        // Check if the output coordinates match local quadratures
        for (unsigned i = 0; i < points.size(); ++i) {
          auto local_point = cell->transform_real_to_unit_cell(points.at(i));
          for (unsigned j = 0; j < Dim; ++j)
            REQUIRE(local_point[j] ==
                    Approx(quadrature(j, i)).epsilon(Tolerance));
        }

        // Assign quadrature 2x2x2
        cell->assign_quadrature(2);

        points = cell->generate_points();
        REQUIRE(points.size() == 8);

        auto quad2 = std::make_shared<mpm::HexahedronQuadrature<Dim, 8>>();
        quadrature = quad2->quadratures();
        // Check if the output coordinates match local quadratures
        for (unsigned i = 0; i < points.size(); ++i) {
          auto local_point = cell->transform_real_to_unit_cell(points.at(i));
          for (unsigned j = 0; j < Dim; ++j)
            REQUIRE(local_point[j] ==
                    Approx(quadrature(j, i)).epsilon(Tolerance));
        }

        // Assign quadrature 3x3x3
        cell->assign_quadrature(3);

        points = cell->generate_points();
        REQUIRE(points.size() == 27);

        auto quad3 = std::make_shared<mpm::HexahedronQuadrature<Dim, 27>>();
        quadrature = quad3->quadratures();
        // Check if the output coordinates match local quadratures
        for (unsigned i = 0; i < points.size(); ++i) {
          auto local_point = cell->transform_real_to_unit_cell(points.at(i));
          for (unsigned j = 0; j < Dim; ++j)
            REQUIRE(local_point[j] ==
                    Approx(quadrature(j, i)).epsilon(Tolerance));
        }

        // Assign quadrature 4x4x4
        cell->assign_quadrature(4);

        points = cell->generate_points();
        REQUIRE(points.size() == 64);

        auto quad4 = std::make_shared<mpm::HexahedronQuadrature<Dim, 64>>();
        quadrature = quad4->quadratures();
        // Check if the output coordinates match local quadratures
        for (unsigned i = 0; i < points.size(); ++i) {
          auto local_point = cell->transform_real_to_unit_cell(points.at(i));
          for (unsigned j = 0; j < Dim; ++j)
            REQUIRE(local_point[j] ==
                    Approx(quadrature(j, i)).epsilon(Tolerance));
        }
      }
    }
  }

  SECTION("Check shape functions") {
    mpm::Index id = 0;
    // Check 8-noded function
    SECTION("Check 8-noded Hexahedron") {
      // 8-noded hexahedron shape functions
      std::shared_ptr<mpm::Element<Dim>> element =
          Factory<mpm::Element<Dim>>::instance()->create("ED3H8");

      auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);
      REQUIRE(cell->nfunctions() == 8);

      // Fail trying to create when number of nodes don't match
      const unsigned nnodes = 9;
      auto cell1 = std::make_shared<mpm::Cell<Dim>>(id, nnodes, element);
    }
    // Check 20-noded function
    SECTION("Check 20-noded Hexahedron") {
      // 20-noded hexahedron shape functions
      const unsigned nnodes = 20;
      std::shared_ptr<mpm::Element<Dim>> element =
          Factory<mpm::Element<Dim>>::instance()->create("ED3H20");

      auto cell = std::make_shared<mpm::Cell<Dim>>(id, nnodes, element);
      REQUIRE(cell->nfunctions() == 20);
    }
  }

  SECTION("Test particle addition deletion") {
    mpm::Index pid = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(0, Nnodes, element);
    REQUIRE(cell->status() == false);
    REQUIRE(cell->nparticles() == 0);
    REQUIRE(cell->add_particle_id(pid) == true);
    REQUIRE(cell->status() == true);
    REQUIRE(cell->nparticles() == 1);
    REQUIRE(cell->particles().size() == 1);
    REQUIRE(cell->particles().at(0) == pid);
    cell->remove_particle_id(pid);
    REQUIRE(cell->status() == false);
    REQUIRE(cell->nparticles() == 0);
    REQUIRE(cell->particles().size() == 0);
  }

  SECTION("Test node status") {
    mpm::Index pid = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(0, Nnodes, element);
    cell->add_node(0, node0);
    cell->add_node(1, node1);
    cell->add_node(2, node2);
    cell->add_node(3, node3);
    cell->add_node(4, node4);
    cell->add_node(5, node5);
    cell->add_node(6, node6);
    cell->add_node(7, node7);

    std::vector<std::shared_ptr<mpm::NodeBase<Dim>>> nodes;
    nodes.emplace_back(node0);
    nodes.emplace_back(node1);
    nodes.emplace_back(node2);
    nodes.emplace_back(node3);
    nodes.emplace_back(node4);
    nodes.emplace_back(node5);
    nodes.emplace_back(node6);
    nodes.emplace_back(node7);

    // When number of particles are zero
    REQUIRE(cell->nparticles() == 0);
    REQUIRE(cell->status() == false);
    cell->activate_nodes();
    for (const auto& node : nodes) REQUIRE(node->status() == false);

    // Add a particle
    REQUIRE(cell->add_particle_id(pid) == true);
    REQUIRE(cell->status() == true);
    REQUIRE(cell->nparticles() == 1);
    cell->activate_nodes();
    for (const auto& node : nodes) REQUIRE(node->status() == true);

    // Remove a particle
    cell->remove_particle_id(pid);
    REQUIRE(cell->status() == false);
    REQUIRE(cell->nparticles() == 0);
    // initialise nodes
    for (const auto& node : nodes) node->initialise();
    cell->activate_nodes();
    for (const auto& node : nodes) REQUIRE(node->status() == false);
  }

  // Check if a point is in an oblique cell
  SECTION("Point in a distorted cell") {
    // Coordinates
    Eigen::Vector3d coords;
    coords << 812481.0000000000, 815877.0000000001, 158.0900000000;
    std::shared_ptr<mpm::NodeBase<Dim>> node0 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

    coords << 812482.9999999999, 815877.0000000001, 159.1500000000;
    std::shared_ptr<mpm::NodeBase<Dim>> node1 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);

    coords << 812482.9999999999, 815879.0000000000, 160.3900000000;
    std::shared_ptr<mpm::NodeBase<Dim>> node2 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);

    coords << 812481.0000000000, 815879.0000000000, 159.4800000000;
    std::shared_ptr<mpm::NodeBase<Dim>> node3 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);

    coords << 812481.0000000000, 815877.0000000001, 160.0900000000;
    std::shared_ptr<mpm::NodeBase<Dim>> node4 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(4, coords);

    coords << 812482.9999999999, 815877.0000000001, 161.1500000000;
    std::shared_ptr<mpm::NodeBase<Dim>> node5 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(5, coords);

    coords << 812482.9999999999, 815879.0000000000, 162.3900000000;
    std::shared_ptr<mpm::NodeBase<Dim>> node6 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(6, coords);

    coords << 812481.0000000000, 815879.0000000000, 161.4800000000;
    std::shared_ptr<mpm::NodeBase<Dim>> node7 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(7, coords);

    // 8-noded hexahedron shape functions
    std::shared_ptr<mpm::Element<Dim>> element =
        Factory<mpm::Element<Dim>>::instance()->create("ED3H8");

    //! Check Cell IDs
    mpm::Index id = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

    // Check if cell is initialised, before addition of nodes
    REQUIRE(cell->is_initialised() == false);

    cell->add_node(0, node0);
    cell->add_node(1, node1);
    cell->add_node(2, node2);
    cell->add_node(3, node3);
    cell->add_node(4, node4);
    cell->add_node(5, node5);
    cell->add_node(6, node6);
    cell->add_node(7, node7);
    REQUIRE(cell->nnodes() == 8);

    REQUIRE(cell->nfunctions() == 8);

    REQUIRE(cell->initialise() == true);

    Eigen::Vector3d xi;
    // Check point in cell
    Eigen::Vector3d point;
    point << 812482.5000000000, 815878.5000000000, 160.0825000000;
    REQUIRE(cell->is_point_in_cell(point, &xi) == true);
  }

  // Check if a point is in an oblique cell
  SECTION("Point in a cell with Affine transformation and Newton Raphson") {

    // Check point in cell
    Eigen::Vector3d point;
    point << 812498.5000000001, 815872.5000000000, 167.459;

    // Coordinates
    Eigen::Vector3d coords;

    coords << 812496.9999999999, 815870.9999999999, 165.4800000000;
    std::shared_ptr<mpm::NodeBase<Dim>> node0 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

    coords << 812499.0000000000, 815870.9999999999, 166.7800000000;
    std::shared_ptr<mpm::NodeBase<Dim>> node1 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);

    coords << 812499.0000000000, 815873.0000000000, 168.0800000000;
    std::shared_ptr<mpm::NodeBase<Dim>> node2 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);

    coords << 812496.9999999999, 815873.0000000000, 166.8400000000;
    std::shared_ptr<mpm::NodeBase<Dim>> node3 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);

    coords << 812496.9999999999, 815870.9999999999, 167.4800000000;
    std::shared_ptr<mpm::NodeBase<Dim>> node4 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(4, coords);

    coords << 812499.0000000000, 815870.9999999999, 168.7800000000;
    std::shared_ptr<mpm::NodeBase<Dim>> node5 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(5, coords);

    coords << 812499.0000000000, 815873.0000000000, 170.0800000000;
    std::shared_ptr<mpm::NodeBase<Dim>> node6 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(6, coords);

    coords << 812496.9999999999, 815873.0000000000, 168.8400000000;
    std::shared_ptr<mpm::NodeBase<Dim>> node7 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(7, coords);

    // 8-noded hexahedron shape functions
    std::shared_ptr<mpm::Element<Dim>> element =
        Factory<mpm::Element<Dim>>::instance()->create("ED3H8");

    //! Check Cell IDs
    mpm::Index id = 0;
    auto cell1 = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

    // Check if cell is initialised, before addition of nodes
    REQUIRE(cell1->is_initialised() == false);

    cell1->add_node(0, node0);
    cell1->add_node(1, node1);
    cell1->add_node(2, node2);
    cell1->add_node(3, node3);
    cell1->add_node(4, node4);
    cell1->add_node(5, node5);
    cell1->add_node(6, node6);
    cell1->add_node(7, node7);
    REQUIRE(cell1->nnodes() == 8);

    REQUIRE(cell1->nfunctions() == 8);

    // Initialise cell
    REQUIRE(cell1->initialise() == true);

    Eigen::Vector3d xi;
    xi = cell1->transform_real_to_unit_cell(point);

    // Check using unit cell with affine transformation / Newton-Raphson
    REQUIRE(cell1->is_point_in_cell(point, &xi) == true);

    // Cell 2
    mpm::Index id2 = 0;
    auto cell2 = std::make_shared<mpm::Cell<Dim>>(id2, Nnodes, element);

    // Check if cell is initialised, before addition of nodes
    REQUIRE(cell2->is_initialised() == false);

    // Element 1
    coords << 812496.9999999999, 815870.9999999999, 163.4800000000;
    std::shared_ptr<mpm::NodeBase<Dim>> node10 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

    coords << 812499.0000000000, 815870.9999999999, 164.7800000000;
    std::shared_ptr<mpm::NodeBase<Dim>> node11 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);

    coords << 812499.0000000000, 815873.0000000000, 166.0800000000;
    std::shared_ptr<mpm::NodeBase<Dim>> node12 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);

    coords << 812496.9999999999, 815873.0000000000, 164.8400000000;
    std::shared_ptr<mpm::NodeBase<Dim>> node13 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);

    coords << 812496.9999999999, 815870.9999999999, 165.4800000000;
    std::shared_ptr<mpm::NodeBase<Dim>> node14 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(4, coords);

    coords << 812499.0000000000, 815870.9999999999, 166.7800000000;
    std::shared_ptr<mpm::NodeBase<Dim>> node15 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(5, coords);

    coords << 812499.0000000000, 815873.0000000000, 168.0800000000;
    std::shared_ptr<mpm::NodeBase<Dim>> node16 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(6, coords);

    coords << 812496.9999999999, 815873.0000000000, 166.8400000000;
    std::shared_ptr<mpm::NodeBase<Dim>> node17 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(7, coords);

    cell2->add_node(0, node10);
    cell2->add_node(1, node11);
    cell2->add_node(2, node12);
    cell2->add_node(3, node13);
    cell2->add_node(4, node14);
    cell2->add_node(5, node15);
    cell2->add_node(6, node16);
    cell2->add_node(7, node17);
    REQUIRE(cell2->nnodes() == 8);
    REQUIRE(cell2->nfunctions() == 8);

    // Initialise cell
    REQUIRE(cell2->initialise() == true);

    // Point in cell fails to capture
    REQUIRE(cell2->point_in_cartesian_cell(point) == false);

    // Use affine transformation
    xi = cell2->transform_real_to_unit_cell(point);

    // Check using unit cell with affine transformation / Newton-Raphson
    REQUIRE(cell2->is_point_in_cell(point, &xi) == false);
  }

  SECTION("Test nglobal particles") {
    mpm::Index pid = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(0, Nnodes, element);
    REQUIRE(cell->nglobal_particles() == 0);
    cell->nglobal_particles(5);
    REQUIRE(cell->nglobal_particles() == 5);
  }
}

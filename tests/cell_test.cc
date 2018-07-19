#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

#include "cell.h"
#include "factory.h"
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

  // Shape function
  // 4-noded quadrilateral shape functions
  std::shared_ptr<mpm::ShapeFn<Dim>> shapefn =
      Factory<mpm::ShapeFn<Dim>>::instance()->create("SFQ4");

  Eigen::Vector2d coords;
  coords.setZero();

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

  //! Check Cell IDs
  SECTION("Check cell ids") {
    //! Check for id = 0
    SECTION("Cell id is zero") {
      mpm::Index id = 0;
      auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, shapefn);
      REQUIRE(cell->id() == 0);
    }

    SECTION("Cell id is positive") {
      //! Check for id is a positive value
      mpm::Index id = std::numeric_limits<mpm::Index>::max();
      auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, shapefn);
      REQUIRE(cell->id() == std::numeric_limits<mpm::Index>::max());
    }
  }

  SECTION("Add nodes") {
    mpm::Index id = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, shapefn);
    // Check if cell is initialised, before addition of nodes
    REQUIRE(cell->is_initialised() == false);

    cell->add_node(0, node0);
    cell->add_node(1, node1);
    cell->add_node(2, node2);
    cell->add_node(3, node3);
    REQUIRE(cell->nnodes() == 4);

    // Check if cell is initialised, after addition of nodes
    REQUIRE(cell->is_initialised() == false);

    // Check cell length calculation
    SECTION("Compute mean length of cell") {
      // Check cell length at initialisation
      REQUIRE(cell->mean_length() == std::numeric_limits<double>::max());

      cell->compute_mean_length();

      // Length of the cell
      const double length = 2.0;

      // Check cell length calculations
      REQUIRE(cell->mean_length() == Approx(length).epsilon(Tolerance));
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

      // Check if cell is initialised, after addition of shapefn
      REQUIRE(cell->is_initialised() == false);

      // Compute volume
      REQUIRE(cell->volume() ==
              Approx(std::numeric_limits<double>::max()).epsilon(Tolerance));
      cell->compute_volume();
      REQUIRE(cell->volume() == Approx(4.0).epsilon(Tolerance));

      // Check if cell is initialised, after volume calculation
      REQUIRE(cell->is_initialised() == true);

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

      SECTION("Transform real to unit cell, xi is (0.5, 0.5)") {
        // Number of nodes in cell
        const unsigned Nnodes = 4;

        // Coordinates
        Eigen::Vector2d coords;

        coords << 2, 1;
        std::shared_ptr<mpm::NodeBase<Dim>> node0 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

        coords << 4, 2;
        std::shared_ptr<mpm::NodeBase<Dim>> node1 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);

        coords << 2, 4;
        std::shared_ptr<mpm::NodeBase<Dim>> node2 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);

        coords << 1, 3;
        std::shared_ptr<mpm::NodeBase<Dim>> node3 =
            std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);

        // 4-noded quadrilateral shape functions
        std::shared_ptr<mpm::ShapeFn<Dim>> shapefn =
            Factory<mpm::ShapeFn<Dim>>::instance()->create("SFQ4");

        mpm::Index id = 0;
        auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, shapefn);

        cell->add_node(0, node0);
        cell->add_node(1, node1);
        cell->add_node(2, node2);
        cell->add_node(3, node3);
        REQUIRE(cell->nnodes() == 4);

        // Coordinates of a point in real cell
        Eigen::Vector2d point;
        point << 2.1875, 3.25;

        // Coordinates of the point in an unit cell
        Eigen::Matrix<double, 2, 1> point_unit_cell;
        point_unit_cell << 0.5, 0.5;

        // Use Newton-raphson iteration to find local coordinates
        auto local_point = cell->transform_real_to_unit_cell(point);
        for (unsigned i = 0; i < local_point.size(); ++i)
          REQUIRE(local_point[i] ==
                  Approx(point_unit_cell[i]).epsilon(Tolerance));
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
        std::shared_ptr<mpm::ShapeFn<Dim>> shapefn =
            Factory<mpm::ShapeFn<Dim>>::instance()->create("SFQ4");

        mpm::Index id = 0;
        auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, shapefn);

        cell->add_node(0, node0);
        cell->add_node(1, node1);
        cell->add_node(2, node2);
        cell->add_node(3, node3);
        REQUIRE(cell->nnodes() == 4);

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
    }
  }

  SECTION("Add neighbours") {
    auto cell = std::make_shared<mpm::Cell<Dim>>(0, Nnodes, shapefn);
    auto neighbourcell = std::make_shared<mpm::Cell<Dim>>(1, Nnodes, shapefn);
    REQUIRE(cell->nneighbours() == 0);
    cell->add_neighbour(0, neighbourcell);
    REQUIRE(cell->nneighbours() == 1);
  }

  SECTION("Check shape functions") {
    mpm::Index id = 0;
    // Check 4-noded function
    SECTION("Check 4-noded Quadrilateral") {
      // 4-noded quadrilateral shape functions
      std::shared_ptr<mpm::ShapeFn<Dim>> shapefn =
          Factory<mpm::ShapeFn<Dim>>::instance()->create("SFQ4");
      auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, shapefn);
      REQUIRE(cell->nfunctions() == 4);
    }
    // Check 8-noded function
    SECTION("Check 8-noded Quadrilateral") {
      // 8-noded quadrilateral shape functions
      std::shared_ptr<mpm::ShapeFn<Dim>> shapefn =
          Factory<mpm::ShapeFn<Dim>>::instance()->create("SFQ8");
      auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, shapefn);
      REQUIRE(cell->nfunctions() == 8);
    }
    // Check 9-noded function
    SECTION("Check 9-noded Quadrilateral") {
      // 9-noded quadrilateral shape functions
      std::shared_ptr<mpm::ShapeFn<Dim>> shapefn =
          Factory<mpm::ShapeFn<Dim>>::instance()->create("SFQ9");
      auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, shapefn);
      REQUIRE(cell->nfunctions() == 9);
    }
  }

  SECTION("Test particle addition deletion") {
    mpm::Index pid = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(0, Nnodes, shapefn);
    REQUIRE(cell->status() == false);
    bool status = cell->add_particle_id(pid);
    REQUIRE(status == true);
    REQUIRE(cell->status() == true);
    cell->remove_particle_id(pid);
    REQUIRE(cell->status() == false);
  }

  SECTION("Test particle information mapping") {
    mpm::Index id = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, shapefn);
    cell->add_node(0, node0);
    cell->add_node(1, node1);
    cell->add_node(2, node2);
    cell->add_node(3, node3);

    // Create a vector of node pointers
    std::vector<std::shared_ptr<mpm::NodeBase<Dim>>> nodes{node0, node1, node2,
                                                           node3};

    // 4-noded quadrilateral shape functions
    std::shared_ptr<mpm::ShapeFn<Dim>> shapefn =
        Factory<mpm::ShapeFn<Dim>>::instance()->create("SFQ4");

    // Local coordinate of a particle
    Eigen::Vector2d xi = Eigen::Vector2d::Zero();
    // Particle mass
    double pmass = 4.;
    // Particle volume
    double pvolume = 8.;
    // Particle velocity
    Eigen::Vector2d pvelocity;
    pvelocity << 1., 1.;
    // Particle gravity
    Eigen::Vector2d pgravity;
    pgravity << 0., 9.814;
    // Phase
    unsigned phase = 0;

    SECTION("Check particle mass mapping") {
      cell->map_particle_mass_to_nodes(xi, phase, pmass);
      for (const auto& node : nodes)
        REQUIRE(node->mass(phase) == Approx(1.0).epsilon(Tolerance));
    }

    SECTION("Check particle volume mapping") {
      cell->map_particle_volume_to_nodes(xi, phase, pvolume);
      for (const auto& node : nodes)
        REQUIRE(node->volume(phase) == Approx(2.0).epsilon(Tolerance));
    }

    SECTION("Check particle momentum mapping") {
      cell->map_particle_mass_to_nodes(xi, phase, pmass);
      cell->compute_nodal_momentum(xi, phase, pmass, pvelocity);
      for (const auto& node : nodes) {
        for (unsigned i = 0; i < pvelocity.size(); ++i)
          REQUIRE(node->momentum(phase)(i) == Approx(1.0).epsilon(Tolerance));
      }
    }

    SECTION("Check particle body force mapping") {
      // Calculate body force at nodes
      cell->compute_nodal_body_force(xi, phase, pmass, pgravity);
      Eigen::Vector2d bodyforce;
      bodyforce << 0., 9.814;
      for (const auto& node : nodes) {
        for (unsigned i = 0; i < bodyforce.size(); ++i)
          REQUIRE(node->external_force(phase)(i) ==
                  Approx(bodyforce(i)).epsilon(Tolerance));
      }
    }

    SECTION("Check particle internal force mapping") {
      // Assign internal force to nodes
      const double pvolume = 0.5;
      Eigen::Matrix<double, 6, 1> pinternal_stress;
      pinternal_stress << 0.5, 0.5, 0.5, 0.5, 0.5, 0.5;

      cell->compute_nodal_internal_force(phase, pvolume, xi, pinternal_stress);

      // Check internal force
      std::vector<Eigen::Vector2d> internal_forces;
      Eigen::Vector2d intforce;
      // Node 0
      intforce << -0.125, -0.125;
      internal_forces.push_back(intforce);
      // Node 1
      intforce << 0., 0.;
      internal_forces.push_back(intforce);
      // Node 2
      intforce << 0.125, 0.125;
      internal_forces.push_back(intforce);
      // Node 3
      intforce << 0, 0;
      internal_forces.push_back(intforce);

      unsigned j = 0;
      for (const auto& node : nodes) {
        for (unsigned i = 0; i < intforce.size(); ++i)
          REQUIRE(node->internal_force(phase)(i) ==
                  Approx(internal_forces.at(j)(i)).epsilon(Tolerance));
        ++j;
      }
    }

    SECTION("Check interpolate velocity") {
      // Assign mass to 100
      const double mass = 100.;

      // Apply momentum
      Eigen::Matrix<double, Dim, 1> momentum;
      unsigned j = 1;
      for (const auto& node : nodes) {
        // Apply momentum
        for (unsigned i = 0; i < momentum.size(); ++i)
          momentum(i) = 10. * static_cast<double>(j);

        // Nodal mass
        node->update_mass(false, phase, mass);
        REQUIRE(node->mass(phase) == Approx(100.0).epsilon(Tolerance));

        // Nodal momentum
        node->update_momentum(false, phase, momentum);
        for (unsigned i = 0; i < momentum.size(); ++i)
          REQUIRE(node->momentum(phase)(i) ==
                  Approx(10. * static_cast<double>(j)).epsilon(Tolerance));

        for (unsigned i = 0; i < momentum.size(); ++i)
          REQUIRE(node->momentum(phase)(i) ==
                  Approx(10. * static_cast<double>(j)).epsilon(Tolerance));

        // Compute and check velocity
        node->compute_velocity();
        for (unsigned i = 0; i < Dim; ++i)
          REQUIRE(node->velocity(phase)(i) ==
                  Approx(0.1 * static_cast<double>(j)).epsilon(Tolerance));
        // Increment j
        ++j;
      }
      // Check interpolate velocity (0, 0)
      Eigen::Vector2d velocity = cell->interpolate_nodal_velocity(xi, phase);

      Eigen::Vector2d interpolated_velocity;
      interpolated_velocity << 0.25, 0.25;
      for (unsigned i = 0; i < velocity.size(); ++i)
        REQUIRE(velocity(i) ==
                Approx(interpolated_velocity(i)).epsilon(Tolerance));

      // Check interpolate velocity (0.5, 0.5)
      xi << 0.5, 0.5;
      velocity = cell->interpolate_nodal_velocity(xi, phase);

      interpolated_velocity << 0.2875, 0.2875;
      for (unsigned i = 0; i < velocity.size(); ++i)
        REQUIRE(velocity(i) ==
                Approx(interpolated_velocity(i)).epsilon(Tolerance));

      // Check interpolate velocity (-0.5, -0.5)
      xi << -0.5, -0.5;
      velocity = cell->interpolate_nodal_velocity(xi, phase);

      interpolated_velocity << 0.1875, 0.1875;
      for (unsigned i = 0; i < velocity.size(); ++i)
        REQUIRE(velocity(i) ==
                Approx(interpolated_velocity(i)).epsilon(Tolerance));
    }

    SECTION("Check interpolate acceleration") {

      // Apply acceleration
      Eigen::Matrix<double, Dim, 1> acceleration;
      unsigned j = 1;
      for (const auto& node : nodes) {
        // Apply acceleration
        for (unsigned i = 0; i < acceleration.size(); ++i)
          acceleration(i) = 10. * static_cast<double>(j);

        // Nodal acceleration
        node->update_acceleration(false, phase, acceleration);
        for (unsigned i = 0; i < acceleration.size(); ++i)
          REQUIRE(node->acceleration(phase)(i) ==
                  Approx(10. * static_cast<double>(j)).epsilon(Tolerance));
        // Increment j
        ++j;
      }
      // Check interpolate acceleration (0, 0)
      Eigen::Vector2d check_acceleration =
          cell->interpolate_nodal_acceleration(xi, phase);

      Eigen::Vector2d interpolated_acceleration;
      interpolated_acceleration << 25., 25.;
      for (unsigned i = 0; i < check_acceleration.size(); ++i)
        REQUIRE(check_acceleration(i) ==
                Approx(interpolated_acceleration(i)).epsilon(Tolerance));

      // Check interpolate acceleration (0.5, 0.5)
      xi << 0.5, 0.5;
      check_acceleration = cell->interpolate_nodal_acceleration(xi, phase);

      interpolated_acceleration << 28.75, 28.75;
      for (unsigned i = 0; i < check_acceleration.size(); ++i)
        REQUIRE(check_acceleration(i) ==
                Approx(interpolated_acceleration(i)).epsilon(Tolerance));

      // Check interpolate acceleration (-0.5, -0.5)
      xi << -0.5, -0.5;
      check_acceleration = cell->interpolate_nodal_acceleration(xi, phase);

      interpolated_acceleration << 18.75, 18.75;
      for (unsigned i = 0; i < check_acceleration.size(); ++i)
        REQUIRE(check_acceleration(i) ==
                Approx(interpolated_acceleration(i)).epsilon(Tolerance));
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
  std::shared_ptr<mpm::ShapeFn<Dim>> shapefn =
      Factory<mpm::ShapeFn<Dim>>::instance()->create("SFH8");

  //! Check Cell IDs
  SECTION("Check cell ids") {
    //! Check for id = 0
    SECTION("Cell id is zero") {
      mpm::Index id = 0;
      auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, shapefn);
      REQUIRE(cell->id() == 0);
    }

    SECTION("Cell id is positive") {
      //! Check for id is a positive value
      mpm::Index id = std::numeric_limits<mpm::Index>::max();
      auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, shapefn);
      REQUIRE(cell->id() == std::numeric_limits<mpm::Index>::max());
    }
  }

  // Check node additions
  SECTION("Add nodes") {
    mpm::Index id = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, shapefn);

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

    // Check if cell is initialised, after addition of nodes
    REQUIRE(cell->is_initialised() == false);

    // Check cell length calculation
    SECTION("Compute mean length of cell") {
      // Check cell length at initialisation
      REQUIRE(cell->mean_length() == std::numeric_limits<double>::max());

      cell->compute_mean_length();

      // Length of the cell
      const double length = 2.0;

      // Check cell length calculations
      REQUIRE(cell->mean_length() == Approx(length).epsilon(Tolerance));
    }

    // Check centroid calculation
    SECTION("Compute centroid of a cell") {
      REQUIRE(cell->nfunctions() == 8);

      // Compute centroid
      cell->compute_centroid();
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

    // Check cell volume calculation
    SECTION("Compute volume of a cell") {
      REQUIRE(cell->nfunctions() == 8);

      // Check if cell is initialised, after addition of shapefn
      REQUIRE(cell->is_initialised() == false);

      // Compute volume
      REQUIRE(cell->volume() ==
              Approx(std::numeric_limits<double>::max()).epsilon(Tolerance));
      cell->compute_volume();
      REQUIRE(cell->volume() == Approx(8.0).epsilon(Tolerance));

      // Check if cell is initialised, after volume computation
      REQUIRE(cell->is_initialised() == true);

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
        point << 2.5, 2.5, 2.5;
        REQUIRE(cell->point_in_cell(point) == false);
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

      SECTION("Transform real to unit cell, xi is (0.5, 0.5, 0.)") {
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
        std::shared_ptr<mpm::ShapeFn<Dim>> shapefn =
            Factory<mpm::ShapeFn<Dim>>::instance()->create("SFH8");

        mpm::Index id = 0;
        auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, shapefn);

        cell->add_node(0, node0);
        cell->add_node(1, node1);
        cell->add_node(2, node2);
        cell->add_node(3, node3);
        cell->add_node(4, node4);
        cell->add_node(5, node5);
        cell->add_node(6, node6);
        cell->add_node(7, node7);
        REQUIRE(cell->nnodes() == 8);

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
        std::shared_ptr<mpm::ShapeFn<Dim>> shapefn =
            Factory<mpm::ShapeFn<Dim>>::instance()->create("SFH8");

        // Cell 1
        mpm::Index id1 = 0;
        auto cell1 = std::make_shared<mpm::Cell<Dim>>(id1, Nnodes, shapefn);

        cell1->add_node(0, node0);
        cell1->add_node(1, node1);
        cell1->add_node(2, node2);
        cell1->add_node(3, node3);
        cell1->add_node(4, node4);
        cell1->add_node(5, node5);
        cell1->add_node(6, node6);
        cell1->add_node(7, node7);
        REQUIRE(cell1->nnodes() == 8);

        // Cell 2
        mpm::Index id2 = 0;
        auto cell2 = std::make_shared<mpm::Cell<Dim>>(id2, Nnodes, shapefn);

        cell2->add_node(0, node3);
        cell2->add_node(1, node2);
        cell2->add_node(2, node8);
        cell2->add_node(3, node9);
        cell2->add_node(4, node7);
        cell2->add_node(5, node6);
        cell2->add_node(6, node10);
        cell2->add_node(7, node11);
        REQUIRE(cell2->nnodes() == 8);

        // Check point 1
        // Coordinates of a point in real cell
        Eigen::Vector3d point1;
        point1 << 2.25, 0.375, 1.5;

        // Coordinates of the point in an unit cell
        Eigen::Matrix<double, 3, 1> point_unit_cell1;
        point_unit_cell1 << 0.5, 0.5, 0.5;

        // Check if point is in cell
        REQUIRE(cell1->point_in_cell(point1) == true);
        REQUIRE(cell2->point_in_cell(point1) == false);

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
        REQUIRE(cell1->point_in_cell(point2) == false);
        REQUIRE(cell2->point_in_cell(point2) == true);

        // Use Newton-raphson iteration to find local coordinates
        auto local_point2 = cell2->transform_real_to_unit_cell(point2);
        for (unsigned i = 0; i < local_point2.size(); ++i)
          REQUIRE(local_point2[i] ==
                  Approx(point_unit_cell2[i]).epsilon(Tolerance));
      }
    }
  }

  SECTION("Add neighbours") {
    auto cell = std::make_shared<mpm::Cell<Dim>>(0, Nnodes, shapefn);
    auto neighbourcell = std::make_shared<mpm::Cell<Dim>>(1, Nnodes, shapefn);
    REQUIRE(cell->nneighbours() == 0);
    cell->add_neighbour(0, neighbourcell);
    REQUIRE(cell->nneighbours() == 1);
  }

  SECTION("Check shape functions") {
    mpm::Index id = 0;
    // Check 8-noded function
    SECTION("Check 8-noded Hexahedron") {
      // 8-noded hexahedron shape functions
      std::shared_ptr<mpm::ShapeFn<Dim>> shapefn =
          Factory<mpm::ShapeFn<Dim>>::instance()->create("SFH8");

      auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, shapefn);
      REQUIRE(cell->nfunctions() == 8);
    }
    // Check 20-noded function
    SECTION("Check 20-noded Hexahedron") {
      // 20-noded hexahedron shape functions
      std::shared_ptr<mpm::ShapeFn<Dim>> shapefn =
          Factory<mpm::ShapeFn<Dim>>::instance()->create("SFH20");

      auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, shapefn);
      REQUIRE(cell->nfunctions() == 20);
    }
  }

  SECTION("Test particle addition deletion") {
    mpm::Index pid = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(0, Nnodes, shapefn);
    REQUIRE(cell->status() == false);
    bool status = cell->add_particle_id(pid);
    REQUIRE(status == true);
    REQUIRE(cell->status() == true);
    cell->remove_particle_id(pid);
    REQUIRE(cell->status() == false);
  }

  SECTION("Test particle information mapping") {
    mpm::Index id = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, shapefn);
    cell->add_node(0, node0);
    cell->add_node(1, node1);
    cell->add_node(2, node2);
    cell->add_node(3, node3);
    cell->add_node(4, node4);
    cell->add_node(5, node5);
    cell->add_node(6, node6);
    cell->add_node(7, node7);
    REQUIRE(cell->nnodes() == 8);

    // Create a vector of node pointers
    std::vector<std::shared_ptr<mpm::NodeBase<Dim>>> nodes{
        node0, node1, node2, node3, node4, node5, node6, node7};

    // Local coordinate of a particle
    Eigen::Vector3d xi = Eigen::Vector3d::Zero();

    // Particle mass
    double pmass = 4.;
    // Particle volume
    double pvolume = 8.;
    // Particle velocity
    Eigen::Vector3d pvelocity;
    pvelocity << 1., 1., 1.;
    // Particle gravity
    Eigen::Vector3d pgravity;
    pgravity << 0., 0., 9.814;
    // Phase
    unsigned phase = 0;

    SECTION("Check particle mass mapping") {
      cell->map_particle_mass_to_nodes(xi, phase, pmass);
      for (const auto& node : nodes)
        REQUIRE(node->mass(phase) == Approx(0.5).epsilon(Tolerance));
    }

    SECTION("Check particle volume mapping") {
      cell->map_particle_volume_to_nodes(xi, phase, pvolume);
      REQUIRE(nodes.size() == 8);
      for (const auto& node : nodes)
        REQUIRE(node->volume(phase) == Approx(1.0).epsilon(Tolerance));
    }

    SECTION("Check particle momentum mapping") {
      cell->map_particle_mass_to_nodes(xi, phase, pmass);
      cell->compute_nodal_momentum(xi, phase, pmass, pvelocity);
      for (const auto& node : nodes) {
        for (unsigned i = 0; i < pvelocity.size(); ++i)
          REQUIRE(node->momentum(phase)(i) == Approx(0.5).epsilon(Tolerance));
      }
    }

    SECTION("Check particle body force mapping") {
      // Compute body force at nodes
      cell->compute_nodal_body_force(xi, phase, pmass, pgravity);
      Eigen::Vector3d bodyforce;
      bodyforce << 0., 0., 0.5 * 9.814;
      for (const auto& node : nodes) {
        for (unsigned i = 0; i < bodyforce.size(); ++i)
          REQUIRE(node->external_force(phase)(i) ==
                  Approx(bodyforce(i)).epsilon(Tolerance));
      }
    }

    SECTION("Check particle internal force mapping") {
      // Assign internal force to nodes
      const double pvolume = 0.5;
      Eigen::Matrix<double, 6, 1> pinternal_stress;
      pinternal_stress << 0.5, 0.5, 0.5, 0.5, 0.5, 0.5;

      cell->compute_nodal_internal_force(phase, pvolume, xi, pinternal_stress);

      // Check internal force
      std::vector<Eigen::Vector3d> internal_forces;
      Eigen::Vector3d intforce;
      // Node 0
      intforce << -0.09375, -0.09375, -0.09375;
      internal_forces.push_back(intforce);
      // Node 1
      intforce << -0.03125, -0.03125, -0.03125;
      internal_forces.push_back(intforce);
      // Node 2
      intforce << 0.03125, 0.03125, 0.03125;
      internal_forces.push_back(intforce);
      // Node 3
      intforce << -0.03125, -0.03125, -0.03125;
      internal_forces.push_back(intforce);
      // Node 4
      intforce << -0.03125, -0.03125, -0.03125;
      internal_forces.push_back(intforce);
      // Node 5
      intforce << 0.03125, 0.03125, 0.03125;
      internal_forces.push_back(intforce);
      // Node 6
      intforce << 0.09375, 0.09375, 0.09375;
      internal_forces.push_back(intforce);
      // Node 7
      intforce << 0.03125, 0.03125, 0.03125;
      internal_forces.push_back(intforce);

      unsigned j = 0;
      for (const auto& node : nodes) {
        for (unsigned i = 0; i < intforce.size(); ++i)
          REQUIRE(node->internal_force(phase)(i) ==
                  Approx(internal_forces.at(j)(i)).epsilon(Tolerance));
        ++j;
      }
    }

    SECTION("Check interpolate velocity") {
      // Assign mass to 100
      const double mass = 100.;

      // Apply momentum
      Eigen::Matrix<double, Dim, 1> momentum;
      unsigned j = 1;
      for (const auto& node : nodes) {
        // Apply momentum
        for (unsigned i = 0; i < momentum.size(); ++i)
          momentum(i) = 10. * static_cast<double>(j);

        // Nodal mass
        node->update_mass(false, phase, mass);
        REQUIRE(node->mass(phase) == Approx(100.0).epsilon(Tolerance));

        // Nodal momentum
        node->update_momentum(false, phase, momentum);
        for (unsigned i = 0; i < momentum.size(); ++i)
          REQUIRE(node->momentum(phase)(i) ==
                  Approx(10. * static_cast<double>(j)).epsilon(Tolerance));

        for (unsigned i = 0; i < momentum.size(); ++i)
          REQUIRE(node->momentum(phase)(i) ==
                  Approx(10. * static_cast<double>(j)).epsilon(Tolerance));

        // Compute and check velocity
        node->compute_velocity();
        for (unsigned i = 0; i < Dim; ++i)
          REQUIRE(node->velocity(phase)(i) ==
                  Approx(0.1 * static_cast<double>(j)).epsilon(Tolerance));
        // Increment j
        ++j;
      }
      // Check interpolate velocity (0, 0)
      xi.setZero();
      Eigen::Vector3d velocity = cell->interpolate_nodal_velocity(xi, phase);

      Eigen::Vector3d interpolated_velocity;
      interpolated_velocity << 0.45, 0.45, 0.45;
      for (unsigned i = 0; i < velocity.size(); ++i)
        REQUIRE(velocity(i) ==
                Approx(interpolated_velocity(i)).epsilon(Tolerance));

      // Check interpolate velocity (0.5, 0.5)
      xi << 0.5, 0.5, 0.5;
      velocity = cell->interpolate_nodal_velocity(xi, phase);

      interpolated_velocity << 0.5875, 0.5875, 0.5875;
      for (unsigned i = 0; i < velocity.size(); ++i)
        REQUIRE(velocity(i) ==
                Approx(interpolated_velocity(i)).epsilon(Tolerance));

      // Check interpolate velocity (-0.5, -0.5)
      xi << -0.5, -0.5, -0.5;
      velocity = cell->interpolate_nodal_velocity(xi, phase);

      interpolated_velocity << 0.2875, 0.2875, 0.2875;
      for (unsigned i = 0; i < velocity.size(); ++i)
        REQUIRE(velocity(i) ==
                Approx(interpolated_velocity(i)).epsilon(Tolerance));
    }

    SECTION("Check interpolate acceleration") {
      // Apply acceleration
      Eigen::Matrix<double, Dim, 1> acceleration;
      unsigned j = 1;
      for (const auto& node : nodes) {
        // Apply acceleration
        for (unsigned i = 0; i < acceleration.size(); ++i)
          acceleration(i) = 10. * static_cast<double>(j);

        // Nodal acceleration
        node->update_acceleration(false, phase, acceleration);
        for (unsigned i = 0; i < acceleration.size(); ++i)
          REQUIRE(node->acceleration(phase)(i) ==
                  Approx(10. * static_cast<double>(j)).epsilon(Tolerance));

        // Increment j
        ++j;
      }
      // Set coordinates as zero
      xi.setZero();
      // Check interpolate acceleration (0, 0, 0)
      Eigen::Vector3d check_acceleration =
          cell->interpolate_nodal_acceleration(xi, phase);

      Eigen::Vector3d interpolated_acceleration;
      interpolated_acceleration << 45.0, 45.0, 45.0;
      for (unsigned i = 0; i < check_acceleration.size(); ++i)
        REQUIRE(check_acceleration(i) ==
                Approx(interpolated_acceleration(i)).epsilon(Tolerance));

      // Check interpolate acceleration (0.5, 0.5, 0.5)
      xi << 0.5, 0.5, 0.5;
      check_acceleration = cell->interpolate_nodal_acceleration(xi, phase);

      interpolated_acceleration << 58.75, 58.75, 58.75;
      for (unsigned i = 0; i < check_acceleration.size(); ++i)
        REQUIRE(check_acceleration(i) ==
                Approx(interpolated_acceleration(i)).epsilon(Tolerance));

      // Check interpolate acceleration (-0.5, -0.5, -0.5)
      xi << -0.5, -0.5, -0.5;
      check_acceleration = cell->interpolate_nodal_acceleration(xi, phase);

      interpolated_acceleration << 28.75, 28.75, 28.75;
      for (unsigned i = 0; i < check_acceleration.size(); ++i)
        REQUIRE(check_acceleration(i) ==
                Approx(interpolated_acceleration(i)).epsilon(Tolerance));
    }
  }
}

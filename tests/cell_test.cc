#include <cmath>
#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

#include "cell.h"
#include "element.h"
#include "factory.h"
#include "hexahedron_element.h"
#include "node.h"
#include "quadrilateral_element.h"

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

  SECTION("Add nodes") {
    mpm::Index id = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);
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
            Approx(std::numeric_limits<double>::max()).epsilon(Tolerance));

    // Initialise cell
    REQUIRE(cell->initialise() == true);
    // Check if cell is initialised, after addition of nodes
    REQUIRE(cell->is_initialised() == true);

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

        // Test if point is in cell
        REQUIRE(cell->point_in_cell(point) == true);

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
        REQUIRE(cell->point_in_cell(point1) == true);

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
    }
  }

  SECTION("Add neighbours") {
    auto cell = std::make_shared<mpm::Cell<Dim>>(0, Nnodes, element);
    auto neighbourcell = std::make_shared<mpm::Cell<Dim>>(1, Nnodes, element);
    REQUIRE(cell->nneighbours() == 0);
    cell->add_neighbour(0, neighbourcell);
    REQUIRE(cell->nneighbours() == 1);
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
    REQUIRE(cell->nparticles() == 0);
    REQUIRE(cell->status() == false);
    bool status = cell->add_particle_id(pid);
    REQUIRE(status == true);
    REQUIRE(cell->status() == true);
    REQUIRE(cell->nparticles() == 1);
    cell->remove_particle_id(pid);
    REQUIRE(cell->status() == false);
    REQUIRE(cell->nparticles() == 0);
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

  SECTION("Test particle information mapping") {
    mpm::Index id = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);
    cell->add_node(0, node0);
    cell->add_node(1, node1);
    cell->add_node(2, node2);
    cell->add_node(3, node3);

    // Initialise cell
    REQUIRE(cell->initialise() == true);

    // Create a vector of node pointers
    std::vector<std::shared_ptr<mpm::NodeBase<Dim>>> nodes{node0, node1, node2,
                                                           node3};

    // 4-noded quadrilateral shape functions
    std::shared_ptr<mpm::Element<Dim>> element =
        Factory<mpm::Element<Dim>>::instance()->create("ED2Q4");

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
    const auto shapefns_xi = element->shapefn(xi);
    const auto bmatrix = element->bmatrix(xi);

    SECTION("Check particle mass mapping") {
      cell->map_particle_mass_to_nodes(shapefns_xi, phase, pmass);
      for (const auto& node : nodes)
        REQUIRE(node->mass(phase) == Approx(1.0).epsilon(Tolerance));
    }

    SECTION("Check particle volume mapping") {
      cell->map_particle_volume_to_nodes(shapefns_xi, phase, pvolume);
      for (const auto& node : nodes)
        REQUIRE(node->volume(phase) == Approx(2.0).epsilon(Tolerance));
    }

    SECTION("Check particle momentum mapping") {
      // Map particle mass to nodes
      cell->map_particle_mass_to_nodes(shapefns_xi, phase, pmass);
      for (const auto& node : nodes)
        REQUIRE(node->mass(phase) == Approx(1.0).epsilon(Tolerance));

      // Map momentum to nodes
      cell->compute_nodal_momentum(shapefns_xi, phase, pmass, pvelocity);
      for (const auto& node : nodes) {
        for (unsigned i = 0; i < pvelocity.size(); ++i)
          REQUIRE(node->momentum(phase)(i) == Approx(1.0).epsilon(Tolerance));
      }

      // Update mass and momentum
      cell->map_mass_momentum_to_nodes(shapefns_xi, phase, pmass, pvelocity);
      for (const auto& node : nodes)
        REQUIRE(node->mass(phase) == Approx(2.0).epsilon(Tolerance));
      for (const auto& node : nodes) {
        for (unsigned i = 0; i < pvelocity.size(); ++i)
          REQUIRE(node->momentum(phase)(i) == Approx(2.0).epsilon(Tolerance));
      }
    }

    SECTION("Check particle strain") {
      // Particle mass
      pmass = 40.;

      // Map particle mass to nodes
      cell->map_particle_mass_to_nodes(shapefns_xi, phase, pmass);
      for (const auto& node : nodes)
        REQUIRE(node->mass(phase) == Approx(10.0).epsilon(Tolerance));

      // Map momentum to nodes
      cell->compute_nodal_momentum(shapefns_xi, phase, pmass, pvelocity);
      for (const auto& node : nodes) {
        for (unsigned i = 0; i < pvelocity.size(); ++i)
          REQUIRE(node->momentum(phase)(i) == Approx(10.0).epsilon(Tolerance));
      }

      // Update mass and momentum
      cell->map_mass_momentum_to_nodes(shapefns_xi, phase, pmass, pvelocity);
      for (const auto& node : nodes)
        REQUIRE(node->mass(phase) == Approx(20.0).epsilon(Tolerance));
      for (const auto& node : nodes) {
        for (unsigned i = 0; i < pvelocity.size(); ++i)
          REQUIRE(node->momentum(phase)(i) == Approx(20.0).epsilon(Tolerance));
      }

      // Update particle mass
      pmass = 80;
      cell->map_particle_mass_to_nodes(shapefns_xi, phase, pmass);

      // Compute nodal velocity
      for (const auto& node : nodes) {
        node->compute_velocity();
        for (unsigned i = 0; i < pvelocity.size(); ++i) {
          REQUIRE(node->momentum(phase)(i) == Approx(20.0).epsilon(Tolerance));
          REQUIRE(node->mass(phase) == Approx(40.0).epsilon(Tolerance));
          REQUIRE(node->velocity(phase)(i) == Approx(0.5).epsilon(Tolerance));
        }
      }

      Eigen::VectorXd strain_rate = cell->compute_strain_rate(bmatrix, phase);
      REQUIRE(strain_rate.size() == 3);
      for (unsigned i = 0; i < strain_rate.size(); ++i)
        REQUIRE(strain_rate(i) == Approx(0.).epsilon(Tolerance));

      Eigen::VectorXd strain_rate_centroid =
          cell->compute_strain_rate_centroid(phase);
      REQUIRE(strain_rate_centroid.size() == 3);
      for (unsigned i = 0; i < strain_rate_centroid.size(); ++i)
        REQUIRE(strain_rate_centroid(i) == Approx(0.).epsilon(Tolerance));
    }

    SECTION("Check particle body force mapping") {
      // Calculate body force at nodes
      cell->compute_nodal_body_force(shapefns_xi, phase, pmass, pgravity);
      Eigen::Vector2d bodyforce;
      bodyforce << 0., 9.814;
      for (const auto& node : nodes) {
        for (unsigned i = 0; i < bodyforce.size(); ++i)
          REQUIRE(node->external_force(phase)(i) ==
                  Approx(bodyforce(i)).epsilon(Tolerance));
      }
    }

    SECTION("Check particle traction force mapping") {
      // Check external force at nodes
      for (const auto& node : nodes)
        for (unsigned i = 0; i < Dim; ++i)
          REQUIRE(node->external_force(phase)(i) ==
                  Approx(0.).epsilon(Tolerance));

      // Apply traction force
      Eigen::Vector2d tractionforce;
      tractionforce << 1.5, 2.5;
      // Calculate traction force at nodes
      cell->compute_nodal_traction_force(shapefns_xi, phase, tractionforce);

      // Check traction force
      tractionforce *= 0.25;  // traction force * shapefn value (0.25)
      for (const auto& node : nodes)
        for (unsigned i = 0; i < tractionforce.size(); ++i)
          REQUIRE(node->external_force(phase)(i) ==
                  Approx(tractionforce(i)).epsilon(Tolerance));
    }

    SECTION("Check particle internal force mapping") {
      // Assign internal force to nodes
      const double pvolume = 0.5;
      Eigen::Matrix<double, 6, 1> pinternal_stress;
      pinternal_stress << 0.5, 0.5, 0.5, 0.5, 0.5, 0.5;

      cell->compute_nodal_internal_force(bmatrix, phase, pvolume,
                                         pinternal_stress);

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
      Eigen::Vector2d velocity =
          cell->interpolate_nodal_velocity(shapefns_xi, phase);

      Eigen::Vector2d interpolated_velocity;
      interpolated_velocity << 0.25, 0.25;
      for (unsigned i = 0; i < velocity.size(); ++i)
        REQUIRE(velocity(i) ==
                Approx(interpolated_velocity(i)).epsilon(Tolerance));

      // Check interpolate velocity (0.5, 0.5)
      xi << 0.5, 0.5;
      auto shapefn_xi = element->shapefn(xi);
      velocity = cell->interpolate_nodal_velocity(shapefn_xi, phase);

      interpolated_velocity << 0.2875, 0.2875;
      for (unsigned i = 0; i < velocity.size(); ++i)
        REQUIRE(velocity(i) ==
                Approx(interpolated_velocity(i)).epsilon(Tolerance));

      // Check interpolate velocity (-0.5, -0.5)
      xi << -0.5, -0.5;
      shapefn_xi = element->shapefn(xi);
      velocity = cell->interpolate_nodal_velocity(shapefn_xi, phase);

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
          cell->interpolate_nodal_acceleration(shapefns_xi, phase);

      Eigen::Vector2d interpolated_acceleration;
      interpolated_acceleration << 25., 25.;
      for (unsigned i = 0; i < check_acceleration.size(); ++i)
        REQUIRE(check_acceleration(i) ==
                Approx(interpolated_acceleration(i)).epsilon(Tolerance));

      // Check interpolate acceleration (0.5, 0.5)
      xi << 0.5, 0.5;
      auto shapefn_xi = element->shapefn(xi);
      check_acceleration =
          cell->interpolate_nodal_acceleration(shapefn_xi, phase);

      interpolated_acceleration << 28.75, 28.75;
      for (unsigned i = 0; i < check_acceleration.size(); ++i)
        REQUIRE(check_acceleration(i) ==
                Approx(interpolated_acceleration(i)).epsilon(Tolerance));

      // Check interpolate acceleration (-0.5, -0.5)
      xi << -0.5, -0.5;
      shapefn_xi = element->shapefn(xi);
      check_acceleration =
          cell->interpolate_nodal_acceleration(shapefn_xi, phase);

      interpolated_acceleration << 18.75, 18.75;
      for (unsigned i = 0; i < check_acceleration.size(); ++i)
        REQUIRE(check_acceleration(i) ==
                Approx(interpolated_acceleration(i)).epsilon(Tolerance));
    }
  }

  SECTION("Check velocity boundary conditions") {

    SECTION("Check normal vector") {

      Eigen::Vector2d coords;
      coords.setZero();

      coords << 0., 0.;
      std::shared_ptr<mpm::NodeBase<Dim>> node0 =
          std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

      coords << 4., 1.;
      std::shared_ptr<mpm::NodeBase<Dim>> node1 =
          std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);

      coords << 2., 2.;
      std::shared_ptr<mpm::NodeBase<Dim>> node2 =
          std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);

      coords << 0., 2.;
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

      cell->assign_cell_velocity_constraint(0, 0, 0.5);
      cell->assign_cell_velocity_constraint(0, 1, 0.6);
      cell->assign_cell_velocity_constraint(1, 1, 0.7);

      cell->compute_normals();

      REQUIRE(cell->nnormal() == 2);

      // Check normal vector when id is out of range range
      Eigen::Vector2d zero_normal_vector;
      zero_normal_vector << 0., 0.;
      auto check_zero_normal_vector = cell->normal(100);
      for (unsigned i = 0; i < Dim; i++) {
        REQUIRE(check_zero_normal_vector(i) == zero_normal_vector(i));
      }

      // Check normal vector
      Eigen::Matrix<double, 2, Dim> normal_vector;

      // clang-format off
      normal_vector << 0.242535625036333, -0.970142500145332,
                       0.447213595499958,  0.894427190999916;
      // clang-format on

      for (unsigned i = 0; i < cell->nnormal(); ++i) {
        REQUIRE(cell->normal(i).size() == Dim);
        for (unsigned j = 0; j < cell->normal(i).size(); ++j) {
          REQUIRE(cell->normal(i)(j) ==
                  Approx(normal_vector(i, j)).epsilon(Tolerance));
        }
      }
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

  // Check node additions
  SECTION("Add nodes") {
    mpm::Index id = 0;
    auto cell = std::make_shared<mpm::Cell<Dim>>(id, Nnodes, element);

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
            Approx(std::numeric_limits<double>::max()).epsilon(Tolerance));

    // Initialise cell
    REQUIRE(cell->initialise() == true);
    // Check if cell is initialised, after addition of nodes
    REQUIRE(cell->is_initialised() == true);

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

        // Check point 3
        // Coordinates of a point in real cell
        Eigen::Vector3d point3;
        point3 << 0., 0., 0.;

        // Coordinates of the point in an unit cell
        Eigen::Matrix<double, 3, 1> point_unit_cell3;
        point_unit_cell3 << -1., -1., -1.;

        // Check if point is in cell
        REQUIRE(cell1->point_in_cell(point3) == true);
        REQUIRE(cell2->point_in_cell(point3) == false);

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
        REQUIRE(cell1->point_in_cell(point4) == false);
        REQUIRE(cell2->point_in_cell(point4) == true);

        // Use Newton-raphson iteration to find local coordinates
        auto local_point4 = cell2->transform_real_to_unit_cell(point4);
        for (unsigned i = 0; i < local_point4.size(); ++i)
          REQUIRE(local_point4[i] ==
                  Approx(point_unit_cell4[i]).epsilon(Tolerance));
      }
    }
  }

  SECTION("Add neighbours") {
    auto cell = std::make_shared<mpm::Cell<Dim>>(0, Nnodes, element);
    auto neighbourcell = std::make_shared<mpm::Cell<Dim>>(1, Nnodes, element);
    REQUIRE(cell->nneighbours() == 0);
    cell->add_neighbour(0, neighbourcell);
    REQUIRE(cell->nneighbours() == 1);
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
    bool status = cell->add_particle_id(pid);
    REQUIRE(status == true);
    REQUIRE(cell->status() == true);
    REQUIRE(cell->nparticles() == 1);
    cell->remove_particle_id(pid);
    REQUIRE(cell->status() == false);
    REQUIRE(cell->nparticles() == 0);
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

    // Check point in cell
    Eigen::Vector3d point;
    point << 812482.5000000000, 815878.5000000000, 160.0825000000;
    REQUIRE(cell->point_in_cell(point) == true);
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

    // Point in cell1 fails to capture
    REQUIRE(cell1->point_in_cell(point) == false);

    Eigen::Vector3d xi;
    xi = cell1->transform_real_to_unit_cell(point);

    // Check using unit cell with affine transformation / Newton-Raphson
    REQUIRE(cell1->is_point_in_cell(point) == true);

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

    // Point in cell fails to capture
    REQUIRE(cell2->point_in_cell(point) == false);

    // Use affine transformation
    xi = cell2->transform_real_to_unit_cell(point);

    // Check using unit cell with affine transformation / Newton-Raphson
    REQUIRE(cell2->is_point_in_cell(point) == false);
  }

  SECTION("Test particle information mapping") {
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

    const auto shapefns_xi = element->shapefn(xi);
    const auto bmatrix = element->bmatrix(xi);

    SECTION("Check particle mass mapping") {
      cell->map_particle_mass_to_nodes(shapefns_xi, phase, pmass);
      for (const auto& node : nodes)
        REQUIRE(node->mass(phase) == Approx(0.5).epsilon(Tolerance));
    }

    SECTION("Check particle volume mapping") {
      cell->map_particle_volume_to_nodes(shapefns_xi, phase, pvolume);
      REQUIRE(nodes.size() == 8);
      for (const auto& node : nodes)
        REQUIRE(node->volume(phase) == Approx(1.0).epsilon(Tolerance));
    }

    SECTION("Check particle momentum mapping") {
      // Map particle mass to nodes
      cell->map_particle_mass_to_nodes(shapefns_xi, phase, pmass);
      for (const auto& node : nodes)
        REQUIRE(node->mass(phase) == Approx(0.5).epsilon(Tolerance));

      // Map momentum to nodes
      cell->compute_nodal_momentum(shapefns_xi, phase, pmass, pvelocity);
      for (const auto& node : nodes) {
        for (unsigned i = 0; i < pvelocity.size(); ++i)
          REQUIRE(node->momentum(phase)(i) == Approx(0.5).epsilon(Tolerance));
      }

      // Update mass and momentum
      cell->map_mass_momentum_to_nodes(shapefns_xi, phase, pmass, pvelocity);
      for (const auto& node : nodes)
        REQUIRE(node->mass(phase) == Approx(1.0).epsilon(Tolerance));
      for (const auto& node : nodes) {
        for (unsigned i = 0; i < pvelocity.size(); ++i)
          REQUIRE(node->momentum(phase)(i) == Approx(1.0).epsilon(Tolerance));
      }
    }

    SECTION("Check particle strain") {
      // Particle mass
      pmass = 40.;

      // Map particle mass to nodes
      cell->map_particle_mass_to_nodes(shapefns_xi, phase, pmass);
      for (const auto& node : nodes)
        REQUIRE(node->mass(phase) == Approx(5.).epsilon(Tolerance));

      // Map momentum to nodes
      cell->compute_nodal_momentum(shapefns_xi, phase, pmass, pvelocity);
      for (const auto& node : nodes) {
        for (unsigned i = 0; i < pvelocity.size(); ++i)
          REQUIRE(node->momentum(phase)(i) == Approx(5.).epsilon(Tolerance));
      }

      // Update mass and momentum
      cell->map_mass_momentum_to_nodes(shapefns_xi, phase, pmass, pvelocity);
      for (const auto& node : nodes)
        REQUIRE(node->mass(phase) == Approx(10.).epsilon(Tolerance));
      for (const auto& node : nodes) {
        for (unsigned i = 0; i < pvelocity.size(); ++i)
          REQUIRE(node->momentum(phase)(i) == Approx(10.).epsilon(Tolerance));
      }

      // Update particle mass
      pmass = 80.;
      cell->map_particle_mass_to_nodes(shapefns_xi, phase, pmass);

      // Compute nodal velocity
      for (const auto& node : nodes) {
        node->compute_velocity();
        for (unsigned i = 0; i < pvelocity.size(); ++i) {
          REQUIRE(node->momentum(phase)(i) == Approx(10.0).epsilon(Tolerance));
          REQUIRE(node->mass(phase) == Approx(20.0).epsilon(Tolerance));
          REQUIRE(node->velocity(phase)(i) == Approx(0.5).epsilon(Tolerance));
        }
      }

      Eigen::VectorXd strain_rate = cell->compute_strain_rate(bmatrix, phase);
      REQUIRE(strain_rate.size() == 6);
      for (unsigned i = 0; i < strain_rate.size(); ++i)
        REQUIRE(strain_rate(i) == Approx(0.).epsilon(Tolerance));

      Eigen::VectorXd strain_rate_centroid =
          cell->compute_strain_rate_centroid(phase);
      REQUIRE(strain_rate_centroid.size() == 6);
      for (unsigned i = 0; i < strain_rate_centroid.size(); ++i)
        REQUIRE(strain_rate_centroid(i) == Approx(0.).epsilon(Tolerance));
    }

    SECTION("Check particle body force mapping") {
      // Compute body force at nodes
      cell->compute_nodal_body_force(shapefns_xi, phase, pmass, pgravity);
      Eigen::Vector3d bodyforce;
      bodyforce << 0., 0., 0.5 * 9.814;
      for (const auto& node : nodes) {
        for (unsigned i = 0; i < bodyforce.size(); ++i)
          REQUIRE(node->external_force(phase)(i) ==
                  Approx(bodyforce(i)).epsilon(Tolerance));
      }
    }

    SECTION("Check particle traction force mapping") {
      // Check external force at nodes
      for (const auto& node : nodes)
        for (unsigned i = 0; i < Dim; ++i)
          REQUIRE(node->external_force(phase)(i) ==
                  Approx(0.).epsilon(Tolerance));

      // Apply traction force
      Eigen::Vector3d tractionforce;
      tractionforce << 1.5, 2.5, 3.7;
      // Calculate traction force at nodes
      cell->compute_nodal_traction_force(shapefns_xi, phase, tractionforce);

      // Check traction force
      tractionforce *= 0.125;  // traction force * shapefn value (0.25)
      for (const auto& node : nodes)
        for (unsigned i = 0; i < tractionforce.size(); ++i)
          REQUIRE(node->external_force(phase)(i) ==
                  Approx(tractionforce(i)).epsilon(Tolerance));
    }

    SECTION("Check particle internal force mapping") {
      // Assign internal force to nodes
      const double pvolume = 0.5;
      Eigen::Matrix<double, 6, 1> pinternal_stress;
      pinternal_stress << 0.5, 0.5, 0.5, 0.5, 0.5, 0.5;

      cell->compute_nodal_internal_force(bmatrix, phase, pvolume,
                                         pinternal_stress);

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
      Eigen::Vector3d velocity =
          cell->interpolate_nodal_velocity(shapefns_xi, phase);

      Eigen::Vector3d interpolated_velocity;
      interpolated_velocity << 0.45, 0.45, 0.45;
      for (unsigned i = 0; i < velocity.size(); ++i)
        REQUIRE(velocity(i) ==
                Approx(interpolated_velocity(i)).epsilon(Tolerance));

      // Check interpolate velocity (0.5, 0.5)
      xi << 0.5, 0.5, 0.5;
      auto shapefn_xi = element->shapefn(xi);
      velocity = cell->interpolate_nodal_velocity(shapefn_xi, phase);

      interpolated_velocity << 0.5875, 0.5875, 0.5875;
      for (unsigned i = 0; i < velocity.size(); ++i)
        REQUIRE(velocity(i) ==
                Approx(interpolated_velocity(i)).epsilon(Tolerance));

      // Check interpolate velocity (-0.5, -0.5)
      xi << -0.5, -0.5, -0.5;
      shapefn_xi = element->shapefn(xi);
      velocity = cell->interpolate_nodal_velocity(shapefn_xi, phase);

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
          cell->interpolate_nodal_acceleration(shapefns_xi, phase);

      Eigen::Vector3d interpolated_acceleration;
      interpolated_acceleration << 45.0, 45.0, 45.0;
      for (unsigned i = 0; i < check_acceleration.size(); ++i)
        REQUIRE(check_acceleration(i) ==
                Approx(interpolated_acceleration(i)).epsilon(Tolerance));

      // Check interpolate acceleration (0.5, 0.5, 0.5)
      xi << 0.5, 0.5, 0.5;
      auto shapefn_xi = element->shapefn(xi);
      check_acceleration =
          cell->interpolate_nodal_acceleration(shapefn_xi, phase);

      interpolated_acceleration << 58.75, 58.75, 58.75;
      for (unsigned i = 0; i < check_acceleration.size(); ++i)
        REQUIRE(check_acceleration(i) ==
                Approx(interpolated_acceleration(i)).epsilon(Tolerance));

      // Check interpolate acceleration (-0.5, -0.5, -0.5)
      xi << -0.5, -0.5, -0.5;
      shapefn_xi = element->shapefn(xi);
      check_acceleration =
          cell->interpolate_nodal_acceleration(shapefn_xi, phase);

      interpolated_acceleration << 28.75, 28.75, 28.75;
      for (unsigned i = 0; i < check_acceleration.size(); ++i)
        REQUIRE(check_acceleration(i) ==
                Approx(interpolated_acceleration(i)).epsilon(Tolerance));
    }
  }

  SECTION("Check velocity boundary conditions") {

    SECTION("Check normal vector") {
      // Coordinates
      Eigen::Vector3d coords;

      coords << 0, 0, 0;
      std::shared_ptr<mpm::NodeBase<Dim>> node0 =
          std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

      coords << 4, 1, 1;
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

      cell->assign_cell_velocity_constraint(0, 0, 0.5);
      cell->assign_cell_velocity_constraint(0, 1, 0.6);
      cell->assign_cell_velocity_constraint(1, 2, 0.7);

      cell->compute_normals();

      REQUIRE(cell->nnormal() == 2);

      // Check normal vector when id is out of range range
      Eigen::Vector3d zero_normal_vector;
      zero_normal_vector << 0., 0., 0.;
      auto check_zero_normal_vector = cell->normal(100);
      for (unsigned i = 0; i < Dim; i++) {
        REQUIRE(check_zero_normal_vector(i) == zero_normal_vector(i));
      }

      // Check normal vector
      Eigen::Matrix<double, 2, Dim> normal_vector;

      // clang-format off
      normal_vector <<  0.242535625036333, -0.970142500145332,  0.,
                       -0.301511344577764,  0.904534033733291,  0.301511344577764;
      // clang-format on

      for (unsigned i = 0; i < cell->nnormal(); ++i) {
        REQUIRE(cell->normal(i).size() == Dim);
        for (unsigned j = 0; j < cell->normal(i).size(); ++j) {
          REQUIRE(cell->normal(i)(j) ==
                  Approx(normal_vector(i, j)).epsilon(Tolerance));
        }
      }
    }
  }
}

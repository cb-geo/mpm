#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

#include "hex_shapefn.h"
#include "mesh.h"
#include "node.h"
#include "quad_shapefn.h"
#include "shapefn.h"

//! \brief Check mesh class for 2D case
TEST_CASE("Mesh is checked for 2D case", "[mesh][2D]") {
  // Dimension
  const unsigned Dim = 2;
  // Degrees of freedom
  const unsigned Dof = 2;
  // Number of nodes per cell
  const unsigned Nnodes = 4;

  //! Check Mesh IDs
  SECTION("Check mesh ids") {
    //! Check for id = 0
    SECTION("Mesh id is zero") {
      unsigned id = 0;
      auto mesh = std::make_shared<mpm::Mesh<Dim>>(id);
      REQUIRE(mesh->id() == 0);
    }

    SECTION("Mesh id is positive") {
      //! Check for id is a positive value
      unsigned id = std::numeric_limits<unsigned>::max();
      auto mesh = std::make_shared<mpm::Mesh<Dim>>(id);
      REQUIRE(mesh->id() == std::numeric_limits<unsigned>::max());
    }
  }

  SECTION("Add neighbours") {
    auto mesh = std::make_shared<mpm::Mesh<Dim>>(0);
    auto neighbourmesh = std::make_shared<mpm::Mesh<Dim>>(1);
    REQUIRE(mesh->nneighbours() == 0);
    mesh->add_neighbour(0, neighbourmesh);
    REQUIRE(mesh->nneighbours() == 1);
  }

  // Check add / remove particle
  SECTION("Check add / remove particle functionality") {
    // Particle 1
    mpm::Index id1 = 0;
    Eigen::Vector2d coords;
    coords.setZero();
    auto particle1 = std::make_shared<mpm::Particle<Dim>>(id1, coords);

    // Particle 2
    mpm::Index id2 = 1;
    auto particle2 = std::make_shared<mpm::Particle<Dim>>(id2, coords);

    auto mesh = std::make_shared<mpm::Mesh<Dim>>(0);
    // Check mesh is active
    REQUIRE(mesh->status() == false);

    // Add particle 1 and check status
    bool status1 = mesh->add_particle(particle1);
    REQUIRE(status1 == true);
    // Add particle 2 and check status
    bool status2 = mesh->add_particle(particle2);
    REQUIRE(status2 == true);
    // Add particle 2 again and check status
    bool status3 = mesh->add_particle(particle2);
    REQUIRE(status3 == false);

    // Check mesh is active
    REQUIRE(mesh->status() == true);
    // Check number of particles in mesh
    REQUIRE(mesh->nparticles() == 2);

    // Remove particle 2 and check status
    bool remove_status = mesh->remove_particle(particle2);
    REQUIRE(remove_status == true);
    // Check number of particles in mesh
    REQUIRE(mesh->nparticles() == 1);
  }

  // Check add / remove node
  SECTION("Check add / remove node functionality") {
    // Node 1
    mpm::Index id1 = 0;
    Eigen::Vector2d coords;
    coords.setZero();
    auto node1 = std::make_shared<mpm::Node<Dim>>(id1, coords, Dof);

    // Node 2
    mpm::Index id2 = 1;
    auto node2 = std::make_shared<mpm::Node<Dim>>(id2, coords, Dof);

    auto mesh = std::make_shared<mpm::Mesh<Dim>>(0);
    // Check mesh is active
    REQUIRE(mesh->status() == false);

    // Add node 1 and check status
    bool status1 = mesh->add_node(node1);
    REQUIRE(status1 == true);
    // Add node 2 and check status
    bool status2 = mesh->add_node(node2);
    REQUIRE(status2 == true);
    // Add node 2 again and check status
    bool status3 = mesh->add_node(node2);
    REQUIRE(status3 == false);

    // Check number of nodes in mesh
    REQUIRE(mesh->nnodes() == 2);

    // Remove node 2 and check status
    bool remove_status = mesh->remove_node(node2);
    REQUIRE(remove_status == true);
    // Check number of nodes in mesh
    REQUIRE(mesh->nnodes() == 1);
  }

  // Check add / remove cell
  SECTION("Check add / remove cell functionality") {
    // Cell 1
    mpm::Index id1 = 0;
    Eigen::Vector2d coords;
    coords.setZero();
    auto cell1 = std::make_shared<mpm::Cell<Dim>>(id1, Nnodes);

    // Cell 2
    mpm::Index id2 = 1;
    auto cell2 = std::make_shared<mpm::Cell<Dim>>(id2, Nnodes);

    auto mesh = std::make_shared<mpm::Mesh<Dim>>(0);
    // Check mesh is active
    REQUIRE(mesh->status() == false);

    // Add cell 1 and check status
    bool status1 = mesh->add_cell(cell1);
    REQUIRE(status1 == true);
    // Add cell 2 and check status
    bool status2 = mesh->add_cell(cell2);
    REQUIRE(status2 == true);
    // Add cell 2 again and check status
    bool status3 = mesh->add_cell(cell2);
    REQUIRE(status3 == false);

    // Check number of cells in mesh
    REQUIRE(mesh->ncells() == 2);

    // Remove cell 2 and check status
    bool remove_status = mesh->remove_cell(cell2);
    REQUIRE(remove_status == true);
    // Check number of cells in mesh
    REQUIRE(mesh->ncells() == 1);
  }
}

//! \brief Check mesh class for 3D case
TEST_CASE("Mesh is checked for 3D case", "[mesh][3D]") {
  // Dimension
  const unsigned Dim = 3;
  // Degrees of freedom
  const unsigned Dof = 6;
  // Number of nodes per cell
  const unsigned Nnodes = 8;
  
  //! Check Mesh IDs
  SECTION("Check mesh ids") {
    //! Check for id = 0
    SECTION("Mesh id is zero") {
      unsigned id = 0;
      auto mesh = std::make_shared<mpm::Mesh<Dim>>(id);
      REQUIRE(mesh->id() == 0);
    }

    SECTION("Mesh id is positive") {
      //! Check for id is a positive value
      unsigned id = std::numeric_limits<unsigned>::max();
      auto mesh = std::make_shared<mpm::Mesh<Dim>>(id);
      REQUIRE(mesh->id() == std::numeric_limits<unsigned>::max());
    }
  }

  SECTION("Add neighbours") {
    auto mesh = std::make_shared<mpm::Mesh<Dim>>(0);
    auto neighbourmesh = std::make_shared<mpm::Mesh<Dim>>(1);
    REQUIRE(mesh->nneighbours() == 0);
    mesh->add_neighbour(0, neighbourmesh);
    REQUIRE(mesh->nneighbours() == 1);
  }

  // Check add / remove particle
  SECTION("Check add / remove particle functionality") {
    // Particle 1
    mpm::Index id1 = 0;
    Eigen::Vector3d coords;
    coords.setZero();
    auto particle1 = std::make_shared<mpm::Particle<Dim>>(id1, coords);

    // Particle 2
    mpm::Index id2 = 1;
    auto particle2 = std::make_shared<mpm::Particle<Dim>>(id2, coords);

    auto mesh = std::make_shared<mpm::Mesh<Dim>>(0);
    // Check mesh is active
    REQUIRE(mesh->status() == false);

    // Add particle 1 and check status
    bool status1 = mesh->add_particle(particle1);
    REQUIRE(status1 == true);
    // Add particle 2 and check status
    bool status2 = mesh->add_particle(particle2);
    REQUIRE(status2 == true);
    // Add particle 2 again and check status
    bool status3 = mesh->add_particle(particle2);
    REQUIRE(status3 == false);

    // Check mesh is active
    REQUIRE(mesh->status() == true);
    // Check number of particles in mesh
    REQUIRE(mesh->nparticles() == 2);

    // Remove particle 2 and check status
    bool remove_status = mesh->remove_particle(particle2);
    REQUIRE(remove_status == true);
    // Check number of particles in mesh
    REQUIRE(mesh->nparticles() == 1);
  }

  // Check add / remove node
  SECTION("Check add / remove node functionality") {
    // Node 1
    mpm::Index id1 = 0;
    Eigen::Vector3d coords;
    coords.setZero();
    auto node1 = std::make_shared<mpm::Node<Dim>>(id1, coords, Dof);

    // Node 2
    mpm::Index id2 = 1;
    auto node2 = std::make_shared<mpm::Node<Dim>>(id2, coords, Dof);

    auto mesh = std::make_shared<mpm::Mesh<Dim>>(0);
    // Check mesh is active
    REQUIRE(mesh->status() == false);

    // Add node 1 and check status
    bool status1 = mesh->add_node(node1);
    REQUIRE(status1 == true);
    // Add node 2 and check status
    bool status2 = mesh->add_node(node2);
    REQUIRE(status2 == true);
    // Add node 2 again and check status
    bool status3 = mesh->add_node(node2);
    REQUIRE(status3 == false);

    // Check number of nodes in mesh
    REQUIRE(mesh->nnodes() == 2);

    // Remove node 2 and check status
    bool remove_status = mesh->remove_node(node2);
    REQUIRE(remove_status == true);
    // Check number of nodes in mesh
    REQUIRE(mesh->nnodes() == 1);
  }


  // Check add / remove cell
  SECTION("Check add / remove cell functionality") {
    // Cell 1
    mpm::Index id1 = 0;
    Eigen::Vector3d coords;
    coords.setZero();
    auto cell1 = std::make_shared<mpm::Cell<Dim>>(id1, Nnodes);

    // Cell 2
    mpm::Index id2 = 1;
    auto cell2 = std::make_shared<mpm::Cell<Dim>>(id2, Nnodes);

    auto mesh = std::make_shared<mpm::Mesh<Dim>>(0);
    // Check mesh is active
    REQUIRE(mesh->status() == false);

    // Add cell 1 and check status
    bool status1 = mesh->add_cell(cell1);
    REQUIRE(status1 == true);
    // Add cell 2 and check status
    bool status2 = mesh->add_cell(cell2);
    REQUIRE(status2 == true);
    // Add cell 2 again and check status
    bool status3 = mesh->add_cell(cell2);
    REQUIRE(status3 == false);

    // Check number of cells in mesh
    REQUIRE(mesh->ncells() == 2);

    // Remove cell 2 and check status
    bool remove_status = mesh->remove_cell(cell2);
    REQUIRE(remove_status == true);
    // Check number of cells in mesh
    REQUIRE(mesh->ncells() == 1);
  }
}

#include <limits>

#include "catch.hpp"

#include "data_types.h"
#include "element.h"
#include "hdf5_particle.h"
#include "hexahedron_element.h"
#include "material/material.h"
#include "mesh.h"
#include "mpi_datatypes.h"
#include "node.h"
#include "particle.h"
#include "quadrilateral_element.h"

//! Check transfer of particles across MPI tasks
TEST_CASE("MPI transfer particle is checked in 2D",
          "[particle][mpi][transfer][2D]") {
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

    // Create cell1
    coords.setZero();
    auto cell1 = std::make_shared<mpm::Cell<Dim>>(id1, Nnodes, element);

    // Add nodes to cell
    cell1->add_node(0, node0);
    cell1->add_node(1, node1);
    cell1->add_node(2, node2);
    cell1->add_node(3, node3);

    // Initialize cell
    REQUIRE(cell1->initialise() == true);

    // Initialize material models
    mesh->initialise_material_models(materials);

    // Add cell 1 and check
    REQUIRE(mesh->add_cell(cell1) == true);

    // Add 1 particle to rank 0
    if (mpi_rank == 0) {
      // Particle 1
      coords << 1.0, 1.0;
      std::shared_ptr<mpm::ParticleBase<Dim>> particle1 =
          std::make_shared<mpm::Particle<Dim>>(0, coords);
      particle1->assign_material(material);

      // Add particle 1 and check
      REQUIRE(mesh->add_particle(particle1) == true);

      // Check mesh is active
      REQUIRE(mesh->status() == true);

      // Locate particles in a mesh
      auto particles = mesh->locate_particles_mesh();

      // Should find all particles in mesh
      REQUIRE(particles.size() == 0);

      // Check location of particle 1
      REQUIRE(particle1->cell_id() == 0);

      // Number of particles in cell 1 is 2
      REQUIRE(cell1->nparticles() == 1);

      if (mpi_size == 2) {
        // Particle 2
        coords << 1.5, 1.5;
        std::shared_ptr<mpm::ParticleBase<Dim>> particle2 =
            std::make_shared<mpm::Particle<Dim>>(1, coords);
        particle2->assign_material(material);

        // Add particle 2 and check
        REQUIRE(mesh->add_particle(particle2) == true);

        // Check mesh is active
        REQUIRE(mesh->status() == true);

        // Locate particles in a mesh
        auto particles = mesh->locate_particles_mesh();

        // Should find all particles in mesh
        REQUIRE(particles.size() == 0);

        // Check location of particle 2
        REQUIRE(particle2->cell_id() == 0);

        // Number of particles in cell 1 is 2
        REQUIRE(cell1->nparticles() == 2);
      }
    }

    // Add 1 particle to rank 2
    if (mpi_rank == 2) {
      // Particle 2
      coords << 1.5, 1.5;
      std::shared_ptr<mpm::ParticleBase<Dim>> particle2 =
          std::make_shared<mpm::Particle<Dim>>(1, coords);
      particle2->assign_material(material);

      // Add particle 2 and check
      REQUIRE(mesh->add_particle(particle2) == true);

      // Check mesh is active
      REQUIRE(mesh->status() == true);

      // Locate particles in a mesh
      auto particles = mesh->locate_particles_mesh();

      // Should find all particles in mesh
      REQUIRE(particles.size() == 0);

      // Check location of particle 2
      REQUIRE(particle2->cell_id() == 0);

      // Number of particles in cell 1 is 2
      REQUIRE(cell1->nparticles() == 1);
    }

    if (mpi_size > 1) {
      // Assign a MPI rank of 1 to cell in all MPI ranks
      cell1->rank(1);

      // Transfer particle to the correct MPI rank
      mesh->transfer_nonrank_particles();
      // Check all non receiver ranks
      if (mpi_rank != receiver) {
        REQUIRE(cell1->nparticles() == 0);
        REQUIRE(mesh->nparticles() == 0);
      }
      // Number of particles in cell 1 is 0 for rank 2
      if (mpi_rank == receiver) {
        auto particles = mesh->locate_particles_mesh();
        REQUIRE(mesh->nparticles() == 2);
        REQUIRE(cell1->nparticles() == 2);
      }
    }
    SECTION("Check node MPI ranks") {
      if (mpi_size > 1) {
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

        auto cell2 = std::make_shared<mpm::Cell<Dim>>(2, Nnodes, element);

        // Add nodes to cell
        cell2->add_node(0, node1);
        cell2->add_node(1, node4);
        cell2->add_node(2, node5);
        cell2->add_node(3, node2);

        // Initialize cell
        REQUIRE(cell2->initialise() == true);

        // Initialize material models
        mesh->initialise_material_models(materials);

        // Add cell 1 and check
        REQUIRE(mesh->add_cell(cell2) == true);

        // Assign MPI ranks
        cell1->rank(0);
        cell2->rank(1);

        // Transfer particle to the correct MPI rank
        mesh->identify_domain_shared_nodes();

        REQUIRE(node0->mpi_ranks().size() == 1);
        REQUIRE(node1->mpi_ranks().size() == 2);
        REQUIRE(node2->mpi_ranks().size() == 2);
        REQUIRE(node3->mpi_ranks().size() == 1);
        REQUIRE(node4->mpi_ranks().size() == 1);
        REQUIRE(node5->mpi_ranks().size() == 1);
      }
    }
  }
}

//! Check transfer of particles across MPI tasks
TEST_CASE("MPI Transfer Particle is checked in 3D",
          "[particle][mpi][transfer][3D]") {
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

    // Create cell1
    auto cell1 = std::make_shared<mpm::Cell<Dim>>(id1, Nnodes, element);

    // Add nodes to cell
    cell1->add_node(0, node0);
    cell1->add_node(1, node1);
    cell1->add_node(2, node2);
    cell1->add_node(3, node3);
    cell1->add_node(4, node4);
    cell1->add_node(5, node5);
    cell1->add_node(6, node6);
    cell1->add_node(7, node7);

    REQUIRE(cell1->nnodes() == 8);

    // Initialize cell
    REQUIRE(cell1->initialise() == true);

    // Initialize material models
    mesh->initialise_material_models(materials);

    // Add cell 1 and check
    REQUIRE(mesh->add_cell(cell1) == true);

    if (mpi_rank == 0) {
      // Particle 1
      coords << 1.0, 1.0, 1.0;
      std::shared_ptr<mpm::ParticleBase<Dim>> particle1 =
          std::make_shared<mpm::Particle<Dim>>(0, coords);
      particle1->assign_material(material);

      // Particle 2
      coords << 1.5, 1.5, 1.5;
      std::shared_ptr<mpm::ParticleBase<Dim>> particle2 =
          std::make_shared<mpm::Particle<Dim>>(1, coords);
      particle2->assign_material(material);

      // Add particle 1 and check
      REQUIRE(mesh->add_particle(particle1) == true);
      // Add particle 2 and check
      REQUIRE(mesh->add_particle(particle2) == true);

      // Check mesh is active
      REQUIRE(mesh->status() == true);

      // Locate particles in a mesh
      auto particles = mesh->locate_particles_mesh();

      // Should find all particles in mesh
      REQUIRE(particles.size() == 0);

      // Check location of particle 1
      REQUIRE(particle1->cell_id() == 0);
      // Check location of particle 2
      REQUIRE(particle2->cell_id() == 0);

      // Number of particles in cell 1 is 2
      REQUIRE(cell1->nparticles() == 2);
    }

    // Assign a MPI rank of 1 to cell in all MPI ranks
    cell1->rank(1);

    if (mpi_size > 1) {
      // Transfer particle to the correct MPI rank
      mesh->transfer_nonrank_particles();
      // Check all non receiver ranks
      if (mpi_rank != receiver) {
        REQUIRE(cell1->nparticles() == 0);
        REQUIRE(mesh->nparticles() == 0);
      }
      // Number of particles in cell 1 is 0 for rank 2
      if (mpi_rank == receiver) {
        auto particles = mesh->locate_particles_mesh();
        REQUIRE(mesh->nparticles() == 2);
        REQUIRE(cell1->nparticles() == 2);
      }
    }
    SECTION("Check node MPI ranks") {
      if (mpi_size > 1) {
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

        auto cell2 = std::make_shared<mpm::Cell<Dim>>(2, Nnodes, element);

        // Add nodes to cell
        cell2->add_node(0, node4);
        cell2->add_node(1, node5);
        cell2->add_node(2, node6);
        cell2->add_node(3, node7);
        cell2->add_node(4, node8);
        cell2->add_node(5, node9);
        cell2->add_node(6, node10);
        cell2->add_node(7, node11);

        // Initialize cell
        REQUIRE(cell2->initialise() == true);

        // Initialize material models
        mesh->initialise_material_models(materials);

        // Add cell 1 and check
        REQUIRE(mesh->add_cell(cell2) == true);

        // Assign MPI ranks
        cell1->rank(0);
        cell2->rank(1);

        // Transfer particle to the correct MPI rank
        mesh->identify_domain_shared_nodes();

        REQUIRE(node0->mpi_ranks().size() == 1);
        REQUIRE(node1->mpi_ranks().size() == 1);
        REQUIRE(node2->mpi_ranks().size() == 1);
        REQUIRE(node3->mpi_ranks().size() == 1);
        REQUIRE(node4->mpi_ranks().size() == 2);
        REQUIRE(node5->mpi_ranks().size() == 2);
        REQUIRE(node6->mpi_ranks().size() == 2);
        REQUIRE(node7->mpi_ranks().size() == 2);
        REQUIRE(node8->mpi_ranks().size() == 1);
        REQUIRE(node9->mpi_ranks().size() == 1);
        REQUIRE(node10->mpi_ranks().size() == 1);
        REQUIRE(node11->mpi_ranks().size() == 1);
      }
    }
  }
}

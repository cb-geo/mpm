#include <limits>
#include <memory>

#include "catch.hpp"

#include "data_types.h"
#include "element.h"
#include "graph.h"
#include "hdf5_particle.h"
#include "hexahedron_element.h"
#include "material.h"
#include "mesh.h"
#include "mpi_datatypes.h"
#include "node.h"
#include "particle.h"
#include "quadrilateral_element.h"

#ifdef USE_MPI
#ifdef USE_KAHIP
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

  SECTION("Transfer particles in mesh domain decomposition of MPI cells") {

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
            "LinearElastic2D", std::move(1), jmaterial);

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

    // Add 1 particle to rank 1
    if (mpi_rank == 1) {
      // Particle 1
      coords << 0.5, 0.5;
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

      // Number of particles in cell 0 is 1
      REQUIRE(cell0->nparticles() == 1);
    }

    if (mpi_rank == 2) {
      // Particle 2
      coords << 1.5, 1.5;
      std::shared_ptr<mpm::ParticleBase<Dim>> particle2 =
          std::make_shared<mpm::Particle<Dim>>(2, coords);
      particle2->assign_material(material);

      // Add particle 2 and check
      REQUIRE(mesh->add_particle(particle2) == true);

      // Check mesh is active
      REQUIRE(mesh->status() == true);

      // Locate particles in a mesh
      auto lparticles = mesh->locate_particles_mesh();

      // Should find all particles in mesh
      REQUIRE(lparticles.size() == 0);

      // Check location of particle 2
      REQUIRE(particle2->cell_id() == 0);

      // Number of particles in cell 0 is 2
      REQUIRE(cell0->nparticles() == 1);
    }

    if (mpi_rank == 3) {
      // Particle 2
      coords << 0.5, 1.5;
      std::shared_ptr<mpm::ParticleBase<Dim>> particle3 =
          std::make_shared<mpm::Particle<Dim>>(3, coords);
      particle3->assign_material(material);

      // Add particle 3 and check
      REQUIRE(mesh->add_particle(particle3) == true);

      // Check mesh is active
      REQUIRE(mesh->status() == true);

      // Locate particles in a mesh
      auto lparticles = mesh->locate_particles_mesh();

      // Should find all particles in mesh
      REQUIRE(lparticles.size() == 0);

      // Check location of particle 3
      REQUIRE(particle3->cell_id() == 0);

      // Number of particles in cell 0 is 1
      REQUIRE(cell0->nparticles() == 1);
    }

    if (mpi_size == 4) {
      // Assign a MPI rank of 1 to cell in all MPI ranks
      cell0->rank(0);
      cell1->rank(1);
      cell2->rank(2);
      cell3->rank(3);
      // Identify ghost boundary cells
      mesh->find_cell_neighbours();
      mesh->find_ghost_boundary_cells();
      REQUIRE_NOTHROW(mesh->find_nglobal_particles_cells());
      // Transfer particle to the correct MPI rank
      mesh->transfer_halo_particles();
      // Check sender ranks
      if (mpi_rank != 0) {
        REQUIRE(cell0->nparticles() == 0);
        REQUIRE(mesh->nparticles() == 0);
      }
      // Number of particles in receiver rank
      if (mpi_rank == 0) {
        auto particles = mesh->locate_particles_mesh();
        REQUIRE(mesh->nparticles() == 3);
        REQUIRE(cell0->nparticles() == 3);
      }

      SECTION("Check node MPI ranks") {
        // Transfer particle to the correct MPI rank
        mesh->find_domain_shared_nodes();

        REQUIRE(node0->mpi_ranks().size() == 1);
        REQUIRE(node1->mpi_ranks().size() == 2);
        REQUIRE(node2->mpi_ranks().size() == 4);
        REQUIRE(node3->mpi_ranks().size() == 2);
        REQUIRE(node4->mpi_ranks().size() == 1);
        REQUIRE(node5->mpi_ranks().size() == 2);
        REQUIRE(node6->mpi_ranks().size() == 1);
        REQUIRE(node7->mpi_ranks().size() == 2);
        REQUIRE(node8->mpi_ranks().size() == 1);
      }
    }

    // Transfer non-rank particles
    if (mpi_size == 4) {
      // Initialize MPI
      MPI_Comm comm;
      MPI_Comm_dup(MPI_COMM_WORLD, &comm);

      // Create graph object if empty
      auto graph = std::make_shared<mpm::Graph<Dim>>(mesh->cells());

      // Find global nparticles
      REQUIRE_NOTHROW(mesh->find_nglobal_particles_cells());

      // Construct a weighted DAG
      REQUIRE_NOTHROW(graph->construct_graph(mpi_size, mpi_rank));

      // Graph partitioning mode
      int mode = 4;  // FAST
      // Create graph partition
      bool graphpartition = graph->create_partitions(&comm, mode);
      // Collect the partitions
      auto exchange_cells =
          graph->collect_partitions(mpi_size, mpi_rank, &comm);
      // Transfer particle to the correct MPI rank
      REQUIRE_NOTHROW(mesh->transfer_nonrank_particles(exchange_cells));
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
            "LinearElastic3D", std::move(1), jmaterial);

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

    if (mpi_rank == 1) {
      // Particle 1
      coords << 1.0, 1.0, 1.0;
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

      // Number of particles in cell 0 is 1
      REQUIRE(cell0->nparticles() == 1);
    }

    if (mpi_rank == 2) {
      // Particle 2
      coords << 1.5, 1.5, 1.5;
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

      // Number of particles in cell 0 is 1
      REQUIRE(cell0->nparticles() == 1);
    }

    if (mpi_rank == 3) {
      // Particle 3
      coords << 0.5, 0.5, 0.5;
      std::shared_ptr<mpm::ParticleBase<Dim>> particle3 =
          std::make_shared<mpm::Particle<Dim>>(3, coords);
      particle3->assign_material(material);

      // Add particle 3 and check
      REQUIRE(mesh->add_particle(particle3) == true);

      // Check mesh is active
      REQUIRE(mesh->status() == true);

      // Locate particles in a mesh
      auto particles = mesh->locate_particles_mesh();

      // Should find all particles in mesh
      REQUIRE(particles.size() == 0);

      // Check location of particle 3
      REQUIRE(particle3->cell_id() == 0);

      // Number of particles in cell 0 is 1
      REQUIRE(cell0->nparticles() == 1);
    }

    if (mpi_size == 4) {
      // Assign a MPI ranks to cells
      cell0->rank(0);
      cell1->rank(1);
      cell2->rank(2);
      cell3->rank(3);
      // Identify ghost boundary cells
      mesh->find_cell_neighbours();
      mesh->find_ghost_boundary_cells();
      // Test find nglobal particles in each cell
      REQUIRE_NOTHROW(mesh->find_nglobal_particles_cells());

      // Transfer particle to the correct MPI rank
      mesh->transfer_halo_particles();
      // Check all non receiver ranks
      if (mpi_rank != 0) {
        REQUIRE(cell0->nparticles() == 0);
        REQUIRE(mesh->nparticles() == 0);
      }
      // Number of particles in cell 0 is 3 in receiver
      if (mpi_rank == 0) {
        auto particles = mesh->locate_particles_mesh();
        REQUIRE(mesh->nparticles() == 3);
        REQUIRE(cell0->nparticles() == 3);
      }

      SECTION("Check node MPI ranks") {
        // Transfer particle to the correct MPI rank
        mesh->find_domain_shared_nodes();

        REQUIRE(node0->mpi_ranks().size() == 1);
        REQUIRE(node1->mpi_ranks().size() == 2);
        REQUIRE(node2->mpi_ranks().size() == 2);
        REQUIRE(node3->mpi_ranks().size() == 1);
        REQUIRE(node4->mpi_ranks().size() == 2);
        REQUIRE(node5->mpi_ranks().size() == 4);
        REQUIRE(node6->mpi_ranks().size() == 4);
        REQUIRE(node7->mpi_ranks().size() == 2);
        REQUIRE(node8->mpi_ranks().size() == 1);
        REQUIRE(node9->mpi_ranks().size() == 2);
        REQUIRE(node10->mpi_ranks().size() == 2);
        REQUIRE(node11->mpi_ranks().size() == 1);
        REQUIRE(node12->mpi_ranks().size() == 1);
        REQUIRE(node13->mpi_ranks().size() == 1);
        REQUIRE(node14->mpi_ranks().size() == 2);
        REQUIRE(node15->mpi_ranks().size() == 2);
        REQUIRE(node16->mpi_ranks().size() == 1);
        REQUIRE(node17->mpi_ranks().size() == 1);
      }
    }

    // Transfer non-rank particles
    if (mpi_size == 4) {
      // Initialize MPI
      MPI_Comm comm;
      MPI_Comm_dup(MPI_COMM_WORLD, &comm);

      // Create graph object if empty
      auto graph = std::make_shared<mpm::Graph<Dim>>(mesh->cells());

      // Find global nparticles
      REQUIRE_NOTHROW(mesh->find_nglobal_particles_cells());

      // Construct a weighted DAG
      REQUIRE_NOTHROW(graph->construct_graph(mpi_size, mpi_rank));

      // Graph partitioning mode
      int mode = 4;  // FAST
      // Create graph partition
      bool graphpartition = graph->create_partitions(&comm, mode);
      // Collect the partitions
      auto exchange_cells =
          graph->collect_partitions(mpi_size, mpi_rank, &comm);
      // Transfer particle to the correct MPI rank
      REQUIRE_NOTHROW(mesh->transfer_nonrank_particles(exchange_cells));
    }
  }
}
#endif
#endif

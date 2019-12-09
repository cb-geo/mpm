#include <iostream>
#include <limits>

#include "catch.hpp"

#include <iostream>

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

#ifdef USE_MPI
//! \brief Check transfer of particles across MPI tasks
TEST_CASE("MPI Transfer Particle is checked", "[particle][mpi][rank]") {
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

    int sender = 0;
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

    if (mpi_rank == 0) {
      // Particle 1
      coords << 1.0, 1.0;
      std::shared_ptr<mpm::ParticleBase<Dim>> particle1 =
          std::make_shared<mpm::Particle<Dim>>(0, coords);
      particle1->assign_material(material);

      // Particle 2
      coords << 1.5, 1.5;
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
    // Assign a MPI rank of 1 to cell
    cell1->rank(1);

    std::cout << "MPI size: " << mpi_size << std::endl;
    if (mpi_size > 1) {
      // Transfer particle to the correct MPI rank
      mesh->transfer_nonrank_particles(mpi_rank);

      if (mpi_rank == sender) {
        std::cout << "MPI sender rank: " << mpi_rank << std::endl;
        REQUIRE(cell1->nparticles() == 0);
        REQUIRE(mesh->nparticles() == 0);
      }
      if (mpi_rank == receiver) {
        // Transfer particle to the correct MPI rank
        std::cout << "MPI receiver rank: " << mpi_rank << std::endl;
        // Number of particles in cell 1 is 0 for rank 2
        REQUIRE(mesh->nparticles() == 2);
        std::cout << "Required number of particles: " << mesh->nparticles()
                  << " rank: " << mpi_rank << std::endl;
      }
    }
  }
}
#endif  // MPI

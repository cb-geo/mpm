#include <cmath>
#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"
// MPI
#ifdef USE_MPI
#include "mpi.h"
#endif

#include "element.h"
#include "hexahedron_element.h"
#include "mesh.h"
#include "node.h"
#include "quadrilateral_element.h"

//! \brief Check mesh class for 2D case
TEST_CASE("Mesh is checked for 2D case", "[mesh][2D]") {
  // Dimension
  const unsigned Dim = 2;
  // Degrees of freedom
  const unsigned Dof = 2;
  // Number of phases
  const unsigned Nphases = 1;
  // Number of nodes per cell
  const unsigned Nnodes = 4;
  // Tolerance
  const double Tolerance = 1.E-9;

  // 4-noded quadrilateral element
  std::shared_ptr<mpm::Element<Dim>> element =
      Factory<mpm::Element<Dim>>::instance()->create("ED2Q4");

  //! Check Mesh IDs
  SECTION("Check mesh ids") {
    //! Check for id = 0
    SECTION("Mesh id is zero") {
      unsigned id = 0;
      auto mesh = std::make_shared<mpm::Mesh<Dim>>(id);
      REQUIRE(mesh->id() == 0);
      REQUIRE(mesh->is_isoparametric() == true);
    }

    SECTION("Mesh id is zero and cartesian") {
      unsigned id = 0;
      auto mesh = std::make_shared<mpm::Mesh<Dim>>(id, false);
      REQUIRE(mesh->id() == 0);
      REQUIRE(mesh->is_isoparametric() == false);
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
    std::shared_ptr<mpm::ParticleBase<Dim>> particle1 =
        std::make_shared<mpm::Particle<Dim, Nphases>>(id1, coords);

    // Particle 2
    mpm::Index id2 = 1;
    coords << 2., 2.;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle2 =
        std::make_shared<mpm::Particle<Dim, Nphases>>(id2, coords);

    auto mesh = std::make_shared<mpm::Mesh<Dim>>(0);
    // Check mesh is active
    REQUIRE(mesh->status() == false);

    // Check nodal coordinates size
    REQUIRE(mesh->nodal_coordinates().size() == 0);
    // Check node pairs size
    REQUIRE(mesh->node_pairs().size() == 0);

    // Define nodes
    coords << 0., 0.;
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

    // Add cell 1 and check
    REQUIRE(mesh->add_cell(cell1) == true);

    // Check nodal coordinates size
    REQUIRE(mesh->nodal_coordinates().size() == 4);
    // Check node pairs size
    REQUIRE(mesh->node_pairs().size() == 4);

    // Add particle 1 and check
    REQUIRE(mesh->add_particle(particle1) == true);
    // Add particle 2 and check
    REQUIRE(mesh->add_particle(particle2) == true);
    // Add particle 2 again and check
    REQUIRE(mesh->add_particle(particle2) == false);

    // Check mesh is active
    REQUIRE(mesh->status() == true);
    // Check number of particles in mesh
    REQUIRE(mesh->nparticles() == 2);

    // Update coordinates
    Eigen::Vector2d coordinates;
    coordinates << 1., 1.;
    // Check iterate over functionality
    mesh->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Dim>::assign_coordinates,
                  std::placeholders::_1, coordinates));

    // Particle 1
    {
      // Check if nodal coordinate update has gone through
      auto check_coords = particle1->coordinates();
      // Check if coordinates for each particle is zero
      for (unsigned i = 0; i < check_coords.size(); ++i)
        REQUIRE(check_coords[i] == Approx(1.).epsilon(Tolerance));
    }
    // Particle 2
    {
      // Check if nodal coordinate update has gone through
      auto check_coords = particle2->coordinates();
      // Check if coordinates for each particle is zero
      for (unsigned i = 0; i < check_coords.size(); ++i)
        REQUIRE(check_coords[i] == Approx(1.).epsilon(Tolerance));
    }

    // Remove particle 2 and check
    REQUIRE(mesh->remove_particle(particle2) == true);
    // Check number of particles in mesh
    REQUIRE(mesh->nparticles() == 1);
  }

  // Check add / remove node
  SECTION("Check add / remove node functionality") {
    // Node 1
    mpm::Index id1 = 0;
    Eigen::Vector2d coords;
    coords.setZero();
    std::shared_ptr<mpm::NodeBase<Dim>> node1 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id1, coords);

    // Node 2
    mpm::Index id2 = 1;
    std::shared_ptr<mpm::NodeBase<Dim>> node2 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id2, coords);

    // Create a mesh
    auto mesh = std::make_shared<mpm::Mesh<Dim>>(0);
    // Check mesh is active
    REQUIRE(mesh->status() == false);

    // Check nodal coordinates size
    REQUIRE(mesh->nodal_coordinates().size() == 0);
    // Check node pairs size
    REQUIRE(mesh->node_pairs().size() == 0);

    // Add node 1 and check
    REQUIRE(mesh->add_node(node1) == true);
    // Add node 2 and check
    REQUIRE(mesh->add_node(node2) == true);
    // Add node 2 again and check
    REQUIRE(mesh->add_node(node2) == false);

    // Check number of nodes in mesh
    REQUIRE(mesh->nnodes() == 2);

    // Check nodal coordinates size
    REQUIRE(mesh->nodal_coordinates().size() == 2);

    // Update coordinates
    Eigen::Vector2d coordinates;
    coordinates << 1., 1.;

    // Set only node2 to be active
    node2->assign_status(true);

    // Check iterate over functionality if nodes are active
    mesh->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Dim>::assign_coordinates,
                  std::placeholders::_1, coordinates),
        std::bind(&mpm::NodeBase<Dim>::status, std::placeholders::_1));

    // Node 1
    {
      // Check if nodal coordinate update has gone through
      auto check_coords = node1->coordinates();
      // Check if coordinates for each node is zero
      for (unsigned i = 0; i < check_coords.size(); ++i)
        REQUIRE(check_coords[i] == Approx(0.).epsilon(Tolerance));
    }
    // Node 2
    {
      // Check if nodal coordinate update has gone through
      auto check_coords = node2->coordinates();
      // Check if coordinates for each node is zero
      for (unsigned i = 0; i < check_coords.size(); ++i)
        REQUIRE(check_coords[i] == Approx(1.).epsilon(Tolerance));
    }

    coordinates.setZero();
    node1->assign_coordinates(coordinates);
    node2->assign_coordinates(coordinates);
    for (unsigned i = 0; i < coordinates.size(); ++i) {
      REQUIRE(node1->coordinates()[i] == Approx(0.).epsilon(Tolerance));
      REQUIRE(node2->coordinates()[i] == Approx(0.).epsilon(Tolerance));
    }

    REQUIRE(node1->status() == false);
    REQUIRE(node2->status() == true);

    mesh->find_active_nodes();

    // Check iterate over functionality if nodes are active
    coordinates.fill(5.3);

    mesh->iterate_over_active_nodes(
        std::bind(&mpm::NodeBase<Dim>::assign_coordinates,
                  std::placeholders::_1, coordinates));

    // Node 1
    {
      // Check if nodal coordinate update has gone through
      auto check_coords = node1->coordinates();
      for (unsigned i = 0; i < check_coords.size(); ++i)
        REQUIRE(check_coords[i] == Approx(0.).epsilon(Tolerance));
    }
    // Node 2
    {
      // Check if nodal coordinate update has gone through
      auto check_coords = node2->coordinates();
      for (unsigned i = 0; i < check_coords.size(); ++i)
        REQUIRE(check_coords[i] == Approx(5.3).epsilon(Tolerance));
    }

    coordinates.fill(1.0);
    // Check iterate over functionality
    mesh->iterate_over_nodes(std::bind(&mpm::NodeBase<Dim>::assign_coordinates,
                                       std::placeholders::_1, coordinates));

    // Node 1
    {
      // Check if nodal coordinate update has gone through
      auto check_coords = node1->coordinates();
      // Check if coordinates for each node is zero
      for (unsigned i = 0; i < check_coords.size(); ++i)
        REQUIRE(check_coords[i] == Approx(1.).epsilon(Tolerance));
    }
    // Node 2
    {
      // Check if nodal coordinate update has gone through
      auto check_coords = node2->coordinates();
      // Check if coordinates for each node is zero
      for (unsigned i = 0; i < check_coords.size(); ++i)
        REQUIRE(check_coords[i] == Approx(1.).epsilon(Tolerance));
    }

    // Remove node 2 and check
    REQUIRE(mesh->remove_node(node2) == true);
    // Check number of nodes in mesh
    REQUIRE(mesh->nnodes() == 1);
  }

  // Check add / remove cell
  SECTION("Check add / remove cell functionality") {

    // Cell 1
    mpm::Index id1 = 0;
    Eigen::Vector2d coords;
    coords.setZero();
    auto cell1 = std::make_shared<mpm::Cell<Dim>>(id1, Nnodes, element);

    // Cell 2
    mpm::Index id2 = 1;
    auto cell2 = std::make_shared<mpm::Cell<Dim>>(id2, Nnodes, element);

    auto mesh = std::make_shared<mpm::Mesh<Dim>>(0);
    // Check mesh is active
    REQUIRE(mesh->status() == false);

    // Add cell 1 and check
    REQUIRE(mesh->add_cell(cell1) == true);
    // Add cell 2 and check
    REQUIRE(mesh->add_cell(cell2) == true);
    // Add cell 2 again and check
    REQUIRE(mesh->add_cell(cell2) == false);

    // Check number of cells in mesh
    REQUIRE(mesh->ncells() == 2);

    // Check iterate over functionality
    mesh->iterate_over_cells(
        std::bind(&mpm::Cell<Dim>::nnodes, std::placeholders::_1));

    // Remove cell 2 and check
    REQUIRE(mesh->remove_cell(cell2) == true);
    // Check number of cells in mesh
    REQUIRE(mesh->ncells() == 1);
  }

  SECTION("Check particle is in cell") {
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

    // Generate material points in cell
    auto points = mesh->generate_material_points(1);
    REQUIRE(points.size() == 0);

    // Add cell 1 and check
    REQUIRE(mesh->add_cell(cell1) == true);

    // Generate material points in cell
    points = mesh->generate_material_points(1);
    REQUIRE(points.size() == 1);

    points = mesh->generate_material_points(2);
    REQUIRE(points.size() == 4);

    points = mesh->generate_material_points(3);
    REQUIRE(points.size() == 9);

    // Particle 1
    coords << 1.0, 1.0;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle1 =
        std::make_shared<mpm::Particle<Dim, Nphases>>(0, coords);

    // Particle 2
    coords << 1.5, 1.5;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle2 =
        std::make_shared<mpm::Particle<Dim, Nphases>>(1, coords);

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
  }

  //! Check create nodes and cells in a mesh
  SECTION("Check create nodes and cells") {
    // Vector of nodal coordinates
    std::vector<Eigen::Matrix<double, Dim, 1>> coordinates;

    // Nodal coordinates
    Eigen::Matrix<double, Dim, 1> node;

    // Cell 0
    // Node 0
    node << 0., 0.;
    coordinates.emplace_back(node);
    // Node 1
    node << 0.5, 0.;
    coordinates.emplace_back(node);
    // Node 2
    node << 0.5, 0.5;
    coordinates.emplace_back(node);
    // Node 3
    node << 0., 0.5;
    coordinates.emplace_back(node);

    // Cell 1
    // Node 4
    node << 1.0, 0.;
    coordinates.emplace_back(node);
    // Node 5
    node << 1.0, 0.5;
    coordinates.emplace_back(node);

    // Create a new mesh
    unsigned meshid = 0;
    auto mesh = std::make_shared<mpm::Mesh<Dim>>(meshid);

    SECTION("Check creation of nodes") {
      // Node type 2D
      const std::string node_type = "N2D";
      // Global node index
      mpm::Index gnid = 0;
      mesh->create_nodes(gnid, node_type, coordinates, false);
      // Check if mesh has added nodes
      REQUIRE(mesh->nnodes() == coordinates.size());
      // Try again this shouldn't add more coordinates
      mesh->create_nodes(gnid, node_type, coordinates);
      // Check if mesh has added nodes
      REQUIRE(mesh->nnodes() == coordinates.size());
      // Clear coordinates and try creating a list of nodes with an empty list
      unsigned nnodes = coordinates.size();
      coordinates.clear();
      // This fails with empty list error in node creation
      mesh->create_nodes(gnid, node_type, coordinates);
      REQUIRE(mesh->nnodes() == nnodes);

      SECTION("Check creation of cells") {
        // Cell with node ids
        std::vector<std::vector<mpm::Index>> cells{// cell #0
                                                   {0, 1, 2, 3},
                                                   // cell #1
                                                   {1, 4, 5, 2}};
        // Assign 4-noded quadrilateral element to cell
        std::shared_ptr<mpm::Element<Dim>> element =
            Factory<mpm::Element<Dim>>::instance()->create("ED2Q4");

        // Global cell index
        mpm::Index gcid = 0;
        mesh->create_cells(gcid, element, cells, false);
        // Check if mesh has added cells
        REQUIRE(mesh->ncells() == cells.size());
        // Try again this shouldn't add more cells
        mesh->create_cells(gcid, element, cells);
        // Check if mesh has added cells
        REQUIRE(mesh->ncells() == cells.size());
        // Clear cells and try creating a list of empty cells
        unsigned ncells = cells.size();
        cells.clear();
        // This fails with empty list error in node creation
        gcid = 10;
        mesh->create_cells(gcid, element, cells);
        REQUIRE(mesh->ncells() == ncells);

        // Try with invalid node ids
        cells = {// cell #0
                 {10, 11, 12, 13},
                 // cell #1
                 {11, 14, 15, 12}};
        gcid = 20;
        mesh->create_cells(gcid, element, cells);
        REQUIRE(mesh->ncells() == ncells);

        SECTION("Check creation of particles") {
          // Vector of particle coordinates
          std::vector<Eigen::Matrix<double, Dim, 1>> coordinates;
          coordinates.clear();

          // Particle coordinates
          Eigen::Matrix<double, Dim, 1> particle;

          // Cell 0
          // Particle 0
          particle << 0.125, 0.125;
          coordinates.emplace_back(particle);
          // Particle 1
          particle << 0.25, 0.125;
          coordinates.emplace_back(particle);
          // Particle 2
          particle << 0.25, 0.25;
          coordinates.emplace_back(particle);
          // Particle 3
          particle << 0.125, 0.25;
          coordinates.emplace_back(particle);

          // Cell 1
          // Particle 4
          particle << 0.675, 0.125;
          coordinates.emplace_back(particle);
          // Particle 5
          particle << 0.85, 0.125;
          coordinates.emplace_back(particle);
          // Particle 6
          particle << 0.85, 0.25;
          coordinates.emplace_back(particle);
          // Particle 7
          particle << 0.675, 0.25;
          coordinates.emplace_back(particle);

          SECTION("Check addition of particles to mesh") {
            // Particle type 2D
            const std::string particle_type = "P2D";
            // Global particle index
            std::vector<mpm::Index> gpid(coordinates.size());
            std::iota(gpid.begin(), gpid.end(), 0);
            mesh->create_particles(gpid, particle_type, coordinates, false);
            // Check if mesh has added particles
            REQUIRE(mesh->nparticles() == coordinates.size());
            // Try again this shouldn't add more coordinates
            mesh->create_particles(gpid, particle_type, coordinates);
            // Check if mesh has added particles
            REQUIRE(mesh->nparticles() == coordinates.size());

            // Test assign particles cells
            SECTION("Check assign particles cells") {
              // Vector of particle cells
              std::vector<std::array<mpm::Index, 2>> particles_cells;
              // Particle cells
              particles_cells.emplace_back(std::array<mpm::Index, 2>({0, 0}));
              particles_cells.emplace_back(std::array<mpm::Index, 2>({1, 0}));
              particles_cells.emplace_back(std::array<mpm::Index, 2>({2, 0}));
              particles_cells.emplace_back(std::array<mpm::Index, 2>({3, 0}));

              particles_cells.emplace_back(std::array<mpm::Index, 2>({4, 1}));
              particles_cells.emplace_back(std::array<mpm::Index, 2>({5, 1}));
              particles_cells.emplace_back(std::array<mpm::Index, 2>({6, 1}));
              particles_cells.emplace_back(std::array<mpm::Index, 2>({7, 1}));

              REQUIRE(mesh->assign_particles_cells(particles_cells) == true);

              // Locate particles
              auto missing_particles = mesh->locate_particles_mesh();
              REQUIRE(missing_particles.size() == 0);

              auto check_particles_cells = mesh->particles_cells();

              REQUIRE(check_particles_cells.size() == mesh->nparticles());

              for (unsigned i = 0; i < particles_cells.size(); ++i)
                for (unsigned j = 0; j < 2; ++j)
                  REQUIRE(check_particles_cells.at(i).at(j) ==
                          particles_cells.at(i).at(j));
            }

            // Locate particles
            auto missing_particles = mesh->locate_particles_mesh();
            REQUIRE(missing_particles.size() == 0);

            // Test assign particles cells again should fail
            SECTION("Check assign particles cells") {
              // Vector of particle cells
              std::vector<std::array<mpm::Index, 2>> particles_cells;
              // Particle cells
              particles_cells.emplace_back(std::array<mpm::Index, 2>({0, 0}));
              particles_cells.emplace_back(std::array<mpm::Index, 2>({1, 0}));
              particles_cells.emplace_back(std::array<mpm::Index, 2>({2, 0}));
              particles_cells.emplace_back(std::array<mpm::Index, 2>(
                  {3, std::numeric_limits<mpm::Index>::max()}));
              particles_cells.emplace_back(std::array<mpm::Index, 2>({50, 0}));

              REQUIRE(mesh->assign_particles_cells(particles_cells) == false);
            }

            // Clear coordinates and try creating a list of particles with
            // an empty list
            unsigned nparticles = coordinates.size();
            coordinates.clear();
            // This fails with empty list error in particle creation
            mesh->create_particles(gpid, particle_type, coordinates);
            REQUIRE(mesh->nparticles() == nparticles);

            const unsigned phase = 0;
            // Particles coordinates
            REQUIRE(mesh->particle_coordinates().size() == mesh->nparticles());
            // Particle stresses
            std::string attribute = "stresses";
            REQUIRE(mesh->particles_vector_data(attribute, phase).size() ==
                    mesh->nparticles());
            // Particle strains
            attribute = "strains";
            REQUIRE(mesh->particles_vector_data(attribute, phase).size() ==
                    mesh->nparticles());
            // Particle velocities
            attribute = "velocities";
            REQUIRE(mesh->particles_vector_data(attribute, phase).size() ==
                    mesh->nparticles());
            // Particle invalid data
            attribute = "invalid";
            REQUIRE(mesh->particles_vector_data(attribute, phase).size() == 0);

            // Locate particles in mesh
            SECTION("Locate particles in mesh") {
              // Locate particles in a mesh
              auto particles = mesh->locate_particles_mesh();

              // Should find all particles in mesh
              REQUIRE(particles.size() == 0);

              // Create particle 100
              Eigen::Vector2d coords;
              coords << 100., 100.;

              mpm::Index pid = 100;
              std::shared_ptr<mpm::ParticleBase<Dim>> particle100 =
                  std::make_shared<mpm::Particle<Dim, Nphases>>(pid, coords);

              // Add particle100 and check
              REQUIRE(mesh->add_particle(particle100) == false);

              // Locate particles in a mesh
              particles = mesh->locate_particles_mesh();
              // Should miss particle100
              REQUIRE(particles.size() == 0);
            }

            // Test HDF5
            SECTION("Write particles HDF5") {
              REQUIRE(mesh->write_particles_hdf5(0, "particles-2d.h5") == true);
            }

            // Test assign particles volumes
            SECTION("Check assign particles volumes") {
              // Vector of particle coordinates
              std::vector<std::tuple<mpm::Index, double>> particles_volumes;
              // Volumes
              particles_volumes.emplace_back(std::make_tuple(0, 10.5));
              particles_volumes.emplace_back(std::make_tuple(1, 10.5));

              REQUIRE(mesh->nparticles() == 8);

              REQUIRE(mesh->assign_particles_volumes(particles_volumes) ==
                      true);

              // When volume assignment fails
              particles_volumes.emplace_back(std::make_tuple(2, 0.0));
              particles_volumes.emplace_back(std::make_tuple(3, -10.0));

              REQUIRE(mesh->assign_particles_volumes(particles_volumes) ==
                      false);
            }

            // Test assign particles tractions
            SECTION("Check assign particles tractions") {
              // Vector of particle coordinates
              std::vector<std::tuple<mpm::Index, unsigned, double>>
                  particles_tractions;
              // Tractions
              particles_tractions.emplace_back(std::make_tuple(0, 0, 10.5));
              particles_tractions.emplace_back(std::make_tuple(1, 1, -10.5));
              particles_tractions.emplace_back(std::make_tuple(2, 0, -12.5));
              particles_tractions.emplace_back(std::make_tuple(3, 1, 0.0));

              REQUIRE(mesh->nparticles() == 8);

              REQUIRE(mesh->assign_particles_tractions(particles_tractions) ==
                      false);
              // Compute volume
              mesh->iterate_over_particles(
                  std::bind(&mpm::ParticleBase<Dim>::compute_volume,
                            std::placeholders::_1, phase));

              REQUIRE(mesh->assign_particles_tractions(particles_tractions) ==
                      true);
              // When tractions fail
              particles_tractions.emplace_back(std::make_tuple(3, 2, 0.0));
              REQUIRE(mesh->assign_particles_tractions(particles_tractions) ==
                      false);
              particles_tractions.emplace_back(std::make_tuple(300, 0, 0.0));
              REQUIRE(mesh->assign_particles_tractions(particles_tractions) ==
                      false);
            }

            // Test assign nodes tractions
            SECTION("Check assign nodes tractions") {
              // Vector of node coordinates
              std::vector<std::tuple<mpm::Index, unsigned, double>>
                  nodes_tractions;
              // Tractions
              nodes_tractions.emplace_back(std::make_tuple(0, 0, 10.5));
              nodes_tractions.emplace_back(std::make_tuple(1, 1, -10.5));
              nodes_tractions.emplace_back(std::make_tuple(2, 0, -12.5));
              nodes_tractions.emplace_back(std::make_tuple(3, 1, 0.0));

              REQUIRE(mesh->nnodes() == 6);

              REQUIRE(mesh->assign_nodal_tractions(nodes_tractions) == true);
              // When tractions fail
              nodes_tractions.emplace_back(std::make_tuple(3, 2, 0.0));
              REQUIRE(mesh->assign_nodal_tractions(nodes_tractions) == false);
              nodes_tractions.emplace_back(std::make_tuple(300, 0, 0.0));
              REQUIRE(mesh->assign_nodal_tractions(nodes_tractions) == false);
            }

            // Test assign particles stresses
            SECTION("Check assign particles stresses") {
              // Vector of particle stresses
              std::vector<Eigen::Matrix<double, 6, 1>> particles_stresses;

              REQUIRE(mesh->nparticles() == 8);

              // Stresses
              particles_stresses.emplace_back(
                  Eigen::Matrix<double, 6, 1>::Constant(0.0));
              particles_stresses.emplace_back(
                  Eigen::Matrix<double, 6, 1>::Constant(0.1));
              particles_stresses.emplace_back(
                  Eigen::Matrix<double, 6, 1>::Constant(-0.2));
              particles_stresses.emplace_back(
                  Eigen::Matrix<double, 6, 1>::Constant(0.3));
              particles_stresses.emplace_back(
                  Eigen::Matrix<double, 6, 1>::Constant(-0.4));
              particles_stresses.emplace_back(
                  Eigen::Matrix<double, 6, 1>::Constant(0.5));
              particles_stresses.emplace_back(
                  Eigen::Matrix<double, 6, 1>::Constant(-0.6));
              particles_stresses.emplace_back(
                  Eigen::Matrix<double, 6, 1>::Constant(-0.7));

              REQUIRE(mesh->assign_particles_stresses(particles_stresses) ==
                      true);
              // When stresses fail
              particles_stresses.emplace_back(
                  Eigen::Matrix<double, 6, 1>::Constant(0.8));
              REQUIRE(mesh->assign_particles_stresses(particles_stresses) ==
                      false);
              unsigned id = 1;
              auto mesh_fail = std::make_shared<mpm::Mesh<Dim>>(id);
              REQUIRE(mesh_fail->assign_particles_stresses(
                          particles_stresses) == false);
            }

            // Test assign particles velocity constraints
            SECTION("Check assign particles velocity constraints") {
              // Vector of particle coordinates
              std::vector<std::tuple<mpm::Index, unsigned, double>>
                  particles_velocity_constraints;
              // Constraint
              particles_velocity_constraints.emplace_back(
                  std::make_tuple(0, 0, 10.5));
              particles_velocity_constraints.emplace_back(
                  std::make_tuple(1, 1, -10.5));
              particles_velocity_constraints.emplace_back(
                  std::make_tuple(2, 0, -12.5));
              particles_velocity_constraints.emplace_back(
                  std::make_tuple(3, 1, 0.0));

              REQUIRE(mesh->assign_particles_velocity_constraints(
                          particles_velocity_constraints) == true);
              // When constraints fail
              particles_velocity_constraints.emplace_back(
                  std::make_tuple(3, 2, 0.0));
              REQUIRE(mesh->assign_particles_velocity_constraints(
                          particles_velocity_constraints) == false);
            }
          }
        }
        // Test assign velocity constraints to nodes
        SECTION("Check assign velocity constraints to nodes") {
          // Vector of particle coordinates
          std::vector<std::tuple<mpm::Index, unsigned, double>>
              velocity_constraints;
          // Constraint
          velocity_constraints.emplace_back(std::make_tuple(0, 0, 10.5));
          velocity_constraints.emplace_back(std::make_tuple(1, 1, -10.5));
          velocity_constraints.emplace_back(std::make_tuple(2, 0, -12.5));
          velocity_constraints.emplace_back(std::make_tuple(3, 1, 0.0));

          REQUIRE(mesh->assign_velocity_constraints(velocity_constraints) ==
                  true);
          // When constraints fail
          velocity_constraints.emplace_back(std::make_tuple(3, 2, 0.0));
          REQUIRE(mesh->assign_velocity_constraints(velocity_constraints) ==
                  false);
        }

        // Test assign rotation matrices to nodes
        SECTION("Check assign rotation matrices to nodes") {
          // Map of nodal id and euler angles
          std::map<mpm::Index, Eigen::Matrix<double, Dim, 1>> euler_angles;
          // Node 0 with Euler angles of 10 and 20 deg
          euler_angles.emplace(std::make_pair(
              0, (Eigen::Matrix<double, Dim, 1>() << 10. * M_PI / 180,
                  20. * M_PI / 180)
                     .finished()));
          // Node 1 with Euler angles of 30 and 40 deg
          euler_angles.emplace(std::make_pair(
              1, (Eigen::Matrix<double, Dim, 1>() << 30. * M_PI / 180,
                  40. * M_PI / 180)
                     .finished()));
          // Node 2 with Euler angles of 50 and 60 deg
          euler_angles.emplace(std::make_pair(
              2, (Eigen::Matrix<double, Dim, 1>() << 50. * M_PI / 180,
                  60. * M_PI / 180)
                     .finished()));
          // Node 3 with Euler angles of 70 and 80 deg
          euler_angles.emplace(std::make_pair(
              3, (Eigen::Matrix<double, Dim, 1>() << 70. * M_PI / 180,
                  80. * M_PI / 180)
                     .finished()));

          // Check compute and assign rotation matrix
          REQUIRE(mesh->compute_nodal_rotation_matrices(euler_angles) == true);

          // Check for failure when missing node id
          // Node 100 (non-existent) with Euler angles of 90 and 100 deg
          euler_angles.emplace(std::make_pair(
              100, (Eigen::Matrix<double, Dim, 1>() << 90. * M_PI / 180,
                    100. * M_PI / 180)
                       .finished()));
          REQUIRE(mesh->compute_nodal_rotation_matrices(euler_angles) == false);

          // Check for failure of empty input
          std::map<mpm::Index, Eigen::Matrix<double, Dim, 1>>
              empty_euler_angles;
          REQUIRE(mesh->compute_nodal_rotation_matrices(empty_euler_angles) ==
                  false);

          // Check for failure when no nodes are assigned
          auto mesh_fail = std::make_shared<mpm::Mesh<Dim>>(1);
          REQUIRE(mesh_fail->compute_nodal_rotation_matrices(euler_angles) ==
                  false);
        }

        // Test assign velocity constraints to cells
        SECTION("Check assign velocity constraints to cells") {
          // Vector of particle coordinates
          std::vector<std::tuple<mpm::Index, unsigned, unsigned, double>>
              velocity_constraints;
          // Constraint
          velocity_constraints.emplace_back(std::make_tuple(0, 3, 0, 10.5));
          velocity_constraints.emplace_back(std::make_tuple(1, 2, 1, -10.5));

          REQUIRE(mesh->assign_cell_velocity_constraints(
                      velocity_constraints) == true);
          // When constraints fail
          velocity_constraints.emplace_back(std::make_tuple(1, 10, 1, 0.0));
          REQUIRE(mesh->assign_cell_velocity_constraints(
                      velocity_constraints) == false);
        }
      }
    }
  }
}

//! \brief Check mesh class for 3D case
TEST_CASE("Mesh is checked for 3D case", "[mesh][3D]") {
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

  // 8-noded hexahedron element
  std::shared_ptr<mpm::Element<Dim>> element =
      Factory<mpm::Element<Dim>>::instance()->create("ED3H8");

  //! Check Mesh IDs
  SECTION("Check mesh ids") {
    //! Check for id = 0
    SECTION("Mesh id is zero") {
      unsigned id = 0;
      auto mesh = std::make_shared<mpm::Mesh<Dim>>(id);
      REQUIRE(mesh->id() == 0);
      REQUIRE(mesh->is_isoparametric() == true);
    }

    SECTION("Mesh id is zero and cartesian") {
      unsigned id = 0;
      auto mesh = std::make_shared<mpm::Mesh<Dim>>(id, false);
      REQUIRE(mesh->id() == 0);
      REQUIRE(mesh->is_isoparametric() == false);
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
    std::shared_ptr<mpm::ParticleBase<Dim>> particle1 =
        std::make_shared<mpm::Particle<Dim, Nphases>>(id1, coords);

    // Particle 2
    mpm::Index id2 = 1;
    coords << 2., 2., 2.;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle2 =
        std::make_shared<mpm::Particle<Dim, Nphases>>(id2, coords);

    auto mesh = std::make_shared<mpm::Mesh<Dim>>(0);
    // Check mesh is active
    REQUIRE(mesh->status() == false);

    // Check nodal coordinates size
    REQUIRE(mesh->nodal_coordinates().size() == 0);
    // Check node pairs size
    REQUIRE(mesh->node_pairs().size() == 0);

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

    REQUIRE(mesh->add_cell(cell1) == true);

    REQUIRE(cell1->initialise() == true);

    // Check nodal coordinates size
    REQUIRE(mesh->nodal_coordinates().size() == 8);
    // Check node pairs size
    REQUIRE(mesh->node_pairs().size() == 12);

    // Add particle 1 and check
    REQUIRE(mesh->add_particle(particle1) == true);
    // Add particle 2 and check
    REQUIRE(mesh->add_particle(particle2) == true);
    // Add particle 2 again and check
    REQUIRE(mesh->add_particle(particle2) == false);

    // Check mesh is active
    REQUIRE(mesh->status() == true);
    // Check number of particles in mesh
    REQUIRE(mesh->nparticles() == 2);

    // Remove particle 2 and check
    REQUIRE(mesh->remove_particle(particle2) == true);
    // Check number of particles in mesh
    REQUIRE(mesh->nparticles() == 1);
  }

  // Check add / remove node
  SECTION("Check add / remove node functionality") {
    // Node 1
    mpm::Index id1 = 0;
    Eigen::Vector3d coords;
    coords.setZero();
    std::shared_ptr<mpm::NodeBase<Dim>> node1 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id1, coords);

    // Node 2
    mpm::Index id2 = 1;
    std::shared_ptr<mpm::NodeBase<Dim>> node2 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id2, coords);

    auto mesh = std::make_shared<mpm::Mesh<Dim>>(0);
    // Check mesh is active
    REQUIRE(mesh->status() == false);

    // Check nodal coordinates size
    REQUIRE(mesh->nodal_coordinates().size() == 0);
    // Check node pairs size
    REQUIRE(mesh->node_pairs().size() == 0);

    // Add node 1 and check
    REQUIRE(mesh->add_node(node1) == true);
    // Add node 2 and check
    REQUIRE(mesh->add_node(node2) == true);
    // Add node 2 again and check
    REQUIRE(mesh->add_node(node2) == false);

    // Check nodal coordinates size
    REQUIRE(mesh->nodal_coordinates().size() == 2);

    // Check number of nodes in mesh
    REQUIRE(mesh->nnodes() == 2);

    // Update coordinates
    Eigen::Vector3d coordinates;
    coordinates << 7., 7., 7.;

    // Set only node2 to be active
    node2->assign_status(true);

    // Check iterate over functionality if nodes are active
    mesh->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Dim>::assign_coordinates,
                  std::placeholders::_1, coordinates),
        std::bind(&mpm::NodeBase<Dim>::status, std::placeholders::_1));

    // Node 1
    {
      // Check if nodal coordinate update has gone through
      auto check_coords = node1->coordinates();
      // Check if coordinates for each node is zero
      for (unsigned i = 0; i < check_coords.size(); ++i)
        REQUIRE(check_coords[i] == Approx(0.).epsilon(Tolerance));
    }
    // Node 2
    {
      // Check if nodal coordinate update has gone through
      auto check_coords = node2->coordinates();
      // Check if coordinates for each node is zero
      for (unsigned i = 0; i < check_coords.size(); ++i)
        REQUIRE(check_coords[i] == Approx(7.).epsilon(Tolerance));
    }

    coordinates.setZero();
    node1->assign_coordinates(coordinates);
    node2->assign_coordinates(coordinates);
    for (unsigned i = 0; i < coordinates.size(); ++i) {
      REQUIRE(node1->coordinates()[i] == Approx(0.).epsilon(Tolerance));
      REQUIRE(node2->coordinates()[i] == Approx(0.).epsilon(Tolerance));
    }

    REQUIRE(node1->status() == false);
    REQUIRE(node2->status() == true);

    mesh->find_active_nodes();

    // Check iterate over functionality if nodes are active
    coordinates.fill(5.3);

    mesh->iterate_over_active_nodes(
        std::bind(&mpm::NodeBase<Dim>::assign_coordinates,
                  std::placeholders::_1, coordinates));

    // Node 1
    {
      // Check if nodal coordinate update has gone through
      auto check_coords = node1->coordinates();
      for (unsigned i = 0; i < check_coords.size(); ++i)
        REQUIRE(check_coords[i] == Approx(0.).epsilon(Tolerance));
    }
    // Node 2
    {
      // Check if nodal coordinate update has gone through
      auto check_coords = node2->coordinates();
      for (unsigned i = 0; i < check_coords.size(); ++i)
        REQUIRE(check_coords[i] == Approx(5.3).epsilon(Tolerance));
    }

    coordinates.fill(7.0);

    // Check iterate over functionality
    mesh->iterate_over_nodes(std::bind(&mpm::NodeBase<Dim>::assign_coordinates,
                                       std::placeholders::_1, coordinates));

    // Node 1
    {
      // Check if nodal coordinate update has gone through
      auto check_coords = node1->coordinates();
      // Check if coordinates for each node is zero
      for (unsigned i = 0; i < check_coords.size(); ++i)
        REQUIRE(check_coords[i] == Approx(7.).epsilon(Tolerance));
    }
    // Node 2
    {
      // Check if nodal coordinate update has gone through
      auto check_coords = node2->coordinates();
      // Check if coordinates for each node is zero
      for (unsigned i = 0; i < check_coords.size(); ++i)
        REQUIRE(check_coords[i] == Approx(7.).epsilon(Tolerance));
    }

    mesh->iterate_over_nodes(std::bind(&mpm::NodeBase<Dim>::update_mass,
                                       std::placeholders::_1, false, 0, 10.0));

#ifdef USE_MPI
    // Get number of MPI ranks
    int mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_size == 1) {
      // Run if there is more than a single MPI task
      // MPI all reduce nodal mass
      mesh->allreduce_nodal_scalar_property(
          std::bind(&mpm::NodeBase<Dim>::mass, std::placeholders::_1, 0),
          std::bind(&mpm::NodeBase<Dim>::update_mass, std::placeholders::_1,
                    false, 0, std::placeholders::_2));
      // MPI all reduce nodal momentum
      mesh->allreduce_nodal_vector_property(
          std::bind(&mpm::NodeBase<Dim>::coordinates, std::placeholders::_1),
          std::bind(&mpm::NodeBase<Dim>::assign_coordinates,
                    std::placeholders::_1, std::placeholders::_2));
    }
#endif
    // Node 1
    {
      // Check mass
      REQUIRE(node1->mass(0) == Approx(10.).epsilon(Tolerance));
      // Check if nodal coordinate update has gone through
      auto check_coords = node1->coordinates();
      // Check if coordinates for each node is zero
      for (unsigned i = 0; i < check_coords.size(); ++i)
        REQUIRE(check_coords[i] == Approx(7.).epsilon(Tolerance));
    }
    // Node 2
    {
      // Check mass
      REQUIRE(node2->mass(0) == Approx(10.).epsilon(Tolerance));
      // Check if nodal coordinate update has gone through
      auto check_coords = node2->coordinates();
      // Check if coordinates for each node is zero
      for (unsigned i = 0; i < check_coords.size(); ++i)
        REQUIRE(check_coords[i] == Approx(7.).epsilon(Tolerance));
    }

    // Remove node 2 and check
    REQUIRE(mesh->remove_node(node2) == true);
    // Check number of nodes in mesh
    REQUIRE(mesh->nnodes() == 1);
  }

  // Check add / remove cell
  SECTION("Check add / remove cell functionality") {
    // Cell 1
    mpm::Index id1 = 0;
    Eigen::Vector3d coords;
    coords.setZero();
    auto cell1 = std::make_shared<mpm::Cell<Dim>>(id1, Nnodes, element);

    // Cell 2
    mpm::Index id2 = 1;
    auto cell2 = std::make_shared<mpm::Cell<Dim>>(id2, Nnodes, element);

    auto mesh = std::make_shared<mpm::Mesh<Dim>>(0);
    // Check mesh is active
    REQUIRE(mesh->status() == false);

    // Add cell 1 and check
    REQUIRE(mesh->add_cell(cell1) == true);
    // Add cell 2 and check
    REQUIRE(mesh->add_cell(cell2) == true);
    // Add cell 2 again and check
    REQUIRE(mesh->add_cell(cell2) == false);

    // Check number of cells in mesh
    REQUIRE(mesh->ncells() == 2);

    // Check iterate over functionality
    mesh->iterate_over_cells(
        std::bind(&mpm::Cell<Dim>::nnodes, std::placeholders::_1));

    // Remove cell 2 and check
    REQUIRE(mesh->remove_cell(cell2) == true);
    // Check number of cells in mesh
    REQUIRE(mesh->ncells() == 1);
  }

  SECTION("Check particle is in cell") {
    // Index
    mpm::Index id1 = 0;

    // Mesh
    auto mesh = std::make_shared<mpm::Mesh<Dim>>(0);
    // Check mesh is inactive (false)
    REQUIRE(mesh->status() == false);

    // Coordinates
    Eigen::Vector3d coords;
    coords.setZero();

    // Define nodes
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

    // Create cell1
    coords.setZero();
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

    // Initialise cell and compute volume
    REQUIRE(cell1->initialise() == true);

    // Generate material points in cell
    auto points = mesh->generate_material_points(1);
    REQUIRE(points.size() == 0);

    // Add cell 1 and check
    REQUIRE(mesh->add_cell(cell1) == true);

    // Generate material points in cell
    points = mesh->generate_material_points(1);
    REQUIRE(points.size() == 1);

    points = mesh->generate_material_points(2);
    REQUIRE(points.size() == 8);

    points = mesh->generate_material_points(3);
    REQUIRE(points.size() == 27);

    // Particle 1
    coords << 1.0, 1.0, 1.0;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle1 =
        std::make_shared<mpm::Particle<Dim, Nphases>>(0, coords);

    // Particle 2
    coords << 1.5, 1.5, 1.5;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle2 =
        std::make_shared<mpm::Particle<Dim, Nphases>>(1, coords);

    // Add particle 1 and check
    REQUIRE(mesh->add_particle(particle1) == true);
    // Add particle 2 and check
    REQUIRE(mesh->add_particle(particle2) == true);

    // Check mesh is active
    REQUIRE(mesh->status() == true);

    // Locate particles in a mesh
    mesh->locate_particles_mesh();

    // Check location of particle 1
    REQUIRE(particle1->cell_id() == 0);
    // Check location of particle 2
    REQUIRE(particle2->cell_id() == 0);
  }

  //! Check create nodes and cells in a mesh
  SECTION("Check create nodes and cells") {
    // Vector of nodal coordinates
    std::vector<Eigen::Matrix<double, Dim, 1>> coordinates;

    // Nodal coordinates
    Eigen::Matrix<double, Dim, 1> node;

    // Cell 0
    // Node 0
    node << 0., 0., 0.;
    coordinates.emplace_back(node);
    // Node 1
    node << 0.5, 0., 0.;
    coordinates.emplace_back(node);
    // Node 2
    node << 0.5, 0.5, 0.;
    coordinates.emplace_back(node);
    // Node 3
    node << 0., 0.5, 0.;
    coordinates.emplace_back(node);
    // Node 4
    node << 0., 0., 0.5;
    coordinates.emplace_back(node);
    // Node 5
    node << 0.5, 0., 0.5;
    coordinates.emplace_back(node);
    // Node 6
    node << 0.5, 0.5, 0.5;
    coordinates.emplace_back(node);
    // Node 7
    node << 0., 0.5, 0.5;
    coordinates.emplace_back(node);

    // Cell 1
    // Node 8
    node << 1.0, 0., 0.;
    coordinates.emplace_back(node);
    // Node 9
    node << 1.0, 0.5, 0.;
    coordinates.emplace_back(node);
    // Node 10
    node << 1.0, 0., 0.5;
    coordinates.emplace_back(node);
    // Node 11
    node << 1.0, 0.5, 0.5;
    coordinates.emplace_back(node);

    // Create a new mesh
    unsigned meshid = 0;
    auto mesh = std::make_shared<mpm::Mesh<Dim>>(meshid);

    SECTION("Check creation of nodes") {
      // Node type 3D
      const std::string node_type = "N3D";
      // Global node index
      mpm::Index gnid = 0;
      mesh->create_nodes(gnid, node_type, coordinates);
      // Check if mesh has added nodes
      REQUIRE(mesh->nnodes() == coordinates.size());
      // Try again this shouldn't add more coordinates
      mesh->create_nodes(gnid, node_type, coordinates);
      // Check if mesh has added nodes
      REQUIRE(mesh->nnodes() == coordinates.size());
      // Clear coordinates and try creating a list of nodes with an empty list
      unsigned nnodes = coordinates.size();
      coordinates.clear();
      // This fails with empty list error in node creation
      mesh->create_nodes(gnid, node_type, coordinates);
      REQUIRE(mesh->nnodes() == nnodes);

      SECTION("Check creation of cells") {
        // Cell with node ids
        std::vector<std::vector<mpm::Index>> cells{// cell #0
                                                   {0, 1, 2, 3, 4, 5, 6, 7},
                                                   // cell #1
                                                   {1, 8, 9, 2, 5, 10, 11, 6}};
        // Assign 8-noded hexahedron element to cell
        std::shared_ptr<mpm::Element<Dim>> element =
            Factory<mpm::Element<Dim>>::instance()->create("ED3H8");

        // Global cell index
        mpm::Index gcid = 0;
        mesh->create_cells(gcid, element, cells);
        // Check if mesh has added cells
        REQUIRE(mesh->ncells() == cells.size());
        // Try again this shouldn't add more cells
        mesh->create_cells(gcid, element, cells);
        // Check if mesh has added cells
        REQUIRE(mesh->ncells() == cells.size());
        // Clear cells and try creating a list of empty cells
        unsigned ncells = cells.size();
        cells.clear();
        // This fails with empty list error in node creation
        gcid = 100;
        mesh->create_cells(gcid, element, cells);
        REQUIRE(mesh->ncells() == ncells);

        // Try with invalid node ids
        cells = {// cell #0
                 {90, 91, 92, 93, 94, 95, 96, 97},
                 // cell #1
                 {71, 88, 89, 82, 85, 80, 81, 86}};
        gcid = 200;
        mesh->create_cells(gcid, element, cells);
        REQUIRE(mesh->ncells() == ncells);

        SECTION("Check creation of particles") {
          // Vector of particle coordinates
          std::vector<Eigen::Matrix<double, Dim, 1>> coordinates;
          coordinates.clear();

          // Particle coordinates
          Eigen::Matrix<double, Dim, 1> particle;

          // Cell 0
          // Particle 0
          particle << 0.125, 0.125, 0.125;
          coordinates.emplace_back(particle);
          // Particle 1
          particle << 0.25, 0.125, 0.125;
          coordinates.emplace_back(particle);
          // Particle 2
          particle << 0.25, 0.25, 0.125;
          coordinates.emplace_back(particle);
          // Particle 3
          particle << 0.125, 0.25, 0.125;
          coordinates.emplace_back(particle);
          // Particle 4
          particle << 0.125, 0.125, 0.25;
          coordinates.emplace_back(particle);
          // Particle 5
          particle << 0.25, 0.125, 0.25;
          coordinates.emplace_back(particle);
          // Particle 6
          particle << 0.25, 0.25, 0.25;
          coordinates.emplace_back(particle);
          // Particle 7
          particle << 0.125, 0.25, 0.25;
          coordinates.emplace_back(particle);

          // Cell 1
          // Particle 8
          particle << 0.675, 0.125, 0.125;
          coordinates.emplace_back(particle);
          // Particle 9
          particle << 0.85, 0.125, 0.125;
          coordinates.emplace_back(particle);
          // Particle 10
          particle << 0.85, 0.25, 0.125;
          coordinates.emplace_back(particle);
          // Particle 11
          particle << 0.675, 0.25, 0.125;
          coordinates.emplace_back(particle);
          // Particle 12
          particle << 0.675, 0.125, 0.25;
          coordinates.emplace_back(particle);
          // Particle 13
          particle << 0.85, 0.125, 0.25;
          coordinates.emplace_back(particle);
          // Particle 14
          particle << 0.85, 0.25, 0.25;
          coordinates.emplace_back(particle);
          // Particle 15
          particle << 0.675, 0.25, 0.25;
          coordinates.emplace_back(particle);

          SECTION("Check addition of particles to mesh") {
            // Particle type 3D
            const std::string particle_type = "P3D";
            // Global particle index
            std::vector<mpm::Index> gpid(coordinates.size());
            std::iota(gpid.begin(), gpid.end(), 0);
            mesh->create_particles(gpid, particle_type, coordinates);
            // Check if mesh has added particles
            REQUIRE(mesh->nparticles() == coordinates.size());
            // Try again this shouldn't add more coordinates
            mesh->create_particles(gpid, particle_type, coordinates);
            // Check if mesh has added particles
            REQUIRE(mesh->nparticles() == coordinates.size());
            // Clear coordinates and try creating a list of particles with an
            // empty list
            unsigned nparticles = coordinates.size();
            coordinates.clear();
            // This fails with empty list error in particle creation
            mesh->create_particles(gpid, particle_type, coordinates);
            REQUIRE(mesh->nparticles() == nparticles);

            // Test assign particles cells again should fail
            SECTION("Check assign particles cells") {
              // Vector of particle cells
              std::vector<std::array<mpm::Index, 2>> particles_cells;
              // Particle cells
              particles_cells.emplace_back(std::array<mpm::Index, 2>({0, 0}));
              particles_cells.emplace_back(std::array<mpm::Index, 2>({1, 0}));
              particles_cells.emplace_back(std::array<mpm::Index, 2>({2, 0}));
              particles_cells.emplace_back(std::array<mpm::Index, 2>(
                  {3, std::numeric_limits<mpm::Index>::max()}));
              particles_cells.emplace_back(std::array<mpm::Index, 2>({50, 0}));

              REQUIRE(mesh->assign_particles_cells(particles_cells) == false);
            }

            const unsigned phase = 0;
            // Particles coordinates
            REQUIRE(mesh->particle_coordinates().size() == mesh->nparticles());
            // Particle stresses
            std::string attribute = "stresses";
            REQUIRE(mesh->particles_vector_data(attribute, phase).size() ==
                    mesh->nparticles());
            // Particle strains
            attribute = "strains";
            REQUIRE(mesh->particles_vector_data(attribute, phase).size() ==
                    mesh->nparticles());
            // Particle velocities
            attribute = "velocities";
            REQUIRE(mesh->particles_vector_data(attribute, phase).size() ==
                    mesh->nparticles());
            // Particle invalid data
            attribute = "invalid";
            REQUIRE(mesh->particles_vector_data(attribute, phase).size() == 0);

            SECTION("Locate particles in mesh") {
              // Locate particles in a mesh
              auto particles = mesh->locate_particles_mesh();

              // Should find all particles in mesh
              REQUIRE(particles.size() == 0);
              // Create particle 100
              Eigen::Vector3d coords;
              coords << 100., 100., 100.;

              mpm::Index pid = 100;
              std::shared_ptr<mpm::ParticleBase<Dim>> particle100 =
                  std::make_shared<mpm::Particle<Dim, Nphases>>(pid, coords);

              // Add particle100 and check
              REQUIRE(mesh->add_particle(particle100) == false);

              // Locate particles in a mesh
              particles = mesh->locate_particles_mesh();
              // Should miss particle100
              REQUIRE(particles.size() == 0);

              SECTION("Check return particles cells") {
                // Vector of particle cells
                std::vector<std::array<mpm::Index, 2>> particles_cells;
                // Particle cells
                particles_cells.emplace_back(std::array<mpm::Index, 2>({0, 0}));
                particles_cells.emplace_back(std::array<mpm::Index, 2>({1, 0}));
                particles_cells.emplace_back(std::array<mpm::Index, 2>({2, 0}));
                particles_cells.emplace_back(std::array<mpm::Index, 2>({3, 0}));
                particles_cells.emplace_back(std::array<mpm::Index, 2>({4, 0}));
                particles_cells.emplace_back(std::array<mpm::Index, 2>({5, 0}));
                particles_cells.emplace_back(std::array<mpm::Index, 2>({6, 0}));
                particles_cells.emplace_back(std::array<mpm::Index, 2>({7, 0}));

                particles_cells.emplace_back(std::array<mpm::Index, 2>({8, 1}));
                particles_cells.emplace_back(std::array<mpm::Index, 2>({9, 1}));
                particles_cells.emplace_back(
                    std::array<mpm::Index, 2>({10, 1}));
                particles_cells.emplace_back(
                    std::array<mpm::Index, 2>({11, 1}));
                particles_cells.emplace_back(
                    std::array<mpm::Index, 2>({12, 1}));
                particles_cells.emplace_back(
                    std::array<mpm::Index, 2>({13, 1}));
                particles_cells.emplace_back(
                    std::array<mpm::Index, 2>({14, 1}));
                particles_cells.emplace_back(
                    std::array<mpm::Index, 2>({15, 1}));

                auto check_particles_cells = mesh->particles_cells();

                REQUIRE(check_particles_cells.size() == mesh->nparticles());

                for (unsigned i = 0; i < particles_cells.size(); ++i)
                  for (unsigned j = 0; j < 2; ++j)
                    REQUIRE(check_particles_cells.at(i).at(j) ==
                            particles_cells.at(i).at(j));
              }
            }
            // Test HDF5
            SECTION("Write particles HDF5") {
              REQUIRE(mesh->write_particles_hdf5(0, "particles-3d.h5") == true);
            }

            // Test assign particles volumes
            SECTION("Check assign particles volumes") {
              // Vector of particle coordinates
              std::vector<std::tuple<mpm::Index, double>> particles_volumes;
              // Volumes
              particles_volumes.emplace_back(std::make_tuple(0, 10.5));
              particles_volumes.emplace_back(std::make_tuple(1, 10.5));

              REQUIRE(mesh->nparticles() == 16);

              REQUIRE(mesh->assign_particles_volumes(particles_volumes) ==
                      true);

              // When volume assignment fails
              particles_volumes.emplace_back(std::make_tuple(2, 0.0));
              particles_volumes.emplace_back(std::make_tuple(3, -10.0));

              REQUIRE(mesh->assign_particles_volumes(particles_volumes) ==
                      false);
            }

            // Test assign particles tractions
            SECTION("Check assign particles tractions") {
              // Vector of particle coordinates
              std::vector<std::tuple<mpm::Index, unsigned, double>>
                  particles_tractions;
              // Tractions
              particles_tractions.emplace_back(std::make_tuple(0, 0, 10.5));
              particles_tractions.emplace_back(std::make_tuple(1, 1, -10.5));
              particles_tractions.emplace_back(std::make_tuple(2, 0, -12.5));
              particles_tractions.emplace_back(std::make_tuple(3, 1, 0.0));

              REQUIRE(mesh->assign_particles_tractions(particles_tractions) ==
                      false);
              // Compute volume
              mesh->iterate_over_particles(
                  std::bind(&mpm::ParticleBase<Dim>::compute_volume,
                            std::placeholders::_1, phase));

              REQUIRE(mesh->assign_particles_tractions(particles_tractions) ==
                      true);
              // When tractions fail
              particles_tractions.emplace_back(std::make_tuple(3, 3, 0.0));
              REQUIRE(mesh->assign_particles_tractions(particles_tractions) ==
                      false);
              particles_tractions.emplace_back(std::make_tuple(300, 0, 0.0));
              REQUIRE(mesh->assign_particles_tractions(particles_tractions) ==
                      false);
            }

            // Test assign nodes tractions
            SECTION("Check assign nodes tractions") {
              // Vector of node coordinates
              std::vector<std::tuple<mpm::Index, unsigned, double>>
                  nodes_tractions;
              // Tractions
              nodes_tractions.emplace_back(std::make_tuple(0, 0, 10.5));
              nodes_tractions.emplace_back(std::make_tuple(1, 1, -10.5));
              nodes_tractions.emplace_back(std::make_tuple(2, 0, -12.5));
              nodes_tractions.emplace_back(std::make_tuple(3, 1, 0.0));

              REQUIRE(mesh->nnodes() == 12);

              REQUIRE(mesh->assign_nodal_tractions(nodes_tractions) == true);
              // When tractions fail
              nodes_tractions.emplace_back(std::make_tuple(3, 4, 0.0));
              REQUIRE(mesh->assign_nodal_tractions(nodes_tractions) == false);
              nodes_tractions.emplace_back(std::make_tuple(300, 0, 0.0));
              REQUIRE(mesh->assign_nodal_tractions(nodes_tractions) == false);
            }

            // Test assign particles stresses
            SECTION("Check assign particles stresses") {
              // Vector of particle stresses
              std::vector<Eigen::Matrix<double, 6, 1>> particles_stresses;

              REQUIRE(mesh->nparticles() == 16);

              // Stresses
              particles_stresses.emplace_back(
                  Eigen::Matrix<double, 6, 1>::Constant(0.0));
              particles_stresses.emplace_back(
                  Eigen::Matrix<double, 6, 1>::Constant(0.1));
              particles_stresses.emplace_back(
                  Eigen::Matrix<double, 6, 1>::Constant(-0.2));
              particles_stresses.emplace_back(
                  Eigen::Matrix<double, 6, 1>::Constant(-0.4));
              particles_stresses.emplace_back(
                  Eigen::Matrix<double, 6, 1>::Constant(0.3));
              particles_stresses.emplace_back(
                  Eigen::Matrix<double, 6, 1>::Constant(0.5));
              particles_stresses.emplace_back(
                  Eigen::Matrix<double, 6, 1>::Constant(-0.6));
              particles_stresses.emplace_back(
                  Eigen::Matrix<double, 6, 1>::Constant(-0.7));
              particles_stresses.emplace_back(
                  Eigen::Matrix<double, 6, 1>::Constant(0.8));
              particles_stresses.emplace_back(
                  Eigen::Matrix<double, 6, 1>::Constant(-0.9));
              particles_stresses.emplace_back(
                  Eigen::Matrix<double, 6, 1>::Constant(0.10));
              particles_stresses.emplace_back(
                  Eigen::Matrix<double, 6, 1>::Constant(0.11));
              particles_stresses.emplace_back(
                  Eigen::Matrix<double, 6, 1>::Constant(-0.12));
              particles_stresses.emplace_back(
                  Eigen::Matrix<double, 6, 1>::Constant(0.13));
              particles_stresses.emplace_back(
                  Eigen::Matrix<double, 6, 1>::Constant(-0.14));
              particles_stresses.emplace_back(
                  Eigen::Matrix<double, 6, 1>::Constant(0.15));

              REQUIRE(mesh->assign_particles_stresses(particles_stresses) ==
                      true);
              // When stresses fail
              particles_stresses.emplace_back(
                  Eigen::Matrix<double, 6, 1>::Constant(0.16));
              REQUIRE(mesh->assign_particles_stresses(particles_stresses) ==
                      false);
              unsigned id = 1;
              auto mesh_fail = std::make_shared<mpm::Mesh<Dim>>(id);
              REQUIRE(mesh_fail->assign_particles_stresses(
                          particles_stresses) == false);
            }

            // Test assign particles velocity constraints
            SECTION("Check assign particles velocity constraints") {
              // Vector of particle coordinates
              std::vector<std::tuple<mpm::Index, unsigned, double>>
                  particles_velocity_constraints;
              // Constraint
              particles_velocity_constraints.emplace_back(
                  std::make_tuple(0, 0, 10.5));
              particles_velocity_constraints.emplace_back(
                  std::make_tuple(1, 1, -10.5));
              particles_velocity_constraints.emplace_back(
                  std::make_tuple(2, 2, -12.5));
              particles_velocity_constraints.emplace_back(
                  std::make_tuple(3, 1, 0.0));

              REQUIRE(mesh->assign_particles_velocity_constraints(
                          particles_velocity_constraints) == true);

              // When constraints fail
              particles_velocity_constraints.emplace_back(
                  std::make_tuple(3, 3, 0.0));
              REQUIRE(mesh->assign_particles_velocity_constraints(
                          particles_velocity_constraints) == false);
            }
          }
        }
        // Test assign velocity constraints to nodes
        SECTION("Check assign velocity constraints to nodes") {
          // Vector of particle coordinates
          std::vector<std::tuple<mpm::Index, unsigned, double>>
              velocity_constraints;
          // Constraint
          velocity_constraints.emplace_back(std::make_tuple(0, 0, 10.5));
          velocity_constraints.emplace_back(std::make_tuple(1, 1, -10.5));
          velocity_constraints.emplace_back(std::make_tuple(2, 2, -12.5));
          velocity_constraints.emplace_back(std::make_tuple(3, 1, 0.0));

          REQUIRE(mesh->assign_velocity_constraints(velocity_constraints) ==
                  true);

          // When constraints fail
          velocity_constraints.emplace_back(std::make_tuple(3, 3, 0.0));
          REQUIRE(mesh->assign_velocity_constraints(velocity_constraints) ==
                  false);
        }

        // Test assign rotation matrices to nodes
        SECTION("Check assign rotation matrices to nodes") {
          // Map of nodal id and euler angles
          std::map<mpm::Index, Eigen::Matrix<double, Dim, 1>> euler_angles;
          // Insert euler angles and node id into map
          // Node 0 with Euler angles of 10, 20 and 30 deg
          euler_angles.emplace(std::make_pair(
              0, (Eigen::Matrix<double, Dim, 1>() << 10. * M_PI / 180,
                  20. * M_PI / 180, 30. * M_PI / 180)
                     .finished()));
          // Node 1 with Euler angles of 40, 50 and 60 deg
          euler_angles.emplace(std::make_pair(
              1, (Eigen::Matrix<double, Dim, 1>() << 40. * M_PI / 180,
                  50. * M_PI / 180, 60. * M_PI / 180)
                     .finished()));
          // Node 2 with Euler angles of 70, 80 and 90 deg
          euler_angles.emplace(std::make_pair(
              2, (Eigen::Matrix<double, Dim, 1>() << 70. * M_PI / 180,
                  80. * M_PI / 180, 90. * M_PI / 180)
                     .finished()));
          // Node 3 with Euler angles of 100, 110 and 120 deg
          euler_angles.emplace(std::make_pair(
              3, (Eigen::Matrix<double, Dim, 1>() << 100. * M_PI / 180,
                  110. * M_PI / 180, 120. * M_PI / 180)
                     .finished()));

          // Check compute and assign rotation matrix
          REQUIRE(mesh->compute_nodal_rotation_matrices(euler_angles) == true);

          // Check for failure when missing node id
          // Node 100 (non-existent) with Euler angles of 130, 140 and 150 deg
          euler_angles.emplace(std::make_pair(
              100, (Eigen::Matrix<double, Dim, 1>() << 130. * M_PI / 180,
                    140. * M_PI / 180, 150. * M_PI / 180)
                       .finished()));
          REQUIRE(mesh->compute_nodal_rotation_matrices(euler_angles) == false);

          // Check for failure of empty input
          std::map<mpm::Index, Eigen::Matrix<double, Dim, 1>>
              empty_euler_angles;
          REQUIRE(mesh->compute_nodal_rotation_matrices(empty_euler_angles) ==
                  false);

          // Check for failure when no nodes are assigned
          auto mesh_fail = std::make_shared<mpm::Mesh<Dim>>(1);
          REQUIRE(mesh_fail->compute_nodal_rotation_matrices(euler_angles) ==
                  false);
        }

        // Test assign velocity constraints to cells
        SECTION("Check assign velocity constraints to cells") {
          // Vector of particle coordinates
          std::vector<std::tuple<mpm::Index, unsigned, unsigned, double>>
              velocity_constraints;
          // Constraint
          velocity_constraints.emplace_back(std::make_tuple(0, 3, 0, 10.5));
          velocity_constraints.emplace_back(std::make_tuple(1, 2, 1, -10.5));

          REQUIRE(mesh->assign_cell_velocity_constraints(
                      velocity_constraints) == true);

          // When constraints fail
          velocity_constraints.emplace_back(std::make_tuple(1, 10, 1, -10.5));
          REQUIRE(mesh->assign_cell_velocity_constraints(
                      velocity_constraints) == false);
        }
      }
    }
  }
}

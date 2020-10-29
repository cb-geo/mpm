#include <cmath>
#include <limits>
#include <memory>

#include <Eigen/Dense>
#include <boost/filesystem.hpp>

#include "catch.hpp"
// MPI
#ifdef USE_MPI
#include "mpi.h"
#endif

#include "constraints.h"
#include "element.h"
#include "function_base.h"
#include "hexahedron_element.h"
#include "linear_function.h"
#include "mesh.h"
#include "node.h"
#include "partio_writer.h"
#include "quadrilateral_element.h"

//! Check mesh class for 2D case
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
  // Json property
  Json jfunctionproperties;
  jfunctionproperties["id"] = 0;
  std::vector<double> x_values{{0.0, 0.5, 1.0}};
  std::vector<double> fx_values{{0.0, 1.0, 1.0}};
  jfunctionproperties["xvalues"] = x_values;
  jfunctionproperties["fxvalues"] = fx_values;

  // math function
  std::shared_ptr<mpm::FunctionBase> mfunction =
      std::make_shared<mpm::LinearFunction>(0, jfunctionproperties);

  // 4-noded quadrilateral element
  std::shared_ptr<mpm::Element<Dim>> element =
      Factory<mpm::Element<Dim>>::instance()->create("ED2Q4");

  // Assign material
  unsigned mid = 0;
  std::vector<unsigned> mids(1, mid);
  // Initialise material
  Json jmaterial;
  jmaterial["density"] = 1000.;
  jmaterial["youngs_modulus"] = 1.0E+7;
  jmaterial["poisson_ratio"] = 0.3;

  auto le_material =
      Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
          "LinearElastic2D", std::move(0), jmaterial);

  std::map<unsigned, std::shared_ptr<mpm::Material<Dim>>> materials;
  materials[mid] = le_material;

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
        std::make_shared<mpm::Particle<Dim>>(id1, coords);

    // Particle 2
    mpm::Index id2 = 1;
    coords << 2., 2.;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle2 =
        std::make_shared<mpm::Particle<Dim>>(id2, coords);

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

    int mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (mpi_size == 1) cell1->rank(1);

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
    REQUIRE(mesh->nparticles("P2D") == 2);

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
    REQUIRE(mesh->nparticles("P2D") == 1);

    // Remove all non-rank particles in mesh
    mesh->remove_all_nonrank_particles();
    // Check number of particles in mesh
    REQUIRE(mesh->nparticles() == 0);
    REQUIRE(mesh->nparticles("P2D") == 0);

    // Add and use remove all particles
    REQUIRE(mesh->add_particle(particle1) == true);
    REQUIRE(mesh->add_particle(particle2) == true);

    // Check number of particles in mesh
    REQUIRE(mesh->nparticles() == 2);
    std::vector<mpm::Index> remove_pids = {{0, 1}};
    // Remove all particles
    mesh->remove_particles(remove_pids);
    // Check number of particles in mesh
    REQUIRE(mesh->nparticles() == 0);

    // Test assign node concentrated force
    SECTION("Check assign node concentrated force") {
      unsigned Nphase = 0;
      // Set external force to zero
      Eigen::Matrix<double, Dim, 1> force;
      force.setZero();
      REQUIRE_NOTHROW(node0->update_external_force(false, Nphase, force));
      REQUIRE_NOTHROW(node1->update_external_force(false, Nphase, force));

      const unsigned Direction = 0;
      // Check external force
      for (unsigned i = 0; i < Dim; ++i) {
        REQUIRE(node0->external_force(Nphase)(i) ==
                Approx(0.).epsilon(Tolerance));
        REQUIRE(node1->external_force(Nphase)(i) ==
                Approx(0.).epsilon(Tolerance));
      }

      tsl::robin_map<mpm::Index, std::vector<mpm::Index>> node_sets;
      node_sets[0] = std::vector<mpm::Index>{0, 1};

      REQUIRE(mesh->create_node_sets(node_sets, true) == true);

      REQUIRE(mesh->assign_nodal_concentrated_forces(mfunction, 0, 0, 10.5) ==
              true);
      REQUIRE(mesh->assign_nodal_concentrated_forces(mfunction, -1, 0, 0.5) ==
              true);
      REQUIRE(mesh->assign_nodal_concentrated_forces(mfunction, 5, 0, 0.5) ==
              false);
      REQUIRE(mesh->assign_nodal_concentrated_forces(mfunction, -5, 1, 0.5) ==
              false);

      double current_time = 0.0;
      node0->apply_concentrated_force(Nphase, current_time);
      node1->apply_concentrated_force(Nphase, current_time);
      // Check external force
      for (unsigned i = 0; i < Dim; ++i) {
        REQUIRE(node0->external_force(Nphase)(i) ==
                Approx(0.).epsilon(Tolerance));
        REQUIRE(node1->external_force(Nphase)(i) ==
                Approx(0.).epsilon(Tolerance));
      }

      current_time = 0.25;
      node0->apply_concentrated_force(Nphase, current_time);
      node1->apply_concentrated_force(Nphase, current_time);
      std::vector<double> ext_forces = {0.25, 0., 0.};
      // Check external force
      for (unsigned i = 0; i < Dim; ++i) {
        REQUIRE(node0->external_force(Nphase)(i) ==
                Approx(ext_forces.at(i)).epsilon(Tolerance));
        REQUIRE(node1->external_force(Nphase)(i) ==
                Approx(ext_forces.at(i)).epsilon(Tolerance));
      }

      current_time = 5.0;
      node0->apply_concentrated_force(Nphase, current_time);
      node1->apply_concentrated_force(Nphase, current_time);
      ext_forces = {0.75, 0., 0.};
      // Check external force
      for (unsigned i = 0; i < Dim; ++i) {
        REQUIRE(node0->external_force(Nphase)(i) ==
                Approx(ext_forces.at(i)).epsilon(Tolerance));
        REQUIRE(node1->external_force(Nphase)(i) ==
                Approx(ext_forces.at(i)).epsilon(Tolerance));
      }
    }
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
    // Particle type 2D
    const std::string particle_type = "P2D";

    // Initialise material models
    mesh->initialise_material_models(materials);

    // Generate material points in cell
    REQUIRE(mesh->nparticles() == 0);

    REQUIRE(mesh->generate_material_points(1, particle_type, mids, -1, 0) ==
            false);
    REQUIRE(mesh->nparticles() == 0);

    // Add cell 1 and check
    REQUIRE(mesh->add_cell(cell1) == true);

    SECTION("Check generating 1 particle / cell") {
      // Generate material points in cell
      REQUIRE(mesh->generate_material_points(1, particle_type, mids, -1, 0) ==
              true);
      REQUIRE(mesh->nparticles() == 1);
    }

    SECTION("Check generating 2 particle / cell") {
      REQUIRE(mesh->generate_material_points(2, particle_type, mids, -1, 0) ==
              true);
      REQUIRE(mesh->nparticles() == 4);
    }

    SECTION("Check generating 3 particle / cell") {
      REQUIRE(mesh->generate_material_points(3, particle_type, mids, -1, 0) ==
              true);
      REQUIRE(mesh->nparticles() == 9);
    }

    SECTION("Check material point generation") {
      // Assign argc and argv to nput arguments of MPM
      int argc = 7;
      char* argv[] = {(char*)"./mpm",   (char*)"-f", (char*)"./",
                      (char*)"-p",      (char*)"8",  (char*)"-i",
                      (char*)"mpm.json"};

      // Create an IO object
      auto io = std::make_shared<mpm::IO>(argc, argv);

      tsl::robin_map<mpm::Index, std::vector<mpm::Index>> cell_sets;
      cell_sets[1] = std::vector<mpm::Index>{0};

      REQUIRE(mesh->create_cell_sets(cell_sets, true) == true);

      REQUIRE(mesh->nparticles() == 0);

      SECTION("Gauss point generation") {
        // Gauss point generation
        Json jgen;
        jgen["type"] = "gauss";
        jgen["material_id"] = {0};
        jgen["cset_id"] = 1;
        jgen["particle_type"] = "P2D";
        jgen["check_duplicates"] = false;
        jgen["nparticles_per_dir"] = 2;
        jgen["pset_id"] = 2;

        // Generate
        REQUIRE(mesh->generate_particles(io, jgen) == true);
        // Number of particles
        REQUIRE(mesh->nparticles() == 4);
      }

      SECTION("Inject points") {
        // Gauss point generation
        Json jgen;
        jgen["type"] = "inject";
        jgen["material_id"] = {0};
        jgen["cset_id"] = 1;
        jgen["particle_type"] = "P2D";
        jgen["check_duplicates"] = false;
        jgen["nparticles_per_dir"] = 2;
        jgen["velocity"] = {0., 0.};
        jgen["duration"] = {0.1, 0.2};

        // Generate
        REQUIRE(mesh->generate_particles(io, jgen) == true);
        // Inject particles
        REQUIRE_NOTHROW(mesh->inject_particles(0.05));
        // Number of particles
        REQUIRE(mesh->nparticles() == 0);
        // Inject particles
        REQUIRE_NOTHROW(mesh->inject_particles(0.25));
        // Number of particles
        REQUIRE(mesh->nparticles() == 0);
        // Inject particles
        REQUIRE_NOTHROW(mesh->inject_particles(0.15));
        // Number of particles
        REQUIRE(mesh->nparticles() == 4);
      }
    }

    // Particle 1
    coords << 1.0, 1.0;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle1 =
        std::make_shared<mpm::Particle<Dim>>(100, coords);

    // Particle 2
    coords << 1.5, 1.5;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle2 =
        std::make_shared<mpm::Particle<Dim>>(101, coords);

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

          // Initialise material models in mesh
          mesh->initialise_material_models(materials);

          SECTION("Check addition of particles to mesh") {
            // Particle type 2D
            const std::string particle_type = "P2D";
            // Create particles from file
            mesh->create_particles(particle_type, coordinates, mids, 0, false);
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
            mesh->create_particles(particle_type, coordinates, mids, 0, false);
            REQUIRE(mesh->nparticles() == nparticles);

            const unsigned phase = 0;
            // Particles coordinates
            REQUIRE(mesh->particle_coordinates().size() == mesh->nparticles());
            // Particle stresses
            std::string attribute = "stresses";
            REQUIRE(mesh->template particles_tensor_data<6>(attribute).size() ==
                    mesh->nparticles());
            // Particle strains
            attribute = "strains";
            REQUIRE(mesh->template particles_tensor_data<6>(attribute).size() ==
                    mesh->nparticles());
            // Particle velocities
            attribute = "velocities";
            REQUIRE(mesh->particles_vector_data(attribute).size() ==
                    mesh->nparticles());
            // Particle mass
            attribute = "mass";
            REQUIRE(mesh->particles_scalar_data(attribute).size() ==
                    mesh->nparticles());

            // Particle invalid data for tensor, vector, scalar
            attribute = "invalid";
            const auto& invalid_tensor =
                mesh->template particles_tensor_data<6>(attribute);
            REQUIRE(invalid_tensor.size() == mesh->nparticles());
            const auto& invalid_vector = mesh->particles_vector_data(attribute);
            REQUIRE(invalid_vector.size() == mesh->nparticles());
            const auto& invalid_scalar = mesh->particles_scalar_data(attribute);
            REQUIRE(invalid_scalar.size() == mesh->nparticles());

            // State variable
            attribute = "pdstrain";
            REQUIRE(mesh->particles_statevars_data(attribute).size() ==
                    mesh->nparticles());

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
                  std::make_shared<mpm::Particle<Dim>>(pid, coords);

              // Add particle100 and check
              REQUIRE(mesh->add_particle(particle100) == false);

              // Locate particles in a mesh
              particles = mesh->locate_particles_mesh();
              // Should miss particle100
              REQUIRE(particles.size() == 0);
            }

            // Test HDF5
            SECTION("Write particles HDF5") {
              REQUIRE(mesh->write_particles_hdf5("particles-2d.h5") == true);

              auto phdf5 = mesh->particles_hdf5();
              REQUIRE(phdf5.size() == mesh->nparticles());

#ifdef USE_PARTIO
              REQUIRE_NOTHROW(mpm::partio::write_particles(
                  "partio-2d.bgeo", mesh->particles_hdf5()));
              // Check if .bgeo exists
              REQUIRE(boost::filesystem::exists("./partio-2d.bgeo") == true);
#endif
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
              tsl::robin_map<mpm::Index, std::vector<mpm::Index>> particle_sets;
              particle_sets[0] = std::vector<mpm::Index>{0};
              particle_sets[1] = std::vector<mpm::Index>{1};
              particle_sets[2] = std::vector<mpm::Index>{2};
              particle_sets[3] = std::vector<mpm::Index>{3};

              REQUIRE(mesh->create_particle_sets(particle_sets, true) == true);

              REQUIRE(mesh->nparticles() == 8);

              REQUIRE(mesh->create_particles_tractions(mfunction, 0, 0, 10.5) ==
                      true);
              REQUIRE(mesh->create_particles_tractions(mfunction, 1, 1,
                                                       -10.5) == true);
              REQUIRE(mesh->create_particles_tractions(mfunction, 2, 0,
                                                       -12.5) == true);
              REQUIRE(mesh->create_particles_tractions(mfunction, 3, 1, 0.5) ==
                      true);

              REQUIRE(mesh->create_particles_tractions(mfunction, -1, 1, 0.5) ==
                      true);
              REQUIRE(mesh->create_particles_tractions(mfunction, 5, 0, 0.5) ==
                      false);
              REQUIRE(mesh->create_particles_tractions(mfunction, -5, 1, 0.5) ==
                      false);

              // Locate particles in a mesh
              auto particles = mesh->locate_particles_mesh();
              REQUIRE(particles.size() == 0);
              mesh->iterate_over_particles(
                  std::bind(&mpm::ParticleBase<Dim>::compute_shapefn,
                            std::placeholders::_1));

              // Compute volume
              mesh->iterate_over_particles(
                  std::bind(&mpm::ParticleBase<Dim>::compute_volume,
                            std::placeholders::_1));

              mesh->apply_traction_on_particles(10);
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
              tsl::robin_map<mpm::Index, std::vector<mpm::Index>> particle_sets;
              particle_sets[0] = std::vector<mpm::Index>{0};
              particle_sets[1] = std::vector<mpm::Index>{1};
              particle_sets[2] = std::vector<mpm::Index>{2};
              particle_sets[3] = std::vector<mpm::Index>{3};

              REQUIRE(mesh->create_particle_sets(particle_sets, true) == true);

              REQUIRE(mesh->nparticles() == 8);

              int set_id = 0;
              int dir = 0;
              double constraint = 10.5;
              // Add velocity constraint to mesh
              auto velocity_constraint =
                  std::make_shared<mpm::VelocityConstraint>(set_id, dir,
                                                            constraint);
              REQUIRE(mesh->create_particle_velocity_constraint(
                          set_id, velocity_constraint) == true);

              // Add velocity constraint to all nodes in mesh
              velocity_constraint = std::make_shared<mpm::VelocityConstraint>(
                  -1, dir, constraint);
              REQUIRE(mesh->create_particle_velocity_constraint(
                          set_id, velocity_constraint) == true);

              // When constraints fail
              dir = 2;
              // Add velocity constraint to mesh
              velocity_constraint = std::make_shared<mpm::VelocityConstraint>(
                  set_id, dir, constraint);
              REQUIRE(mesh->create_particle_velocity_constraint(
                          set_id, velocity_constraint) == false);

              mesh->apply_particle_velocity_constraints();
            }
          }
        }
        // Test assign velocity constraints to nodes
        SECTION("Check assign velocity constraints to nodes") {
          tsl::robin_map<mpm::Index, std::vector<mpm::Index>> node_sets;
          node_sets[0] = std::vector<mpm::Index>{0, 2};
          node_sets[1] = std::vector<mpm::Index>{1, 3};

          REQUIRE(mesh->create_node_sets(node_sets, true) == true);

          //! Constraints object
          auto constraints = std::make_shared<mpm::Constraints<Dim>>(mesh);

          int set_id = 0;
          int dir = 0;
          double constraint = 10.5;
          // Add velocity constraint to mesh
          auto velocity_constraint = std::make_shared<mpm::VelocityConstraint>(
              set_id, dir, constraint);
          REQUIRE(constraints->assign_nodal_velocity_constraint(
                      set_id, velocity_constraint) == true);

          set_id = 1;
          dir = 1;
          constraint = -12.5;
          // Add velocity constraint to mesh
          velocity_constraint = std::make_shared<mpm::VelocityConstraint>(
              set_id, dir, constraint);
          REQUIRE(constraints->assign_nodal_velocity_constraint(
                      set_id, velocity_constraint) == true);

          // Add velocity constraint to all nodes in mesh
          velocity_constraint =
              std::make_shared<mpm::VelocityConstraint>(-1, dir, constraint);
          REQUIRE(constraints->assign_nodal_velocity_constraint(
                      set_id, velocity_constraint) == true);

          // When constraints fail
          dir = 2;
          // Add velocity constraint to mesh
          velocity_constraint = std::make_shared<mpm::VelocityConstraint>(
              set_id, dir, constraint);
          REQUIRE(constraints->assign_nodal_velocity_constraint(
                      set_id, velocity_constraint) == false);
        }

        SECTION("Check assign friction constraints to nodes") {
          tsl::robin_map<mpm::Index, std::vector<mpm::Index>> node_sets;
          node_sets[0] = std::vector<mpm::Index>{0, 2};
          node_sets[1] = std::vector<mpm::Index>{1, 3};

          REQUIRE(mesh->create_node_sets(node_sets, true) == true);

          //! Constraints object
          auto constraints = std::make_shared<mpm::Constraints<Dim>>(mesh);

          int set_id = 0;
          int dir = 0;
          int sign_n = 1;
          double friction = 0.5;
          // Add friction constraint to mesh
          auto friction_constraint = std::make_shared<mpm::FrictionConstraint>(
              set_id, dir, sign_n, friction);
          REQUIRE(constraints->assign_nodal_frictional_constraint(
                      set_id, friction_constraint) == true);

          set_id = 1;
          dir = 1;
          sign_n = -1;
          friction = -0.25;
          // Add friction constraint to mesh
          friction_constraint = std::make_shared<mpm::FrictionConstraint>(
              set_id, dir, sign_n, friction);
          REQUIRE(constraints->assign_nodal_frictional_constraint(
                      set_id, friction_constraint) == true);

          // Add friction constraint to all nodes in mesh
          friction_constraint = std::make_shared<mpm::FrictionConstraint>(
              -1, dir, sign_n, friction);
          REQUIRE(constraints->assign_nodal_frictional_constraint(
                      set_id, friction_constraint) == true);

          // When constraints fail
          dir = 2;
          // Add friction constraint to mesh
          friction_constraint = std::make_shared<mpm::FrictionConstraint>(
              set_id, dir, sign_n, friction);
          REQUIRE(constraints->assign_nodal_frictional_constraint(
                      set_id, friction_constraint) == false);
        }

        SECTION("Check assign pressure constraints to nodes") {
          tsl::robin_map<mpm::Index, std::vector<mpm::Index>> node_sets;
          node_sets[0] = std::vector<mpm::Index>{0, 2};
          node_sets[1] = std::vector<mpm::Index>{1, 3};

          REQUIRE(mesh->create_node_sets(node_sets, true) == true);

          //! Constraints object
          auto constraints = std::make_shared<mpm::Constraints<Dim>>(mesh);

          int set_id = 0;
          double pressure = 500.2;
          // Add pressure constraint to mesh
          REQUIRE(constraints->assign_nodal_pressure_constraint(
                      mfunction, set_id, 0, pressure) == true);
          REQUIRE(constraints->assign_nodal_pressure_constraint(
                      mfunction, set_id, 1, pressure) == true);

          // Add pressure constraint to all nodes in mesh
          REQUIRE(constraints->assign_nodal_pressure_constraint(
                      mfunction, -1, 0, pressure) == true);
          REQUIRE(constraints->assign_nodal_pressure_constraint(
                      mfunction, -1, 1, pressure) == true);
        }

        // Test assign velocity constraints to nodes
        SECTION("Check assign velocity constraints to nodes") {
          // Vector of particle coordinates
          std::vector<std::tuple<mpm::Index, unsigned, double>>
              velocity_constraints;
          //! Constraints object
          auto constraints = std::make_shared<mpm::Constraints<Dim>>(mesh);
          // Constraint
          velocity_constraints.emplace_back(std::make_tuple(0, 0, 10.5));
          velocity_constraints.emplace_back(std::make_tuple(1, 1, -10.5));
          velocity_constraints.emplace_back(std::make_tuple(2, 0, -12.5));
          velocity_constraints.emplace_back(std::make_tuple(3, 1, 0.0));

          REQUIRE(constraints->assign_nodal_velocity_constraints(
                      velocity_constraints) == true);
          // When constraints fail
          velocity_constraints.emplace_back(std::make_tuple(3, 2, 0.0));
          REQUIRE(constraints->assign_nodal_velocity_constraints(
                      velocity_constraints) == false);
        }

        // Test assign friction constraints to nodes
        SECTION("Check assign friction constraints to nodes") {
          // Vector of particle coordinates
          std::vector<std::tuple<mpm::Index, unsigned, int, double>>
              friction_constraints;
          //! Constraints object
          auto constraints = std::make_shared<mpm::Constraints<Dim>>(mesh);
          // Constraint
          friction_constraints.emplace_back(std::make_tuple(0, 0, 1, 0.5));
          friction_constraints.emplace_back(std::make_tuple(1, 1, -1, 0.5));
          friction_constraints.emplace_back(std::make_tuple(2, 0, 1, 0.25));
          friction_constraints.emplace_back(std::make_tuple(3, 1, -1, 0.0));

          REQUIRE(constraints->assign_nodal_friction_constraints(
                      friction_constraints) == true);
          // When constraints fail
          friction_constraints.emplace_back(std::make_tuple(3, 2, -1, 0.0));
          REQUIRE(constraints->assign_nodal_friction_constraints(
                      friction_constraints) == false);
        }

        // Test assign pressure constraints to nodes
        SECTION("Check assign pressure constraints to nodes") {
          // Vector of pressure constraints
          std::vector<std::tuple<mpm::Index, double>> pressure_constraints;
          //! Constraints object
          auto constraints = std::make_shared<mpm::Constraints<Dim>>(mesh);
          // Constraint
          pressure_constraints.emplace_back(std::make_tuple(0, 500.5));
          pressure_constraints.emplace_back(std::make_tuple(1, 210.5));
          pressure_constraints.emplace_back(std::make_tuple(2, 320.2));
          pressure_constraints.emplace_back(std::make_tuple(3, 0.0));

          REQUIRE(constraints->assign_nodal_pressure_constraints(
                      0, pressure_constraints) == true);
          REQUIRE(constraints->assign_nodal_pressure_constraints(
                      1, pressure_constraints) == true);
          // When constraints fail
          REQUIRE(constraints->assign_nodal_pressure_constraints(
                      4, pressure_constraints) == false);

          pressure_constraints.emplace_back(std::make_tuple(100, 0.0));
          REQUIRE(constraints->assign_nodal_pressure_constraints(
                      0, pressure_constraints) == false);
        }

        // Test assign nodes concentrated_forces
        SECTION("Check assign nodes concentrated_forces") {
          // Vector of node coordinates
          std::vector<std::tuple<mpm::Index, unsigned, double>>
              nodes_concentrated_forces;
          // Concentrated_Forces
          nodes_concentrated_forces.emplace_back(std::make_tuple(0, 0, 10.5));
          nodes_concentrated_forces.emplace_back(std::make_tuple(1, 1, -10.5));
          nodes_concentrated_forces.emplace_back(std::make_tuple(2, 0, -12.5));
          nodes_concentrated_forces.emplace_back(std::make_tuple(3, 1, 0.0));

          REQUIRE(mesh->nnodes() == 6);

          REQUIRE(mesh->assign_nodal_concentrated_forces(
                      nodes_concentrated_forces) == true);
          // When concentrated_forces fail
          nodes_concentrated_forces.emplace_back(std::make_tuple(3, 2, 0.0));
          REQUIRE(mesh->assign_nodal_concentrated_forces(
                      nodes_concentrated_forces) == false);
          nodes_concentrated_forces.emplace_back(std::make_tuple(300, 0, 0.0));
          REQUIRE(mesh->assign_nodal_concentrated_forces(
                      nodes_concentrated_forces) == false);
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
      }
    }
  }

  //! Check if nodal properties is initialised
  SECTION("Check nodal properties initialisation") {
    // Create the different meshes
    std::shared_ptr<mpm::Mesh<Dim>> mesh = std::make_shared<mpm::Mesh<Dim>>(0);

    // Define nodes
    Eigen::Vector2d coords;
    coords << 0., 0.;
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

    // Add nodes 0 to 3 to the mesh
    REQUIRE(mesh->add_node(node0) == true);
    REQUIRE(mesh->add_node(node1) == true);
    REQUIRE(mesh->add_node(node2) == true);
    REQUIRE(mesh->add_node(node3) == true);

    // Initialise material models
    mesh->initialise_material_models(materials);

    // Check nodal properties creation
    REQUIRE_NOTHROW(mesh->create_nodal_properties());

    // Check nodal properties initialisation
    REQUIRE_NOTHROW(mesh->initialise_nodal_properties());
  }
}

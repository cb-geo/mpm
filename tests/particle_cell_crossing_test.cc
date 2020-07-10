#include <limits>

#include "catch.hpp"

#include "cell.h"
#include "element.h"
#include "hdf5_particle.h"
#include "hexahedron_element.h"
#include "material.h"
#include "mesh.h"
#include "node.h"
#include "particle.h"
#include "quadrilateral_element.h"

//! \brief Check particle cell crossing for 2D case
TEST_CASE("Particle cell crossing is checked for 2D case",
          "[particle][2D][cellcrossing]") {
  // Dimension
  const unsigned Dim = 2;
  // Degree of freedom
  const unsigned Dof = 2;
  // Number of phases
  const unsigned Nphases = 1;
  // Phase
  const unsigned Phase = 0;
  // Number of nodes per cell
  const unsigned Nnodes = 4;
  // Tolerance
  const double Tolerance = 1.E-7;
  // Time step
  const double dt = 0.0004;

  // Element
  std::shared_ptr<mpm::Element<Dim>> element =
      std::make_shared<mpm::QuadrilateralElement<Dim, 4>>();

  // Create cell
  auto cell0 = std::make_shared<mpm::Cell<Dim>>(0, Nnodes, element);
  // Create cell
  auto cell1 = std::make_shared<mpm::Cell<Dim>>(1, Nnodes, element);

  Eigen::Matrix<double, Dim, 1> coords;

  // Add nodes to cell
  coords << 0.0, 0.0;
  std::shared_ptr<mpm::NodeBase<Dim>> node0 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

  coords << 1.0, 0.0;
  std::shared_ptr<mpm::NodeBase<Dim>> node1 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);

  coords << 1.0, 1.0;
  std::shared_ptr<mpm::NodeBase<Dim>> node2 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);

  coords << 0.0, 1.0;
  std::shared_ptr<mpm::NodeBase<Dim>> node3 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);

  coords << 2.0, 0.0;
  std::shared_ptr<mpm::NodeBase<Dim>> node4 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(4, coords);

  coords << 2.0, 1.0;
  std::shared_ptr<mpm::NodeBase<Dim>> node5 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(5, coords);

  cell0->add_node(0, node0);
  cell0->add_node(1, node1);
  cell0->add_node(2, node2);
  cell0->add_node(3, node3);
  REQUIRE(cell0->nnodes() == 4);

  cell1->add_node(0, node1);
  cell1->add_node(1, node4);
  cell1->add_node(2, node5);
  cell1->add_node(3, node2);
  REQUIRE(cell1->nnodes() == 4);

  // Initialise cell properties
  cell0->initialise();
  cell1->initialise();

  // Check if cell is initialised
  REQUIRE(cell0->is_initialised() == true);
  REQUIRE(cell1->is_initialised() == true);

  // Create mesh
  mpm::Index id = 0;
  auto mesh = std::make_shared<mpm::Mesh<Dim>>(id);
  REQUIRE(mesh->id() == 0);

  // Add nodes and check
  REQUIRE(mesh->add_node(node0) == true);
  REQUIRE(mesh->add_node(node1) == true);
  REQUIRE(mesh->add_node(node2) == true);
  REQUIRE(mesh->add_node(node3) == true);
  REQUIRE(mesh->add_node(node4) == true);
  REQUIRE(mesh->add_node(node5) == true);

  // Add cells and check status
  REQUIRE(mesh->add_cell(cell0) == true);
  REQUIRE(mesh->add_cell(cell1) == true);

  // Add particle
  id = 0;
  coords << 0.25, 0.25;
  auto particle0 = std::make_shared<mpm::Particle<Dim>>(id, coords);

  // Add particle
  id = 1;
  coords << 0.75, 0.25;
  auto particle1 = std::make_shared<mpm::Particle<Dim>>(id, coords);

  // Add particle
  id = 2;
  coords << 0.75, 0.75;
  auto particle2 = std::make_shared<mpm::Particle<Dim>>(id, coords);

  // Add particle
  id = 3;
  coords << 0.25, 0.75;
  auto particle3 = std::make_shared<mpm::Particle<Dim>>(id, coords);

  // Add particles and check status
  REQUIRE(mesh->add_particle(particle0) == true);
  REQUIRE(mesh->add_particle(particle1) == true);
  REQUIRE(mesh->add_particle(particle2) == true);
  REQUIRE(mesh->add_particle(particle3) == true);

  // Locate particles in a mesh
  auto particles = mesh->locate_particles_mesh();

  // Should find all particles in mesh
  REQUIRE(particles.size() == 0);

  // Number of particles in each cell
  REQUIRE(cell0->nparticles() == 4);
  REQUIRE(cell1->nparticles() == 0);

  // Activate nodes
  mesh->iterate_over_cells(
      std::bind(&mpm::Cell<Dim>::activate_nodes, std::placeholders::_1));

  // Assign material
  unsigned mid = 0;
  // Initialise material
  Json jmaterial;
  jmaterial["density"] = 1000.;
  jmaterial["youngs_modulus"] = 1.0E+7;
  jmaterial["poisson_ratio"] = 0.3;

  auto material =
      Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
          "LinearElastic2D", std::move(mid), jmaterial);

  // Iterate over each particle to initialise
  mesh->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Dim>::initialise, std::placeholders::_1));

  // Iterate over each particle to assign material
  mesh->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Dim>::assign_material, std::placeholders::_1,
                material));

  // Compute volume
  mesh->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Dim>::compute_volume, std::placeholders::_1));

  // Compute mass
  mesh->iterate_over_particles([](std::shared_ptr<mpm::ParticleBase<Dim>> ptr) {
    return mpm::particle::compute_mass<Dim>(ptr);
  });

  // Initialise nodes
  mesh->iterate_over_nodes(
      std::bind(&mpm::NodeBase<Dim>::initialise, std::placeholders::_1));

  mesh->iterate_over_cells(
      std::bind(&mpm::Cell<Dim>::activate_nodes, std::placeholders::_1));

  // Apply velocity constraints
  REQUIRE(node0->assign_velocity_constraint(0, 0) == true);
  REQUIRE(node1->assign_velocity_constraint(0, 1000) == true);
  REQUIRE(node2->assign_velocity_constraint(0, 1000) == true);
  REQUIRE(node3->assign_velocity_constraint(0, 0) == true);

  // Iterate over each particle to compute shapefn
  mesh->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Dim>::compute_shapefn, std::placeholders::_1));

  // Assign mass and momentum to nodes
  mesh->iterate_over_particles([](std::shared_ptr<mpm::ParticleBase<Dim>> ptr) {
    return mpm::particle::map_mass_momentum_to_nodes<Dim>(ptr);
  });

  // Iterate over active nodes to compute acceleratation and velocity
  mesh->iterate_over_nodes_predicate(
      std::bind(&mpm::NodeBase<Dim>::compute_acceleration_velocity,
                std::placeholders::_1, Phase, dt),
      std::bind(&mpm::NodeBase<Dim>::status, std::placeholders::_1));

  // Iterate over each particle to compute updated position
  mesh->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Dim>::compute_updated_position,
                std::placeholders::_1, dt, false));

  // Locate particles in a mesh
  particles = mesh->locate_particles_mesh();

  // Should find all particles in mesh
  REQUIRE(particles.size() == 0);

  // Number of particles in each cell
  REQUIRE(cell0->nparticles() == 2);
  REQUIRE(cell1->nparticles() == 2);
}

//! \brief Check particle cell crossing for 3D case
TEST_CASE("Particle cell crossing is checked for 3D case",
          "[particle][3D][cellcrossing]") {
  // Dimension
  const unsigned Dim = 3;
  // Degree of freedom
  const unsigned Dof = 3;
  // Number of phases
  const unsigned Nphases = 1;
  // Phase
  const unsigned Phase = 0;
  // Number of nodes per cell
  const unsigned Nnodes = 8;
  // Tolerance
  const double Tolerance = 1.E-7;
  // Time step
  const double dt = 0.0004;

  // Element
  std::shared_ptr<mpm::Element<Dim>> element =
      std::make_shared<mpm::HexahedronElement<Dim, Nnodes>>();

  // Create cell
  auto cell0 = std::make_shared<mpm::Cell<Dim>>(0, Nnodes, element);
  // Create cell
  auto cell1 = std::make_shared<mpm::Cell<Dim>>(1, Nnodes, element);

  Eigen::Matrix<double, Dim, 1> coords;

  // Add nodes to cell
  // Define nodes
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

  // Cell 1
  coords << 2.0, 0., 0.;
  std::shared_ptr<mpm::NodeBase<Dim>> node8 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(8, coords);

  coords << 2.0, 1.0, 0.;
  std::shared_ptr<mpm::NodeBase<Dim>> node9 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(9, coords);

  coords << 2.0, 0., 1.0;
  std::shared_ptr<mpm::NodeBase<Dim>> node10 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(10, coords);

  coords << 2.0, 1.0, 1.0;
  std::shared_ptr<mpm::NodeBase<Dim>> node11 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(11, coords);

  // Add nodes to cell
  cell1->add_node(0, node1);
  cell1->add_node(1, node8);
  cell1->add_node(2, node9);
  cell1->add_node(3, node2);
  cell1->add_node(4, node5);
  cell1->add_node(5, node10);
  cell1->add_node(6, node11);
  cell1->add_node(7, node6);

  REQUIRE(cell1->nnodes() == 8);

  // Initialise cell properties
  cell0->initialise();
  cell1->initialise();

  // Check if cell is initialised
  REQUIRE(cell0->is_initialised() == true);
  REQUIRE(cell1->is_initialised() == true);

  // Create mesh
  mpm::Index id = 0;
  auto mesh = std::make_shared<mpm::Mesh<Dim>>(id);
  REQUIRE(mesh->id() == 0);

  // Add nodes and check
  REQUIRE(mesh->add_node(node0) == true);
  REQUIRE(mesh->add_node(node1) == true);
  REQUIRE(mesh->add_node(node2) == true);
  REQUIRE(mesh->add_node(node3) == true);
  REQUIRE(mesh->add_node(node4) == true);
  REQUIRE(mesh->add_node(node5) == true);
  REQUIRE(mesh->add_node(node6) == true);
  REQUIRE(mesh->add_node(node7) == true);
  REQUIRE(mesh->add_node(node8) == true);
  REQUIRE(mesh->add_node(node9) == true);
  REQUIRE(mesh->add_node(node10) == true);
  REQUIRE(mesh->add_node(node11) == true);

  // Add cells and check status
  REQUIRE(mesh->add_cell(cell0) == true);
  REQUIRE(mesh->add_cell(cell1) == true);

  // Add particle
  coords << 0.25, 0.25, 0.25;
  auto particle0 = std::make_shared<mpm::Particle<Dim>>(0, coords);

  coords << 0.75, 0.25, 0.25;
  auto particle1 = std::make_shared<mpm::Particle<Dim>>(1, coords);

  coords << 0.25, 0.75, 0.25;
  auto particle2 = std::make_shared<mpm::Particle<Dim>>(2, coords);

  coords << 0.75, 0.75, 0.25;
  auto particle3 = std::make_shared<mpm::Particle<Dim>>(3, coords);

  coords << 0.25, 0.25, 0.75;
  auto particle4 = std::make_shared<mpm::Particle<Dim>>(4, coords);

  coords << 0.25, 0.25, 0.75;
  auto particle5 = std::make_shared<mpm::Particle<Dim>>(5, coords);

  coords << 0.75, 0.75, 0.75;
  auto particle6 = std::make_shared<mpm::Particle<Dim>>(6, coords);

  coords << 0.75, 0.75, 0.75;
  auto particle7 = std::make_shared<mpm::Particle<Dim>>(7, coords);

  // Add particles and check status
  REQUIRE(mesh->add_particle(particle0) == true);
  REQUIRE(mesh->add_particle(particle1) == true);
  REQUIRE(mesh->add_particle(particle2) == true);
  REQUIRE(mesh->add_particle(particle3) == true);
  REQUIRE(mesh->add_particle(particle4) == true);
  REQUIRE(mesh->add_particle(particle5) == true);
  REQUIRE(mesh->add_particle(particle6) == true);
  REQUIRE(mesh->add_particle(particle7) == true);

  // Locate particles in a mesh
  auto particles = mesh->locate_particles_mesh();

  // Should find all particles in mesh
  REQUIRE(particles.size() == 0);

  // Number of particles in each cell
  REQUIRE(cell0->nparticles() == 8);
  REQUIRE(cell1->nparticles() == 0);

  // Activate nodes
  mesh->iterate_over_cells(
      std::bind(&mpm::Cell<Dim>::activate_nodes, std::placeholders::_1));

  // Assign material
  unsigned mid = 0;
  // Initialise material
  Json jmaterial;
  jmaterial["density"] = 1000.;
  jmaterial["youngs_modulus"] = 1.0E+7;
  jmaterial["poisson_ratio"] = 0.3;
  auto material =
      Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
          "LinearElastic3D", std::move(mid), jmaterial);

  // Iterate over each particle to initialise
  mesh->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Dim>::initialise, std::placeholders::_1));

  // Iterate over each particle to assign material
  mesh->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Dim>::assign_material, std::placeholders::_1,
                material));

  // Compute volume
  mesh->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Dim>::compute_volume, std::placeholders::_1));

  // Compute mass
  mesh->iterate_over_particles([](std::shared_ptr<mpm::ParticleBase<Dim>> ptr) {
    return mpm::particle::compute_mass<Dim>(ptr);
  });

  // Initialise nodes
  mesh->iterate_over_nodes(
      std::bind(&mpm::NodeBase<Dim>::initialise, std::placeholders::_1));

  mesh->iterate_over_cells(
      std::bind(&mpm::Cell<Dim>::activate_nodes, std::placeholders::_1));

  // Apply velocity constraints
  REQUIRE(node0->assign_velocity_constraint(0, 0) == true);
  REQUIRE(node1->assign_velocity_constraint(0, 1000) == true);
  REQUIRE(node2->assign_velocity_constraint(0, 1000) == true);
  REQUIRE(node3->assign_velocity_constraint(0, 0) == true);
  REQUIRE(node4->assign_velocity_constraint(0, 0) == true);
  REQUIRE(node5->assign_velocity_constraint(0, 1000) == true);
  REQUIRE(node6->assign_velocity_constraint(0, 1000) == true);
  REQUIRE(node7->assign_velocity_constraint(0, 0) == true);

  // Iterate over each particle to compute shapefn
  mesh->iterate_over_particles(std::bind(
      &mpm::ParticleBase<Dim>::compute_shapefn, std::placeholders::_1));

  // Assign mass and momentum to nodes
  mesh->iterate_over_particles([](std::shared_ptr<mpm::ParticleBase<Dim>> ptr) {
    return mpm::particle::map_mass_momentum_to_nodes<Dim>(ptr);
  });

  // Iterate over active nodes to compute acceleratation and velocity
  mesh->iterate_over_nodes_predicate(
      std::bind(&mpm::NodeBase<Dim>::compute_acceleration_velocity,
                std::placeholders::_1, Phase, dt),
      std::bind(&mpm::NodeBase<Dim>::status, std::placeholders::_1));

  // Iterate over each particle to compute updated position
  mesh->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Dim>::compute_updated_position,
                std::placeholders::_1, dt, false));

  // Locate particles in a mesh
  particles = mesh->locate_particles_mesh();

  // Should find all particles in mesh
  REQUIRE(particles.size() == 0);

  // Number of particles in each cell
  REQUIRE(cell0->nparticles() == 4);
  REQUIRE(cell1->nparticles() == 4);
}

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
#include "contact.h"
#include "contact_friction.h"
#include "element.h"
#include "function_base.h"
#include "hexahedron_element.h"
#include "linear_function.h"
#include "material.h"
#include "mesh.h"
#include "node.h"
#include "partio_writer.h"
#include "quadrilateral_element.h"
#include "stress_update_usf.h"

//! \brief Check stress update 3D case
TEST_CASE("Contact test case", "[contact][friction][3D]") {
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

  // Assign material
  unsigned mid = 0;
  std::vector<unsigned> mids(1, mid);
  // Initialise material
  Json jmaterial;
  jmaterial["density"] = 1000.;
  jmaterial["youngs_modulus"] = 1.0E+7;
  jmaterial["poisson_ratio"] = 0.3;

  auto material =
      Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
          "LinearElastic3D", std::move(0), jmaterial);

  std::map<unsigned, std::shared_ptr<mpm::Material<Dim>>> materials;
  materials[mid] = material;

  // 8-noded hexahedron element
  std::shared_ptr<mpm::Element<Dim>> element =
      Factory<mpm::Element<Dim>>::instance()->create("ED3H8");

  // Particle 1
  mpm::Index id1 = 0;
  Eigen::Vector3d coords;
  coords.setZero();
  std::shared_ptr<mpm::ParticleBase<Dim>> particle1 =
      std::make_shared<mpm::Particle<Dim>>(id1, coords);

  // Particle 2
  mpm::Index id2 = 1;
  coords << 2., 2., 2.;
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

  REQUIRE(particle1->assign_material(material) == true);
  REQUIRE(particle2->assign_material(material) == true);

  // Check mesh is active
  REQUIRE(mesh->status() == true);
  // Check number of particles in mesh
  REQUIRE(mesh->nparticles() == 2);

  REQUIRE_NOTHROW(mesh->locate_particles_mesh());

  SECTION("Check ContactFriction") {

    unsigned phase = 0;

    auto le_material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "LinearElastic3D", std::move(0), jmaterial);

    std::map<unsigned, std::shared_ptr<mpm::Material<Dim>>> materials;
    materials[mid] = le_material;

    auto stress_update =
        std::make_shared<mpm::StressUpdateUSF<Dim>>(mesh, 0.01);

    auto contact = std::make_shared<mpm::ContactFriction<Dim>>(mesh);

    // Initialise material models
    REQUIRE_NOTHROW(mesh->initialise_material_models(materials));

    // Create nodal properties
    REQUIRE_NOTHROW(mesh->create_nodal_properties());
    // Initialise
    REQUIRE_NOTHROW(stress_update->initialise());
    // Contact initialize
    REQUIRE_NOTHROW(contact->initialise());

    // Mass momentum and compute velocity at nodes
    REQUIRE_NOTHROW(stress_update->momentum_nodes(phase));
    // Contact compute forces
    REQUIRE_NOTHROW(contact->compute_contact_forces());
  }
}

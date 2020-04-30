#include <limits>

#include "catch.hpp"

#include "nodal_properties.h"

#include "cell.h"
#include "element.h"
#include "function_base.h"
#include "hdf5_particle.h"
#include "hexahedron_element.h"
#include "linear_function.h"
#include "material.h"
#include "node.h"
#include "particle.h"
#include "quadrilateral_element.h"

//! \brief Check interface functions
TEST_CASE("Interface functions are checked", "[interface]") {
  // Number of materials
  const unsigned Nmaterials = 2;
  // Number of nodes
  const unsigned Nnodes = 4;
  // Number of phases
  const unsigned Nphases = 1;
  // Dimension
  const unsigned Dim = 2;
  // Degrees of freedom
  const unsigned Dof = 2;
  // Tolerance
  const double tolerance = 1.E-7;
  // Initialise a nodal property pool
  auto nodal_properties = std::make_shared<mpm::NodalProperties>();

  // Add masses and momenta to the nodal properties pool
  nodal_properties->create_property("masses", Nnodes, Nmaterials);
  nodal_properties->create_property("momenta", Nnodes * Dim, Nmaterials);

  // Element
  std::shared_ptr<mpm::Element<Dim>> element =
      std::make_shared<mpm::QuadrilateralElement<Dim, 4>>();

  // Create cell
  auto cell = std::make_shared<mpm::Cell<Dim>>(1, Nnodes, element);
  // Add nodes to cell
  Eigen::Vector2d coords;
  coords << 0.0, 0.0;
  std::shared_ptr<mpm::NodeBase<Dim>> node0 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

  coords << 0.5, 0.0;
  std::shared_ptr<mpm::NodeBase<Dim>> node1 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);

  coords << 0.0, 0.5;
  std::shared_ptr<mpm::NodeBase<Dim>> node2 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);

  coords << 0.5, 0.5;
  std::shared_ptr<mpm::NodeBase<Dim>> node3 =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);

  cell->add_node(0, node0);
  cell->add_node(1, node1);
  cell->add_node(2, node3);
  cell->add_node(3, node2);

  // Initialise cell properties
  cell->initialise();

  // Initialise property handle in each node
  REQUIRE_NOTHROW(node0->initialise_property_handle(0, nodal_properties));
  REQUIRE_NOTHROW(node1->initialise_property_handle(1, nodal_properties));
  REQUIRE_NOTHROW(node3->initialise_property_handle(2, nodal_properties));
  REQUIRE_NOTHROW(node2->initialise_property_handle(3, nodal_properties));

  // Create particle
  // Particle coordinates
  coords << 0.1, 0.2;
  auto particle1 = std::make_shared<mpm::Particle<Dim>>(0, coords);
  coords << 0.2, 0.1;
  auto particle2 = std::make_shared<mpm::Particle<Dim>>(1, coords);
  coords << 0.4, 0.4;
  auto particle3 = std::make_shared<mpm::Particle<Dim>>(2, coords);

  // Assign particle to cell
  particle1->assign_cell_id(1);
  particle2->assign_cell_id(1);
  particle3->assign_cell_id(1);
  particle1->assign_cell(cell);
  particle2->assign_cell(cell);
  particle3->assign_cell(cell);

  // Initialise material
  Json jmaterial;
  jmaterial["density"] = 1000.;
  jmaterial["youngs_modulus"] = 1.0E+7;
  jmaterial["poisson_ratio"] = 0.3;

  auto material1 =
      Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
          "LinearElastic2D", std::move(0), jmaterial);
  auto material2 =
      Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
          "LinearElastic2D", std::move(1), jmaterial);

  // Assign materials to particles
  particle1->assign_material(material1);
  particle2->assign_material(material2);
  particle3->assign_material(material2);

  // Assign mass
  particle1->assign_mass(2.0);
  particle2->assign_mass(3.0);
  particle3->assign_mass(0.5);

  // Assign velocity
  Eigen::VectorXd velocity;
  velocity.resize(Dim);
  for (unsigned i = 0; i < velocity.size(); ++i) velocity(i) = i + 1;
  particle1->assign_velocity(velocity);
  particle2->assign_velocity(velocity);
  particle3->assign_velocity(velocity);

  // Compute shape functions
  particle1->compute_shapefn();
  particle2->compute_shapefn();
  particle3->compute_shapefn();

  // Map multimaterial properties from the particles to the nodes
  particle1->map_multimaterial_properties_to_nodes();
  particle2->map_multimaterial_properties_to_nodes();
  particle3->map_multimaterial_properties_to_nodes();

  Eigen::Matrix<double, 4, 2> masses;
  // clang-format off
  masses << 0.96, 1.46,
            0.24, 1.04,
            0.16, 0.56,
            0.64, 0.44;
  // clang-format on
  Eigen::Matrix<double, 8, 2> momenta;
  // clang-format off
  momenta << 0.96,  1.46,
             1.92,  2.92,
             0.24,  1.04,
             0.48,  2.08,
             0.16,  0.56,
             0.32,  1.12,
             0.64,  0.44,
             1.28,  0.88;
  // clang-format on

  for (int i = 0; i < Nnodes; ++i) {
    for (int j = 0; j < Nmaterials; ++j) {
      REQUIRE(nodal_properties->property("masses", i, j, 1)(0, 0) ==
              Approx(masses(i, j)).epsilon(tolerance));
      for (int k = 0; k < Dim; ++k)
        REQUIRE(nodal_properties->property("momenta", i, j, Dim)(k, 0) ==
                Approx(momenta(i * Dim + k, j)).epsilon(tolerance));
    }
  }
}

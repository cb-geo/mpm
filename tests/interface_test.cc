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
  nodal_properties->create_property("change_in_momenta", Nnodes * Dim,
                                    Nmaterials);
  nodal_properties->create_property("displacements", Nnodes * Dim, Nmaterials);
  nodal_properties->create_property("separation_vectors", Nnodes * Dim,
                                    Nmaterials);

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
  Eigen::VectorXd velocity1;
  Eigen::VectorXd velocity2;
  velocity1.resize(Dim);
  velocity2.resize(Dim);
  for (unsigned i = 0; i < velocity1.size(); ++i) {
    velocity1(i) = i + 1;
    velocity2(i) = i - 0.5;
  }
  particle1->assign_velocity(velocity1);
  particle2->assign_velocity(velocity2);
  particle3->assign_velocity(velocity1);

  // Compute shape functions
  particle1->compute_shapefn();
  particle2->compute_shapefn();
  particle3->compute_shapefn();

  // Map multimaterial properties from the particles to the nodes
  particle1->map_multimaterial_mass_momentum_to_nodes();
  particle2->map_multimaterial_mass_momentum_to_nodes();
  particle3->map_multimaterial_mass_momentum_to_nodes();

  // Map masses and momenta from particles to nodes
  particle1->map_mass_momentum_to_nodes();
  particle2->map_mass_momentum_to_nodes();
  particle3->map_mass_momentum_to_nodes();

  // Compute velocities at nodes
  node0->compute_velocity();
  node1->compute_velocity();
  node2->compute_velocity();
  node3->compute_velocity();

  // Append material ids to node
  for (unsigned i = 0; i < 2; ++i) {
    node0->append_material_id(i);
    node1->append_material_id(i);
    node2->append_material_id(i);
    node3->append_material_id(i);
  }

  // Compute multimaterial change in momentum
  node0->compute_multimaterial_change_in_momentum();
  node1->compute_multimaterial_change_in_momentum();
  node2->compute_multimaterial_change_in_momentum();
  node3->compute_multimaterial_change_in_momentum();

  // Compute displacements of next time step with dt = 0.05
  const double dt = 0.05;
  particle1->compute_updated_position(dt, true);
  particle2->compute_updated_position(dt, true);
  particle3->compute_updated_position(dt, true);

  // Map multimaterial displacements to nodes
  particle1->map_multimaterial_displacements_to_nodes();
  particle2->map_multimaterial_displacements_to_nodes();
  particle3->map_multimaterial_displacements_to_nodes();

  // Determine separation vectors
  node0->compute_multimaterial_separation_vector();
  node1->compute_multimaterial_separation_vector();
  node2->compute_multimaterial_separation_vector();
  node3->compute_multimaterial_separation_vector();

  Eigen::Matrix<double, 4, 2> masses;
  // clang-format off
  masses << 0.96, 1.46,
            0.24, 1.04,
            0.16, 0.56,
            0.64, 0.44;
  // clang-format on
  Eigen::Matrix<double, 8, 2> momenta;
  // clang-format off
  momenta << 0.96,  -0.70,
             1.92,   0.76,
             0.24,  -0.40,
             0.48,   0.64,
             0.16,   0.20,
             0.32,   0.76,
             0.64,  -0.10,
             1.28,   0.34;
  // clang-format on

  Eigen::Matrix<double, 8, 2> delta_momenta;
  // clang-format off
  delta_momenta << -0.8568595041322, 0.8568595041322,
                   -0.8568595041322, 0.8568595041322,
                   -0.2700000000000, 0.2700000000000,
                   -0.2700000000000, 0.2700000000000,
                   -0.0800000000000, 0.0800000000000,
                   -0.0800000000000, 0.0800000000000,
                   -0.3200000000000, 0.3200000000000,
                   -0.3200000000000, 0.3200000000000;
  // clang-format on

  Eigen::Matrix<double, 8, 2> displacements;
  // clang-format off
  displacements << 0.01182851239670, 0.00576531189856,
                   0.06182851239669, 0.05576531189856,
                   0.01182851239670, 0.00662746344565,
                   0.06182851239669, 0.05662746344565,
                   0.01182851239670, 0.01337072018890,
                   0.06182851239669, 0.06337072018890,
                   0.01182851239670, 0.00805785123967,
                   0.06182851239669, 0.05805785123967;
  // clang-format on

  Eigen::Matrix<double, 8, 2> separation;
  // clang-format off
  separation << -0.00606320049813,  0.00606320049813,
                -0.00606320049813,  0.00606320049813,
                -0.00520104895105,  0.00520104895105,
                -0.00520104895105,  0.00520104895105,
                 0.00154220779221, -0.00154220779221,
                 0.00154220779221, -0.00154220779221,
                -0.00377066115702,  0.00377066115702,
                -0.00377066115702,  0.00377066115702;
  // clang-format on

  // Check if masses, momenta and velocities were properly mapped and computed
  for (int i = 0; i < Nnodes; ++i) {
    for (int j = 0; j < Nmaterials; ++j) {
      REQUIRE(nodal_properties->property("masses", i, j, 1)(0, 0) ==
              Approx(masses(i, j)).epsilon(tolerance));
      for (int k = 0; k < Dim; ++k) {
        REQUIRE(nodal_properties->property("momenta", i, j, Dim)(k, 0) ==
                Approx(momenta(i * Dim + k, j)).epsilon(tolerance));
        REQUIRE(
            nodal_properties->property("change_in_momenta", i, j, Dim)(k, 0) ==
            Approx(delta_momenta(i * Dim + k, j)).epsilon(tolerance));
        REQUIRE(nodal_properties->property("displacements", i, j, Dim)(k, 0) ==
                Approx(displacements(i * Dim + k, j)).epsilon(tolerance));
        REQUIRE(
            nodal_properties->property("separation_vectors", i, j, Dim)(k, 0) ==
            Approx(separation(i * Dim + k, j)).epsilon(tolerance));
      }
    }
  }
}

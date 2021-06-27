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
  nodal_properties->create_property("domains", Nnodes, Nmaterials);
  nodal_properties->create_property("domain_gradients", Nnodes * Dim,
                                    Nmaterials);
  nodal_properties->create_property("normal_unit_vectors", Nnodes * Dim,
                                    Nmaterials);
  nodal_properties->create_property("internal_forces", Nnodes * Dim,
                                    Nmaterials);
  nodal_properties->create_property("external_forces", Nnodes * Dim,
                                    Nmaterials);

  // Element
  std::shared_ptr<mpm::Element<Dim>> element =
      std::make_shared<mpm::QuadrilateralElement<Dim, 4>>();

  // Create cell
  auto cell = std::make_shared<mpm::Cell<Dim>>(1, Nnodes, element);
  
  // Create nodes
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

  // Create vector of nodes
  std::vector<std::shared_ptr<mpm::NodeBase<Dim>>> nodes;
  nodes.emplace_back(node0);
  nodes.emplace_back(node1);
  nodes.emplace_back(node3);
  nodes.emplace_back(node2);

  // Add nodes to cell
  for (unsigned i = 0; i < Nnodes; ++i) cell->add_node(i, nodes[i]);

  // Initialise cell properties
  cell->initialise();

  // Initialise property handle in each node
  for (unsigned i = 0; i < Nnodes; ++i)
    REQUIRE_NOTHROW(nodes[i]->initialise_property_handle(i, nodal_properties));
  
  // Create particle
  // Particle coordinates
  coords << 0.1, 0.2;
  auto particle1 = std::make_shared<mpm::Particle<Dim>>(0, coords);
  coords << 0.2, 0.1;
  auto particle2 = std::make_shared<mpm::Particle<Dim>>(1, coords);
  coords << 0.4, 0.4;
  auto particle3 = std::make_shared<mpm::Particle<Dim>>(2, coords);

  // Create vector of particles
  std::vector<std::shared_ptr<mpm::Particle<Dim>>> particles;
  particles.emplace_back(particle1);
  particles.emplace_back(particle2);
  particles.emplace_back(particle3);
  int Nparticles = particles.size();

  // Assign particle to cell
  for (unsigned i = 0; i < Nparticles; ++i) {
    particles[i]->assign_cell_id(1);
    particles[i]->assign_cell(cell);
  }
  
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
  particles[0]->assign_material(material1);
  particles[1]->assign_material(material2);
  particles[2]->assign_material(material2);

  // Assign mass
  particles[0]->assign_mass(2.0);
  particles[1]->assign_mass(3.0);
  particles[2]->assign_mass(0.5);

  // Assign volume
  particles[0]->assign_volume(4.0);
  particles[1]->assign_volume(3.0);
  particles[2]->assign_volume(0.5);

  // Assign velocity
  Eigen::VectorXd velocity1;
  Eigen::VectorXd velocity2;
  velocity1.resize(Dim);
  velocity2.resize(Dim);
  for (unsigned i = 0; i < velocity1.size(); ++i) {
    velocity1(i) = i + 1;
    velocity2(i) = i - 0.5;
  }
  particles[0]->assign_velocity(velocity1);
  particles[1]->assign_velocity(velocity2);
  particles[2]->assign_velocity(velocity1);

  // Compute shape functions
  for (unsigned i = 0; i < Nparticles; ++i) particles[i]->compute_shapefn();

  // Append material ids to node
  for (unsigned i = 0; i < 2; ++i)
    for (unsigned j = 0; j < Nnodes; ++j)
      REQUIRE_NOTHROW(nodes[j]->append_material_id(i));

  // Check computation of nodal mass and momentum
  SECTION("Check mass and momentum at nodes") {
    for (unsigned i = 0; i < Nparticles; ++i) {
      // Map masses and momenta from particles to nodes
      particles[i]->map_mass_momentum_to_nodes();
      
      // Map multimaterial properties from the particles to the nodes
      REQUIRE_NOTHROW(particles[i]->map_multimaterial_mass_momentum_to_nodes());
    }

    for (unsigned i = 0; i < Nnodes; ++i) {
      // Compute velocities at nodes
      nodes[i]->compute_velocity();

      // Compute multimaterial change in momentum
      REQUIRE_NOTHROW(nodes[i]->compute_multimaterial_change_in_momentum());
    }

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

    // Check values of mass and momentum at each node
    for (int i = 0; i < Nnodes; ++i) {
      for (int j = 0; j < Nmaterials; ++j) {
        REQUIRE(nodal_properties->property("masses", i, j, 1)(0, 0) ==
                Approx(masses(i, j)).epsilon(tolerance));
        for (int k = 0; k < Dim; ++k) {
          REQUIRE(nodal_properties->property("momenta", i, j, Dim)(k, 0) ==
                  Approx(momenta(i * Dim + k, j)).epsilon(tolerance));
          REQUIRE(nodal_properties->property("change_in_momenta", i, j, Dim)(
                      k, 0) ==
                  Approx(delta_momenta(i * Dim + k, j)).epsilon(tolerance));
        }
      }
    }
  }

  SECTION("Check normal unit vector") {
    for (unsigned i = 0; i < Nparticles; ++i) {
      // Map multimaterial domains (volumes) to nodes
      REQUIRE_NOTHROW(particles[i]->map_multimaterial_domain_to_nodes());

      // Map multimaterial domain gradients to nodes
      REQUIRE_NOTHROW(
          particles[i]->map_multimaterial_domain_gradients_to_nodes());
    }

    // Compute normal unit vectors at nodes
    for (unsigned i = 0; i < Nnodes; ++i)
      REQUIRE_NOTHROW(nodes[i]->compute_multimaterial_normal_unit_vector());

    Eigen::Matrix<double, 4, 2> domains;
    // clang-format off
    domains << 1.92, 1.46,
              0.48, 1.04,
              0.32, 0.56,
              1.28, 0.44;
    // clang-format on

    Eigen::Matrix<double, 8, 2> gradients;
    // clang-format off
    gradients << -4.8, -5.0,
                 -6.4, -3.8,
                  4.8,  5.0,
                 -1.6, -3.2,
                  3.2,  2.0,
                  1.6,  3.2,
                 -3.2, -2.0,
                  6.4,  3.8;
    // clang-format on

    Eigen::Matrix<double, 8, 2> normal;
    // clang-format off
    normal << -0.60000000000000, -0.79616219412310,
              -0.80000000000000, -0.60508326753356,
               0.94868329805051,  0.84227140066151,
              -0.31622776601684, -0.53905369642337,
               0.89442719099992,  0.52999894000318,
               0.44721359549996,  0.84799830400509,
              -0.44721359549996, -0.46574643283262,
               0.89442719099992,  0.88491822238198;
    // clang-format on

    // Check if nodal properties were properly mapped and computed
    for (int i = 0; i < Nnodes; ++i) {
      for (int j = 0; j < Nmaterials; ++j) {
        REQUIRE(nodal_properties->property("domains", i, j, 1)(0, 0) ==
                Approx(domains(i, j)).epsilon(tolerance));
        for (int k = 0; k < Dim; ++k) {
          REQUIRE(
              nodal_properties->property("domain_gradients", i, j, Dim)(k, 0) ==
              Approx(gradients(i * Dim + k, j)).epsilon(tolerance));
          REQUIRE(nodal_properties->property("normal_unit_vectors", i, j, Dim)(
                      k, 0) ==
                  Approx(normal(i * Dim + k, j)).epsilon(tolerance));
        }
        // Check if normal vector are also unit vectors
        REQUIRE(nodal_properties->property("normal_unit_vectors", i, j, Dim)
                    .norm() == Approx(1.0).epsilon(tolerance));
      }
    }
  }

  // Check internal forces
  SECTION("Check internal forces") {}

  // Check external forces
  SECTION("Check external forces") {}
}

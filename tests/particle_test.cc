#include <iostream>
#include <limits>

#include "catch.hpp"

#include "cell.h"
#include "hex_shapefn.h"
#include "material/material.h"
#include "node.h"
#include "particle.h"
#include "quad_shapefn.h"
#include "shapefn.h"

//! \brief Check particle class for 1D case
TEST_CASE("Particle is checked for 1D case", "[particle][1D]") {
  // Dimension
  const unsigned Dim = 1;
  // Dimension
  const unsigned Dof = 1;
  // Phases
  const unsigned Nphases = 1;
  // Phase
  const unsigned Phase = 0;

  // Coordinates
  Eigen::Matrix<double, 1, 1> coords;
  coords.setZero();

  //! Check for id = 0
  SECTION("Particle id is zero") {
    mpm::Index id = 0;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);
    REQUIRE(particle->id() == 0);
    REQUIRE(particle->status() == true);
  }

  SECTION("Particle id is positive") {
    //! Check for id is a positive value
    mpm::Index id = std::numeric_limits<mpm::Index>::max();
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);
    REQUIRE(particle->id() == std::numeric_limits<mpm::Index>::max());
    REQUIRE(particle->status() == true);
  }

  //! Construct with id, coordinates and status
  SECTION("Particle with id, coordinates, and status") {
    mpm::Index id = 0;
    bool status = true;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords, status);
    REQUIRE(particle->id() == 0);
    REQUIRE(particle->status() == true);
    particle->assign_status(false);
    REQUIRE(particle->status() == false);
  }

  //! Test coordinates function
  SECTION("coordinates function is checked") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;

    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);

    // Check for coordinates being zero
    auto coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));
    REQUIRE(coordinates.size() == Dim);

    // Check for negative value of coordinates
    for (unsigned i = 0; i < coordinates.size(); ++i)
      coords(i) = -1. * std::numeric_limits<double>::max();
    particle->assign_coordinates(coords);
    coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    REQUIRE(coordinates.size() == Dim);

    // Check for positive value of coordinates
    for (unsigned i = 0; i < coordinates.size(); ++i)
      coords(i) = std::numeric_limits<double>::max();
    particle->assign_coordinates(coords);
    coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    REQUIRE(coordinates.size() == Dim);
  }

  SECTION("Check particle properties") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);

    // Check mass
    REQUIRE(particle->mass(Phase) == Approx(0.0).epsilon(Tolerance));
    double mass = 100.5;
    particle->assign_mass(Phase, mass);
    REQUIRE(particle->mass(Phase) == Approx(100.5).epsilon(Tolerance));

    // Check stress
    Eigen::Matrix<double, 6, 1> stress;
    for (unsigned i = 0; i < stress.size(); ++i) stress(i) = 17.51;

    for (unsigned i = 0; i < stress.size(); ++i)
      REQUIRE(particle->stress(Phase)(i) == Approx(0.).epsilon(Tolerance));

    particle->assign_stress(Phase, stress);
    for (unsigned i = 0; i < stress.size(); ++i)
      REQUIRE(particle->stress(Phase)(i) == Approx(17.51).epsilon(Tolerance));

    // Check velocity
    Eigen::VectorXd velocity;
    velocity.resize(Dim);
    for (unsigned i = 0; i < velocity.size(); ++i) velocity(i) = 17.51;

    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(particle->velocity(Phase)(i) == Approx(0.).epsilon(Tolerance));

    bool status = particle->assign_velocity(Phase, velocity);
    REQUIRE(status == true);
    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(particle->velocity(Phase)(i) == Approx(17.51).epsilon(Tolerance));
    // Check for incorrect dimension of velocity
    velocity.resize(Dim * 2);
    for (unsigned i = 0; i < velocity.size(); ++i) velocity(i) = 17.51;
    status = particle->assign_velocity(Phase, velocity);
    REQUIRE(status == false);

    // Check momentum
    Eigen::VectorXd momentum;
    momentum.resize(Dim);
    for (unsigned i = 0; i < momentum.size(); ++i) momentum(i) = 17.51;

    for (unsigned i = 0; i < momentum.size(); ++i)
      REQUIRE(particle->momentum(Phase)(i) == Approx(0.).epsilon(Tolerance));

    status = particle->assign_momentum(Phase, momentum);
    REQUIRE(status == true);
    for (unsigned i = 0; i < momentum.size(); ++i)
      REQUIRE(particle->momentum(Phase)(i) == Approx(17.51).epsilon(Tolerance));
    // Check for incorrect dimension of momentum
    momentum.resize(Dim * 2);
    for (unsigned i = 0; i < momentum.size(); ++i) momentum(i) = 17.51;
    status = particle->assign_momentum(Phase, momentum);
    REQUIRE(status == false);

    // Check acceleration
    Eigen::VectorXd acceleration;
    acceleration.resize(Dim);
    for (unsigned i = 0; i < acceleration.size(); ++i) acceleration(i) = 17.51;

    for (unsigned i = 0; i < acceleration.size(); ++i)
      REQUIRE(particle->acceleration(Phase)(i) ==
              Approx(0.).epsilon(Tolerance));

    status = particle->assign_acceleration(Phase, acceleration);
    REQUIRE(status == true);
    for (unsigned i = 0; i < acceleration.size(); ++i)
      REQUIRE(particle->acceleration(Phase)(i) ==
              Approx(17.51).epsilon(Tolerance));
    // Check for incorrect dimension of acceleration
    acceleration.resize(Dim * 2);
    for (unsigned i = 0; i < acceleration.size(); ++i) acceleration(i) = 17.51;
    status = particle->assign_acceleration(Phase, acceleration);
    REQUIRE(status == false);
  }
}

//! \brief Check particle class for 2D case
TEST_CASE("Particle is checked for 2D case", "[particle][2D]") {
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
  // Coordinates
  Eigen::Vector2d coords;
  coords.setZero();

  //! Check for id = 0
  SECTION("Particle id is zero") {
    mpm::Index id = 0;
    auto particle = std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);
    REQUIRE(particle->id() == 0);
  }

  SECTION("Particle id is positive") {
    //! Check for id is a positive value
    mpm::Index id = std::numeric_limits<mpm::Index>::max();
    auto particle = std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);
    REQUIRE(particle->id() == std::numeric_limits<mpm::Index>::max());
  }

  //! Test coordinates function
  SECTION("coordinates function is checked") {
    mpm::Index id = 0;
    auto particle = std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);

    //! Check for coordinates being zero
    auto coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));
    REQUIRE(coordinates.size() == Dim);

    //! Check for negative value of coordinates
    for (unsigned i = 0; i < coordinates.size(); ++i)
      coords(i) = -1. * std::numeric_limits<double>::max();
    particle->assign_coordinates(coords);
    coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    REQUIRE(coordinates.size() == Dim);

    //! Check for positive value of coordinates
    for (unsigned i = 0; i < coordinates.size(); ++i)
      coords(i) = std::numeric_limits<double>::max();
    particle->assign_coordinates(coords);
    coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    REQUIRE(coordinates.size() == Dim);
  }

  SECTION("Add a pointer to a cell to particle") {
    // Add particle
    mpm::Index id = 0;
    coords << 0.75, 0.75;
    auto particle = std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);

    // Check particle coordinates
    auto coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    // Shape function
    std::shared_ptr<mpm::ShapeFn<Dim>> shapefn =
        std::make_shared<mpm::QuadrilateralShapeFn<Dim, 4>>();

    // Create cell
    auto cell = std::make_shared<mpm::Cell<Dim>>(10, Nnodes, shapefn);
    // Add nodes to cell
    coords << 0.5, 0.5;
    std::shared_ptr<mpm::NodeBase<Dim>> node0 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

    coords << 1.5, 0.5;
    std::shared_ptr<mpm::NodeBase<Dim>> node1 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);

    coords << 1.5, 1.5;
    std::shared_ptr<mpm::NodeBase<Dim>> node2 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);

    coords << 0.5, 1.5;
    std::shared_ptr<mpm::NodeBase<Dim>> node3 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);
    cell->add_node(0, node0);
    cell->add_node(1, node1);
    cell->add_node(2, node2);
    cell->add_node(3, node3);
    REQUIRE(cell->nnodes() == 4);

    std::vector<std::shared_ptr<mpm::NodeBase<Dim>>> nodes;
    nodes.emplace_back(node0);
    nodes.emplace_back(node1);
    nodes.emplace_back(node2);
    nodes.emplace_back(node3);

    // Initialise cell properties
    cell->initialise();

    // Add cell to particle
    REQUIRE(cell->status() == false);
    // Check compute shape functions of a particle
    REQUIRE(particle->compute_shapefn() == false);
    // Compute reference location should throw
    REQUIRE(particle->compute_reference_location() == false);
    // Compute volume
    REQUIRE(particle->compute_volume() == false);

    particle->assign_cell(cell);
    REQUIRE(cell->status() == true);
    REQUIRE(particle->cell_id() == 10);

    // Check if cell is initialised
    REQUIRE(cell->is_initialised() == true);

    // Check compute shape functions of a particle
    REQUIRE(particle->compute_shapefn() == true);

    // Assign volume
    particle->assign_volume(2.0);
    // Check volume
    REQUIRE(particle->volume() == Approx(2.0).epsilon(Tolerance));
    // Compute volume
    REQUIRE(particle->compute_volume() == true);
    // Check volume
    REQUIRE(particle->volume() == Approx(1.0).epsilon(Tolerance));

    // Check reference location
    coords << -0.5, -0.5;
    REQUIRE(particle->compute_reference_location() == true);
    auto ref_coordinates = particle->reference_location();
    for (unsigned i = 0; i < ref_coordinates.size(); ++i)
      REQUIRE(ref_coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    // Assign material
    unsigned mid = 0;
    auto material = Factory<mpm::Material, unsigned>::instance()->create(
        "LinearElastic", std::move(mid));

    // Initialise material
    Json jmaterial;
    jmaterial["density"] = 1000.;
    jmaterial["youngs_modulus"] = 1.0E+7;
    jmaterial["poisson_ratio"] = 0.3;

    unsigned phase = 0;
    // Check compute mass before material and volume
    REQUIRE(particle->compute_mass(phase) == false);

    // Assign material properties
    material->properties(jmaterial);
    REQUIRE(particle->assign_material(material) == true);

    // Compute volume
    REQUIRE(particle->compute_volume() == true);

    // Compute mass
    REQUIRE(particle->compute_mass(phase) == true);
    // Mass
    REQUIRE(particle->mass(phase) == Approx(1000.).epsilon(Tolerance));

    // Map particle mass to nodes
    particle->assign_mass(phase, std::numeric_limits<double>::max());
    REQUIRE(particle->map_mass_momentum_to_nodes(phase) == false);

    // Assign mass to nodes
    REQUIRE(particle->compute_reference_location() == true);
    REQUIRE(particle->compute_shapefn() == true);

    // Check velocity
    Eigen::VectorXd velocity;
    velocity.resize(Dim);
    for (unsigned i = 0; i < velocity.size(); ++i) velocity(i) = i;
    REQUIRE(particle->assign_velocity(Phase, velocity) == true);
    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(particle->velocity(Phase)(i) == Approx(i).epsilon(Tolerance));

    REQUIRE(particle->compute_mass(phase) == true);
    REQUIRE(particle->map_mass_momentum_to_nodes(phase) == true);

    // Values of nodal mass
    std::array<double, 4> nodal_mass{562.5, 187.5, 62.5, 187.5};
    // Check nodal mass
    for (unsigned i = 0; i < nodes.size(); ++i)
      REQUIRE(nodes.at(i)->mass(phase) ==
              Approx(nodal_mass.at(i)).epsilon(Tolerance));

    // Compute nodal velocity
    for (const auto node : nodes) node->compute_velocity();

    // Values of nodal momentum
    Eigen::Matrix<double, 4, 2> nodal_momentum;
    // clang-format off
    nodal_momentum << 0., 562.5,
                      0., 187.5,
                      0., 62.5,
                      0., 187.5;
    // clang-format on
    // Check nodal momentum
    for (unsigned i = 0; i < nodal_momentum.rows(); ++i)
      for (unsigned j = 0; j < nodal_momentum.cols(); ++j)
        REQUIRE(nodes.at(i)->momentum(phase)(j) ==
                Approx(nodal_momentum(i, j)).epsilon(Tolerance));

    // Values of nodal velocity
    Eigen::Matrix<double, 4, 2> nodal_velocity;
    // clang-format off
    nodal_velocity << 0., 1.,
                      0., 1.,
                      0., 1.,
                      0., 1.;
    // clang-format on
    // Check nodal velocity
    for (unsigned i = 0; i < nodal_velocity.rows(); ++i)
      for (unsigned j = 0; j < nodal_velocity.cols(); ++j)
        REQUIRE(nodes.at(i)->velocity(phase)(j) ==
                Approx(nodal_velocity(i, j)).epsilon(Tolerance));

    // Set momentum to get non-zero strain
    // clang-format off
    nodal_momentum << 0., 562.5 * 1.,
                      0., 187.5 * 2.,
                      0.,  62.5 * 3.,
                      0., 187.5 * 4.;
    // clang-format on
    for (unsigned i = 0; i < nodes.size(); ++i)
      nodes.at(i)->update_momentum(false, phase, nodal_momentum.row(i));

    // nodal velocity
    // clang-format off
    nodal_velocity << 0., 1.,
                      0., 2.,
                      0., 3.,
                      0., 4.;
    // clang-format on
    // Compute nodal velocity
    for (const auto node : nodes) node->compute_velocity();
    // Check nodal velocity
    for (unsigned i = 0; i < nodal_velocity.rows(); ++i)
      for (unsigned j = 0; j < nodal_velocity.cols(); ++j)
        REQUIRE(nodes.at(i)->velocity(phase)(j) ==
                Approx(nodal_velocity(i, j)).epsilon(Tolerance));

    // Compute strain
    particle->compute_strain(phase, 0.1);
    // Strain
    Eigen::Matrix<double, 6, 1> strain;
    strain << 0., 0.125, 0., 0.025, 0., 0.;
    // Check strains
    for (unsigned i = 0; i < strain.rows(); ++i)
      REQUIRE(particle->strain(phase)(i) ==
              Approx(strain(i)).epsilon(Tolerance));
  }

  SECTION("Check assign material to particle") {
    // Add particle
    mpm::Index id = 0;
    coords << 0.75, 0.75;
    auto particle = std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);

    unsigned mid = 0;
    auto material = Factory<mpm::Material, unsigned>::instance()->create(
        "LinearElastic", std::move(mid));
    REQUIRE(material->id() == 0);

    // Initialise material
    Json jmaterial;
    jmaterial["density"] = 1000.;
    jmaterial["youngs_modulus"] = 1.0E+7;
    jmaterial["poisson_ratio"] = 0.3;

    // Check material status before assigning material property
    REQUIRE(material->status() == false);

    // Check if particle can be assigned a material without properties
    REQUIRE(particle->assign_material(material) == false);
    // Assign material properties
    material->properties(jmaterial);

    // Assign material to particle
    REQUIRE(particle->assign_material(material) == true);
  }

  SECTION("Check particle properties") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);

    // Check mass
    REQUIRE(particle->mass(Phase) == Approx(0.0).epsilon(Tolerance));
    double mass = 100.5;
    particle->assign_mass(Phase, mass);
    REQUIRE(particle->mass(Phase) == Approx(100.5).epsilon(Tolerance));

    // Check stress
    Eigen::Matrix<double, 6, 1> stress;
    for (unsigned i = 0; i < stress.size(); ++i) stress(i) = 17.52;

    for (unsigned i = 0; i < stress.size(); ++i)
      REQUIRE(particle->stress(Phase)(i) == Approx(0.).epsilon(Tolerance));

    particle->assign_stress(Phase, stress);
    for (unsigned i = 0; i < stress.size(); ++i)
      REQUIRE(particle->stress(Phase)(i) == Approx(17.52).epsilon(Tolerance));

    // Check velocity
    Eigen::VectorXd velocity;
    velocity.resize(Dim);
    for (unsigned i = 0; i < velocity.size(); ++i) velocity(i) = 19.745;

    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(particle->velocity(Phase)(i) == Approx(0.).epsilon(Tolerance));

    bool status = particle->assign_velocity(Phase, velocity);
    REQUIRE(status == true);
    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(particle->velocity(Phase)(i) ==
              Approx(19.745).epsilon(Tolerance));
    // Check for incorrect dimension of velocity
    velocity.resize(Dim * 2);
    for (unsigned i = 0; i < velocity.size(); ++i) velocity(i) = 19.745;
    status = particle->assign_velocity(Phase, velocity);
    REQUIRE(status == false);

    // Check momentum
    Eigen::VectorXd momentum;
    momentum.resize(Dim);
    for (unsigned i = 0; i < momentum.size(); ++i) momentum(i) = 19.745;

    for (unsigned i = 0; i < momentum.size(); ++i)
      REQUIRE(particle->momentum(Phase)(i) == Approx(0.).epsilon(Tolerance));

    status = particle->assign_momentum(Phase, momentum);
    REQUIRE(status == true);
    for (unsigned i = 0; i < momentum.size(); ++i)
      REQUIRE(particle->momentum(Phase)(i) ==
              Approx(19.745).epsilon(Tolerance));
    // Check for incorrect dimension of momentum
    momentum.resize(Dim * 2);
    for (unsigned i = 0; i < momentum.size(); ++i) momentum(i) = 19.745;
    status = particle->assign_momentum(Phase, momentum);
    REQUIRE(status == false);

    // Check acceleration
    Eigen::VectorXd acceleration;
    acceleration.resize(Dim);
    for (unsigned i = 0; i < acceleration.size(); ++i) acceleration(i) = 19.745;

    for (unsigned i = 0; i < acceleration.size(); ++i)
      REQUIRE(particle->acceleration(Phase)(i) ==
              Approx(0.).epsilon(Tolerance));

    status = particle->assign_acceleration(Phase, acceleration);
    REQUIRE(status == true);
    for (unsigned i = 0; i < acceleration.size(); ++i)
      REQUIRE(particle->acceleration(Phase)(i) ==
              Approx(19.745).epsilon(Tolerance));
    // Check for incorrect dimension of acceleration
    acceleration.resize(Dim * 2);
    for (unsigned i = 0; i < acceleration.size(); ++i) acceleration(i) = 19.745;
    status = particle->assign_acceleration(Phase, acceleration);
    REQUIRE(status == false);
  }
}

//! \brief Check particle class for 3D case
TEST_CASE("Particle is checked for 3D case", "[particle][3D]") {
  // Dimension
  const unsigned Dim = 3;
  // Dimension
  const unsigned Dof = 6;
  // Nnumber of phases
  const unsigned Nphases = 1;
  // Phase
  const unsigned Phase = 0;
  // Number of nodes per cell
  const unsigned Nnodes = 8;
  // Tolerance
  const double Tolerance = 1.E-7;

  // Coordinates
  Eigen::Vector3d coords;
  coords.setZero();

  //! Check for id = 0
  SECTION("Particle id is zero") {
    mpm::Index id = 0;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);
    REQUIRE(particle->id() == 0);
    REQUIRE(particle->status() == true);
  }

  SECTION("Particle id is positive") {
    //! Check for id is a positive value
    mpm::Index id = std::numeric_limits<mpm::Index>::max();
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);
    REQUIRE(particle->id() == std::numeric_limits<mpm::Index>::max());
    REQUIRE(particle->status() == true);
  }

  //! Construct with id, coordinates and status
  SECTION("Particle with id, coordinates, and status") {
    mpm::Index id = 0;
    bool status = true;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords, status);
    REQUIRE(particle->id() == 0);
    REQUIRE(particle->status() == true);
    particle->assign_status(false);
    REQUIRE(particle->status() == false);
  }

  //! Test coordinates function
  SECTION("coordinates function is checked") {
    mpm::Index id = 0;
    // Create particle
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);

    //! Check for coordinates being zero
    auto coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));
    REQUIRE(coordinates.size() == Dim);

    //! Check for negative value of coordinates
    for (unsigned i = 0; i < coordinates.size(); ++i)
      coords(i) = -1. * std::numeric_limits<double>::max();
    particle->assign_coordinates(coords);
    coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    REQUIRE(coordinates.size() == Dim);

    //! Check for positive value of coordinates
    for (unsigned i = 0; i < coordinates.size(); ++i)
      coords(i) = std::numeric_limits<double>::max();
    particle->assign_coordinates(coords);
    coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    REQUIRE(coordinates.size() == Dim);
  }

  //! Test assign cell pointer to particle
  SECTION("Add a pointer to a cell to particle") {
    // Add particle
    mpm::Index id = 0;
    coords << 1.5, 1.5, 1.5;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);

    // Check particle coordinates
    auto coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    // Assign hexahedron shape function
    std::shared_ptr<mpm::ShapeFn<Dim>> shapefn =
        std::make_shared<mpm::HexahedronShapeFn<Dim, 8>>();

    // Create cell
    auto cell = std::make_shared<mpm::Cell<Dim>>(10, Nnodes, shapefn);
    // Add nodes
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

    std::vector<std::shared_ptr<mpm::NodeBase<Dim>>> nodes;
    nodes.emplace_back(node0);
    nodes.emplace_back(node1);
    nodes.emplace_back(node2);
    nodes.emplace_back(node3);
    nodes.emplace_back(node4);
    nodes.emplace_back(node5);
    nodes.emplace_back(node6);
    nodes.emplace_back(node7);

    cell->add_node(0, node0);
    cell->add_node(1, node1);
    cell->add_node(2, node2);
    cell->add_node(3, node3);
    cell->add_node(4, node4);
    cell->add_node(5, node5);
    cell->add_node(6, node6);
    cell->add_node(7, node7);
    REQUIRE(cell->nnodes() == 8);

    // Initialise cell properties
    cell->initialise();

    // Check if cell is initialised
    REQUIRE(cell->is_initialised() == true);

    // Add cell to particle
    REQUIRE(cell->status() == false);
    // Check compute shape functions of a particle
    REQUIRE(particle->compute_shapefn() == false);
    // Compute reference location should throw
    REQUIRE(particle->compute_reference_location() == false);
    // Compute volume
    REQUIRE(particle->compute_volume() == false);

    particle->assign_cell(cell);
    REQUIRE(cell->status() == true);
    REQUIRE(particle->cell_id() == 10);

    // Check if cell is initialised
    REQUIRE(cell->is_initialised() == true);

    // Check compute shape functions of a particle
    REQUIRE(particle->compute_shapefn() == true);

    // Assign volume
    particle->assign_volume(2.0);
    // Check volume
    REQUIRE(particle->volume() == Approx(2.0).epsilon(Tolerance));
    // Compute volume
    REQUIRE(particle->compute_volume() == true);
    // Check volume
    REQUIRE(particle->volume() == Approx(8.0).epsilon(Tolerance));

    // Check reference location
    coords << 0.5, 0.5, 0.5;
    REQUIRE(particle->compute_reference_location() == true);
    auto ref_coordinates = particle->reference_location();
    for (unsigned i = 0; i < ref_coordinates.size(); ++i)
      REQUIRE(ref_coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    // Assign material
    unsigned mid = 0;
    auto material = Factory<mpm::Material, unsigned>::instance()->create(
        "LinearElastic", std::move(mid));

    // Initialise material
    Json jmaterial;
    jmaterial["density"] = 1000.;
    jmaterial["youngs_modulus"] = 1.0E+7;
    jmaterial["poisson_ratio"] = 0.3;

    unsigned phase = 0;
    // Check compute mass before material and volume
    REQUIRE(particle->compute_mass(phase) == false);

    // Assign material properties
    material->properties(jmaterial);
    REQUIRE(particle->assign_material(material) == true);

    // Compute volume
    REQUIRE(particle->compute_volume() == true);

    // Compute mass
    REQUIRE(particle->compute_mass(phase) == true);
    // Mass
    REQUIRE(particle->mass(phase) == Approx(8000.).epsilon(Tolerance));

    // Map particle mass to nodes
    particle->assign_mass(phase, std::numeric_limits<double>::max());
    REQUIRE(particle->map_mass_momentum_to_nodes(phase) == false);

    // Assign mass to nodes
    REQUIRE(particle->compute_reference_location() == true);
    REQUIRE(particle->compute_shapefn() == true);

    REQUIRE(particle->compute_mass(phase) == true);
    REQUIRE(particle->map_mass_momentum_to_nodes(phase) == true);

    // Values of nodal mass
    std::array<double, 8> nodal_mass{125., 375.,  1125., 375.,
                                     375., 1125., 3375., 1125.};
    // Check nodal mass
    for (unsigned i = 0; i < nodes.size(); ++i)
      REQUIRE(nodes.at(i)->mass(phase) ==
              Approx(nodal_mass.at(i)).epsilon(Tolerance));
  }

  SECTION("Check assign material to particle") {
    // Add particle
    mpm::Index id = 0;
    coords << 0.75, 0.75, 0.75;
    auto particle = std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);

    unsigned mid = 0;
    auto material = Factory<mpm::Material, unsigned>::instance()->create(
        "LinearElastic", std::move(mid));
    REQUIRE(material->id() == 0);

    // Initialise material
    Json jmaterial;
    jmaterial["density"] = 1000.;
    jmaterial["youngs_modulus"] = 1.0E+7;
    jmaterial["poisson_ratio"] = 0.3;

    // Check material status before assigning material property
    REQUIRE(material->status() == false);

    // Check if particle can be assigned a material without properties
    REQUIRE(particle->assign_material(material) == false);
    // Assign material properties
    material->properties(jmaterial);

    // Assign material to particle
    REQUIRE(particle->assign_material(material) == true);
  }

  SECTION("Check particle properties") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim, Nphases>>(id, coords);

    // Check mass
    REQUIRE(particle->mass(Phase) == Approx(0.0).epsilon(Tolerance));
    double mass = 100.5;
    particle->assign_mass(Phase, mass);
    REQUIRE(particle->mass(Phase) == Approx(100.5).epsilon(Tolerance));

    // Check stress
    Eigen::Matrix<double, 6, 1> stress;
    for (unsigned i = 0; i < stress.size(); ++i) stress(i) = 1.;

    for (unsigned i = 0; i < stress.size(); ++i)
      REQUIRE(particle->stress(Phase)(i) == Approx(0.).epsilon(Tolerance));

    particle->assign_stress(Phase, stress);
    for (unsigned i = 0; i < stress.size(); ++i)
      REQUIRE(particle->stress(Phase)(i) == Approx(1.).epsilon(Tolerance));

    // Check velocity
    Eigen::VectorXd velocity;
    velocity.resize(Dim);
    for (unsigned i = 0; i < velocity.size(); ++i) velocity(i) = 17.51;

    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(particle->velocity(Phase)(i) == Approx(0.).epsilon(Tolerance));

    bool status = particle->assign_velocity(Phase, velocity);
    REQUIRE(status == true);
    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(particle->velocity(Phase)(i) == Approx(17.51).epsilon(Tolerance));
    // Check for incorrect dimension of velocity
    velocity.resize(Dim * 2);
    for (unsigned i = 0; i < velocity.size(); ++i) velocity(i) = 17.51;
    status = particle->assign_velocity(Phase, velocity);
    REQUIRE(status == false);

    // Check momentum
    Eigen::VectorXd momentum;
    momentum.resize(Dim);
    for (unsigned i = 0; i < momentum.size(); ++i) momentum(i) = 17.51;

    for (unsigned i = 0; i < momentum.size(); ++i)
      REQUIRE(particle->momentum(Phase)(i) == Approx(0.).epsilon(Tolerance));

    status = particle->assign_momentum(Phase, momentum);
    REQUIRE(status == true);
    for (unsigned i = 0; i < momentum.size(); ++i)
      REQUIRE(particle->momentum(Phase)(i) == Approx(17.51).epsilon(Tolerance));
    // Check for incorrect dimension of momentum
    momentum.resize(Dim * 2);
    for (unsigned i = 0; i < momentum.size(); ++i) momentum(i) = 17.51;
    status = particle->assign_momentum(Phase, momentum);
    REQUIRE(status == false);

    // Check acceleration
    Eigen::VectorXd acceleration;
    acceleration.resize(Dim);
    for (unsigned i = 0; i < acceleration.size(); ++i) acceleration(i) = 17.51;

    for (unsigned i = 0; i < acceleration.size(); ++i)
      REQUIRE(particle->acceleration(Phase)(i) ==
              Approx(0.).epsilon(Tolerance));

    status = particle->assign_acceleration(Phase, acceleration);
    REQUIRE(status == true);
    for (unsigned i = 0; i < acceleration.size(); ++i)
      REQUIRE(particle->acceleration(Phase)(i) ==
              Approx(17.51).epsilon(Tolerance));
    // Check for incorrect dimension of acceleration
    acceleration.resize(Dim * 2);
    for (unsigned i = 0; i < acceleration.size(); ++i) acceleration(i) = 17.51;
    status = particle->assign_acceleration(Phase, acceleration);
    REQUIRE(status == false);
  }
}

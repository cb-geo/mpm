#include <limits>

#include "catch.hpp"

#include "cell.h"
#include "element.h"
#include "hdf5_particle.h"
#include "hexahedron_element.h"
#include "material/material.h"
#include "node.h"
#include "particle.h"
#include "quadrilateral_element.h"

//! \brief Check particle class for 1D case
TEST_CASE("Particle is checked for 1D case", "[particle][1D]") {
  // Dimension
  const unsigned Dim = 1;
  // Dimension
  const unsigned Dof = 1;
  // Number of phases
  const unsigned Nphases = 1;
  // Phase
  const unsigned phase = 0;

  // Coordinates
  Eigen::Matrix<double, 1, 1> coords;
  coords.setZero();

  //! Check for id = 0
  SECTION("Particle id is zero") {
    mpm::Index id = 0;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim>>(id, coords);
    REQUIRE(particle->id() == 0);
    REQUIRE(particle->status() == true);
  }

  SECTION("Particle id is positive") {
    //! Check for id is a positive value
    mpm::Index id = std::numeric_limits<mpm::Index>::max();
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim>>(id, coords);
    REQUIRE(particle->id() == std::numeric_limits<mpm::Index>::max());
    REQUIRE(particle->status() == true);
  }

  //! Construct with id, coordinates and status
  SECTION("Particle with id, coordinates, and status") {
    mpm::Index id = 0;
    bool status = true;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim>>(id, coords, status);
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
        std::make_shared<mpm::Particle<Dim>>(id, coords);

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

  //! Test initialise particle stresses
  SECTION("Particle with initial stress") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    bool status = true;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim>>(id, coords, status);
    Eigen::Matrix<double, 6, 1> stress =
        Eigen::Matrix<double, 6, 1>::Constant(5.7);
    particle->initial_stress(stress);
    REQUIRE(particle->stress().size() == stress.size());
    auto pstress = particle->stress();
    for (unsigned i = 0; i < pstress.size(); ++i)
      REQUIRE(pstress[i] == Approx(stress[i]).epsilon(Tolerance));

    auto pstress_data = particle->vector_data("stresses");
    for (unsigned i = 0; i < pstress_data.size(); ++i)
      REQUIRE(pstress_data[i] == Approx(stress[i]).epsilon(Tolerance));
  }

  //! Test particles velocity constraints
  SECTION("Particle with velocity constraints") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    bool status = true;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim>>(id, coords, status);
    // Apply particles velocity constraints
    REQUIRE(particle->assign_particle_velocity_constraint(0, 10.5) == true);
    // Check out of bounds condition
    REQUIRE(particle->assign_particle_velocity_constraint(1, 0) == false);

    // Apply constraints
    particle->apply_particle_velocity_constraints();

    // Check apply constraints
    REQUIRE(particle->velocity()(0) == Approx(10.5).epsilon(Tolerance));
  }

  SECTION("Check particle properties") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim>>(id, coords);

    // Check mass
    REQUIRE(particle->mass() == Approx(0.0).epsilon(Tolerance));
    double mass = 100.5;
    particle->assign_mass(mass);
    REQUIRE(particle->mass() == Approx(100.5).epsilon(Tolerance));

    // Check stress
    Eigen::Matrix<double, 6, 1> stress;
    for (unsigned i = 0; i < stress.size(); ++i) stress(i) = 17.51;

    for (unsigned i = 0; i < stress.size(); ++i)
      REQUIRE(particle->stress()(i) == Approx(0.).epsilon(Tolerance));

    // Check velocity
    Eigen::VectorXd velocity;
    velocity.resize(Dim);
    for (unsigned i = 0; i < velocity.size(); ++i) velocity(i) = 17.51;

    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(particle->velocity()(i) == Approx(0.).epsilon(Tolerance));

    REQUIRE(particle->assign_velocity(velocity) == true);
    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(particle->velocity()(i) == Approx(17.51).epsilon(Tolerance));

    // Assign volume
    REQUIRE(particle->assign_volume(0.0) == false);
    REQUIRE(particle->assign_volume(-5.0) == false);
    REQUIRE(particle->assign_volume(2.0) == true);
    // Check volume
    REQUIRE(particle->volume() == Approx(2.0).epsilon(Tolerance));
    // Traction
    double traction = 65.32;
    const unsigned Direction = 0;
    // Check traction
    for (unsigned i = 0; i < Dim; ++i)
      REQUIRE(particle->traction()(i) == Approx(0.).epsilon(Tolerance));

    REQUIRE(particle->assign_traction(Direction, traction) == true);

    for (unsigned i = 0; i < Dim; ++i) {
      if (i == Direction)
        REQUIRE(particle->traction()(i) == Approx(traction).epsilon(Tolerance));
      else
        REQUIRE(particle->traction()(i) == Approx(0.).epsilon(Tolerance));
    }

    // Check for incorrect direction
    const unsigned wrong_dir = 4;
    REQUIRE(particle->assign_traction(wrong_dir, traction) == false);

    // Check again to ensure value hasn't been updated
    for (unsigned i = 0; i < Dim; ++i) {
      if (i == Direction)
        REQUIRE(particle->traction()(i) == Approx(traction).epsilon(Tolerance));
      else
        REQUIRE(particle->traction()(i) == Approx(0.).epsilon(Tolerance));
    }
  }

  SECTION("Check initialise particle HDF5") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim>>(id, coords);

    mpm::HDF5Particle h5_particle;
    h5_particle.id = 13;
    h5_particle.mass = 501.5;

    Eigen::Vector3d coords;
    coords << 1., 0., 0.;
    h5_particle.coord_x = coords[0];
    h5_particle.coord_y = coords[1];
    h5_particle.coord_z = coords[2];

    Eigen::Vector3d displacement;
    displacement << 0.01, 0.0, 0.0;
    h5_particle.displacement_x = displacement[0];
    h5_particle.displacement_y = displacement[1];
    h5_particle.displacement_z = displacement[2];

    Eigen::Vector3d lsize;
    lsize << 0.25, 0.0, 0.0;
    h5_particle.nsize_x = lsize[0];
    h5_particle.nsize_y = lsize[1];
    h5_particle.nsize_z = lsize[2];

    Eigen::Vector3d velocity;
    velocity << 1.5, 0., 0.;
    h5_particle.velocity_x = velocity[0];
    h5_particle.velocity_y = velocity[1];
    h5_particle.velocity_z = velocity[2];

    Eigen::Matrix<double, 6, 1> stress;
    stress << 11.5, -12.5, 13.5, 14.5, -15.5, 16.5;
    h5_particle.stress_xx = stress[0];
    h5_particle.stress_yy = stress[1];
    h5_particle.stress_zz = stress[2];
    h5_particle.tau_xy = stress[3];
    h5_particle.tau_yz = stress[4];
    h5_particle.tau_xz = stress[5];

    Eigen::Matrix<double, 6, 1> strain;
    strain << 0.115, -0.125, 0.135, 0.145, -0.155, 0.165;
    h5_particle.strain_xx = strain[0];
    h5_particle.strain_yy = strain[1];
    h5_particle.strain_zz = strain[2];
    h5_particle.gamma_xy = strain[3];
    h5_particle.gamma_yz = strain[4];
    h5_particle.gamma_xz = strain[5];

    h5_particle.epsilon_v = strain.head(Dim).sum();

    h5_particle.status = true;

    h5_particle.cell_id = 1;

    h5_particle.volume = 2.;

    h5_particle.material_id = 1;

    // Reinitialise particle from HDF5 data
    REQUIRE(particle->initialise_particle(h5_particle) == true);

    // Check particle id
    REQUIRE(particle->id() == h5_particle.id);
    // Check particle mass
    REQUIRE(particle->mass() == h5_particle.mass);
    // Check particle volume
    REQUIRE(particle->volume() == h5_particle.volume);
    // Check particle mass density
    REQUIRE(particle->mass_density() == h5_particle.mass / h5_particle.volume);
    // Check particle status
    REQUIRE(particle->status() == h5_particle.status);

    // Check for coordinates
    auto coordinates = particle->coordinates();
    REQUIRE(coordinates.size() == Dim);
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    // Check for displacement
    auto pdisplacement = particle->displacement();
    REQUIRE(pdisplacement.size() == Dim);
    for (unsigned i = 0; i < Dim; ++i)
      REQUIRE(pdisplacement(i) == Approx(displacement(i)).epsilon(Tolerance));

    // Check for size
    auto size = particle->natural_size();
    REQUIRE(size.size() == Dim);
    for (unsigned i = 0; i < size.size(); ++i)
      REQUIRE(size(i) == Approx(lsize(i)).epsilon(Tolerance));

    // Check velocity
    auto pvelocity = particle->velocity();
    REQUIRE(pvelocity.size() == Dim);
    for (unsigned i = 0; i < Dim; ++i)
      REQUIRE(pvelocity(i) == Approx(velocity(i)).epsilon(Tolerance));

    // Check stress
    auto pstress = particle->stress();
    REQUIRE(pstress.size() == stress.size());
    for (unsigned i = 0; i < stress.size(); ++i)
      REQUIRE(pstress(i) == Approx(stress(i)).epsilon(Tolerance));

    // Check strain
    auto pstrain = particle->strain();
    REQUIRE(pstrain.size() == strain.size());
    for (unsigned i = 0; i < strain.size(); ++i)
      REQUIRE(pstrain(i) == Approx(strain(i)).epsilon(Tolerance));

    // Check particle volumetric strain centroid
    REQUIRE(particle->volumetric_strain_centroid() == h5_particle.epsilon_v);

    // Check cell id
    REQUIRE(particle->cell_id() == h5_particle.cell_id);

    // Check material id
    REQUIRE(particle->material_id() == h5_particle.material_id);

    // Write Particle HDF5 data
    const auto h5_test = particle->hdf5();

    REQUIRE(h5_particle.id == h5_test.id);
    REQUIRE(h5_particle.mass == h5_test.mass);

    REQUIRE(h5_particle.coord_x == Approx(h5_test.coord_x).epsilon(Tolerance));
    REQUIRE(h5_particle.coord_y == Approx(h5_test.coord_y).epsilon(Tolerance));
    REQUIRE(h5_particle.coord_z == Approx(h5_test.coord_z).epsilon(Tolerance));

    REQUIRE(h5_particle.displacement_x ==
            Approx(h5_test.displacement_x).epsilon(Tolerance));
    REQUIRE(h5_particle.displacement_y ==
            Approx(h5_test.displacement_y).epsilon(Tolerance));
    REQUIRE(h5_particle.displacement_z ==
            Approx(h5_test.displacement_z).epsilon(Tolerance));

    REQUIRE(h5_particle.nsize_x == h5_test.nsize_x);
    REQUIRE(h5_particle.nsize_y == h5_test.nsize_y);
    REQUIRE(h5_particle.nsize_z == h5_test.nsize_z);

    REQUIRE(h5_particle.velocity_x ==
            Approx(h5_test.velocity_x).epsilon(Tolerance));
    REQUIRE(h5_particle.velocity_y ==
            Approx(h5_test.velocity_y).epsilon(Tolerance));
    REQUIRE(h5_particle.velocity_z ==
            Approx(h5_test.velocity_z).epsilon(Tolerance));

    REQUIRE(h5_particle.stress_xx ==
            Approx(h5_test.stress_xx).epsilon(Tolerance));
    REQUIRE(h5_particle.stress_yy ==
            Approx(h5_test.stress_yy).epsilon(Tolerance));
    REQUIRE(h5_particle.stress_zz ==
            Approx(h5_test.stress_zz).epsilon(Tolerance));
    REQUIRE(h5_particle.tau_xy == Approx(h5_test.tau_xy).epsilon(Tolerance));
    REQUIRE(h5_particle.tau_yz == Approx(h5_test.tau_yz).epsilon(Tolerance));
    REQUIRE(h5_particle.tau_xz == Approx(h5_test.tau_xz).epsilon(Tolerance));

    REQUIRE(h5_particle.strain_xx ==
            Approx(h5_test.strain_xx).epsilon(Tolerance));
    REQUIRE(h5_particle.strain_yy ==
            Approx(h5_test.strain_yy).epsilon(Tolerance));
    REQUIRE(h5_particle.strain_zz ==
            Approx(h5_test.strain_zz).epsilon(Tolerance));
    REQUIRE(h5_particle.gamma_xy ==
            Approx(h5_test.gamma_xy).epsilon(Tolerance));
    REQUIRE(h5_particle.gamma_yz ==
            Approx(h5_test.gamma_yz).epsilon(Tolerance));
    REQUIRE(h5_particle.gamma_xz ==
            Approx(h5_test.gamma_xz).epsilon(Tolerance));

    REQUIRE(h5_particle.epsilon_v ==
            Approx(h5_test.epsilon_v).epsilon(Tolerance));
    REQUIRE(h5_particle.status == h5_test.status);
    REQUIRE(h5_particle.cell_id == h5_test.cell_id);
    REQUIRE(h5_particle.material_id == h5_test.material_id);
  }
}

//! \brief Check particle class for 2D case
TEST_CASE("Particle is checked for 2D case", "[particle][2D]") {
  // Dimension
  const unsigned Dim = 2;
  // Degree of freedom
  const unsigned Dof = 2;
  // Number of nodes per cell
  const unsigned Nnodes = 4;
  // Number of phases
  const unsigned Nphases = 1;
  // Phase
  const unsigned phase = 0;
  // Tolerance
  const double Tolerance = 1.E-7;
  // Coordinates
  Eigen::Vector2d coords;
  coords.setZero();

  //! Check for id = 0
  SECTION("Particle id is zero") {
    mpm::Index id = 0;
    auto particle = std::make_shared<mpm::Particle<Dim>>(id, coords);
    REQUIRE(particle->id() == 0);
  }

  SECTION("Particle id is positive") {
    //! Check for id is a positive value
    mpm::Index id = std::numeric_limits<mpm::Index>::max();
    auto particle = std::make_shared<mpm::Particle<Dim>>(id, coords);
    REQUIRE(particle->id() == std::numeric_limits<mpm::Index>::max());
  }

  //! Test coordinates function
  SECTION("coordinates function is checked") {
    mpm::Index id = 0;
    auto particle = std::make_shared<mpm::Particle<Dim>>(id, coords);

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

  // Test assign cell to a particle
  SECTION("Add a pointer to a cell to particle") {
    // Add particle
    mpm::Index id = 0;
    coords << 0.75, 0.75;
    auto particle = std::make_shared<mpm::Particle<Dim>>(id, coords);

    // Check particle coordinates
    auto coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    // Element
    std::shared_ptr<mpm::Element<Dim>> element =
        std::make_shared<mpm::QuadrilateralElement<Dim, 4>>();

    // Create cell
    auto cell = std::make_shared<mpm::Cell<Dim>>(10, Nnodes, element);
    // Add nodes to cell
    coords << 0.5, 0.5;
    std::shared_ptr<mpm::NodeBase<Dim>> node0 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

    coords << 1.5, 0.5;
    std::shared_ptr<mpm::NodeBase<Dim>> node1 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);

    coords << 0.5, 1.5;
    std::shared_ptr<mpm::NodeBase<Dim>> node2 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);

    coords << 1.5, 1.5;
    std::shared_ptr<mpm::NodeBase<Dim>> node3 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);

    coords << 0.5, 3.0;
    std::shared_ptr<mpm::NodeBase<Dim>> node4 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);

    coords << 1.5, 3.0;
    std::shared_ptr<mpm::NodeBase<Dim>> node5 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);

    cell->add_node(0, node0);
    cell->add_node(1, node1);
    cell->add_node(2, node3);
    cell->add_node(3, node2);
    REQUIRE(cell->nnodes() == 4);

    // Initialise cell properties
    cell->initialise();

    // Check if cell is initialised
    REQUIRE(cell->is_initialised() == true);

    // Add cell to particle
    REQUIRE(cell->status() == false);
    // Check particle cell status
    REQUIRE(particle->cell_ptr() == false);
    // Assign cell id
    REQUIRE(particle->assign_cell_id(10) == true);
    // Require cell id
    REQUIRE(particle->cell_id() == 10);
    // Assign a very large cell id
    REQUIRE(particle->assign_cell_id(std::numeric_limits<mpm::Index>::max()) ==
            false);
    // Require cell id
    REQUIRE(particle->cell_id() == 10);
    // Assign particle to cell
    REQUIRE(particle->assign_cell(cell) == true);
    // Local coordinates
    Eigen::Vector2d xi;
    xi.fill(std::numeric_limits<double>::max());
    // Assign particle to cell
    REQUIRE(particle->assign_cell_xi(cell, xi) == false);
    // Compute reference location
    cell->is_point_in_cell(particle->coordinates(), &xi);
    // Assign particle to cell
    REQUIRE(particle->assign_cell_xi(cell, xi) == true);

    // Assign cell id again
    REQUIRE(particle->assign_cell_id(10) == false);
    // Check particle cell status
    REQUIRE(particle->cell_ptr() == true);
    // Check cell status on addition of particle
    REQUIRE(cell->status() == true);

    // Create cell
    auto cell2 = std::make_shared<mpm::Cell<Dim>>(20, Nnodes, element);

    cell2->add_node(0, node2);
    cell2->add_node(1, node3);
    cell2->add_node(2, node5);
    cell2->add_node(3, node4);
    REQUIRE(cell2->nnodes() == 4);

    // Initialise cell2 properties
    cell2->initialise();

    // Check if cell2 is initialised
    REQUIRE(cell2->is_initialised() == true);

    // Add cell2 to particle
    REQUIRE(cell2->status() == false);
    // Assign particle to cell2
    REQUIRE(particle->assign_cell(cell2) == false);
    // Check cell2 status for failed addition of particle
    REQUIRE(cell2->status() == false);
    // Check cell status because this should not have removed the particle
    REQUIRE(cell->status() == true);

    // Remove assigned cell
    particle->remove_cell();
    REQUIRE(particle->assign_cell(cell) == true);

    // Clear all particle ids
    REQUIRE(cell->nparticles() == 1);
    cell->clear_particle_ids();
    REQUIRE(cell->nparticles() == 0);
  }

  //! Test initialise particle stresses
  SECTION("Particle with initial stress") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    bool status = true;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim>>(id, coords, status);
    Eigen::Matrix<double, 6, 1> stress =
        Eigen::Matrix<double, 6, 1>::Constant(5.7);
    particle->initial_stress(stress);
    REQUIRE(particle->stress().size() == stress.size());
    auto pstress = particle->stress();
    for (unsigned i = 0; i < pstress.size(); ++i)
      REQUIRE(pstress[i] == Approx(stress[i]).epsilon(Tolerance));
  }

  //! Test particles velocity constraints
  SECTION("Particle with velocity constraints") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    bool status = true;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim>>(id, coords, status);

    // Apply particles velocity constraints
    REQUIRE(particle->assign_particle_velocity_constraint(0, 10.5) == true);
    REQUIRE(particle->assign_particle_velocity_constraint(1, -12.5) == true);
    // Check out of bounds condition
    REQUIRE(particle->assign_particle_velocity_constraint(2, 0) == false);

    // Apply constraints
    particle->apply_particle_velocity_constraints();

    // Check apply constraints
    REQUIRE(particle->velocity()(0) == Approx(10.5).epsilon(Tolerance));
    REQUIRE(particle->velocity()(1) == Approx(-12.5).epsilon(Tolerance));
  }

  //! Test particle, cell and node functions
  SECTION("Test particle, cell and node functions") {
    // Add particle
    mpm::Index id = 0;
    coords << 0.75, 0.75;
    auto particle = std::make_shared<mpm::Particle<Dim>>(id, coords);

    // Time-step
    const double dt = 0.1;

    // Check particle coordinates
    auto coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    // Shape function
    std::shared_ptr<mpm::Element<Dim>> element =
        std::make_shared<mpm::QuadrilateralElement<Dim, 4>>();

    // Create cell
    auto cell = std::make_shared<mpm::Cell<Dim>>(10, Nnodes, element);
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
    // Compute updated particle location should fail
    REQUIRE(particle->compute_updated_position(dt) == false);
    // Compute updated particle location from nodal velocity should fail
    REQUIRE(particle->compute_updated_position_velocity(dt) == false);
    // Compute volume
    REQUIRE(particle->compute_volume() == false);
    // Update volume should fail
    REQUIRE(particle->update_volume_strainrate(dt) == false);

    REQUIRE(particle->assign_cell(cell) == true);
    REQUIRE(cell->status() == true);
    REQUIRE(particle->cell_id() == 10);

    // Check if cell is initialised
    REQUIRE(cell->is_initialised() == true);

    // Check compute shape functions of a particle
    REQUIRE(particle->compute_shapefn() == true);

    // Assign volume
    REQUIRE(particle->assign_volume(0.0) == false);
    REQUIRE(particle->assign_volume(-5.0) == false);
    REQUIRE(particle->assign_volume(2.0) == true);
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
    unsigned mid = 1;
    // Initialise material
    Json jmaterial;
    jmaterial["density"] = 1000.;
    jmaterial["youngs_modulus"] = 1.0E+7;
    jmaterial["poisson_ratio"] = 0.3;

    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "LinearElastic2D", std::move(mid), jmaterial);

    // Check compute mass before material and volume
    REQUIRE(particle->compute_mass() == false);

    // Test compute stress before material assignment
    REQUIRE(particle->compute_stress() == false);

    // Test compute internal force before material assignment
    REQUIRE(particle->map_internal_force() == false);

    // Assign material properties
    REQUIRE(particle->assign_material(material) == true);

    // Check material id
    REQUIRE(particle->material_id() == 1);

    // Compute volume
    REQUIRE(particle->compute_volume() == true);

    // Compute mass
    REQUIRE(particle->compute_mass() == true);
    // Mass
    REQUIRE(particle->mass() == Approx(1000.).epsilon(Tolerance));

    // Map particle mass to nodes
    particle->assign_mass(std::numeric_limits<double>::max());
    REQUIRE(particle->map_mass_momentum_to_nodes() == false);

    // Map particle pressure to nodes
    REQUIRE(particle->map_pressure_to_nodes() == false);

    // Assign mass to nodes
    REQUIRE(particle->compute_reference_location() == true);
    REQUIRE(particle->compute_shapefn() == true);

    // Map particle material id to nodes
    particle->append_material_id_to_nodes();
    REQUIRE(node0->material_ids()[0] == 1);
    REQUIRE(node1->material_ids()[0] == 1);
    REQUIRE(node2->material_ids()[0] == 1);
    REQUIRE(node3->material_ids()[0] == 1);

    // Check velocity
    Eigen::VectorXd velocity;
    velocity.resize(Dim);
    for (unsigned i = 0; i < velocity.size(); ++i) velocity(i) = i;
    REQUIRE(particle->assign_velocity(velocity) == true);
    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(particle->velocity()(i) == Approx(i).epsilon(Tolerance));

    REQUIRE(particle->compute_mass() == true);
    REQUIRE(particle->map_mass_momentum_to_nodes() == true);

    REQUIRE(particle->map_pressure_to_nodes() == true);
    REQUIRE(particle->compute_pressure_smoothing() == true);

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

    // Check pressure
    REQUIRE(particle->pressure() == Approx(0.).epsilon(Tolerance));

    // Compute strain
    particle->compute_strain(dt);
    // Strain
    Eigen::Matrix<double, 6, 1> strain;
    strain << 0., 0.25, 0., 0.050, 0., 0.;
    // Check strains
    for (unsigned i = 0; i < strain.rows(); ++i)
      REQUIRE(particle->strain()(i) == Approx(strain(i)).epsilon(Tolerance));

    // Check volumetric strain at centroid
    const double volumetric_strain = 0.2;
    REQUIRE(particle->volumetric_strain_centroid() ==
            Approx(volumetric_strain).epsilon(Tolerance));

    // Check updated pressure
    const double K = 8333333.333333333;
    REQUIRE(particle->pressure() ==
            Approx(-K * volumetric_strain).epsilon(Tolerance));

    // Update volume strain rate
    REQUIRE(particle->volume() == Approx(1.0).epsilon(Tolerance));
    REQUIRE(particle->update_volume_strainrate(dt) == true);
    REQUIRE(particle->volume() == Approx(1.2).epsilon(Tolerance));

    // Compute stress
    REQUIRE(particle->compute_stress() == true);

    Eigen::Matrix<double, 6, 1> stress;
    // clang-format off
    stress <<  721153.8461538460 * 2.,
              1682692.3076923075 * 2.,
               721153.8461538460 * 2.,
                96153.8461538462 * 2.,
                    0.0000000000 * 2.,
                    0.0000000000 * 2.;
    // clang-format on
    // Check stress
    for (unsigned i = 0; i < stress.rows(); ++i)
      REQUIRE(particle->stress()(i) == Approx(stress(i)).epsilon(Tolerance));

    // Check body force
    Eigen::Matrix<double, 2, 1> gravity;
    gravity << 0., -9.81;

    particle->map_body_force(gravity);

    // Body force
    Eigen::Matrix<double, 4, 2> body_force;
    // clang-format off
    body_force << 0., -5518.125,
                  0., -1839.375,
                  0.,  -613.125,
                  0., -1839.375;
    // clang-format on

    // Check nodal body force
    for (unsigned i = 0; i < body_force.rows(); ++i)
      for (unsigned j = 0; j < body_force.cols(); ++j)
        REQUIRE(nodes[i]->external_force(phase)[j] ==
                Approx(body_force(i, j)).epsilon(Tolerance));

    // Check traction force
    double traction = 7.68;
    const unsigned direction = 1;
    // Assign volume
    REQUIRE(particle->assign_volume(0.0) == false);
    REQUIRE(particle->assign_volume(-5.0) == false);
    REQUIRE(particle->assign_volume(2.0) == true);
    // Assign traction to particle
    particle->assign_traction(direction, traction);
    // Map traction force
    particle->map_traction_force();

    // Traction force
    Eigen::Matrix<double, 4, 2> traction_force;
    // shapefn * volume / size_(dir) * traction
    // clang-format off
    traction_force << 0., 0.5625 * 1.414213562 * 7.68,
                      0., 0.1875 * 1.414213562 * 7.68,
                      0., 0.0625 * 1.414213562 * 7.68,
                      0., 0.1875 * 1.414213562 * 7.68;
    // clang-format on
    // Add previous external body force
    traction_force += body_force;

    // Check nodal traction force
    for (unsigned i = 0; i < traction_force.rows(); ++i)
      for (unsigned j = 0; j < traction_force.cols(); ++j)
        REQUIRE(nodes[i]->external_force(phase)[j] ==
                Approx(traction_force(i, j)).epsilon(Tolerance));
    // Reset traction
    particle->assign_traction(direction, -traction);
    // Map traction force
    particle->map_traction_force();
    // Check nodal external force
    for (unsigned i = 0; i < traction_force.rows(); ++i)
      for (unsigned j = 0; j < traction_force.cols(); ++j)
        REQUIRE(nodes[i]->external_force(phase)[j] ==
                Approx(body_force(i, j)).epsilon(Tolerance));

    // Internal force
    Eigen::Matrix<double, 4, 2> internal_force;
    // clang-format off
    internal_force <<  1225961.538461538,  2668269.23076923,
                      -1033653.846153846,  697115.3846153845,
                      -408653.8461538461, -889423.0769230769,
                       216346.1538461538, -2475961.538461538;
    // clang-format on

    // Map particle internal force
    particle->assign_volume(1.0);
    REQUIRE(particle->map_internal_force() == true);

    // Check nodal internal force
    for (unsigned i = 0; i < internal_force.rows(); ++i)
      for (unsigned j = 0; j < internal_force.cols(); ++j)
        REQUIRE(nodes[i]->internal_force(phase)[j] ==
                Approx(internal_force(i, j)).epsilon(Tolerance));

    // Calculate nodal acceleration and velocity
    for (const auto& node : nodes)
      node->compute_acceleration_velocity(phase, dt);

    // Check nodal velocity
    // clang-format off
    nodal_velocity <<  217.9487179487179,  474.3779743589742,
                      -551.2820512820512,  372.8138717948718,
                      -653.8461538461538, -1421.057923076923,
                       115.3846153846153, -1317.49382051282;
    // clang-format on
    // Check nodal velocity
    for (unsigned i = 0; i < nodal_velocity.rows(); ++i)
      for (unsigned j = 0; j < nodal_velocity.cols(); ++j)
        REQUIRE(nodes[i]->velocity(phase)[j] ==
                Approx(nodal_velocity(i, j)).epsilon(Tolerance));

    // Check nodal acceleration
    Eigen::Matrix<double, 4, 2> nodal_acceleration;
    // clang-format off
    nodal_acceleration <<  2179.487179487179, 4733.779743589742,
                          -5512.820512820512, 3708.138717948717,
                          -6538.461538461537, -14240.57923076923,
                           1153.846153846153, -13214.9382051282;
    // clang-format on
    // Check nodal acceleration
    for (unsigned i = 0; i < nodal_acceleration.rows(); ++i)
      for (unsigned j = 0; j < nodal_acceleration.cols(); ++j)
        REQUIRE(nodes[i]->acceleration(phase)[j] ==
                Approx(nodal_acceleration(i, j)).epsilon(Tolerance));
    // Approx(nodal_velocity(i, j) / dt).epsilon(Tolerance));

    // Check original particle coordinates
    coords << 0.75, 0.75;
    coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    // Compute updated particle location
    REQUIRE(particle->compute_updated_position(dt) == true);
    // Check particle velocity
    velocity << 0., 0.019;
    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(particle->velocity()(i) ==
              Approx(velocity(i)).epsilon(Tolerance));

    // Check particle displacement
    Eigen::Vector2d displacement;
    displacement << 0., 0.0894;
    for (unsigned i = 0; i < displacement.size(); ++i)
      REQUIRE(particle->displacement()(i) ==
              Approx(displacement(i)).epsilon(Tolerance));

    // Updated particle coordinate
    coords << 0.75, .8394;
    // Check particle coordinates
    coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    // Compute updated particle location from nodal velocity
    REQUIRE(particle->compute_updated_position_velocity(dt) == true);
    // Check particle velocity
    velocity << 0., 0.894;
    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(particle->velocity()(i) ==
              Approx(velocity(i)).epsilon(Tolerance));

    // Check particle displacement
    displacement << 0., 0.1788;
    for (unsigned i = 0; i < displacement.size(); ++i)
      REQUIRE(particle->displacement()(i) ==
              Approx(displacement(i)).epsilon(Tolerance));

    // Updated particle coordinate
    coords << 0.75, .9288;
    // Check particle coordinates
    coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));
  }

  SECTION("Check assign material to particle") {
    // Add particle
    mpm::Index id = 0;
    coords << 0.75, 0.75;
    auto particle = std::make_shared<mpm::Particle<Dim>>(id, coords);

    unsigned mid = 1;
    // Initialise material
    Json jmaterial;
    jmaterial["density"] = 1000.;
    jmaterial["youngs_modulus"] = 1.0E+7;
    jmaterial["poisson_ratio"] = 0.3;

    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "LinearElastic2D", std::move(mid), jmaterial);
    REQUIRE(material->id() == 1);

    // Check if particle can be assigned a material is null
    REQUIRE(particle->assign_material(nullptr) == false);

    // Check material id
    REQUIRE(particle->material_id() == std::numeric_limits<unsigned>::max());

    // Assign material to particle
    REQUIRE(particle->assign_material(material) == true);

    // Check material id
    REQUIRE(particle->material_id() == 1);
  }

  SECTION("Check particle properties") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim>>(id, coords);

    // Check mass
    REQUIRE(particle->mass() == Approx(0.0).epsilon(Tolerance));
    double mass = 100.5;
    particle->assign_mass(mass);
    REQUIRE(particle->mass() == Approx(100.5).epsilon(Tolerance));

    // Check stress
    Eigen::Matrix<double, 6, 1> stress;
    for (unsigned i = 0; i < stress.size(); ++i) stress(i) = 17.52;

    for (unsigned i = 0; i < stress.size(); ++i)
      REQUIRE(particle->stress()(i) == Approx(0.).epsilon(Tolerance));

    // Check velocity
    Eigen::VectorXd velocity;
    velocity.resize(Dim);
    for (unsigned i = 0; i < velocity.size(); ++i) velocity(i) = 19.745;

    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(particle->velocity()(i) == Approx(0.).epsilon(Tolerance));

    REQUIRE(particle->assign_velocity(velocity) == true);
    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(particle->velocity()(i) == Approx(19.745).epsilon(Tolerance));

    // Assign volume
    REQUIRE(particle->assign_volume(0.0) == false);
    REQUIRE(particle->assign_volume(-5.0) == false);
    REQUIRE(particle->assign_volume(2.0) == true);
    // Check volume
    REQUIRE(particle->volume() == Approx(2.0).epsilon(Tolerance));
    // Traction
    double traction = 65.32;
    const unsigned Direction = 1;
    // Check traction
    for (unsigned i = 0; i < Dim; ++i)
      REQUIRE(particle->traction()(i) == Approx(0.).epsilon(Tolerance));

    REQUIRE(particle->assign_traction(Direction, traction) == true);

    // Calculate traction force = traction * volume / spacing
    traction *= 2.0 / (std::pow(2.0, 1. / Dim));

    for (unsigned i = 0; i < Dim; ++i) {
      if (i == Direction)
        REQUIRE(particle->traction()(i) == Approx(traction).epsilon(Tolerance));
      else
        REQUIRE(particle->traction()(i) == Approx(0.).epsilon(Tolerance));
    }

    // Check for incorrect direction
    const unsigned wrong_dir = 4;
    REQUIRE(particle->assign_traction(wrong_dir, traction) == false);

    // Check again to ensure value hasn't been updated
    for (unsigned i = 0; i < Dim; ++i) {
      if (i == Direction)
        REQUIRE(particle->traction()(i) == Approx(traction).epsilon(Tolerance));
      else
        REQUIRE(particle->traction()(i) == Approx(0.).epsilon(Tolerance));
    }
  }

  // Check initialise particle from HDF5 file
  SECTION("Check initialise particle HDF5") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim>>(id, coords);

    mpm::HDF5Particle h5_particle;
    h5_particle.id = 13;
    h5_particle.mass = 501.5;

    Eigen::Vector3d coords;
    coords << 1., 2., 0.;
    h5_particle.coord_x = coords[0];
    h5_particle.coord_y = coords[1];
    h5_particle.coord_z = coords[2];

    Eigen::Vector3d displacement;
    displacement << 0.01, 0.02, 0.0;
    h5_particle.displacement_x = displacement[0];
    h5_particle.displacement_y = displacement[1];
    h5_particle.displacement_z = displacement[2];

    Eigen::Vector3d lsize;
    lsize << 0.25, 0.5, 0.;
    h5_particle.nsize_x = lsize[0];
    h5_particle.nsize_y = lsize[1];
    h5_particle.nsize_z = lsize[2];

    Eigen::Vector3d velocity;
    velocity << 1.5, 2.5, 0.0;
    h5_particle.velocity_x = velocity[0];
    h5_particle.velocity_y = velocity[1];
    h5_particle.velocity_z = velocity[2];

    Eigen::Matrix<double, 6, 1> stress;
    stress << 11.5, -12.5, 13.5, 14.5, -15.5, 16.5;
    h5_particle.stress_xx = stress[0];
    h5_particle.stress_yy = stress[1];
    h5_particle.stress_zz = stress[2];
    h5_particle.tau_xy = stress[3];
    h5_particle.tau_yz = stress[4];
    h5_particle.tau_xz = stress[5];

    Eigen::Matrix<double, 6, 1> strain;
    strain << 0.115, -0.125, 0.135, 0.145, -0.155, 0.165;
    h5_particle.strain_xx = strain[0];
    h5_particle.strain_yy = strain[1];
    h5_particle.strain_zz = strain[2];
    h5_particle.gamma_xy = strain[3];
    h5_particle.gamma_yz = strain[4];
    h5_particle.gamma_xz = strain[5];

    h5_particle.epsilon_v = strain.head(Dim).sum();

    h5_particle.status = true;

    h5_particle.cell_id = 1;

    h5_particle.volume = 2.;

    h5_particle.material_id = 1;

    // Reinitialise particle from HDF5 data
    REQUIRE(particle->initialise_particle(h5_particle) == true);

    // Check particle id
    REQUIRE(particle->id() == h5_particle.id);
    // Check particle mass
    REQUIRE(particle->mass() == h5_particle.mass);
    // Check particle volume
    REQUIRE(particle->volume() == h5_particle.volume);
    // Check particle mass density
    REQUIRE(particle->mass_density() == h5_particle.mass / h5_particle.volume);
    // Check particle status
    REQUIRE(particle->status() == h5_particle.status);

    // Check for coordinates
    auto coordinates = particle->coordinates();
    REQUIRE(coordinates.size() == Dim);
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    // Check for displacement
    auto pdisplacement = particle->displacement();
    REQUIRE(pdisplacement.size() == Dim);
    for (unsigned i = 0; i < Dim; ++i)
      REQUIRE(pdisplacement(i) == Approx(displacement(i)).epsilon(Tolerance));

    // Check for size
    auto size = particle->natural_size();
    REQUIRE(size.size() == Dim);
    for (unsigned i = 0; i < size.size(); ++i)
      REQUIRE(size(i) == Approx(lsize(i)).epsilon(Tolerance));

    // Check velocity
    auto pvelocity = particle->velocity();
    REQUIRE(pvelocity.size() == Dim);
    for (unsigned i = 0; i < Dim; ++i)
      REQUIRE(pvelocity(i) == Approx(velocity(i)).epsilon(Tolerance));

    // Check stress
    auto pstress = particle->stress();
    REQUIRE(pstress.size() == stress.size());
    for (unsigned i = 0; i < stress.size(); ++i)
      REQUIRE(pstress(i) == Approx(stress(i)).epsilon(Tolerance));

    // Check strain
    auto pstrain = particle->strain();
    REQUIRE(pstrain.size() == strain.size());
    for (unsigned i = 0; i < strain.size(); ++i)
      REQUIRE(pstrain(i) == Approx(strain(i)).epsilon(Tolerance));

    // Check particle volumetric strain centroid
    REQUIRE(particle->volumetric_strain_centroid() == h5_particle.epsilon_v);

    // Check cell id
    REQUIRE(particle->cell_id() == h5_particle.cell_id);

    // Check material id
    REQUIRE(particle->material_id() == h5_particle.material_id);

    // Write Particle HDF5 data
    const auto h5_test = particle->hdf5();

    REQUIRE(h5_particle.id == h5_test.id);
    REQUIRE(h5_particle.mass == h5_test.mass);

    REQUIRE(h5_particle.coord_x == Approx(h5_test.coord_x).epsilon(Tolerance));
    REQUIRE(h5_particle.coord_y == Approx(h5_test.coord_y).epsilon(Tolerance));
    REQUIRE(h5_particle.coord_z == Approx(h5_test.coord_z).epsilon(Tolerance));

    REQUIRE(h5_particle.displacement_x ==
            Approx(h5_test.displacement_x).epsilon(Tolerance));
    REQUIRE(h5_particle.displacement_y ==
            Approx(h5_test.displacement_y).epsilon(Tolerance));
    REQUIRE(h5_particle.displacement_z ==
            Approx(h5_test.displacement_z).epsilon(Tolerance));

    REQUIRE(h5_particle.nsize_x == h5_test.nsize_x);
    REQUIRE(h5_particle.nsize_y == h5_test.nsize_y);
    REQUIRE(h5_particle.nsize_z == h5_test.nsize_z);

    REQUIRE(h5_particle.velocity_x ==
            Approx(h5_test.velocity_x).epsilon(Tolerance));
    REQUIRE(h5_particle.velocity_y ==
            Approx(h5_test.velocity_y).epsilon(Tolerance));
    REQUIRE(h5_particle.velocity_z ==
            Approx(h5_test.velocity_z).epsilon(Tolerance));

    REQUIRE(h5_particle.stress_xx ==
            Approx(h5_test.stress_xx).epsilon(Tolerance));
    REQUIRE(h5_particle.stress_yy ==
            Approx(h5_test.stress_yy).epsilon(Tolerance));
    REQUIRE(h5_particle.stress_zz ==
            Approx(h5_test.stress_zz).epsilon(Tolerance));
    REQUIRE(h5_particle.tau_xy == Approx(h5_test.tau_xy).epsilon(Tolerance));
    REQUIRE(h5_particle.tau_yz == Approx(h5_test.tau_yz).epsilon(Tolerance));
    REQUIRE(h5_particle.tau_xz == Approx(h5_test.tau_xz).epsilon(Tolerance));

    REQUIRE(h5_particle.strain_xx ==
            Approx(h5_test.strain_xx).epsilon(Tolerance));
    REQUIRE(h5_particle.strain_yy ==
            Approx(h5_test.strain_yy).epsilon(Tolerance));
    REQUIRE(h5_particle.strain_zz ==
            Approx(h5_test.strain_zz).epsilon(Tolerance));
    REQUIRE(h5_particle.gamma_xy ==
            Approx(h5_test.gamma_xy).epsilon(Tolerance));
    REQUIRE(h5_particle.gamma_yz ==
            Approx(h5_test.gamma_yz).epsilon(Tolerance));
    REQUIRE(h5_particle.gamma_xz ==
            Approx(h5_test.gamma_xz).epsilon(Tolerance));

    REQUIRE(h5_particle.epsilon_v ==
            Approx(h5_test.epsilon_v).epsilon(Tolerance));
    REQUIRE(h5_particle.status == h5_test.status);
    REQUIRE(h5_particle.cell_id == h5_test.cell_id);
    REQUIRE(h5_particle.material_id == h5_test.material_id);
  }
}

//! \brief Check particle class for 3D case
TEST_CASE("Particle is checked for 3D case", "[particle][3D]") {
  // Dimension
  const unsigned Dim = 3;
  // Dimension
  const unsigned Dof = 6;
  // Number of nodes per cell
  const unsigned Nnodes = 8;
  // Number of phases
  const unsigned Nphases = 1;
  // Phase
  const unsigned phase = 0;
  // Tolerance
  const double Tolerance = 1.E-7;

  // Coordinates
  Eigen::Vector3d coords;
  coords.setZero();

  //! Check for id = 0
  SECTION("Particle id is zero") {
    mpm::Index id = 0;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim>>(id, coords);
    REQUIRE(particle->id() == 0);
    REQUIRE(particle->status() == true);
  }

  SECTION("Particle id is positive") {
    //! Check for id is a positive value
    mpm::Index id = std::numeric_limits<mpm::Index>::max();
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim>>(id, coords);
    REQUIRE(particle->id() == std::numeric_limits<mpm::Index>::max());
    REQUIRE(particle->status() == true);
  }

  //! Construct with id, coordinates and status
  SECTION("Particle with id, coordinates, and status") {
    mpm::Index id = 0;
    bool status = true;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim>>(id, coords, status);
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
        std::make_shared<mpm::Particle<Dim>>(id, coords);

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
        std::make_shared<mpm::Particle<Dim>>(id, coords);

    // Check particle coordinates
    auto coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    // Assign hexahedron shape function
    std::shared_ptr<mpm::Element<Dim>> element =
        std::make_shared<mpm::HexahedronElement<Dim, 8>>();

    // Create cell
    auto cell = std::make_shared<mpm::Cell<Dim>>(10, Nnodes, element);
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

    coords << 0, 0, 4;
    std::shared_ptr<mpm::NodeBase<Dim>> node8 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(4, coords);

    coords << 2, 0, 4;
    std::shared_ptr<mpm::NodeBase<Dim>> node9 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(5, coords);

    coords << 2, 2, 4;
    std::shared_ptr<mpm::NodeBase<Dim>> node10 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(6, coords);

    coords << 0, 2, 4;
    std::shared_ptr<mpm::NodeBase<Dim>> node11 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(7, coords);

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
    // Check particle cell status
    REQUIRE(particle->cell_ptr() == false);
    // Assign cell id
    REQUIRE(particle->assign_cell_id(10) == true);
    // Require cell id
    REQUIRE(particle->cell_id() == 10);
    // Assign a very large cell id
    REQUIRE(particle->assign_cell_id(std::numeric_limits<mpm::Index>::max()) ==
            false);
    // Require cell id
    REQUIRE(particle->cell_id() == 10);
    // Assign particle to cell
    REQUIRE(particle->assign_cell(cell) == true);
    // Local coordinates
    Eigen::Vector3d xi;
    xi.fill(std::numeric_limits<double>::max());
    // Assign particle to cell
    REQUIRE(particle->assign_cell_xi(cell, xi) == false);
    // Compute reference location
    cell->is_point_in_cell(particle->coordinates(), &xi);
    // Assign particle to cell
    REQUIRE(particle->assign_cell_xi(cell, xi) == true);

    // Assign cell id again
    REQUIRE(particle->assign_cell_id(10) == false);
    // Check particle cell status
    REQUIRE(particle->cell_ptr() == true);
    // Check cell status on addition of particle
    REQUIRE(cell->status() == true);

    // Create cell
    auto cell2 = std::make_shared<mpm::Cell<Dim>>(20, Nnodes, element);

    cell2->add_node(0, node4);
    cell2->add_node(1, node5);
    cell2->add_node(2, node6);
    cell2->add_node(3, node7);
    cell2->add_node(4, node8);
    cell2->add_node(5, node9);
    cell2->add_node(6, node10);
    cell2->add_node(7, node11);
    REQUIRE(cell2->nnodes() == 8);

    // Initialise cell2 properties
    cell2->initialise();

    // Check if cell2 is initialised
    REQUIRE(cell2->is_initialised() == true);

    // Add cell2 to particle
    REQUIRE(cell2->status() == false);
    // Assign particle to cell2
    REQUIRE(particle->assign_cell(cell2) == false);
    // Check cell2 status for failed addition of particle
    REQUIRE(cell2->status() == false);
    // Check cell status because this should not have removed the particle
    REQUIRE(cell->status() == true);

    // Remove assigned cell
    particle->remove_cell();
    REQUIRE(particle->assign_cell(cell) == true);

    // Clear all particle ids
    REQUIRE(cell->nparticles() == 1);
    cell->clear_particle_ids();
    REQUIRE(cell->nparticles() == 0);
  }

  //! Test initialise particle stresses
  SECTION("Particle with initial stress") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    bool status = true;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim>>(id, coords, status);
    Eigen::Matrix<double, 6, 1> stress =
        Eigen::Matrix<double, 6, 1>::Constant(5.7);

    particle->initial_stress(stress);
    REQUIRE(particle->stress().size() == stress.size());
    auto pstress = particle->stress();
    for (unsigned i = 0; i < pstress.size(); ++i)
      REQUIRE(pstress[i] == Approx(stress[i]).epsilon(Tolerance));
  }

  //! Test particles velocity constraints
  SECTION("Particle with velocity constraints") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    bool status = true;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim>>(id, coords, status);

    // Apply particles velocity constraints
    REQUIRE(particle->assign_particle_velocity_constraint(0, 10.5) == true);
    REQUIRE(particle->assign_particle_velocity_constraint(1, -12.5) == true);
    REQUIRE(particle->assign_particle_velocity_constraint(2, 14.5) == true);
    // Check out of bounds condition
    REQUIRE(particle->assign_particle_velocity_constraint(3, 0) == false);

    // Apply constraints
    particle->apply_particle_velocity_constraints();

    // Check apply constraints
    REQUIRE(particle->velocity()(0) == Approx(10.5).epsilon(Tolerance));
    REQUIRE(particle->velocity()(1) == Approx(-12.5).epsilon(Tolerance));
    REQUIRE(particle->velocity()(2) == Approx(14.5).epsilon(Tolerance));
  }

  //! Test particle, cell and node functions
  SECTION("Test particle, cell and node functions") {
    // Add particle
    mpm::Index id = 0;
    coords << 1.5, 1.5, 1.5;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim>>(id, coords);

    // Time-step
    const double dt = 0.1;

    // Check particle coordinates
    auto coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    // Assign hexahedron shape function
    std::shared_ptr<mpm::Element<Dim>> element =
        std::make_shared<mpm::HexahedronElement<Dim, 8>>();

    // Create cell
    auto cell = std::make_shared<mpm::Cell<Dim>>(10, Nnodes, element);
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
    // Compute updated particle location should fail
    REQUIRE(particle->compute_updated_position(dt) == false);
    // Compute updated particle location from nodal velocity should fail
    REQUIRE(particle->compute_updated_position_velocity(dt) == false);
    // Compute volume
    REQUIRE(particle->compute_volume() == false);
    // Update volume should fail
    REQUIRE(particle->update_volume_strainrate(dt) == false);

    REQUIRE(particle->assign_cell(cell) == true);
    REQUIRE(cell->status() == true);
    REQUIRE(particle->cell_id() == 10);

    // Check if cell is initialised
    REQUIRE(cell->is_initialised() == true);

    // Check compute shape functions of a particle
    REQUIRE(particle->compute_shapefn() == true);

    // Assign volume
    REQUIRE(particle->assign_volume(0.0) == false);
    REQUIRE(particle->assign_volume(-5.0) == false);
    REQUIRE(particle->assign_volume(2.0) == true);
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
    // Initialise material
    Json jmaterial;
    jmaterial["density"] = 1000.;
    jmaterial["youngs_modulus"] = 1.0E+7;
    jmaterial["poisson_ratio"] = 0.3;

    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "LinearElastic3D", std::move(mid), jmaterial);

    // Check compute mass before material and volume
    REQUIRE(particle->compute_mass() == false);

    // Test compute stress before material assignment
    REQUIRE(particle->compute_stress() == false);

    // Test compute internal force before material assignment
    REQUIRE(particle->map_internal_force() == false);

    // Assign material properties
    REQUIRE(particle->assign_material(material) == true);

    // Check material id from particle
    REQUIRE(particle->material_id() == 0);

    // Compute volume
    REQUIRE(particle->compute_volume() == true);

    // Compute mass
    REQUIRE(particle->compute_mass() == true);
    // Mass
    REQUIRE(particle->mass() == Approx(8000.).epsilon(Tolerance));

    // Map particle mass to nodes
    particle->assign_mass(std::numeric_limits<double>::max());
    REQUIRE(particle->map_mass_momentum_to_nodes() == false);

    // Map particle pressure to nodes
    REQUIRE(particle->map_pressure_to_nodes() == false);

    // Assign mass to nodes
    REQUIRE(particle->compute_reference_location() == true);
    REQUIRE(particle->compute_shapefn() == true);

    // Map particle material id to nodes
    particle->append_material_id_to_nodes();
    REQUIRE(node0->material_ids()[0] == 0);
    REQUIRE(node1->material_ids()[0] == 0);
    REQUIRE(node2->material_ids()[0] == 0);
    REQUIRE(node3->material_ids()[0] == 0);
    REQUIRE(node4->material_ids()[0] == 0);
    REQUIRE(node5->material_ids()[0] == 0);
    REQUIRE(node6->material_ids()[0] == 0);
    REQUIRE(node7->material_ids()[0] == 0);

    // Check velocity
    Eigen::VectorXd velocity;
    velocity.resize(Dim);
    for (unsigned i = 0; i < velocity.size(); ++i) velocity(i) = i;
    REQUIRE(particle->assign_velocity(velocity) == true);
    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(particle->velocity()(i) == Approx(i).epsilon(Tolerance));

    REQUIRE(particle->compute_mass() == true);
    REQUIRE(particle->map_mass_momentum_to_nodes() == true);

    REQUIRE(particle->map_pressure_to_nodes() == true);
    REQUIRE(particle->compute_pressure_smoothing() == true);

    // Values of nodal mass
    std::array<double, 8> nodal_mass{125., 375.,  1125., 375.,
                                     375., 1125., 3375., 1125.};
    // Check nodal mass
    for (unsigned i = 0; i < nodes.size(); ++i)
      REQUIRE(nodes.at(i)->mass(phase) ==
              Approx(nodal_mass.at(i)).epsilon(Tolerance));

    // Compute nodal velocity
    for (const auto node : nodes) node->compute_velocity();

    // Values of nodal momentum
    Eigen::Matrix<double, 8, 3> nodal_momentum;
    // clang-format off
    nodal_momentum << 0.,  125.,  250.,
                      0.,  375.,  750.,
                      0., 1125., 2250.,
                      0.,  375.,  750.,
                      0.,  375.,  750.,
                      0., 1125., 2250.,
                      0., 3375., 6750.,
                      0., 1125., 2250.;
    // clang-format on

    // Check nodal momentum
    for (unsigned i = 0; i < nodal_momentum.rows(); ++i)
      for (unsigned j = 0; j < nodal_momentum.cols(); ++j)
        REQUIRE(nodes.at(i)->momentum(phase)[j] ==
                Approx(nodal_momentum(i, j)).epsilon(Tolerance));

    // Values of nodal velocity
    Eigen::Matrix<double, 8, 3> nodal_velocity;
    // clang-format off
    nodal_velocity << 0., 1., 2.,
                      0., 1., 2.,
                      0., 1., 2.,
                      0., 1., 2.,
                      0., 1., 2.,
                      0., 1., 2.,
                      0., 1., 2.,
                      0., 1., 2.;
    // clang-format on
    // Check nodal velocity
    for (unsigned i = 0; i < nodal_velocity.rows(); ++i)
      for (unsigned j = 0; j < nodal_velocity.cols(); ++j)
        REQUIRE(nodes.at(i)->velocity(phase)(j) ==
                Approx(nodal_velocity(i, j)).epsilon(Tolerance));

    // Set momentum to get non-zero strain
    // clang-format off
    nodal_momentum << 0.,  125. * 1.,  250. * 1.,
                      0.,  375. * 2.,  750. * 2.,
                      0., 1125. * 3., 2250. * 3.,
                      0.,  375. * 4.,  750. * 4.,
                      0.,  375. * 5.,  750. * 5.,
                      0., 1125. * 6., 2250. * 6.,
                      0., 3375. * 7., 6750. * 7.,
                      0., 1125. * 8., 2250. * 8.;
    // clang-format on
    for (unsigned i = 0; i < nodes.size(); ++i)
      REQUIRE(nodes.at(i)->update_momentum(false, phase,
                                           nodal_momentum.row(i)) == true);

    // nodal velocity
    // clang-format off
    nodal_velocity << 0., 1.,  2.,
                      0., 2.,  4.,
                      0., 3.,  6.,
                      0., 4.,  8.,
                      0., 5., 10.,
                      0., 6., 12.,
                      0., 7., 14.,
                      0., 8., 16.;
    // clang-format on
    // Compute nodal velocity
    for (const auto node : nodes) node->compute_velocity();
    // Check nodal velocity
    for (unsigned i = 0; i < nodal_velocity.rows(); ++i)
      for (unsigned j = 0; j < nodal_velocity.cols(); ++j)
        REQUIRE(nodes.at(i)->velocity(phase)(j) ==
                Approx(nodal_velocity(i, j)).epsilon(Tolerance));

    // Check pressure
    REQUIRE(particle->pressure() == Approx(0.).epsilon(Tolerance));

    // Compute strain
    particle->compute_strain(dt);
    // Strain
    Eigen::Matrix<double, 6, 1> strain;
    strain << 0.00000, 0.07500, 0.40000, -0.02500, 0.35000, -0.05000;

    // Check strains
    for (unsigned i = 0; i < strain.rows(); ++i)
      REQUIRE(particle->strain()(i) == Approx(strain(i)).epsilon(Tolerance));

    // Check volumetric strain at centroid
    const double volumetric_strain = 0.5;
    REQUIRE(particle->volumetric_strain_centroid() ==
            Approx(volumetric_strain).epsilon(Tolerance));

    // Check updated pressure
    const double K = 8333333.333333333;
    REQUIRE(particle->pressure() ==
            Approx(-K * volumetric_strain).epsilon(Tolerance));

    // Update volume strain rate
    REQUIRE(particle->volume() == Approx(8.0).epsilon(Tolerance));
    REQUIRE(particle->update_volume_strainrate(dt) == true);
    REQUIRE(particle->volume() == Approx(12.0).epsilon(Tolerance));

    // Compute stress
    REQUIRE(particle->compute_stress() == true);

    Eigen::Matrix<double, 6, 1> stress;
    // clang-format off
    stress << 2740384.6153846150,
              3317307.6923076920,
              5817307.6923076920,
               -96153.8461538463,
              1346153.8461538465,
              -192307.6923076927;
    // clang-format on
    // Check stress
    for (unsigned i = 0; i < stress.rows(); ++i)
      REQUIRE(particle->stress()(i) == Approx(stress(i)).epsilon(Tolerance));

    // Check body force
    Eigen::Matrix<double, 3, 1> gravity;
    gravity << 0., 0., -9.81;

    particle->map_body_force(gravity);

    // Body force
    Eigen::Matrix<double, 8, 3> body_force;
    // clang-format off
    body_force << 0., 0.,  -1226.25,
                  0., 0.,  -3678.75,
                  0., 0., -11036.25,
                  0., 0.,  -3678.75,
                  0., 0.,  -3678.75,
                  0., 0., -11036.25,
                  0., 0., -33108.75,
                  0., 0., -11036.25;
    // clang-format on

    // Check nodal body force
    for (unsigned i = 0; i < body_force.rows(); ++i)
      for (unsigned j = 0; j < body_force.cols(); ++j)
        REQUIRE(nodes[i]->external_force(phase)[j] ==
                Approx(body_force(i, j)).epsilon(Tolerance));

    // Check traction force
    double traction = 7.68;
    const unsigned direction = 2;
    // Assign volume
    REQUIRE(particle->assign_volume(0.0) == false);
    REQUIRE(particle->assign_volume(-5.0) == false);
    REQUIRE(particle->assign_volume(2.0) == true);
    // Assign traction to particle
    particle->assign_traction(direction, traction);
    // Map traction force
    particle->map_traction_force();

    // Traction force
    Eigen::Matrix<double, 8, 3> traction_force;
    // shapefn * volume / size_(dir) * traction
    // clang-format off
    traction_force << 0., 0., 0.015625 * 1.587401052 * 7.68,
                      0., 0., 0.046875 * 1.587401052 * 7.68,
                      0., 0., 0.140625 * 1.587401052 * 7.68,
                      0., 0., 0.046875 * 1.587401052 * 7.68,
                      0., 0., 0.046875 * 1.587401052 * 7.68,
                      0., 0., 0.140625 * 1.587401052 * 7.68,
                      0., 0., 0.421875 * 1.587401052 * 7.68,
                      0., 0., 0.140625 * 1.587401052 * 7.68;
    // clang-format on
    // Add previous external body force
    traction_force += body_force;

    // Check nodal traction force
    for (unsigned i = 0; i < traction_force.rows(); ++i)
      for (unsigned j = 0; j < traction_force.cols(); ++j)
        REQUIRE(nodes[i]->external_force(phase)[j] ==
                Approx(traction_force(i, j)).epsilon(Tolerance));
    // Reset traction
    particle->assign_traction(direction, -traction);
    // Map traction force
    particle->map_traction_force();
    // Check nodal external force
    for (unsigned i = 0; i < traction_force.rows(); ++i)
      for (unsigned j = 0; j < traction_force.cols(); ++j)
        REQUIRE(nodes[i]->external_force(phase)[j] ==
                Approx(body_force(i, j)).epsilon(Tolerance));

    // Internal force
    Eigen::Matrix<double, 8, 3> internal_force;
    // clang-format off
    internal_force <<  612980.7692307689,  1141826.923076923,  1742788.461538461,
                      -901442.3076923079,  3521634.615384615,  5420673.076923076,
                      -2415865.384615385,  612980.7692307703,  12223557.69230769,
                       1935096.153846153,  108173.0769230771,  3882211.538461538,
                                 2031250,  2079326.923076922, -588942.3076923075,
                      -2127403.846153846,  6526442.307692306, -1189903.846153845,
                       -5516826.92307692, -10276442.30769231, -15685096.15384615,
                       6382211.538461537, -3713942.307692308, -5805288.461538462;
    // clang-format on

    // Map particle internal force
    particle->assign_volume(8.0);
    REQUIRE(particle->map_internal_force() == true);

    // Check nodal internal force
    for (unsigned i = 0; i < internal_force.rows(); ++i)
      for (unsigned j = 0; j < internal_force.cols(); ++j)
        REQUIRE(nodes[i]->internal_force(phase)[j] ==
                Approx(internal_force(i, j)).epsilon(Tolerance));

    // Calculate nodal acceleration and velocity
    for (const auto& node : nodes)
      node->compute_acceleration_velocity(phase, dt);

    // Check nodal velocity
    // clang-format off
    nodal_velocity <<  490.3846153846152,  914.4615384615383, 1395.249769230769,
                      -240.3846153846155,  941.1025641025641, 1448.531820512821,
                      -214.7435897435898,   57.4871794871796, 1091.557461538462,
                       516.0256410256410,   32.8461538461539, 1042.275410256410,
                       541.6666666666666,  559.4871794871794, -148.032282051282,
                      -189.1025641025641,  586.1282051282051, -94.75023076923067,
                      -163.4615384615384, -297.4871794871795, -451.7245897435898,
                       567.3076923076923, -322.1282051282053, -501.0066410256412;
    // clang-format on
    // Check nodal velocity
    for (unsigned i = 0; i < nodal_velocity.rows(); ++i)
      for (unsigned j = 0; j < nodal_velocity.cols(); ++j)
        REQUIRE(nodes[i]->velocity(phase)[j] ==
                Approx(nodal_velocity(i, j)).epsilon(Tolerance));

    // Check nodal acceleration
    Eigen::Matrix<double, 8, 3> nodal_acceleration;
    // clang-format off
    nodal_acceleration << 4903.846153846152, 9134.615384615383, 13932.49769230769,
                         -2403.846153846155, 9391.025641025641, 14445.31820512821,
                         -2147.435897435898, 544.8717948717959, 10855.57461538462,
                          5160.256410256409, 288.461538461539,  10342.7541025641,
                          5416.666666666666, 5544.871794871794, -1580.32282051282,
                         -1891.025641025641, 5801.282051282051, -1067.502307692307,
                         -1634.615384615384, -3044.871794871795, -4657.245897435898,
                          5673.076923076923, -3301.282051282052, -5170.066410256411;
    // clang-format on
    // Check nodal acceleration
    for (unsigned i = 0; i < nodal_acceleration.rows(); ++i)
      for (unsigned j = 0; j < nodal_acceleration.cols(); ++j)
        REQUIRE(nodes[i]->acceleration(phase)[j] ==
                Approx(nodal_acceleration(i, j)).epsilon(Tolerance));
    // Approx(nodal_velocity(i, j) / dt).epsilon(Tolerance));

    // Check original particle coordinates
    coords << 1.5, 1.5, 1.5;
    coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    // Compute updated particle location
    REQUIRE(particle->compute_updated_position(dt) == true);
    // Check particle velocity
    velocity << 0., 1., 1.019;
    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(particle->velocity()(i) ==
              Approx(velocity(i)).epsilon(Tolerance));

    // Check particle displacement
    Eigen::Vector3d displacement;
    displacement << 0.0, 0.5875, 1.0769;
    for (unsigned i = 0; i < displacement.size(); ++i)
      REQUIRE(particle->displacement()(i) ==
              Approx(displacement(i)).epsilon(Tolerance));

    // Updated particle coordinate
    coords << 1.5, 2.0875, 2.5769;
    // Check particle coordinates
    coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    // Compute updated particle location based on nodal velocity
    REQUIRE(particle->compute_updated_position_velocity(dt) == true);
    // Check particle velocity
    velocity << 0., 5.875, 10.769;
    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(particle->velocity()(i) ==
              Approx(velocity(i)).epsilon(Tolerance));

    // Check particle displacement
    displacement << 0.0, 1.175, 2.1538;
    for (unsigned i = 0; i < displacement.size(); ++i)
      REQUIRE(particle->displacement()(i) ==
              Approx(displacement(i)).epsilon(Tolerance));

    // Updated particle coordinate
    coords << 1.5, 2.675, 3.6538;
    // Check particle coordinates
    coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));
  }

  SECTION("Check assign material to particle") {
    // Add particle
    mpm::Index id = 0;
    coords << 0.75, 0.75, 0.75;
    auto particle = std::make_shared<mpm::Particle<Dim>>(id, coords);

    unsigned mid = 1;
    // Initialise material
    Json jmaterial;
    jmaterial["density"] = 1000.;
    jmaterial["youngs_modulus"] = 1.0E+7;
    jmaterial["poisson_ratio"] = 0.3;

    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "LinearElastic3D", std::move(mid), jmaterial);
    REQUIRE(material->id() == 1);

    // Check if particle can be assigned a null material
    REQUIRE(particle->assign_material(nullptr) == false);
    // Check material id
    REQUIRE(particle->material_id() == std::numeric_limits<unsigned>::max());

    // Assign material to particle
    REQUIRE(particle->assign_material(material) == true);
    // Check material id
    REQUIRE(particle->material_id() == 1);
  }

  SECTION("Check particle properties") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim>>(id, coords);

    // Check mass
    REQUIRE(particle->mass() == Approx(0.0).epsilon(Tolerance));
    double mass = 100.5;
    particle->assign_mass(mass);
    REQUIRE(particle->mass() == Approx(100.5).epsilon(Tolerance));

    // Check stress
    Eigen::Matrix<double, 6, 1> stress;
    for (unsigned i = 0; i < stress.size(); ++i) stress(i) = 1.;

    for (unsigned i = 0; i < stress.size(); ++i)
      REQUIRE(particle->stress()(i) == Approx(0.).epsilon(Tolerance));

    // Check velocity
    Eigen::VectorXd velocity;
    velocity.resize(Dim);
    for (unsigned i = 0; i < velocity.size(); ++i) velocity(i) = 17.51;

    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(particle->velocity()(i) == Approx(0.).epsilon(Tolerance));

    REQUIRE(particle->assign_velocity(velocity) == true);
    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(particle->velocity()(i) == Approx(17.51).epsilon(Tolerance));

    // Assign volume
    REQUIRE(particle->assign_volume(0.0) == false);
    REQUIRE(particle->assign_volume(-5.0) == false);
    REQUIRE(particle->assign_volume(2.0) == true);
    // Check volume
    REQUIRE(particle->volume() == Approx(2.0).epsilon(Tolerance));
    // Traction
    double traction = 65.32;
    const unsigned Direction = 1;
    // Check traction
    for (unsigned i = 0; i < Dim; ++i)
      REQUIRE(particle->traction()(i) == Approx(0.).epsilon(Tolerance));

    REQUIRE(particle->assign_traction(Direction, traction) == true);

    // Calculate traction force = traction * volume / spacing
    traction *= 2.0 / (std::pow(2.0, 1. / Dim));

    for (unsigned i = 0; i < Dim; ++i) {
      if (i == Direction)
        REQUIRE(particle->traction()(i) == Approx(traction).epsilon(Tolerance));
      else
        REQUIRE(particle->traction()(i) == Approx(0.).epsilon(Tolerance));
    }

    // Check for incorrect direction
    const unsigned wrong_dir = 4;
    REQUIRE(particle->assign_traction(wrong_dir, traction) == false);

    // Check again to ensure value hasn't been updated
    for (unsigned i = 0; i < Dim; ++i) {
      if (i == Direction)
        REQUIRE(particle->traction()(i) == Approx(traction).epsilon(Tolerance));
      else
        REQUIRE(particle->traction()(i) == Approx(0.).epsilon(Tolerance));
    }
  }

  // Check initialise particle from HDF5 file
  SECTION("Check initialise particle HDF5") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim>>(id, coords);

    mpm::HDF5Particle h5_particle;
    h5_particle.id = 13;
    h5_particle.mass = 501.5;

    Eigen::Vector3d coords;
    coords << 1., 2., 3.;
    h5_particle.coord_x = coords[0];
    h5_particle.coord_y = coords[1];
    h5_particle.coord_z = coords[2];

    Eigen::Vector3d displacement;
    displacement << 0.01, 0.02, 0.03;
    h5_particle.displacement_x = displacement[0];
    h5_particle.displacement_y = displacement[1];
    h5_particle.displacement_z = displacement[2];

    Eigen::Vector3d lsize;
    lsize << 0.25, 0.5, 0.75;
    h5_particle.nsize_x = lsize[0];
    h5_particle.nsize_y = lsize[1];
    h5_particle.nsize_z = lsize[2];

    Eigen::Vector3d velocity;
    velocity << 1.5, 2.5, 3.5;
    h5_particle.velocity_x = velocity[0];
    h5_particle.velocity_y = velocity[1];
    h5_particle.velocity_z = velocity[2];

    Eigen::Matrix<double, 6, 1> stress;
    stress << 11.5, -12.5, 13.5, 14.5, -15.5, 16.5;
    h5_particle.stress_xx = stress[0];
    h5_particle.stress_yy = stress[1];
    h5_particle.stress_zz = stress[2];
    h5_particle.tau_xy = stress[3];
    h5_particle.tau_yz = stress[4];
    h5_particle.tau_xz = stress[5];

    Eigen::Matrix<double, 6, 1> strain;
    strain << 0.115, -0.125, 0.135, 0.145, -0.155, 0.165;
    h5_particle.strain_xx = strain[0];
    h5_particle.strain_yy = strain[1];
    h5_particle.strain_zz = strain[2];
    h5_particle.gamma_xy = strain[3];
    h5_particle.gamma_yz = strain[4];
    h5_particle.gamma_xz = strain[5];

    h5_particle.epsilon_v = strain.head(Dim).sum();

    h5_particle.status = true;

    h5_particle.cell_id = 1;

    h5_particle.volume = 2.;

    h5_particle.material_id = 1;

    // Reinitialise particle from HDF5 data
    REQUIRE(particle->initialise_particle(h5_particle) == true);

    // Check particle id
    REQUIRE(particle->id() == h5_particle.id);
    // Check particle mass
    REQUIRE(particle->mass() == h5_particle.mass);
    // Check particle volume
    REQUIRE(particle->volume() == h5_particle.volume);
    // Check particle mass density
    REQUIRE(particle->mass_density() == h5_particle.mass / h5_particle.volume);
    // Check particle status
    REQUIRE(particle->status() == h5_particle.status);

    // Check for coordinates
    auto coordinates = particle->coordinates();
    REQUIRE(coordinates.size() == Dim);
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));
    REQUIRE(coordinates.size() == Dim);

    // Check for displacement
    auto pdisplacement = particle->displacement();
    REQUIRE(pdisplacement.size() == Dim);
    for (unsigned i = 0; i < Dim; ++i)
      REQUIRE(pdisplacement(i) == Approx(displacement(i)).epsilon(Tolerance));

    // Check for size
    auto size = particle->natural_size();
    REQUIRE(size.size() == Dim);
    for (unsigned i = 0; i < size.size(); ++i)
      REQUIRE(size(i) == Approx(lsize(i)).epsilon(Tolerance));

    // Check velocity
    auto pvelocity = particle->velocity();
    REQUIRE(pvelocity.size() == Dim);
    for (unsigned i = 0; i < Dim; ++i)
      REQUIRE(pvelocity(i) == Approx(velocity(i)).epsilon(Tolerance));

    // Check stress
    auto pstress = particle->stress();
    REQUIRE(pstress.size() == stress.size());
    for (unsigned i = 0; i < stress.size(); ++i)
      REQUIRE(pstress(i) == Approx(stress(i)).epsilon(Tolerance));

    // Check strain
    auto pstrain = particle->strain();
    REQUIRE(pstrain.size() == strain.size());
    for (unsigned i = 0; i < strain.size(); ++i)
      REQUIRE(pstrain(i) == Approx(strain(i)).epsilon(Tolerance));

    // Check particle volumetric strain centroid
    REQUIRE(particle->volumetric_strain_centroid() == h5_particle.epsilon_v);

    // Check cell id
    REQUIRE(particle->cell_id() == h5_particle.cell_id);

    // Check material id
    REQUIRE(particle->material_id() == h5_particle.material_id);

    // Write Particle HDF5 data
    const auto h5_test = particle->hdf5();

    REQUIRE(h5_particle.id == h5_test.id);
    REQUIRE(h5_particle.mass == h5_test.mass);

    REQUIRE(h5_particle.coord_x == Approx(h5_test.coord_x).epsilon(Tolerance));
    REQUIRE(h5_particle.coord_y == Approx(h5_test.coord_y).epsilon(Tolerance));
    REQUIRE(h5_particle.coord_z == Approx(h5_test.coord_z).epsilon(Tolerance));

    REQUIRE(h5_particle.displacement_x ==
            Approx(h5_test.displacement_x).epsilon(Tolerance));
    REQUIRE(h5_particle.displacement_y ==
            Approx(h5_test.displacement_y).epsilon(Tolerance));
    REQUIRE(h5_particle.displacement_z ==
            Approx(h5_test.displacement_z).epsilon(Tolerance));

    REQUIRE(h5_particle.nsize_x == h5_test.nsize_x);
    REQUIRE(h5_particle.nsize_y == h5_test.nsize_y);
    REQUIRE(h5_particle.nsize_z == h5_test.nsize_z);

    REQUIRE(h5_particle.velocity_x ==
            Approx(h5_test.velocity_x).epsilon(Tolerance));
    REQUIRE(h5_particle.velocity_y ==
            Approx(h5_test.velocity_y).epsilon(Tolerance));
    REQUIRE(h5_particle.velocity_z ==
            Approx(h5_test.velocity_z).epsilon(Tolerance));

    REQUIRE(h5_particle.stress_xx ==
            Approx(h5_test.stress_xx).epsilon(Tolerance));
    REQUIRE(h5_particle.stress_yy ==
            Approx(h5_test.stress_yy).epsilon(Tolerance));
    REQUIRE(h5_particle.stress_zz ==
            Approx(h5_test.stress_zz).epsilon(Tolerance));
    REQUIRE(h5_particle.tau_xy == Approx(h5_test.tau_xy).epsilon(Tolerance));
    REQUIRE(h5_particle.tau_yz == Approx(h5_test.tau_yz).epsilon(Tolerance));
    REQUIRE(h5_particle.tau_xz == Approx(h5_test.tau_xz).epsilon(Tolerance));

    REQUIRE(h5_particle.strain_xx ==
            Approx(h5_test.strain_xx).epsilon(Tolerance));
    REQUIRE(h5_particle.strain_yy ==
            Approx(h5_test.strain_yy).epsilon(Tolerance));
    REQUIRE(h5_particle.strain_zz ==
            Approx(h5_test.strain_zz).epsilon(Tolerance));
    REQUIRE(h5_particle.gamma_xy ==
            Approx(h5_test.gamma_xy).epsilon(Tolerance));
    REQUIRE(h5_particle.gamma_yz ==
            Approx(h5_test.gamma_yz).epsilon(Tolerance));
    REQUIRE(h5_particle.gamma_xz ==
            Approx(h5_test.gamma_xz).epsilon(Tolerance));

    REQUIRE(h5_particle.epsilon_v ==
            Approx(h5_test.epsilon_v).epsilon(Tolerance));
    REQUIRE(h5_particle.status == h5_test.status);
    REQUIRE(h5_particle.cell_id == h5_test.cell_id);
    REQUIRE(h5_particle.material_id == h5_test.material_id);
  }
}

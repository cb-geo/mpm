#include <limits>

#include "catch.hpp"

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
#include "twophase_particle.h"

//! \brief Check twophase particle class for 1D case
TEST_CASE("TwoPhase Particle is checked for 1D case",
          "[twophase_particle][1D]") {
  // Dimension
  const unsigned Dim = 1;
  // Dimension
  const unsigned Dof = 1;
  // Number of phases
  const unsigned Nphases = 2;
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

  // Coordinates
  Eigen::Matrix<double, 1, 1> coords;
  coords.setZero();

  //! Check for id = 0
  SECTION("TwoPhase Particle id is zero") {
    mpm::Index id = 0;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, coords);
    REQUIRE(particle->id() == 0);
    REQUIRE(particle->status() == true);
  }

  SECTION("TwoPhase Particle id is positive") {
    //! Check for id is a positive value
    mpm::Index id = std::numeric_limits<mpm::Index>::max();
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, coords);
    REQUIRE(particle->id() == std::numeric_limits<mpm::Index>::max());
    REQUIRE(particle->status() == true);
  }

  //! Construct with id, coordinates and status
  SECTION("TwoPhase Particle with id, coordinates, and status") {
    mpm::Index id = 0;
    bool status = true;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, coords, status);
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
        std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, coords);

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
  SECTION("TwoPhase Particle with initial stress") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    bool status = true;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, coords, status);
    //! Test initialise particle stresses
    Eigen::Matrix<double, 6, 1> stress =
        Eigen::Matrix<double, 6, 1>::Constant(5.7);
    particle->initial_stress(stress);
    REQUIRE(particle->stress().size() == stress.size());
    auto pstress = particle->stress();
    for (unsigned i = 0; i < pstress.size(); ++i)
      REQUIRE(pstress[i] == Approx(stress[i]).epsilon(Tolerance));

    auto pstress_data = particle->tensor_data("stresses");
    for (unsigned i = 0; i < pstress_data.size(); ++i)
      REQUIRE(pstress_data[i] == Approx(stress[i]).epsilon(Tolerance));
  }

  //! Test particles velocity constraints
  SECTION("TwoPhase Particle with velocity constraints") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    bool status = true;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, coords, status);
    // Apply constraints for solid phase
    particle->apply_particle_velocity_constraints(0, 10.5);
    particle->apply_particle_velocity_constraints(1, 20.5);

    // Check apply constraints
    REQUIRE(particle->velocity()(0) == Approx(10.5).epsilon(Tolerance));
    REQUIRE(particle->liquid_velocity()(0) == Approx(20.5).epsilon(Tolerance));
  }

  SECTION("Check particle properties") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, coords);

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
    // Try with a null math fuction ptr
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
        std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, coords);

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

//! \brief Check twophase particle class for 2D case
TEST_CASE("TwoPhase Particle is checked for 2D case",
          "[twophase_particle][2D]") {
  // Dimension
  const unsigned Dim = 2;
  // Degree of freedom
  const unsigned Dof = 2;
  // Number of nodes per cell
  const unsigned Nnodes = 4;
  // Number of phases
  const unsigned Nphases = 2;
  // Tolerance
  const double Tolerance = 1.E-7;
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
  // Coordinates
  Eigen::Vector2d coords;
  coords.setZero();

  //! Check for id = 0
  SECTION("TwoPhase Particle id is zero") {
    mpm::Index id = 0;
    auto particle = std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, coords);
    REQUIRE(particle->id() == 0);
  }

  SECTION("TwoPhase Particle id is positive") {
    //! Check for id is a positive value
    mpm::Index id = std::numeric_limits<mpm::Index>::max();
    auto particle = std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, coords);
    REQUIRE(particle->id() == std::numeric_limits<mpm::Index>::max());
  }

  //! Test coordinates function
  SECTION("coordinates function is checked") {
    mpm::Index id = 0;
    auto particle = std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, coords);

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
    auto particle = std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, coords);

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
  SECTION("TwoPhase Particle with initial stress") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    bool status = true;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, coords, status);
    //! Test initialise particle stresses
    Eigen::Matrix<double, 6, 1> stress =
        Eigen::Matrix<double, 6, 1>::Constant(5.7);
    particle->initial_stress(stress);
    REQUIRE(particle->stress().size() == stress.size());
    auto pstress = particle->stress();
    for (unsigned i = 0; i < pstress.size(); ++i)
      REQUIRE(pstress[i] == Approx(stress[i]).epsilon(Tolerance));
  }

  //! Test particles velocity constraints
  SECTION("TwoPhase Particle with velocity constraints") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    bool status = true;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, coords, status);

    // Apply constraints
    particle->apply_particle_velocity_constraints(0, 10.5);
    particle->apply_particle_velocity_constraints(1, -12.5);
    particle->apply_particle_velocity_constraints(2, 20.5);
    particle->apply_particle_velocity_constraints(3, -22.5);

    // Check apply constraints
    REQUIRE(particle->velocity()(0) == Approx(10.5).epsilon(Tolerance));
    REQUIRE(particle->velocity()(1) == Approx(-12.5).epsilon(Tolerance));
    REQUIRE(particle->liquid_velocity()(0) == Approx(20.5).epsilon(Tolerance));
    REQUIRE(particle->liquid_velocity()(1) == Approx(-22.5).epsilon(Tolerance));
  }

  //! Test particle, cell and node functions
  SECTION("Test twophase particle, cell and node functions") {
    // Add particle
    mpm::Index id = 0;
    coords << 0.75, 0.75;
    auto particle = std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, coords);

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
    // TODO Assert: REQUIRE_NOTHROW(particle->compute_shapefn());
    // Compute reference location should throw
    REQUIRE(particle->compute_reference_location() == false);
    // Compute updated particle location should fail
    // TODO Assert:
    // REQUIRE_NOTHROW(particle->compute_updated_position(dt) == false);
    // Compute updated particle location from nodal velocity should fail
    // TODO Assert: REQUIRE_NOTHROW(particle->compute_updated_position(dt,
    // true)); Compute volume
    // TODO Assert: REQUIRE(particle->compute_volume() == false);
    // Update volume should fail
    // TODO Assert: REQUIRE(particle->update_volume() == false);

    REQUIRE(particle->assign_cell(cell) == true);
    REQUIRE(cell->status() == true);
    REQUIRE(particle->cell_id() == 10);

    // Check if cell is initialised
    REQUIRE(cell->is_initialised() == true);

    // Check compute shape functions of a particle
    REQUIRE_NOTHROW(particle->compute_shapefn());

    // Assign volume
    REQUIRE(particle->assign_volume(0.0) == false);
    REQUIRE(particle->assign_volume(-5.0) == false);
    REQUIRE(particle->assign_volume(2.0) == true);
    // Check volume
    REQUIRE(particle->volume() == Approx(2.0).epsilon(Tolerance));
    // Compute volume
    REQUIRE_NOTHROW(particle->compute_volume());
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
    jmaterial["porosity"] = 0.3;
    jmaterial["k_x"] = 0.001;
    jmaterial["k_y"] = 0.001;

    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "LinearElastic2D", std::move(mid), jmaterial);

    // Assign liquid material
    unsigned liquid_mid = 2;
    // Initialise material
    Json jmaterial_liquid;
    jmaterial_liquid["density"] = 1000.;
    jmaterial_liquid["bulk_modulus"] = 1.0E+9;
    jmaterial_liquid["mu"] = 0.3;
    jmaterial_liquid["dynamic_viscosity"] = 0.;

    auto liquid_material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "Newtonian2D", std::move(liquid_mid), jmaterial_liquid);

    // Check compute mass before material and volume
    // TODO Assert: REQUIRE(particle->compute_mass() == false);

    // Test compute stress before material assignment
    // TODO Assert: REQUIRE(particle->compute_stress() == false);

    // Assign material properties
    REQUIRE(particle->assign_material(material) == true);
    REQUIRE(particle->assign_material(liquid_material,
                                      mpm::ParticlePhase::Liquid) == true);

    // Assign porosity
    REQUIRE(particle->assign_porosity() == true);

    // Assign permeability
    REQUIRE(particle->assign_permeability() == true);

    // Check material id
    REQUIRE(particle->material_id() == 1);
    REQUIRE(particle->material_id(mpm::ParticlePhase::Liquid) == 2);

    // Compute volume
    REQUIRE_NOTHROW(particle->compute_volume());

    // Compute mass
    REQUIRE_NOTHROW(particle->compute_mass());
    // Mass
    REQUIRE(particle->mass() == Approx(700.).epsilon(Tolerance));
    REQUIRE(particle->liquid_mass() == Approx(300.).epsilon(Tolerance));

    // Map particle mass to nodes
    particle->assign_mass(std::numeric_limits<double>::max());
    // TODO Assert: REQUIRE_NOTHROW(particle->map_mass_momentum_to_nodes());

    // Map particle pressure to nodes
    // TODO Assert: REQUIRE(particle->map_pressure_to_nodes() == false);

    // Assign mass to nodes
    REQUIRE(particle->compute_reference_location() == true);
    REQUIRE_NOTHROW(particle->compute_shapefn());

    // Check velocity
    Eigen::VectorXd velocity;
    velocity.resize(Dim);
    for (unsigned i = 0; i < velocity.size(); ++i) velocity(i) = i;
    REQUIRE(particle->assign_velocity(velocity) == true);
    REQUIRE(particle->assign_liquid_velocity(velocity) == true);
    for (unsigned i = 0; i < velocity.size(); ++i) {
      REQUIRE(particle->velocity()(i) == Approx(i).epsilon(Tolerance));
      REQUIRE(particle->liquid_velocity()(i) == Approx(i).epsilon(Tolerance));
    }

    REQUIRE_NOTHROW(particle->compute_mass());
    REQUIRE_NOTHROW(particle->map_mass_momentum_to_nodes());

    // TODO Assert: REQUIRE(particle->map_pressure_to_nodes() == false);
    REQUIRE(particle->compute_pressure_smoothing() == false);

    // Values of nodal mass
    std::array<double, 4> nodal_mass{562.5, 187.5, 62.5, 187.5};
    // Check nodal mass
    for (unsigned i = 0; i < nodes.size(); ++i) {
      // Solid phase
      REQUIRE(nodes.at(i)->mass(mpm::NodePhase::nSolid) ==
              Approx(nodal_mass.at(i) * (1 - particle->porosity()))
                  .epsilon(Tolerance));
      // Liquid phase
      REQUIRE(
          nodes.at(i)->mass(mpm::NodePhase::nLiquid) ==
          Approx(nodal_mass.at(i) * particle->porosity()).epsilon(Tolerance));
    }

    // Compute nodal velocity
    for (const auto& node : nodes) node->compute_velocity();

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
      for (unsigned j = 0; j < nodal_momentum.cols(); ++j) {
        // Solid phase
        REQUIRE(nodes.at(i)->momentum(mpm::NodePhase::nSolid)(j) ==
                Approx(nodal_momentum(i, j) * (1 - particle->porosity()))
                    .epsilon(Tolerance));
        // Liquid phase
        REQUIRE(nodes.at(i)->momentum(mpm::NodePhase::nLiquid)(j) ==
                Approx(nodal_momentum(i, j) * particle->porosity())
                    .epsilon(Tolerance));
      }
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
      for (unsigned j = 0; j < nodal_velocity.cols(); ++j) {
        // Solid phase
        REQUIRE(nodes.at(i)->velocity(mpm::NodePhase::nSolid)(j) ==
                Approx(nodal_velocity(i, j)).epsilon(Tolerance));
        // Liquid phase
        REQUIRE(nodes.at(i)->velocity(mpm::NodePhase::nLiquid)(j) ==
                Approx(nodal_velocity(i, j)).epsilon(Tolerance));
      }

    // Set momentum to get non-zero strain
    // clang-format off
    nodal_momentum << 0., 562.5 * 1.,
                      0., 187.5 * 2.,
                      0.,  62.5 * 3.,
                      0., 187.5 * 4.;
    // clang-format on
    for (unsigned i = 0; i < nodes.size(); ++i) {
      // Solid phase
      REQUIRE_NOTHROW(nodes.at(i)->update_momentum(
          false, mpm::NodePhase::nSolid,
          nodal_momentum.row(i) * (1 - particle->porosity())));
      // Liquid phase
      REQUIRE_NOTHROW(nodes.at(i)->update_momentum(
          false, mpm::NodePhase::nLiquid,
          nodal_momentum.row(i) * particle->porosity()));
    }

    // nodal velocity
    // clang-format off
    nodal_velocity << 0., 1.,
                      0., 2.,
                      0., 3.,
                      0., 4.;
    // clang-format on
    // Compute nodal velocity
    for (const auto& node : nodes) node->compute_velocity();
    // Check nodal velocity
    for (unsigned i = 0; i < nodal_velocity.rows(); ++i)
      for (unsigned j = 0; j < nodal_velocity.cols(); ++j) {
        // Solid phase
        REQUIRE(nodes.at(i)->velocity(mpm::NodePhase::nSolid)(j) ==
                Approx(nodal_velocity(i, j)).epsilon(Tolerance));
        // Liquid phase
        REQUIRE(nodes.at(i)->velocity(mpm::NodePhase::nLiquid)(j) ==
                Approx(nodal_velocity(i, j)).epsilon(Tolerance));
      }

    // Check pressure
    REQUIRE(std::isnan(particle->pressure()) == true);

    // Compute strain
    particle->compute_strain(dt);
    // Strain
    Eigen::Matrix<double, 6, 1> strain;
    strain << 0., 0.25, 0., 0.050, 0., 0.;
    // Check strains
    for (unsigned i = 0; i < strain.rows(); ++i)
      REQUIRE(particle->strain()(i) == Approx(strain(i)).epsilon(Tolerance));

    // Check volumetric strain at centroid
    double volumetric_strain = 0.2;
    REQUIRE(particle->volumetric_strain_centroid() ==
            Approx(volumetric_strain).epsilon(Tolerance));

    // Check updated pressure
    const double K = 8333333.333333333;
    REQUIRE(std::isnan(particle->pressure()) == true);

    // Update volume strain rate
    REQUIRE(particle->volume() == Approx(1.0).epsilon(Tolerance));
    particle->compute_strain(dt);
    REQUIRE_NOTHROW(particle->update_volume());
    REQUIRE(particle->volume() == Approx(1.2).epsilon(Tolerance));

    // Compute stress
    REQUIRE_NOTHROW(particle->compute_stress());

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

    // Compute pore_pressure
    REQUIRE_NOTHROW(particle->compute_pore_pressure(dt));
    // Check pore pressure
    REQUIRE(particle->pressure(mpm::ParticlePhase::Liquid) ==
            Approx(-783333333.3333334923).epsilon(Tolerance));

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
        REQUIRE(nodes[i]->external_force(mpm::NodePhase::nSolid)[j] ==
                Approx(body_force(i, j)).epsilon(Tolerance));

    // Check traction force
    double traction = 7.68;
    const unsigned direction = 1;
    // Assign volume
    REQUIRE(particle->assign_volume(0.0) == false);
    REQUIRE(particle->assign_volume(-5.0) == false);
    REQUIRE(particle->assign_volume(2.0) == true);
    // Map traction force
    double current_time = 5.0;
    // Assign traction to particle
    particle->assign_traction(direction,
                              mfunction->value(current_time) * traction);
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
        REQUIRE(nodes[i]->external_force(mpm::NodePhase::nSolid)[j] ==
                Approx(traction_force(i, j)).epsilon(Tolerance));
    // Reset traction
    particle->assign_traction(direction,
                              -traction * mfunction->value(current_time));
    // Map traction force
    particle->map_traction_force();
    // Check nodal external force
    for (unsigned i = 0; i < traction_force.rows(); ++i)
      for (unsigned j = 0; j < traction_force.cols(); ++j)
        REQUIRE(nodes[i]->external_force(mpm::NodePhase::nSolid)[j] ==
                Approx(body_force(i, j)).epsilon(Tolerance));

    // Internal force
    Eigen::Matrix<double, 4, 2> internal_force;
    // clang-format off
    internal_force <<  588725961.538461566,  590168269.2307692766,
                      -588533653.8461539745,  196530448.7179487348,
                      -196241987.1794872284, -196722756.4102564454,
                       196049679.4871795177, -589975961.5384616852;
    // clang-format on

    // Map particle internal force
    particle->assign_volume(1.0);
    particle->map_internal_force();

    // Check nodal internal force
    for (unsigned i = 0; i < internal_force.rows(); ++i)
      for (unsigned j = 0; j < internal_force.cols(); ++j)
        REQUIRE(nodes[i]->internal_force(mpm::NodePhase::nMixture)[j] ==
                Approx(internal_force(i, j)).epsilon(Tolerance));

    // Internal force
    Eigen::Matrix<double, 4, 2> drag_force_coefficient;
    // clang-format off
    drag_force_coefficient << 496631.25,  496631.25,
                              165543.75,  165543.75,
                              55181.25,   55181.25,
                              165543.75,  165543.75;

    // Map drag force coefficient
    particle->map_drag_force_coefficient();

    // Check nodal drag force coefficient
    for (unsigned i = 0; i < drag_force_coefficient.rows(); ++i)
      for (unsigned j = 0; j < drag_force_coefficient.cols(); ++j)
        REQUIRE(nodes[i]->drag_force_coefficient()[j] ==
                Approx(drag_force_coefficient(i, j)).epsilon(Tolerance));

    // Calculate nodal acceleration and velocity
    for (const auto& node : nodes)
      node->compute_acceleration_velocity_twophase_explicit(dt);

    // Check nodal velocity
    Eigen::Matrix<double, 4, 2> nodal_liquid_velocity;
    // clang-format off
    nodal_velocity <<  104755.7997557998,  105122.1191221001,
                      -314120.8791208792,  104976.59897558,
                      -314267.3992673994, -315364.2813663005,
                       104609.2796092796, -315216.7612197804;
    nodal_liquid_velocity <<  104444.4444444445,  104444.4634444445,
                             -313333.3333333334,  104445.4634444445,
                             -313333.3333333334,  -313331.3143333333,
                              104444.4444444445,  -313330.3143333334;
    // clang-format on
    // Check nodal velocity
    for (unsigned i = 0; i < nodal_velocity.rows(); ++i)
      for (unsigned j = 0; j < nodal_velocity.cols(); ++j) {
        // Solid phase
        REQUIRE(nodes[i]->velocity(mpm::NodePhase::nSolid)[j] ==
                Approx(nodal_velocity(i, j)).epsilon(Tolerance));
        // Liquid phase
        REQUIRE(nodes[i]->velocity(mpm::NodePhase::nLiquid)[j] ==
                Approx(nodal_liquid_velocity(i, j)).epsilon(Tolerance));
      }

    // Check nodal acceleration
    Eigen::Matrix<double, 4, 2> nodal_acceleration;
    Eigen::Matrix<double, 4, 2> nodal_liquid_acceleration;
    // clang-format off
    nodal_acceleration <<  1047557.9975579976, 1051211.1912210013,
                           -3141208.791208792, 1049745.9897558,
                           -3142673.9926739936, -3153672.8136630044,
                            1046092.7960927964, -3152207.6121978033;
    nodal_liquid_acceleration <<  1044444.4444444446, 1044434.6344444447,
                                  -3133333.333333334, 1044434.6344444446,
                                  -3133333.333333334, -3133343.1433333335,
                                   1044444.4444444446, -3133343.1433333335;
    // clang-format on
    // Check nodal acceleration
    for (unsigned i = 0; i < nodal_acceleration.rows(); ++i)
      for (unsigned j = 0; j < nodal_acceleration.cols(); ++j) {
        // Solid phase
        REQUIRE(nodes[i]->acceleration(mpm::NodePhase::nSolid)[j] ==
                Approx(nodal_acceleration(i, j)).epsilon(Tolerance));
        // Liquid phase
        REQUIRE(nodes[i]->acceleration(mpm::NodePhase::nLiquid)[j] ==
                Approx(nodal_liquid_acceleration(i, j)).epsilon(Tolerance));
      }
    // Approx(nodal_velocity(i, j) / dt).epsilon(Tolerance));

    // Check original particle coordinates
    coords << 0.75, 0.75;
    coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    // Compute updated particle location
    REQUIRE_NOTHROW(particle->compute_updated_position(dt));
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
    REQUIRE_NOTHROW(particle->compute_updated_position(dt, true));
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

    // Update porosity
    REQUIRE_NOTHROW(particle->update_porosity(dt));

    // Check porosity
    REQUIRE(particle->porosity() == Approx(0.44).epsilon(Tolerance));

    SECTION("TwoPhase Particle assign state variables") {
      SECTION("Assign state variable fail") {
        mid = 0;
        Json jmaterial;
        jmaterial["density"] = 1000.;
        jmaterial["youngs_modulus"] = 1.0E+7;
        jmaterial["poisson_ratio"] = 0.3;
        jmaterial["softening"] = false;
        jmaterial["friction"] = 0.;
        jmaterial["dilation"] = 0.;
        jmaterial["cohesion"] = 2000.;
        jmaterial["residual_friction"] = 0.;
        jmaterial["residual_dilation"] = 0.;
        jmaterial["residual_cohesion"] = 1000.;
        jmaterial["peak_pdstrain"] = 0.;
        jmaterial["residual_pdstrain"] = 0.;
        jmaterial["tension_cutoff"] = 0.;

        auto mc_material =
            Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()
                ->create("MohrCoulomb2D", std::move(id), jmaterial);
        REQUIRE(mc_material->id() == 0);

        mpm::dense_map state_variables =
            mc_material->initialise_state_variables();
        REQUIRE(state_variables.at("phi") ==
                Approx(jmaterial["friction"]).epsilon(Tolerance));
        REQUIRE(state_variables.at("psi") ==
                Approx(jmaterial["dilation"]).epsilon(Tolerance));
        REQUIRE(state_variables.at("cohesion") ==
                Approx(jmaterial["cohesion"]).epsilon(Tolerance));
        REQUIRE(state_variables.at("epsilon") == Approx(0.).epsilon(Tolerance));
        REQUIRE(state_variables.at("rho") == Approx(0.).epsilon(Tolerance));
        REQUIRE(state_variables.at("theta") == Approx(0.).epsilon(Tolerance));
        REQUIRE(state_variables.at("pdstrain") ==
                Approx(0.).epsilon(Tolerance));

        SECTION("Assign state variables") {
          // Assign material properties
          REQUIRE(particle->assign_material(mc_material) == true);
          // Assign state variables
          REQUIRE(particle->assign_material_state_vars(state_variables,
                                                       mc_material) == true);
          // Assign and read a state variable
          REQUIRE_NOTHROW(particle->assign_state_variable("phi", 30.));
          REQUIRE(particle->state_variable("phi") == 30.);
          // Assign and read pressure though MC does not contain pressure
          REQUIRE_THROWS(particle->assign_pressure(1000));
          REQUIRE(std::isnan(particle->pressure()) == true);
        }

        SECTION("Assign state variables fail on state variables size") {
          // Assign material
          unsigned mid1 = 0;
          // Initialise material
          Json jmaterial1;
          jmaterial1["density"] = 1000.;
          jmaterial1["bulk_modulus"] = 8333333.333333333;
          jmaterial1["dynamic_viscosity"] = 8.9E-4;

          auto newtonian_material =
              Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()
                  ->create("Newtonian2D", std::move(mid1), jmaterial1);

          // Assign material properties
          REQUIRE(particle->assign_material(newtonian_material) == true);
          // Assign state variables
          REQUIRE(particle->assign_material_state_vars(state_variables,
                                                       mc_material) == false);
        }

        SECTION("Assign state variables fail on material id") {
          // Assign material
          unsigned mid1 = 1;
          // Initialise material
          Json jmaterial1;
          jmaterial1["density"] = 1000.;
          jmaterial1["bulk_modulus"] = 8333333.333333333;
          jmaterial1["dynamic_viscosity"] = 8.9E-4;

          auto newtonian_material =
              Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()
                  ->create("Newtonian2D", std::move(mid1), jmaterial1);

          // Assign material properties
          REQUIRE(particle->assign_material(newtonian_material) == true);
          // Assign state variables
          REQUIRE(particle->assign_material_state_vars(state_variables,
                                                       mc_material) == false);
        }
      }
    }
  }

  SECTION("Check assign material to particle") {
    // Add particle
    mpm::Index id = 0;
    coords << 0.75, 0.75;
    auto particle = std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, coords);

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

  SECTION("Check twophase particle properties") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, coords);

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
        std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, coords);

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

  // Check twophase particle's material id maping to nodes
  SECTION("Check twophase particle's material id maping to nodes") {
    // Add particle
    mpm::Index id1 = 0;
    coords << 0.75, 0.75;
    auto particle1 = std::make_shared<mpm::TwoPhaseParticle<Dim>>(id1, coords);

    // Add particle
    mpm::Index id2 = 1;
    coords << 0.25, 0.25;
    auto particle2 = std::make_shared<mpm::TwoPhaseParticle<Dim>>(id2, coords);

    // Element
    std::shared_ptr<mpm::Element<Dim>> element =
        std::make_shared<mpm::QuadrilateralElement<Dim, 4>>();

    // Create cell
    auto cell = std::make_shared<mpm::Cell<Dim>>(10, Nnodes, element);
    // Create vector of nodes and add them to cell
    coords << 0., 0.;
    std::shared_ptr<mpm::NodeBase<Dim>> node0 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);

    coords << 1., 0.;
    std::shared_ptr<mpm::NodeBase<Dim>> node1 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);

    coords << 1., 1.;
    std::shared_ptr<mpm::NodeBase<Dim>> node2 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);

    coords << 0., 1.;
    std::shared_ptr<mpm::NodeBase<Dim>> node3 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);
    std::vector<std::shared_ptr<mpm::NodeBase<Dim>>> nodes = {node0, node1,
                                                              node2, node3};

    for (int j = 0; j < nodes.size(); ++j) cell->add_node(j, nodes[j]);

    // Initialise cell properties and assign cell to particle
    cell->initialise();
    particle1->assign_cell(cell);
    particle2->assign_cell(cell);

    // Assign material 1
    unsigned mid1 = 0;
    // Initialise material 1
    Json jmaterial1;
    jmaterial1["density"] = 1000.;
    jmaterial1["youngs_modulus"] = 1.0E+7;
    jmaterial1["poisson_ratio"] = 0.3;

    auto material1 =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "LinearElastic2D", std::move(mid1), jmaterial1);

    particle1->assign_material(material1);

    // Assign material 2
    unsigned mid2 = 1;
    // Initialise material 2
    Json jmaterial2;
    jmaterial2["density"] = 2000.;
    jmaterial2["youngs_modulus"] = 2.0E+7;
    jmaterial2["poisson_ratio"] = 0.25;

    auto material2 =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "LinearElastic2D", std::move(mid2), jmaterial2);

    particle2->assign_material(material2);

    // Append particle's material id to nodes in cell
    particle1->append_material_id_to_nodes();
    particle2->append_material_id_to_nodes();

    // check if the correct amount of material ids were added to node and if
    // their indexes are correct
    std::vector<unsigned> material_ids = {0, 1};
    for (const auto& node : nodes) {
      REQUIRE(node->material_ids().size() == 2);
      auto mat_ids = node->material_ids();
      unsigned i = 0;
      for (auto mitr = mat_ids.begin(); mitr != mat_ids.end(); ++mitr, ++i)
        REQUIRE(*mitr == material_ids.at(i));
    }
  }
}

//! \brief Check twophase particle class for 3D case
TEST_CASE("TwoPhase Particle is checked for 3D case",
          "[twophase_particle][3D]") {
  // Dimension
  const unsigned Dim = 3;
  // Dimension
  const unsigned Dof = 6;
  // Number of nodes per cell
  const unsigned Nnodes = 8;
  // Number of phases
  const unsigned Nphases = 2;
  // Tolerance
  const double Tolerance = 1.E-7;
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
  // Current time for traction force
  double current_time = 10.0;

  // Coordinates
  Eigen::Vector3d coords;
  coords.setZero();

  //! Check for id = 0
  SECTION("TwoPhase Particle id is zero") {
    mpm::Index id = 0;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, coords);
    REQUIRE(particle->id() == 0);
    REQUIRE(particle->status() == true);
  }

  SECTION("TwoPhase Particle id is positive") {
    //! Check for id is a positive value
    mpm::Index id = std::numeric_limits<mpm::Index>::max();
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, coords);
    REQUIRE(particle->id() == std::numeric_limits<mpm::Index>::max());
    REQUIRE(particle->status() == true);
  }

  //! Construct with id, coordinates and status
  SECTION("TwoPhase Particle with id, coordinates, and status") {
    mpm::Index id = 0;
    bool status = true;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, coords, status);
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
        std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, coords);

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
        std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, coords);

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

  //! Test initialise twophase particle stresses
  SECTION("TwoPhase Particle with initial stress") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    bool status = true;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, coords, status);
    //! Test initialise particle stresses
    Eigen::Matrix<double, 6, 1> stress =
        Eigen::Matrix<double, 6, 1>::Constant(5.7);
    particle->initial_stress(stress);
    REQUIRE(particle->stress().size() == stress.size());
    auto pstress = particle->stress();
    for (unsigned i = 0; i < pstress.size(); ++i)
      REQUIRE(pstress[i] == Approx(stress[i]).epsilon(Tolerance));
  }

  //! Test twophase particles velocity constraints
  SECTION("TwoPhase Particle with velocity constraints") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    bool status = true;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, coords, status);
    // Apply constraints
    particle->apply_particle_velocity_constraints(0, 10.5);
    particle->apply_particle_velocity_constraints(1, -12.5);
    particle->apply_particle_velocity_constraints(2, 14.5);
    particle->apply_particle_velocity_constraints(3, 20.5);
    particle->apply_particle_velocity_constraints(4, -22.5);
    particle->apply_particle_velocity_constraints(5, 24.5);

    // Check apply constraints
    REQUIRE(particle->velocity()(0) == Approx(10.5).epsilon(Tolerance));
    REQUIRE(particle->velocity()(1) == Approx(-12.5).epsilon(Tolerance));
    REQUIRE(particle->velocity()(2) == Approx(14.5).epsilon(Tolerance));
    REQUIRE(particle->liquid_velocity()(0) == Approx(20.5).epsilon(Tolerance));
    REQUIRE(particle->liquid_velocity()(1) == Approx(-22.5).epsilon(Tolerance));
    REQUIRE(particle->liquid_velocity()(2) == Approx(24.5).epsilon(Tolerance));
  }

  //! Test particle, cell and node functions
  SECTION("Test particle, cell and node functions") {
    // Add particle
    mpm::Index id = 0;
    coords << 1.5, 1.5, 1.5;
    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, coords);

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
    // TODO Assert: REQUIRE(particle->compute_shapefn() == false);
    // Compute reference location should throw
    REQUIRE(particle->compute_reference_location() == false);
    // Compute updated particle location should fail
    // TODO Assert: REQUIRE(particle->compute_updated_position(dt) == false);
    // Compute updated particle location from nodal velocity should fail
    // TODO Assert: REQUIRE_NOTHROW(particle->compute_updated_position(dt,
    // true));

    // Compute volume
    // TODO Assert: REQUIRE(particle->compute_volume() == false);
    // Update volume should fail
    // TODO Assert: REQUIRE(particle->update_volume() == false);

    REQUIRE(particle->assign_cell(cell) == true);
    REQUIRE(cell->status() == true);
    REQUIRE(particle->cell_id() == 10);

    // Check if cell is initialised
    REQUIRE(cell->is_initialised() == true);

    // Check compute shape functions of a particle
    REQUIRE_NOTHROW(particle->compute_shapefn());

    // Assign volume
    REQUIRE(particle->assign_volume(0.0) == false);
    REQUIRE(particle->assign_volume(-5.0) == false);
    REQUIRE(particle->assign_volume(2.0) == true);
    // Check volume
    REQUIRE(particle->volume() == Approx(2.0).epsilon(Tolerance));
    // Compute volume
    REQUIRE_NOTHROW(particle->compute_volume());
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
    jmaterial["porosity"] = 0.3;
    jmaterial["k_x"] = 0.001;
    jmaterial["k_y"] = 0.001;
    jmaterial["k_z"] = 0.001;

    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "LinearElastic3D", std::move(mid), jmaterial);

    // Assign liquid material
    unsigned liquid_mid = 1;
    // Initialise material
    Json jmaterial_liquid;
    jmaterial_liquid["density"] = 1000.;
    jmaterial_liquid["bulk_modulus"] = 1.0E+9;
    jmaterial_liquid["mu"] = 0.3;
    jmaterial_liquid["dynamic_viscosity"] = 0.;

    auto liquid_material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "Newtonian3D", std::move(liquid_mid), jmaterial_liquid);

    // Check compute mass before material and volume
    // TODO Assert: REQUIRE(particle->compute_mass() == false);

    // Test compute stress before material assignment
    // TODO Assert: REQUIRE(particle->compute_stress() == false);

    // Assign material properties
    REQUIRE(particle->assign_material(material) == true);
    REQUIRE(particle->assign_material(liquid_material,
                                      mpm::ParticlePhase::Liquid) == true);

    // Check material id from particle
    REQUIRE(particle->material_id() == 0);
    REQUIRE(particle->material_id(mpm::ParticlePhase::Liquid) == 1);

    // Assign porosity
    REQUIRE(particle->assign_porosity() == true);

    // Assign permeability
    REQUIRE(particle->assign_permeability() == true);

    // Compute volume
    REQUIRE_NOTHROW(particle->compute_volume());

    // Compute mass
    REQUIRE_NOTHROW(particle->compute_mass());
    // Mass
    REQUIRE(particle->mass() == Approx(5600.).epsilon(Tolerance));
    REQUIRE(particle->liquid_mass() == Approx(2400.).epsilon(Tolerance));

    // Map particle mass to nodes
    particle->assign_mass(std::numeric_limits<double>::max());
    // TODO Assert: REQUIRE(particle->map_mass_momentum_to_nodes() == false);

    // Map particle pressure to nodes
    // TODO Assert: REQUIRE(particle->map_pressure_to_nodes() == false);

    // Assign mass to nodes
    REQUIRE(particle->compute_reference_location() == true);
    REQUIRE_NOTHROW(particle->compute_shapefn());

    // Check velocity
    Eigen::VectorXd velocity;
    velocity.resize(Dim);
    for (unsigned i = 0; i < velocity.size(); ++i) velocity(i) = i;
    REQUIRE(particle->assign_velocity(velocity) == true);
    REQUIRE(particle->assign_liquid_velocity(velocity) == true);
    for (unsigned i = 0; i < velocity.size(); ++i) {
      REQUIRE(particle->velocity()(i) == Approx(i).epsilon(Tolerance));
      REQUIRE(particle->liquid_velocity()(i) == Approx(i).epsilon(Tolerance));
    }

    REQUIRE_NOTHROW(particle->compute_mass());
    REQUIRE_NOTHROW(particle->map_mass_momentum_to_nodes());

    // TODO Assert: REQUIRE(particle->map_pressure_to_nodes() == false);
    REQUIRE(particle->compute_pressure_smoothing() == false);

    // Values of nodal mass
    std::array<double, 8> nodal_mass{125., 375.,  1125., 375.,
                                     375., 1125., 3375., 1125.};
    // Check nodal mass
    for (unsigned i = 0; i < nodes.size(); ++i) {
      // Solid phase
      REQUIRE(nodes.at(i)->mass(mpm::NodePhase::nSolid) ==
              Approx(nodal_mass.at(i) * (1 - particle->porosity())));
      // Liquid phase
      REQUIRE(
          nodes.at(i)->mass(mpm::NodePhase::nLiquid) ==
          Approx(nodal_mass.at(i) * particle->porosity()).epsilon(Tolerance));
    }

    // Compute nodal velocity
    for (const auto& node : nodes) node->compute_velocity();

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
      for (unsigned j = 0; j < nodal_momentum.cols(); ++j) {
        // Solid phase
        REQUIRE(nodes.at(i)->momentum(mpm::NodePhase::nSolid)(j) ==
                Approx(nodal_momentum(i, j) * (1 - particle->porosity()))
                    .epsilon(Tolerance));
        // Liquid phase
        REQUIRE(nodes.at(i)->momentum(mpm::NodePhase::nLiquid)(j) ==
                Approx(nodal_momentum(i, j) * particle->porosity())
                    .epsilon(Tolerance));
      }

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
      for (unsigned j = 0; j < nodal_velocity.cols(); ++j) {
        // Solid phase
        REQUIRE(nodes.at(i)->velocity(mpm::NodePhase::nSolid)(j) ==
                Approx(nodal_velocity(i, j)).epsilon(Tolerance));
        // Liquid phase
        REQUIRE(nodes.at(i)->velocity(mpm::NodePhase::nLiquid)(j) ==
                Approx(nodal_velocity(i, j)).epsilon(Tolerance));
      }

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
    for (unsigned i = 0; i < nodes.size(); ++i) {
      // Solid phase
      REQUIRE_NOTHROW(nodes.at(i)->update_momentum(
          false, mpm::NodePhase::nSolid,
          nodal_momentum.row(i) * (1 - particle->porosity())));
      // Liquid phase
      REQUIRE_NOTHROW(nodes.at(i)->update_momentum(
          false, mpm::NodePhase::nLiquid,
          nodal_momentum.row(i) * particle->porosity()));
    }

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
    for (const auto& node : nodes) node->compute_velocity();
    // Check nodal velocity
    for (unsigned i = 0; i < nodal_velocity.rows(); ++i)
      for (unsigned j = 0; j < nodal_velocity.cols(); ++j) {
        // Solid phase
        REQUIRE(nodes.at(i)->velocity(mpm::NodePhase::nSolid)(j) ==
                Approx(nodal_velocity(i, j)).epsilon(Tolerance));
        // Liquid phase
        REQUIRE(nodes.at(i)->velocity(mpm::NodePhase::nLiquid)(j) ==
                Approx(nodal_velocity(i, j)).epsilon(Tolerance));
      }

    // Check pressure
    REQUIRE(std::isnan(particle->pressure()) == true);

    // Compute strain
    particle->compute_strain(dt);
    // Strain
    Eigen::Matrix<double, 6, 1> strain;
    strain << 0.00000, 0.07500, 0.40000, -0.02500, 0.35000, -0.05000;

    // Check strains
    for (unsigned i = 0; i < strain.rows(); ++i)
      REQUIRE(particle->strain()(i) == Approx(strain(i)).epsilon(Tolerance));

    // Check volumetric strain at centroid
    double volumetric_strain = 0.5;
    REQUIRE(particle->volumetric_strain_centroid() ==
            Approx(volumetric_strain).epsilon(Tolerance));

    // Check updated pressure
    REQUIRE(std::isnan(particle->pressure()) == true);

    // Update volume strain rate
    REQUIRE(particle->volume() == Approx(8.0).epsilon(Tolerance));
    particle->compute_strain(dt);
    REQUIRE_NOTHROW(particle->update_volume());
    REQUIRE(particle->volume() == Approx(12.0).epsilon(Tolerance));

    // Compute stress
    REQUIRE_NOTHROW(particle->compute_stress());

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

    // Compute pore_pressure
    REQUIRE_NOTHROW(particle->compute_pore_pressure(dt));
    // Check pore pressure
    REQUIRE(particle->pressure(mpm::ParticlePhase::Liquid) ==
            Approx(-1608333333.3333332539).epsilon(Tolerance));

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
        REQUIRE(nodes[i]->external_force(mpm::NodePhase::nSolid)[j] ==
                Approx(body_force(i, j)).epsilon(Tolerance));

    // Check traction force
    double traction = 7.68;
    const unsigned direction = 2;
    // Assign volume
    REQUIRE(particle->assign_volume(0.0) == false);
    REQUIRE(particle->assign_volume(-5.0) == false);
    REQUIRE(particle->assign_volume(2.0) == true);
    // Assign traction to particle
    particle->assign_traction(direction,
                              mfunction->value(current_time) * traction);
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
        REQUIRE(nodes[i]->external_force(mpm::NodePhase::nSolid)[j] ==
                Approx(traction_force(i, j)).epsilon(Tolerance));
    // Reset traction
    particle->assign_traction(direction,
                              mfunction->value(current_time) * -traction);
    // Map traction force
    particle->map_traction_force();
    // Check nodal external force
    for (unsigned i = 0; i < traction_force.rows(); ++i)
      for (unsigned j = 0; j < traction_force.cols(); ++j)
        REQUIRE(nodes[i]->external_force(mpm::NodePhase::nSolid)[j] ==
                Approx(body_force(i, j)).epsilon(Tolerance));

    // Internal force
    Eigen::Matrix<double, 8, 3> internal_force;
    // clang-format off
    internal_force << 402696314.1025640965,  403225160.2564102411, 403826121.7948717475,
                      -402984775.6410256028, 1209771634.6153848171, 1211670673.0769231319,
                     -1208665865.3846151829, -1205637019.2307691574, 3630973557.6923074722,
                      1208185096.1538460255,  -401975160.2564102411, 1210132211.5384614468,
                      1208281249.9999997616, 1208329326.9230768681, -402672275.64102566242,
                      -1208377403.8461537361, 3625276442.3076920509, -1207439903.8461537361,
                      -3624266826.9230761528, -3629026442.3076920509, -3634435096.153845787,
                      3625132211.5384612083,  -1209963942.3076925278, -1212055288.4615385532;
    // clang-format on

    // Map particle internal force
    particle->assign_volume(8.0);
    particle->map_internal_force();

    // Check nodal internal force
    for (unsigned i = 0; i < internal_force.rows(); ++i)
      for (unsigned j = 0; j < internal_force.cols(); ++j)
        REQUIRE(nodes[i]->internal_force(mpm::NodePhase::nMixture)[j] ==
                Approx(internal_force(i, j)).epsilon(Tolerance));

    // Calculate nodal acceleration and velocity
    for (const auto& node : nodes)
      node->compute_acceleration_velocity(mpm::NodePhase::nSolid, dt);

    // Check nodal velocity
    // clang-format off
    nodal_velocity <<  460224.35897435899824,  460829.75457875471329, 461516.16633699642261,
                      -153518.00976800979697,  460867.38461538479896, 461591.42641025653575,
                      -153481.37973137974041, -153093.76434676436475, 461080.60589743591845,
                       460260.98901098914212, -153129.39438339442131, 461009.34582417592173,
                       460297.61904761905316,  460320.93406593421241, -153390.36357753362972,
                      -153444.74969474971294,  460358.56410256412346, -153315.10350427351659,
                      -153408.11965811965638, -153602.58485958489473, -153825.92401709404658,
                       460334.24908424913883, -153638.2148962149513, -153897.1840903541306;
    // clang-format on
    // Check nodal velocity
    for (unsigned i = 0; i < nodal_velocity.rows(); ++i)
      for (unsigned j = 0; j < nodal_velocity.cols(); ++j)
        REQUIRE(nodes[i]->velocity(mpm::NodePhase::nSolid)[j] ==
                Approx(nodal_velocity(i, j)).epsilon(Tolerance));

    // Check nodal acceleration
    Eigen::Matrix<double, 8, 3> nodal_acceleration;
    // clang-format off
    nodal_acceleration << 4602243.5897435899824, 4608287.5457875467837, 4615141.6633699638769,
                         -1535180.0976800979115, 4608653.8461538478732, 4615874.2641025651246,
                         -1534813.7973137972876, -1530967.643467643531, 4610746.0589743591845,
                          4602609.8901098910719, -1531333.9438339441549, 4610013.4582417588681,
                          4602976.1904761902988, 4603159.3406593417749, -1534003.6357753362972,
                         -1534447.4969474971294, 4603525.6410256410018, -1533271.0350427350495,
                         -1534081.1965811965056, -1536095.8485958487727, -1538399.2401709402911,
                          4603342.4908424913883, -1536462.1489621493965, -1539131.840903541306;
    // clang-format on
    // Check nodal acceleration
    for (unsigned i = 0; i < nodal_acceleration.rows(); ++i)
      for (unsigned j = 0; j < nodal_acceleration.cols(); ++j)
        REQUIRE(nodes[i]->acceleration(mpm::NodePhase::nSolid)[j] ==
                Approx(nodal_acceleration(i, j)).epsilon(Tolerance));

    // Approx(nodal_velocity(i, j) / dt).epsilon(Tolerance));

    // Check original particle coordinates
    coords << 1.5, 1.5, 1.5;
    coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    SECTION("Particle pressure smoothing") {
      // Assign material
      unsigned mid1 = 0;
      // Initialise material
      Json jmaterial1;
      jmaterial1["density"] = 1000.;
      jmaterial1["bulk_modulus"] = 8333333.333333333;
      jmaterial1["dynamic_viscosity"] = 8.9E-4;

      auto material1 =
          Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()
              ->create("Newtonian3D", std::move(mid1), jmaterial1);

      // Assign material properties
      REQUIRE(particle->assign_material(material1) == true);

      // Compute volume
      REQUIRE_NOTHROW(particle->compute_volume());

      // Compute mass
      REQUIRE_NOTHROW(particle->compute_mass());
      // Mass
      REQUIRE(particle->mass() == Approx(5600.).epsilon(Tolerance));
      REQUIRE(particle->liquid_mass() == Approx(2400.).epsilon(Tolerance));

      // Map particle mass to nodes
      particle->assign_mass(std::numeric_limits<double>::max());
      // TODO Assert: REQUIRE(particle->map_mass_momentum_to_nodes() == false);

      // Map particle pressure to nodes
      // TODO Assert: REQUIRE(particle->map_pressure_to_nodes() == false);

      // Assign mass to nodes
      REQUIRE(particle->compute_reference_location() == true);
      REQUIRE_NOTHROW(particle->compute_shapefn());

      // Check velocity
      velocity.resize(Dim);
      for (unsigned i = 0; i < velocity.size(); ++i) velocity(i) = i;
      REQUIRE(particle->assign_velocity(velocity) == true);
      for (unsigned i = 0; i < velocity.size(); ++i)
        REQUIRE(particle->velocity()(i) == Approx(i).epsilon(Tolerance));

      REQUIRE_NOTHROW(particle->compute_mass());
      REQUIRE_NOTHROW(particle->map_mass_momentum_to_nodes());

      // Check volumetric strain at centroid
      volumetric_strain = 0.5;
      REQUIRE(particle->dvolumetric_strain() ==
              Approx(volumetric_strain).epsilon(Tolerance));

      // Compute stress
      REQUIRE_NOTHROW(particle->compute_stress());

      REQUIRE(
          particle->pressure() ==
          Approx(-8333333.333333333 * volumetric_strain).epsilon(Tolerance));

      REQUIRE_NOTHROW(particle->map_pressure_to_nodes());
      REQUIRE(particle->compute_pressure_smoothing() == true);
    }

    SECTION("Particle assign state variables") {
      SECTION("Assign state variable fail") {
        mid = 0;
        Json jmaterial;
        jmaterial["density"] = 1000.;
        jmaterial["youngs_modulus"] = 1.0E+7;
        jmaterial["poisson_ratio"] = 0.3;
        jmaterial["softening"] = false;
        jmaterial["friction"] = 0.;
        jmaterial["dilation"] = 0.;
        jmaterial["cohesion"] = 2000.;
        jmaterial["residual_friction"] = 0.;
        jmaterial["residual_dilation"] = 0.;
        jmaterial["residual_cohesion"] = 1000.;
        jmaterial["peak_pdstrain"] = 0.;
        jmaterial["residual_pdstrain"] = 0.;
        jmaterial["tension_cutoff"] = 0.;

        auto mc_material =
            Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()
                ->create("MohrCoulomb3D", std::move(id), jmaterial);
        REQUIRE(mc_material->id() == 0);

        mpm::dense_map state_variables =
            mc_material->initialise_state_variables();
        REQUIRE(state_variables.at("phi") ==
                Approx(jmaterial["friction"]).epsilon(Tolerance));
        REQUIRE(state_variables.at("psi") ==
                Approx(jmaterial["dilation"]).epsilon(Tolerance));
        REQUIRE(state_variables.at("cohesion") ==
                Approx(jmaterial["cohesion"]).epsilon(Tolerance));
        REQUIRE(state_variables.at("epsilon") == Approx(0.).epsilon(Tolerance));
        REQUIRE(state_variables.at("rho") == Approx(0.).epsilon(Tolerance));
        REQUIRE(state_variables.at("theta") == Approx(0.).epsilon(Tolerance));
        REQUIRE(state_variables.at("pdstrain") ==
                Approx(0.).epsilon(Tolerance));

        SECTION("Assign state variables") {
          // Assign material properties
          REQUIRE(particle->assign_material(mc_material) == true);
          // Assign state variables
          REQUIRE(particle->assign_material_state_vars(state_variables,
                                                       mc_material) == true);
          // Assign and read a state variable
          REQUIRE_NOTHROW(particle->assign_state_variable("phi", 30.));
          REQUIRE(particle->state_variable("phi") == 30.);
          // Assign and read pressure though MC does not contain pressure
          REQUIRE_THROWS(particle->assign_pressure(1000));
          REQUIRE(std::isnan(particle->pressure()) == true);
        }

        SECTION("Assign state variables fail on state variables size") {
          // Assign material
          unsigned mid1 = 0;
          // Initialise material
          Json jmaterial1;
          jmaterial1["density"] = 1000.;
          jmaterial1["bulk_modulus"] = 8333333.333333333;
          jmaterial1["dynamic_viscosity"] = 8.9E-4;

          auto newtonian_material =
              Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()
                  ->create("Newtonian3D", std::move(mid1), jmaterial1);

          // Assign material properties
          REQUIRE(particle->assign_material(newtonian_material) == true);
          // Assign state variables
          REQUIRE(particle->assign_material_state_vars(state_variables,
                                                       mc_material) == false);
        }

        SECTION("Assign state variables fail on material id") {
          // Assign material
          unsigned mid1 = 1;
          // Initialise material
          Json jmaterial1;
          jmaterial1["density"] = 1000.;
          jmaterial1["bulk_modulus"] = 8333333.333333333;
          jmaterial1["dynamic_viscosity"] = 8.9E-4;

          auto newtonian_material =
              Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()
                  ->create("Newtonian3D", std::move(mid1), jmaterial1);

          // Assign material properties
          REQUIRE(particle->assign_material(newtonian_material) == true);
          // Assign state variables
          REQUIRE(particle->assign_material_state_vars(state_variables,
                                                       mc_material) == false);
        }
      }
    }

    // Compute updated particle location
    REQUIRE_NOTHROW(particle->compute_updated_position(dt));
    // Check particle velocity
    velocity << 0., 1., 0.5985714286;
    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(particle->velocity()(i) ==
              Approx(velocity(i)).epsilon(Tolerance));

    // Check particle displacement
    Eigen::Vector3d displacement;
    displacement << 0.0, 0.5875, 1.0348571429;
    for (unsigned i = 0; i < displacement.size(); ++i)
      REQUIRE(particle->displacement()(i) ==
              Approx(displacement(i)).epsilon(Tolerance));

    // Updated particle coordinate
    coords << 1.5, 2.0875, 2.5348571429;
    // Check particle coordinates
    coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    // Compute updated particle location based on nodal velocity
    REQUIRE_NOTHROW(particle->compute_updated_position(dt, true));
    // Check particle velocity
    velocity << 0., 5.875, 10.3485714286;
    for (unsigned i = 0; i < velocity.size(); ++i)
      REQUIRE(particle->velocity()(i) ==
              Approx(velocity(i)).epsilon(Tolerance));

    // Check particle displacement
    displacement << 0.0, 1.175, 2.0697142857;
    for (unsigned i = 0; i < displacement.size(); ++i)
      REQUIRE(particle->displacement()(i) ==
              Approx(displacement(i)).epsilon(Tolerance));

    // Updated particle coordinate
    coords << 1.5, 2.675, 3.5697142857;
    // Check particle coordinates
    coordinates = particle->coordinates();
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));
  }

  SECTION("Check assign material to particle") {
    // Add particle
    mpm::Index id = 0;
    coords << 0.75, 0.75, 0.75;
    auto particle = std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, coords);

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
        std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, coords);

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
        std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, coords);

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

  // Check particle's material id maping to nodes
  SECTION("Check particle's material id maping to nodes") {
    // Add particle
    mpm::Index id1 = 0;
    coords << 1.5, 1.5, 1.5;
    auto particle1 = std::make_shared<mpm::TwoPhaseParticle<Dim>>(id1, coords);

    // Add particle
    mpm::Index id2 = 1;
    coords << 0.5, 0.5, 0.5;
    auto particle2 = std::make_shared<mpm::TwoPhaseParticle<Dim>>(id2, coords);

    // Element
    std::shared_ptr<mpm::Element<Dim>> element =
        std::make_shared<mpm::HexahedronElement<Dim, 8>>();

    // Create cell
    auto cell = std::make_shared<mpm::Cell<Dim>>(10, Nnodes, element);
    // Create vector of nodes and add them to cell
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
    std::vector<std::shared_ptr<mpm::NodeBase<Dim>>> nodes = {
        node0, node1, node2, node3, node4, node5, node6, node7};

    for (int j = 0; j < nodes.size(); ++j) cell->add_node(j, nodes[j]);

    // Initialise cell properties and assign cell to particle
    cell->initialise();
    particle1->assign_cell(cell);
    particle2->assign_cell(cell);

    // Assign material 1
    unsigned mid1 = 0;
    // Initialise material 1
    Json jmaterial1;
    jmaterial1["density"] = 1000.;
    jmaterial1["youngs_modulus"] = 1.0E+7;
    jmaterial1["poisson_ratio"] = 0.3;

    auto material1 =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "LinearElastic3D", std::move(mid1), jmaterial1);

    particle1->assign_material(material1);

    // Assign material 2
    unsigned mid2 = 1;
    // Initialise material 2
    Json jmaterial2;
    jmaterial2["density"] = 2000.;
    jmaterial2["youngs_modulus"] = 2.0E+7;
    jmaterial2["poisson_ratio"] = 0.25;

    auto material2 =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "LinearElastic3D", std::move(mid2), jmaterial2);

    particle2->assign_material(material2);

    // Append particle's material id to nodes in cell
    particle1->append_material_id_to_nodes();
    particle2->append_material_id_to_nodes();

    // check if the correct amount of material ids were added to node and if
    // their indexes are correct
    std::vector<unsigned> material_ids = {0, 1};
    for (const auto& node : nodes) {
      REQUIRE(node->material_ids().size() == 2);
      auto mat_ids = node->material_ids();
      unsigned i = 0;
      for (auto mitr = mat_ids.begin(); mitr != mat_ids.end(); ++mitr, ++i)
        REQUIRE(*mitr == material_ids.at(i));
    }
  }
}
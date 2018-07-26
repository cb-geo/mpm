#include <iostream>
#include <limits>

#include "Eigen/Dense"
#include "catch.hpp"
#include "json.hpp"

#include "material/material.h"
// #include "cell.h"
#include "factory.h"
#include "hex_shapefn.h"
#include "node.h"
#include "shapefn.h"

//! \brief Check Bingham class
TEST_CASE("Bingham is checked in 3D", "[material][bingham][3D]") {
  // Tolerance
  const double Tolerance = 1.E-7;

  const unsigned Dim = 3;

  //! Check for id = 0
  SECTION("Bingham id is zero") {
    unsigned id = 0;
    auto material = Factory<mpm::Material<Dim>, unsigned>::instance()->create(
        "Bingham3D", std::move(id));
    REQUIRE(material->id() == 0);
  }

  SECTION("Bingham id is positive") {
    //! Check for id is a positive value
    unsigned id = std::numeric_limits<unsigned>::max();
    auto material = Factory<mpm::Material<Dim>, unsigned>::instance()->create(
        "Bingham3D", std::move(id));
    REQUIRE(material->id() == std::numeric_limits<unsigned>::max());
  }

  //! Read material properties
  SECTION("Bingham check properties") {
    unsigned id = 0;
    auto material = Factory<mpm::Material<Dim>, unsigned>::instance()->create(
        "Bingham3D", std::move(id));
    REQUIRE(material->id() == 0);

    // Initialise material
    Json jmaterial;
    jmaterial["density"] = 1000.;
    jmaterial["youngs_modulus"] = 1.0E+7;
    jmaterial["poisson_ratio"] = 0.3;
    jmaterial["tau0"] = 771.8;
    jmaterial["mu"] = 0.0451;
    jmaterial["strain_cutoff"] = 0.2;

    // Check material status before assigning material property
    REQUIRE(material->status() == false);

    // Get material properties
    REQUIRE(material->property("density") ==
            Approx(std::numeric_limits<double>::max()).epsilon(Tolerance));

    // Check for property that does not exist
    REQUIRE(material->property("noproperty") ==
            Approx(std::numeric_limits<double>::max()).epsilon(Tolerance));

    material->properties(jmaterial);

    // Check material status after assigning material property
    REQUIRE(material->status() == true);

    // Get material properties
    REQUIRE(material->property("density") ==
            Approx(jmaterial["density"]).epsilon(Tolerance));
  }

  SECTION("Bingham check stresses with no strain rate") {
    unsigned id = 0;
    auto material = Factory<mpm::Material<Dim>, unsigned>::instance()->create(
        "Bingham3D", std::move(id));
    REQUIRE(material->id() == 0);

    // Initialise material
    Json jmaterial;
    jmaterial["density"] = 1000.;
    jmaterial["youngs_modulus"] = 1.0E+7;
    jmaterial["poisson_ratio"] = 0.3;
    jmaterial["tau0"] = 771.8;
    jmaterial["mu"] = 0.0451;
    jmaterial["strain_cutoff"] = 0.2;

    material->properties(jmaterial);

    // Check material status after assigning material property
    REQUIRE(material->property_handle() == true);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    REQUIRE(stress(0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.).epsilon(Tolerance));

    // Initialise dstrain
    mpm::Material<Dim>::Vector6d dstrain;
    dstrain.setZero();
    dstrain(0) = 0.0010000;
    dstrain(1) = 0.0005000;
    dstrain(2) = 0.0005000;
    dstrain(3) = 0.0000000;
    dstrain(4) = 0.0000000;
    dstrain(5) = 0.0000000;

    // Add particle
    mpm::Index pid = 0;
    Eigen::Matrix<double, Dim, 1> coords;
    coords << 0.5, 0.5, 0.5;
    auto particle = std::make_shared<mpm::Particle<Dim, 1>>(pid, coords);

    // Coordinates of nodes for the cell
    mpm::Index cell_id = 0;
    const unsigned Dim = 3;
    const unsigned Dof = 6;
    const unsigned Nphases = 1;
    const unsigned Nnodes = 8;

    coords << -2, 2, -2;
    std::shared_ptr<mpm::NodeBase<Dim>> node0 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);
    coords << 2, 2, -2.;
    std::shared_ptr<mpm::NodeBase<Dim>> node1 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);
    coords << 2, 2, 2.;
    std::shared_ptr<mpm::NodeBase<Dim>> node2 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);
    coords << -2, 2, 2.;
    std::shared_ptr<mpm::NodeBase<Dim>> node3 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);
    coords << -2, -2, -2;
    std::shared_ptr<mpm::NodeBase<Dim>> node4 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(4, coords);
    coords << 2, -2, -2;
    std::shared_ptr<mpm::NodeBase<Dim>> node5 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(5, coords);
    coords << 2, -2, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node6 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(6, coords);
    coords << -2, -2, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node7 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(7, coords);

    std::shared_ptr<mpm::ShapeFn<Dim>> shapefn =
        Factory<mpm::ShapeFn<Dim>>::instance()->create("SFH8");

    auto cell = std::make_shared<mpm::Cell<Dim>>(cell_id, Nnodes, shapefn);

    cell->add_node(0, node0);
    cell->add_node(1, node1);
    cell->add_node(2, node2);
    cell->add_node(3, node3);
    cell->add_node(4, node4);
    cell->add_node(5, node5);
    cell->add_node(6, node6);
    cell->add_node(7, node7);

    particle->assign_cell(cell);
    auto check = particle->compute_shapefn();

    // Reset stress
    stress.setZero();

    // Compute updated stress
    stress = material->compute_stress(stress, dstrain, particle.get());

    // Check stressees
    REQUIRE(stress(0) == Approx(16666.66666667).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(16666.66666667).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(16666.66666667).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.000000e+00).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.000000e+00).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.000000e+00).epsilon(Tolerance));
  }

  SECTION("Bingham check stresses with some strain rate") {
    unsigned id = 0;
    auto material = Factory<mpm::Material<Dim>, unsigned>::instance()->create(
        "Bingham3D", std::move(id));
    REQUIRE(material->id() == 0);

    // Initialise material
    Json jmaterial;
    jmaterial["density"] = 1000.;
    jmaterial["youngs_modulus"] = 1.0E+7;
    jmaterial["poisson_ratio"] = 0.3;
    jmaterial["tau0"] = 771.8;
    jmaterial["mu"] = 0.0451;
    jmaterial["strain_cutoff"] = 0.2;

    material->properties(jmaterial);

    // Check material status after assigning material property
    REQUIRE(material->property_handle() == true);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    REQUIRE(stress(0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.).epsilon(Tolerance));

    // Initialise dstrain
    mpm::Material<Dim>::Vector6d dstrain;
    dstrain.setZero();
    dstrain(0) = 0.0010000;
    dstrain(1) = 0.0005000;
    dstrain(2) = 0.0005000;
    dstrain(3) = 0.0000000;
    dstrain(4) = 0.0000000;
    dstrain(5) = 0.0000000;

    // Add particle
    mpm::Index pid = 0;
    Eigen::Matrix<double, Dim, 1> coords;
    coords << 0.5, 0.5, 0.5;
    auto particle = std::make_shared<mpm::Particle<Dim, 1>>(pid, coords);

    // Coordinates of nodes for the cell
    mpm::Index cell_id = 0;
    const unsigned Dim = 3;
    const unsigned Dof = 6;
    const unsigned Nphases = 1;
    const unsigned Nnodes = 8;

    coords << -2, 2, -2;
    std::shared_ptr<mpm::NodeBase<Dim>> node0 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(0, coords);
    coords << 2, 2, -2.;
    std::shared_ptr<mpm::NodeBase<Dim>> node1 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(1, coords);
    coords << 2, 2, 2.;
    std::shared_ptr<mpm::NodeBase<Dim>> node2 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(2, coords);
    coords << -2, 2, 2.;
    std::shared_ptr<mpm::NodeBase<Dim>> node3 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(3, coords);
    coords << -2, -2, -2;
    std::shared_ptr<mpm::NodeBase<Dim>> node4 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(4, coords);
    coords << 2, -2, -2;
    std::shared_ptr<mpm::NodeBase<Dim>> node5 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(5, coords);
    coords << 2, -2, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node6 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(6, coords);
    coords << -2, -2, 2;
    std::shared_ptr<mpm::NodeBase<Dim>> node7 =
        std::make_shared<mpm::Node<Dim, Dof, Nphases>>(7, coords);

    std::shared_ptr<mpm::ShapeFn<Dim>> shapefn =
        Factory<mpm::ShapeFn<Dim>>::instance()->create("SFH8");

    auto cell = std::make_shared<mpm::Cell<Dim>>(cell_id, Nnodes, shapefn);

    cell->add_node(0, node0);
    cell->add_node(1, node1);
    cell->add_node(2, node2);
    cell->add_node(3, node3);
    cell->add_node(4, node4);
    cell->add_node(5, node5);
    cell->add_node(6, node6);
    cell->add_node(7, node7);

    particle->assign_cell(cell);
    auto check = particle->compute_shapefn();

    // Yet to implement: add velocity constraint to have certain strain_rate

    // Reset stress
    stress.setZero();

    // Compute updated stress
    stress = material->compute_stress(stress, dstrain, particle.get());

    // Check stressees
    REQUIRE(stress(0) == Approx(16666.66666667).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(16666.66666667).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(16666.66666667).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.000000e+00).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.000000e+00).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.000000e+00).epsilon(Tolerance));
  }
}

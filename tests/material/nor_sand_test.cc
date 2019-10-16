#include <limits>
#include <iostream>
#include <fstream>

#include "Eigen/Dense"
#include "catch.hpp"
#include "json.hpp"

#include "cell.h"
#include "material/material.h"
#include "node.h"
#include "particle.h"

//! Check NorSand class in 3D
TEST_CASE("NorSand is checked in 3D", "[material][NorSand][3D]") {
  // Tolerance
  const double Tolerance = 1.E-7;

  const unsigned Dim = 3;

  // Add particle
  mpm::Index pid = 0;
  Eigen::Matrix<double, Dim, 1> coords;
  coords.setZero();
  auto particle = std::make_shared<mpm::Particle<Dim, 1>>(pid, coords);

  // Initialise material
  Json jmaterial;
  jmaterial["density"] = 1800.;
  jmaterial["youngs_modulus"] = 1.0E+7;
  jmaterial["poisson_ratio"] = 0.3;
  jmaterial["shear_modulus_constant"] = 2500000.;
  jmaterial["shear_modulus_exponent"] = 0.5;
  jmaterial["reference_pressure"] = 1000.;
  jmaterial["M"] = 1.33;
  jmaterial["N"] = 0.3;
  jmaterial["e_min"] = 0.542;
  jmaterial["e_max"] = 1.000;
  jmaterial["crushing_pressure"] = 10000000.0;
  jmaterial["chi"] = 3.5;
  jmaterial["hardening_modulus"] = 200.0;
  jmaterial["void_ratio_initial"] = 0.85;
  jmaterial["p_image_initial"] = 8701.46;

  //! Check for id = 0
  SECTION("NorSand id is zero") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "NorSand3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);
  }

  SECTION("NorSand id is positive") {
    //! Check for id is a positive value
    unsigned id = std::numeric_limits<unsigned>::max();
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "NorSand3D", std::move(id), jmaterial);
    REQUIRE(material->id() == std::numeric_limits<unsigned>::max());
  }

  //! Check failed initialisation
  SECTION("NorSand failed initialisation") {
    unsigned id = 0;
    // Initialise material
    Json jmaterial;
    jmaterial["density"] = 1000.;
    jmaterial["poisson_ratio"] = 0.3;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "NorSand3D", std::move(id), jmaterial);
  }

  //! Check material properties
  SECTION("NorSand check material properties") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "NorSand3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Get material properties
    REQUIRE(material->property("density") ==
            Approx(jmaterial["density"]).epsilon(Tolerance));
    REQUIRE(material->property("youngs_modulus") ==
            Approx(jmaterial["youngs_modulus"]).epsilon(Tolerance));
    REQUIRE(material->property("poisson_ratio") ==
            Approx(jmaterial["poisson_ratio"]).epsilon(Tolerance));
    REQUIRE(material->property("shear_modulus_constant") ==
            Approx(jmaterial["shear_modulus_constant"]).epsilon(Tolerance));
    REQUIRE(material->property("shear_modulus_exponent") ==
            Approx(jmaterial["shear_modulus_exponent"]).epsilon(Tolerance));
    REQUIRE(material->property("reference_pressure") ==
            Approx(jmaterial["reference_pressure"]).epsilon(Tolerance));
    REQUIRE(material->property("M") ==
            Approx(jmaterial["M"]).epsilon(Tolerance));
    REQUIRE(material->property("N") ==
            Approx(jmaterial["N"]).epsilon(Tolerance));
    REQUIRE(material->property("e_min") ==
            Approx(jmaterial["e_min"]).epsilon(Tolerance));
    REQUIRE(material->property("e_max") ==
            Approx(jmaterial["e_max"]).epsilon(Tolerance));
    REQUIRE(material->property("crushing_pressure") ==
            Approx(jmaterial["crushing_pressure"]).epsilon(Tolerance));
    REQUIRE(material->property("chi") ==
            Approx(jmaterial["chi"]).epsilon(Tolerance));
    REQUIRE(material->property("hardening_modulus") ==
            Approx(jmaterial["hardening_modulus"]).epsilon(Tolerance));
    REQUIRE(material->property("void_ratio_initial") ==
            Approx(jmaterial["void_ratio_initial"]).epsilon(Tolerance));
    REQUIRE(material->property("p_image_initial") ==
            Approx(jmaterial["p_image_initial"]).epsilon(Tolerance));

    // Check if state variable is initialised
    SECTION("State variable is initialised") {
      mpm::dense_map state_variables = material->initialise_state_variables();
      REQUIRE(state_variables.empty() == false);
      REQUIRE(state_variables.at("p") == Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("q") == Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("void_ratio") == Approx(jmaterial["void_ratio_initial"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("e_image") == Approx(0.9350064171).epsilon(Tolerance));
      REQUIRE(state_variables.at("p_image") == Approx(8701.46).epsilon(Tolerance));
      REQUIRE(state_variables.at("psi_image") == Approx(-0.0850064171).epsilon(Tolerance));

    }
  }

  SECTION("NorSand check stresses") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "NorSand3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    REQUIRE(stress(0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.).epsilon(Tolerance));

    stress(0) = -20000;
    stress(1) = -20000;
    stress(2) = -20000;
    stress(3) = 0;
    stress(4) = 0;
    stress(5) = 0;

    // Initialise dstrain
    mpm::Material<Dim>::Vector6d dstrain;
    dstrain.setZero();
    dstrain(0) = 0.00050000;
    dstrain(1) = -0.0010000;
    dstrain(2) = 0.0005000;
    dstrain(3) = 0.0000000;
    dstrain(4) = 0.0000000;
    dstrain(5) = 0.0000000;

    // Compute updated stress
    mpm::dense_map state_vars = material->initialise_state_variables();
    std::ofstream myfile;
    myfile.open("norsand_stress.txt");
    for (unsigned i = 0; i < 400; ++i) {
      stress =
          material->compute_stress(stress, dstrain, particle.get(), &state_vars);
      myfile << stress(0) << '\t' << stress(1) << '\t' << stress(2) << '\t' << stress(3) << '\t' << stress(4) << '\t' << stress(5) << '\n';      
    }
    myfile.close();

    // // Check stressees
    // REQUIRE(stress(0) == Approx(1.63461538461538e+04).epsilon(Tolerance));
    // REQUIRE(stress(1) == Approx(1.25000000000000e+04).epsilon(Tolerance));
    // REQUIRE(stress(2) == Approx(0.86538461538462e+04).epsilon(Tolerance));
    // REQUIRE(stress(3) == Approx(0.000000e+00).epsilon(Tolerance));
    // REQUIRE(stress(4) == Approx(0.000000e+00).epsilon(Tolerance));
    // REQUIRE(stress(5) == Approx(0.000000e+00).epsilon(Tolerance));

    // // Initialise strain
    // strain(0) = 0.0010000;
    // strain(1) = 0.0005000;
    // strain(2) = 0.0000000;
    // strain(3) = 0.0000100;
    // strain(4) = 0.0000000;
    // strain(5) = 0.0000000;

    // // Reset stress
    // stress.setZero();

    // // Compute updated stress
    // stress =
    //     material->compute_stress(stress, strain, particle.get(), &state_vars);

    // // Check stressees
    // REQUIRE(stress(0) == Approx(1.63461538461538e+04).epsilon(Tolerance));
    // REQUIRE(stress(1) == Approx(1.25000000000000e+04).epsilon(Tolerance));
    // REQUIRE(stress(2) == Approx(0.86538461538462e+04).epsilon(Tolerance));
    // REQUIRE(stress(3) == Approx(3.84615384615385e+01).epsilon(Tolerance));
    // REQUIRE(stress(4) == Approx(0.00000000000000e+00).epsilon(Tolerance));
    // REQUIRE(stress(5) == Approx(0.00000000000000e+00).epsilon(Tolerance));
  }
}



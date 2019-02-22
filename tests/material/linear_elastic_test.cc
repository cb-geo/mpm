#include <limits>

#include "Eigen/Dense"
#include "catch.hpp"
#include "json.hpp"

#include "cell.h"
#include "material/material.h"
#include "node.h"
#include "particle.h"

//! Check linearelastic class in 2D
TEST_CASE("LinearElastic is checked in 2D", "[material][linear_elastic][2D]") {
  // Tolerance
  const double Tolerance = 1.E-7;

  const unsigned Dim = 2;

  // Add particle
  mpm::Index pid = 0;
  Eigen::Matrix<double, Dim, 1> coords;
  coords.setZero();
  auto particle = std::make_shared<mpm::Particle<Dim, 1>>(pid, coords);

  // Initialise material
  Json jmaterial;
  jmaterial["density"] = 1000.;
  jmaterial["youngs_modulus"] = 1.0E+7;
  jmaterial["poisson_ratio"] = 0.3;

  //! Check for id = 0
  SECTION("LinearElastic id is zero") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "LinearElastic2D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);
  }

  SECTION("LinearElastic id is positive") {
    //! Check for id is a positive value
    unsigned id = std::numeric_limits<unsigned>::max();
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "LinearElastic2D", std::move(id), jmaterial);
    REQUIRE(material->id() == std::numeric_limits<unsigned>::max());
  }

  //! Check failed initialisation
  SECTION("LinearElastic failed initialisation") {
    unsigned id = 0;
    // Initialise material
    Json jmaterial;
    jmaterial["density"] = 1000.;
    jmaterial["poisson_ratio"] = 0.3;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "LinearElastic2D", std::move(id), jmaterial);
  }

  //! Check material properties
  SECTION("LinearElastic check material properties") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "LinearElastic2D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Get material properties
    REQUIRE(material->property("density") ==
            Approx(jmaterial["density"]).epsilon(Tolerance));
    REQUIRE(material->property("youngs_modulus") ==
            Approx(jmaterial["youngs_modulus"]).epsilon(Tolerance));
    REQUIRE(material->property("poisson_ratio") ==
            Approx(jmaterial["poisson_ratio"]).epsilon(Tolerance));

    // Check if state variable is initialised
    SECTION("State variable is initialised") {
      mpm::dense_map state_variables = material->initialise_state_variables();
      REQUIRE(state_variables.empty() == true);
    }
  }

  //! Check thermodynamic pressure
  SECTION("LinearElastic check thermodynamic pressure") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "LinearElastic2D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Get material properties
    REQUIRE(material->property("density") ==
            Approx(jmaterial["density"]).epsilon(Tolerance));

    // Calculate modulus values
    const double K = 8333333.333333333;
    const double G = 3846153.846153846;
    const double a1 = 13461538.461566667;
    const double a2 = 5769230.769166667;

    // Calculate pressure
    const double volumetric_strain = 1.0E-5;
    REQUIRE(material->thermodynamic_pressure(volumetric_strain) ==
            Approx(-K * volumetric_strain).epsilon(Tolerance));
  }

  SECTION("LinearElastic check stresses") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "LinearElastic2D", std::move(id), jmaterial);
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

    // Initialise strain
    mpm::Material<Dim>::Vector6d strain;
    strain.setZero();
    strain(0) = 0.0010000;
    strain(1) = 0.0005000;
    strain(2) = 0.0000000;
    strain(3) = 0.0000000;
    strain(4) = 0.0000000;
    strain(5) = 0.0000000;

    // Compute updated stress
    mpm::dense_map state_vars = material->initialise_state_variables();
    stress =
        material->compute_stress(stress, strain, particle.get(), &state_vars);

    // Check stressees
    REQUIRE(stress(0) == Approx(1.63461538461538e+04).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(1.25000000000000e+04).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(0.86538461538462e+04).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.000000e+00).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.000000e+00).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.000000e+00).epsilon(Tolerance));

    // Initialise strain
    strain(0) = 0.0010000;
    strain(1) = 0.0005000;
    strain(2) = 0.0000000;
    strain(3) = 0.0000100;
    strain(4) = 0.0000000;
    strain(5) = 0.0000000;

    // Reset stress
    stress.setZero();

    // Compute updated stress
    stress =
        material->compute_stress(stress, strain, particle.get(), &state_vars);

    // Check stressees
    REQUIRE(stress(0) == Approx(1.63461538461538e+04).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(1.25000000000000e+04).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(0.86538461538462e+04).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(3.84615384615385e+01).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.00000000000000e+00).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.00000000000000e+00).epsilon(Tolerance));
  }
}

//! Check linearelastic class in 3D
TEST_CASE("LinearElastic is checked in 3D", "[material][linear_elastic][3D]") {
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
  jmaterial["density"] = 1000.;
  jmaterial["youngs_modulus"] = 1.0E+7;
  jmaterial["poisson_ratio"] = 0.3;

  //! Check for id = 0
  SECTION("LinearElastic id is zero") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "LinearElastic3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);
  }

  SECTION("LinearElastic id is positive") {
    //! Check for id is a positive value
    unsigned id = std::numeric_limits<unsigned>::max();
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "LinearElastic3D", std::move(id), jmaterial);
    REQUIRE(material->id() == std::numeric_limits<unsigned>::max());
  }

  //! Check failed initialisation
  SECTION("LinearElastic failed initialisation") {
    unsigned id = 0;
    // Initialise material
    Json jmaterial;
    jmaterial["density"] = 1000.;
    jmaterial["poisson_ratio"] = 0.3;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "LinearElastic3D", std::move(id), jmaterial);
  }

  //! Check material properties
  SECTION("LinearElastic check material properties") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "LinearElastic3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Get material properties
    REQUIRE(material->property("density") ==
            Approx(jmaterial["density"]).epsilon(Tolerance));
    REQUIRE(material->property("youngs_modulus") ==
            Approx(jmaterial["youngs_modulus"]).epsilon(Tolerance));
    REQUIRE(material->property("poisson_ratio") ==
            Approx(jmaterial["poisson_ratio"]).epsilon(Tolerance));

    // Check if state variable is initialised
    SECTION("State variable is initialised") {
      mpm::dense_map state_variables = material->initialise_state_variables();
      REQUIRE(state_variables.empty() == true);
    }
  }

  //! Check thermodynamic pressure
  SECTION("LinearElastic check thermodynamic pressure") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "LinearElastic3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Get material properties
    REQUIRE(material->property("density") ==
            Approx(jmaterial["density"]).epsilon(Tolerance));

    // Calculate modulus values
    const double K = 8333333.333333333;
    const double G = 3846153.846153846;
    const double a1 = 13461538.461566667;
    const double a2 = 5769230.769166667;

    // Calculate pressure
    const double volumetric_strain = 1.0E-5;
    REQUIRE(material->thermodynamic_pressure(volumetric_strain) ==
            Approx(-K * volumetric_strain).epsilon(Tolerance));
  }

  SECTION("LinearElastic check stresses") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "LinearElastic3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    //    mpm::Material<Dim>::Matrix6x6 de = material->elastic_tensor();

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    REQUIRE(stress(0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.).epsilon(Tolerance));

    // Initialise strain
    mpm::Material<Dim>::Vector6d strain;
    strain.setZero();
    strain(0) = 0.0010000;
    strain(1) = 0.0005000;
    strain(2) = 0.0005000;
    strain(3) = 0.0000000;
    strain(4) = 0.0000000;
    strain(5) = 0.0000000;

    // Compute updated stress
    mpm::dense_map state_vars = material->initialise_state_variables();
    stress =
        material->compute_stress(stress, strain, particle.get(), &state_vars);

    // Check stressees
    REQUIRE(stress(0) == Approx(1.92307692307333e+04).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(1.53846153845333e+04).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(1.53846153845333e+04).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.000000e+00).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.000000e+00).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.000000e+00).epsilon(Tolerance));

    // Initialise strain
    strain(0) = 0.0010000;
    strain(1) = 0.0005000;
    strain(2) = 0.0005000;
    strain(3) = 0.0000100;
    strain(4) = 0.0000200;
    strain(5) = 0.0000300;

    // Reset stress
    stress.setZero();

    // Compute updated stress
    stress =
        material->compute_stress(stress, strain, particle.get(), &state_vars);

    // Check stressees
    REQUIRE(stress(0) == Approx(1.92307692307333e+04).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(1.53846153845333e+04).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(1.53846153845333e+04).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(3.84615384615385e+01).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(7.69230769230769e+01).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(1.15384615384615e+02).epsilon(Tolerance));
  }
}

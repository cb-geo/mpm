#include <limits>

#include "Eigen/Dense"
#include "catch.hpp"
#include "json.hpp"

#include "cell.h"
#include "material/material.h"
#include "material/mohr_coulomb.h"
#include "node.h"
#include "particle.h"

//! Check MohrCoulomb (without softening) class in 2D
TEST_CASE("MohrCoulomb (without softening) is checked in 2D",
          "[material][mohr_coulomb][2D]") {
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
  jmaterial["softening"] = false;
  jmaterial["friction"] = 0.;
  jmaterial["dilation"] = 0.;
  jmaterial["cohesion"] = 2000.;
  jmaterial["residual_friction"] = 0.;
  jmaterial["residual_dilation"] = 0.;
  jmaterial["residual_cohesion"] = 1000.;
  jmaterial["peak_epds"] = 0.;
  jmaterial["critical_epds"] = 0.;
  jmaterial["tension_cutoff"] = 0.;

  //! Check for id = 0
  SECTION("MohrCoulomb id is zero") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "MohrCoulomb2D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);
  }

  SECTION("MohrCoulomb id is positive") {
    //! Check for id is a positive value
    unsigned id = std::numeric_limits<unsigned>::max();
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "MohrCoulomb2D", std::move(id), jmaterial);
    REQUIRE(material->id() == std::numeric_limits<unsigned>::max());
  }

  //! Check failed initialisation
  SECTION("MohrCoulomb failed initialisation") {
    unsigned id = 0;
    // Initialise material
    Json jmaterial;
    jmaterial["density"] = 1000.;
    jmaterial["poisson_ratio"] = 0.3;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "MohrCoulomb2D", std::move(id), jmaterial);
  }

  //! Check material properties
  SECTION("MohrCoulomb check material properties") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "MohrCoulomb2D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Get material properties
    REQUIRE(material->property("density") ==
            Approx(jmaterial["density"]).epsilon(Tolerance));
    REQUIRE(material->property("youngs_modulus") ==
            Approx(jmaterial["youngs_modulus"]).epsilon(Tolerance));
    REQUIRE(material->property("poisson_ratio") ==
            Approx(jmaterial["poisson_ratio"]).epsilon(Tolerance));
    REQUIRE(material->property("friction") ==
            Approx(jmaterial["friction"]).epsilon(Tolerance));
    REQUIRE(material->property("dilation") ==
            Approx(jmaterial["dilation"]).epsilon(Tolerance));
    REQUIRE(material->property("cohesion") ==
            Approx(jmaterial["cohesion"]).epsilon(Tolerance));
    REQUIRE(material->property("tension_cutoff") ==
            Approx(jmaterial["tension_cutoff"]).epsilon(Tolerance));

    // Check if state variable is initialised
    SECTION("State variable is initialised") {
      mpm::dense_map state_variables = material->initialise_state_variables();
      REQUIRE(state_variables.at("phi") ==
              Approx(jmaterial["friction"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("psi") ==
              Approx(jmaterial["dilation"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("cohesion") ==
              Approx(jmaterial["cohesion"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("j2") == Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("j3") == Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("epsilon") == Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("rho") == Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("theta") == Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("epds") == Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain0") ==
              Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain1") ==
              Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain2") ==
              Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain3") ==
              Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain4") ==
              Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain5") ==
              Approx(0.).epsilon(Tolerance));
    }
  }

  //! Check thermodynamic pressure
  SECTION("MohrCoulomb check thermodynamic pressure") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "MohrCoulomb2D", std::move(id), jmaterial);
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
            Approx(0.).epsilon(Tolerance));
  }

  //! Check stress invariants
  SECTION("MohrCoulomb check stress invariants and yield state") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "MohrCoulomb2D", std::move(id), jmaterial);

    auto mohr_coulomb = std::make_shared<mpm::MohrCoulomb<Dim>>(id, jmaterial);

    REQUIRE(material->id() == 0);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    stress(0) = -5000.;
    stress(1) = -6000.;
    stress(2) = -7000.;
    stress(3) = -1000.;

    // Initialise state variables
    mpm::dense_map state_variables = material->initialise_state_variables();
    // Check if stress invariants is computed correctly
    REQUIRE(mohr_coulomb->compute_stress_invariants(stress, &state_variables) ==
            true);
    REQUIRE(state_variables.at("phi") ==
            Approx(jmaterial["friction"]).epsilon(Tolerance));
    REQUIRE(state_variables.at("psi") ==
            Approx(jmaterial["dilation"]).epsilon(Tolerance));
    REQUIRE(state_variables.at("cohesion") ==
            Approx(jmaterial["cohesion"]).epsilon(Tolerance));
    REQUIRE(state_variables.at("j2") == Approx(2000000.).epsilon(Tolerance));
    REQUIRE(state_variables.at("j3") == Approx(1000000000.).epsilon(Tolerance));
    REQUIRE(state_variables.at("epsilon") ==
            Approx(-10392.30484541).epsilon(Tolerance));
    REQUIRE(state_variables.at("rho") == Approx(2000.).epsilon(Tolerance));
    REQUIRE(state_variables.at("theta") ==
            Approx(0.13545926).epsilon(Tolerance));
    REQUIRE(state_variables.at("epds") == Approx(0.).epsilon(Tolerance));
    REQUIRE(state_variables.at("plastic_strain0") ==
            Approx(0.).epsilon(Tolerance));
    REQUIRE(state_variables.at("plastic_strain1") ==
            Approx(0.).epsilon(Tolerance));
    REQUIRE(state_variables.at("plastic_strain2") ==
            Approx(0.).epsilon(Tolerance));
    REQUIRE(state_variables.at("plastic_strain3") ==
            Approx(0.).epsilon(Tolerance));
    REQUIRE(state_variables.at("plastic_strain4") ==
            Approx(0.).epsilon(Tolerance));
    REQUIRE(state_variables.at("plastic_strain5") ==
            Approx(0.).epsilon(Tolerance));

    // Initialise values of yield functions
    Eigen::Matrix<double, 2, 1> yield_function;
    auto yield_type =
        mohr_coulomb->compute_yield_state(&yield_function, &state_variables);
    // Check if yield function and yield state is computed correctly
    REQUIRE(yield_function(0) == Approx(-4381.96601125).epsilon(Tolerance));
    REQUIRE(yield_function(1) == Approx(-690.98300563).epsilon(Tolerance));
    REQUIRE(yield_type == mohr_coulomb->FailureState::Elastic);
  }
}

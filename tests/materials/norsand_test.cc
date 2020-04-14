#include <limits>

#include <cmath>

#include "Eigen/Dense"
#include "catch.hpp"
#include "json.hpp"

#include "cell.h"
#include "material.h"
#include "node.h"
#include "particle.h"

//! Check NorSand class in 3D non-bonded model
TEST_CASE("NorSand is checked in 3D non-bonded model",
          "[material][NorSand][3D]") {
  // Tolerance
  const double Tolerance = 1.E-7;

  const unsigned Dim = 3;

  // Add particle
  mpm::Index pid = 0;
  Eigen::Matrix<double, Dim, 1> coords;
  coords.setZero();
  auto particle = std::make_shared<mpm::Particle<Dim>>(pid, coords);

  // Initialise material
  Json jmaterial;
  jmaterial["density"] = 1800.;
  jmaterial["poisson_ratio"] = 0.3;
  jmaterial["reference_pressure"] = 1000.;
  jmaterial["friction_cs"] = 30;
  jmaterial["N"] = 0.3;
  jmaterial["lambda"] = 0.1;
  jmaterial["kappa"] = 0.03;
  jmaterial["gamma"] = 1.3;
  jmaterial["chi"] = 3.5;
  jmaterial["hardening_modulus"] = 200.0;
  jmaterial["void_ratio_initial"] = 0.85;
  jmaterial["p_image_initial"] = 87014.6;
  jmaterial["bond_model"] = false;
  jmaterial["p_cohesion_initial"] = 0.0;
  jmaterial["p_dilation_initial"] = 0.0;
  jmaterial["m_cohesion"] = 0.0;
  jmaterial["m_dilation"] = 0.0;
  jmaterial["m_modulus"] = 0.0;

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

  //! Check material properties
  SECTION("NorSand check material properties") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "NorSand3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Get material properties
    REQUIRE(material->template property<double>("density") ==
            Approx(jmaterial.at("density")).epsilon(Tolerance));
    REQUIRE(material->template property<double>("poisson_ratio") ==
            Approx(jmaterial.at("poisson_ratio")).epsilon(Tolerance));
    REQUIRE(material->template property<double>("reference_pressure") ==
            Approx(jmaterial.at("reference_pressure")).epsilon(Tolerance));
    REQUIRE(material->template property<double>("friction_cs") ==
            Approx(jmaterial.at("friction_cs")).epsilon(Tolerance));
    REQUIRE(material->template property<double>("N") ==
            Approx(jmaterial.at("N")).epsilon(Tolerance));
    REQUIRE(material->template property<double>("lambda") ==
            Approx(jmaterial.at("lambda")).epsilon(Tolerance));
    REQUIRE(material->template property<double>("kappa") ==
            Approx(jmaterial.at("kappa")).epsilon(Tolerance));
    REQUIRE(material->template property<double>("gamma") ==
            Approx(jmaterial.at("gamma")).epsilon(Tolerance));
    REQUIRE(material->template property<double>("chi") ==
            Approx(jmaterial.at("chi")).epsilon(Tolerance));
    REQUIRE(material->template property<double>("hardening_modulus") ==
            Approx(jmaterial.at("hardening_modulus")).epsilon(Tolerance));
    REQUIRE(material->template property<double>("void_ratio_initial") ==
            Approx(jmaterial.at("void_ratio_initial")).epsilon(Tolerance));
    REQUIRE(material->template property<double>("p_image_initial") ==
            Approx(jmaterial.at("p_image_initial")).epsilon(Tolerance));
    REQUIRE(material->template property<double>("p_cohesion_initial") ==
            Approx(jmaterial.at("p_cohesion_initial")).epsilon(Tolerance));
    REQUIRE(material->template property<double>("p_dilation_initial") ==
            Approx(jmaterial.at("p_dilation_initial")).epsilon(Tolerance));
    REQUIRE(material->template property<double>("m_cohesion") ==
            Approx(jmaterial.at("m_cohesion")).epsilon(Tolerance));
    REQUIRE(material->template property<double>("m_dilation") ==
            Approx(jmaterial.at("m_dilation")).epsilon(Tolerance));
    REQUIRE(material->template property<double>("m_modulus") ==
            Approx(jmaterial.at("m_modulus")).epsilon(Tolerance));
    REQUIRE(material->template property<bool>("bond_model") ==
            jmaterial.at("bond_model"));

    // Check if state variable is initialised
    SECTION("State variable is initialised") {
      const double sin_friction_cs =
          sin(jmaterial.at("friction_cs").template get<double>() * M_PI / 180.);

      mpm::dense_map state_variables = material->initialise_state_variables();
      REQUIRE(state_variables.empty() == false);
      REQUIRE(state_variables.at("M_theta") ==
              Approx((6 * sin_friction_cs) / (3 - sin_friction_cs))
                  .epsilon(Tolerance));
      REQUIRE(state_variables.at("void_ratio") ==
              Approx(jmaterial["void_ratio_initial"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("e_image") ==
              Approx(0.8533924079).epsilon(Tolerance));
      REQUIRE(state_variables.at("p_image") ==
              Approx(87014.6).epsilon(Tolerance));
      REQUIRE(state_variables.at("p_cohesion") ==
              Approx(jmaterial["p_cohesion_initial"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("p_dilation") ==
              Approx(jmaterial["p_dilation_initial"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("pdstrain") == Approx(0.0).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain0") ==
              Approx(0.0).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain1") ==
              Approx(0.0).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain2") ==
              Approx(0.0).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain3") ==
              Approx(0.0).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain4") ==
              Approx(0.0).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain5") ==
              Approx(0.0).epsilon(Tolerance));
    }
  }

  SECTION("NorSand check elastic stresses") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "NorSand3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    stress(0) = -200000.;
    stress(1) = -200000.;
    stress(2) = -200000.;
    REQUIRE(stress(0) == Approx(-200000.).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(-200000.).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(-200000.).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.).epsilon(Tolerance));

    // Initialise dstrain
    mpm::Material<Dim>::Vector6d dstrain;
    dstrain.setZero();
    dstrain(0) = -0.00010000;
    dstrain(1) = 0.00005000;
    dstrain(2) = 0.00005000;
    dstrain(3) = -0.00000000;
    dstrain(4) = -0.00000000;
    dstrain(5) = -0.00000000;

    // Compute updated stress
    mpm::dense_map state_vars = material->initialise_state_variables();
    stress =
        material->compute_stress(stress, dstrain, particle.get(), &state_vars);

    // Check stresses
    REQUIRE(stress(0) == Approx(-2.011384615384615E+05).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(-1.994307692307692E+05).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(-1.994307692307692E+05).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.000000).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.000000).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.000000).epsilon(Tolerance));

    // Initialise dstrain
    dstrain(0) = -0.00010000;
    dstrain(1) = -0.00005000;
    dstrain(2) = -0.00005000;
    dstrain(3) = -0.00000100;
    dstrain(4) = -0.00000200;
    dstrain(5) = -0.00000300;

    // Reset stress
    stress.setZero();
    stress(0) = -200000.;
    stress(1) = -200000.;
    stress(2) = -200000.;
    REQUIRE(stress(0) == Approx(-200000.).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(-200000.).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(-200000.).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.).epsilon(Tolerance));

    // Compute updated stress
    stress =
        material->compute_stress(stress, dstrain, particle.get(), &state_vars);

    // Check stresses
    REQUIRE(stress(0) == Approx(-2.003794871794872E+05).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(-1.998102564102564E+05).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(-1.998102564102564E+05).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(-0.000056923076923E+05).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(-0.000113846153846E+05).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(-0.000170769230769E+05).epsilon(Tolerance));
  }

  SECTION("NorSand check undrained stresses after yielding with bonded model") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "NorSand3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    stress(0) = -200000;
    stress(1) = -200000;
    stress(2) = -200000;
    stress(3) = 0;
    stress(4) = 0;
    stress(5) = 0;
    REQUIRE(stress(0) == Approx(-200000.).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(-200000.).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(-200000.).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.).epsilon(Tolerance));

    // Initialise dstrain
    mpm::Material<Dim>::Vector6d dstrain;
    dstrain.setZero();
    dstrain(0) = -0.010000;
    dstrain(1) = 0.005000;
    dstrain(2) = 0.005000;
    dstrain(3) = 0.000000;
    dstrain(4) = 0.000000;
    dstrain(5) = 0.000000;

    // Compute updated stress two times to get yielding
    mpm::dense_map state_vars = material->initialise_state_variables();
    stress =
        material->compute_stress(stress, dstrain, particle.get(), &state_vars);
    stress =
        material->compute_stress(stress, dstrain, particle.get(), &state_vars);

    // Check stresses
    REQUIRE(stress(0) == Approx(-3.13851220054474E+05).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(-1.42967854828796E+05).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(-1.42967854828796E+05).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.000000).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.000000).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.000000).epsilon(Tolerance));

    // Check state variables
    REQUIRE(state_vars.empty() == false);
    REQUIRE(state_vars.at("M_theta") == Approx(1.2).epsilon(Tolerance));
    REQUIRE(state_vars.at("void_ratio") ==
            Approx(jmaterial["void_ratio_initial"]).epsilon(Tolerance));
    REQUIRE(state_vars.at("e_image") ==
            Approx(0.7973365716).epsilon(Tolerance));
    REQUIRE(state_vars.at("p_image") ==
            Approx(152419.1488475165).epsilon(Tolerance));
    REQUIRE(state_vars.at("p_cohesion") == Approx(0.000).epsilon(Tolerance));
    REQUIRE(state_vars.at("p_dilation") == Approx(0.000).epsilon(Tolerance));
    REQUIRE(state_vars.at("pdstrain") ==
            Approx(0.0099933165).epsilon(Tolerance));
    REQUIRE(state_vars.at("plastic_strain0") ==
            Approx(0.009995236).epsilon(Tolerance));
    REQUIRE(state_vars.at("plastic_strain1") ==
            Approx(-0.0049947387).epsilon(Tolerance));
    REQUIRE(state_vars.at("plastic_strain2") ==
            Approx(-0.0049947387).epsilon(Tolerance));
    REQUIRE(state_vars.at("plastic_strain3") == Approx(0.0).epsilon(Tolerance));
    REQUIRE(state_vars.at("plastic_strain4") == Approx(0.0).epsilon(Tolerance));
    REQUIRE(state_vars.at("plastic_strain5") == Approx(0.0).epsilon(Tolerance));
  }
}

//! Check NorSand class in 3D bonded model
TEST_CASE("NorSand is checked in 3D bonded model", "[material][NorSand][3D]") {
  // Tolerance
  const double Tolerance = 1.E-7;

  const unsigned Dim = 3;

  // Add particle
  mpm::Index pid = 0;
  Eigen::Matrix<double, Dim, 1> coords;
  coords.setZero();
  auto particle = std::make_shared<mpm::Particle<Dim>>(pid, coords);

  // Initialise material
  Json jmaterial;
  jmaterial["density"] = 1800.;
  jmaterial["poisson_ratio"] = 0.3;
  jmaterial["reference_pressure"] = 1000.;
  jmaterial["friction_cs"] = 30;
  jmaterial["N"] = 0.3;
  jmaterial["lambda"] = 0.1;
  jmaterial["kappa"] = 0.03;
  jmaterial["gamma"] = 1.3;
  jmaterial["chi"] = 3.5;
  jmaterial["hardening_modulus"] = 200.0;
  jmaterial["void_ratio_initial"] = 0.85;
  jmaterial["p_image_initial"] = 87014.6;
  jmaterial["bond_model"] = true;
  jmaterial["p_cohesion_initial"] = 10000.0;
  jmaterial["p_dilation_initial"] = 20000.0;
  jmaterial["m_cohesion"] = 20.0;
  jmaterial["m_dilation"] = 5.0;
  jmaterial["m_modulus"] = 10;

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

  //! Check material properties
  SECTION("NorSand check material properties") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "NorSand3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Get material properties
    REQUIRE(material->template property<double>("density") ==
            Approx(jmaterial["density"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("poisson_ratio") ==
            Approx(jmaterial["poisson_ratio"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("reference_pressure") ==
            Approx(jmaterial["reference_pressure"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("friction_cs") ==
            Approx(jmaterial["friction_cs"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("N") ==
            Approx(jmaterial["N"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("lambda") ==
            Approx(jmaterial["lambda"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("kappa") ==
            Approx(jmaterial["kappa"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("gamma") ==
            Approx(jmaterial["gamma"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("chi") ==
            Approx(jmaterial["chi"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("hardening_modulus") ==
            Approx(jmaterial["hardening_modulus"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("void_ratio_initial") ==
            Approx(jmaterial["void_ratio_initial"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("p_image_initial") ==
            Approx(jmaterial["p_image_initial"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("p_cohesion_initial") ==
            Approx(jmaterial["p_cohesion_initial"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("p_dilation_initial") ==
            Approx(jmaterial["p_dilation_initial"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("m_cohesion") ==
            Approx(jmaterial["m_cohesion"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("m_dilation") ==
            Approx(jmaterial["m_dilation"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("m_modulus") ==
            Approx(jmaterial["m_modulus"]).epsilon(Tolerance));
    REQUIRE(material->template property<bool>("bond_model") ==
            jmaterial["bond_model"]);

    // Check if state variable is initialised
    SECTION("State variable is initialised") {
      mpm::dense_map state_variables = material->initialise_state_variables();
      REQUIRE(state_variables.empty() == false);
      REQUIRE(state_variables.at("void_ratio") ==
              Approx(jmaterial["void_ratio_initial"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("e_image") ==
              Approx(0.8533924079).epsilon(Tolerance));
      REQUIRE(state_variables.at("p_image") ==
              Approx(87014.6).epsilon(Tolerance));
      REQUIRE(state_variables.at("p_cohesion") ==
              Approx(jmaterial["p_cohesion_initial"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("p_dilation") ==
              Approx(jmaterial["p_dilation_initial"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("pdstrain") == Approx(0.0).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain0") ==
              Approx(0.0).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain1") ==
              Approx(0.0).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain2") ==
              Approx(0.0).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain3") ==
              Approx(0.0).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain4") ==
              Approx(0.0).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain5") ==
              Approx(0.0).epsilon(Tolerance));
    }
  }

  SECTION("NorSand check elastic stresses") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "NorSand3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    stress(0) = -200000.;
    stress(1) = -200000.;
    stress(2) = -200000.;
    REQUIRE(stress(0) == Approx(-200000.).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(-200000.).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(-200000.).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.).epsilon(Tolerance));

    // Initialise dstrain
    mpm::Material<Dim>::Vector6d dstrain;
    dstrain.setZero();
    dstrain(0) = -0.00010000;
    dstrain(1) = 0.00005000;
    dstrain(2) = 0.00005000;
    dstrain(3) = -0.00000000;
    dstrain(4) = -0.00000000;
    dstrain(5) = -0.00000000;

    // Compute updated stress
    mpm::dense_map state_vars = material->initialise_state_variables();
    stress =
        material->compute_stress(stress, dstrain, particle.get(), &state_vars);

    // Check stresses
    REQUIRE(stress(0) == Approx(-2.011661538461538E+05).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(-1.994169230769231E+05).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(-1.994169230769231E+05).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.000000).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.000000).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.000000).epsilon(Tolerance));

    // Initialise dstrain
    dstrain(0) = -0.00010000;
    dstrain(1) = -0.00005000;
    dstrain(2) = -0.00005000;
    dstrain(3) = -0.00000100;
    dstrain(4) = -0.00000200;
    dstrain(5) = -0.00000300;

    // Reset stress
    stress.setZero();
    stress(0) = -200000.;
    stress(1) = -200000.;
    stress(2) = -200000.;
    REQUIRE(stress(0) == Approx(-200000.).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(-200000.).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(-200000.).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.).epsilon(Tolerance));

    // Compute updated stress
    stress =
        material->compute_stress(stress, dstrain, particle.get(), &state_vars);

    // Check stresses
    REQUIRE(stress(0) == Approx(-2.029153846153846E+05).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(-2.023323076923077E+05).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(-2.023323076923077E+05).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(-0.000058307692308E+05).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(-0.000116615384615E+05).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(-0.000174923076923E+05).epsilon(Tolerance));
  }

  SECTION("NorSand check undrained stresses after yielding with bonded model") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "NorSand3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    stress(0) = -200000;
    stress(1) = -200000;
    stress(2) = -200000;
    stress(3) = 0;
    stress(4) = 0;
    stress(5) = 0;
    REQUIRE(stress(0) == Approx(-200000.).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(-200000.).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(-200000.).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.).epsilon(Tolerance));

    // Initialise dstrain
    mpm::Material<Dim>::Vector6d dstrain;
    dstrain.setZero();
    dstrain(0) = -0.010000;
    dstrain(1) = 0.005000;
    dstrain(2) = 0.005000;
    dstrain(3) = 0.000000;
    dstrain(4) = 0.000000;
    dstrain(5) = 0.000000;

    // Compute updated stress two times to get yielding
    mpm::dense_map state_vars = material->initialise_state_variables();
    stress =
        material->compute_stress(stress, dstrain, particle.get(), &state_vars);
    stress =
        material->compute_stress(stress, dstrain, particle.get(), &state_vars);

    // Check stresses
    REQUIRE(stress(0) == Approx(-3.16690670676614E+05).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(-1.41479312850677E+05).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(-1.41479312850677E+05).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.000000).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.000000).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.000000).epsilon(Tolerance));

    // Check state variables
    REQUIRE(state_vars.empty() == false);
    REQUIRE(state_vars.at("M_theta") == Approx(1.2).epsilon(Tolerance));
    REQUIRE(state_vars.at("void_ratio") ==
            Approx(jmaterial["void_ratio_initial"]).epsilon(Tolerance));
    REQUIRE(state_vars.at("e_image") ==
            Approx(0.8150562893).epsilon(Tolerance));
    REQUIRE(state_vars.at("p_image") ==
            Approx(127668.5058878436).epsilon(Tolerance));
    REQUIRE(state_vars.at("p_cohesion") ==
            Approx(8.1900065839521E+03).epsilon(Tolerance));
    REQUIRE(state_vars.at("p_dilation") ==
            Approx(1.9026156222418E+04).epsilon(Tolerance));
    REQUIRE(state_vars.at("pdstrain") ==
            Approx(0.0099835196).epsilon(Tolerance));
    REQUIRE(state_vars.at("plastic_strain0") ==
            Approx(0.009986604).epsilon(Tolerance));
    REQUIRE(state_vars.at("plastic_strain1") ==
            Approx(-0.0049886753).epsilon(Tolerance));
    REQUIRE(state_vars.at("plastic_strain2") ==
            Approx(-0.0049886753).epsilon(Tolerance));
    REQUIRE(state_vars.at("plastic_strain3") == Approx(0.0).epsilon(Tolerance));
    REQUIRE(state_vars.at("plastic_strain4") == Approx(0.0).epsilon(Tolerance));
    REQUIRE(state_vars.at("plastic_strain5") == Approx(0.0).epsilon(Tolerance));
  }
}
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

#include <cmath>

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
  jmaterial["reference_pressure"] = 1000.;
  jmaterial["friction_cs"] = 30;
  jmaterial["N"] = 0.3;
  jmaterial["e_min"] = 0.542;
  jmaterial["e_max"] = 1.000;
  jmaterial["crushing_pressure"] = 10000000.0;
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
    REQUIRE(material->template property<double>("youngs_modulus") ==
            Approx(jmaterial["youngs_modulus"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("poisson_ratio") ==
            Approx(jmaterial["poisson_ratio"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("reference_pressure") ==
            Approx(jmaterial["reference_pressure"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("friction_cs") ==
            Approx(jmaterial["friction_cs"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("N") ==
            Approx(jmaterial["N"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("e_min") ==
            Approx(jmaterial["e_min"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("e_max") ==
            Approx(jmaterial["e_max"]).epsilon(Tolerance));
    REQUIRE(material->template property<double>("crushing_pressure") ==
            Approx(jmaterial["crushing_pressure"]).epsilon(Tolerance));
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
    REQUIRE(material->template property<bool>("bond_model") ==
            jmaterial["bond_model"]);

    // Check if state variable is initialised
    SECTION("State variable is initialised") {
      mpm::dense_map state_variables = material->initialise_state_variables();
      REQUIRE(state_variables.empty() == false);
      REQUIRE(state_variables.at("p") == Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("q") == Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("void_ratio") ==
              Approx(jmaterial["void_ratio_initial"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("e_image") ==
              Approx(0.8970333385).epsilon(Tolerance));
      REQUIRE(state_variables.at("p_image") ==
              Approx(87014.6).epsilon(Tolerance));
      REQUIRE(state_variables.at("psi_image") ==
              Approx(-0.0033924079).epsilon(Tolerance));
      REQUIRE(state_variables.at("p_cohesion") ==
              Approx(jmaterial["p_cohesion_initial"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("zeta_cohesion") ==
              Approx(1.0).epsilon(Tolerance));
      REQUIRE(state_variables.at("p_dilation") ==
              Approx(jmaterial["p_dilation_initial"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("zeta_dilation") ==
              Approx(1.0).epsilon(Tolerance));
      REQUIRE(state_variables.at("epds") == Approx(0.0).epsilon(Tolerance));
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
    dstrain(0) = -0.0010000;
    dstrain(1) = 0.0005000;
    dstrain(2) = 0.0005000;
    dstrain(3) = -0.0000000;
    dstrain(4) = -0.0000000;
    dstrain(5) = -0.0000000;

    // Compute updated stress
    mpm::dense_map state_vars = material->initialise_state_variables();
    stress =
        material->compute_stress(stress, dstrain, particle.get(), &state_vars);

    // Check stresses
    REQUIRE(stress(0) == Approx(-2.113846153846154E+05).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(-1.943076923076923E+05).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(-1.943076923076923E+05).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.000000).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.000000).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.000000).epsilon(Tolerance));

    // Initialise dstrain
    dstrain(0) = -0.0010000;
    dstrain(1) = -0.0005000;
    dstrain(2) = -0.0005000;
    dstrain(3) = -0.0000100;
    dstrain(4) = -0.0000200;
    dstrain(5) = -0.0000300;

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
    REQUIRE(stress(0) == Approx(-2.284615384615384E+05).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(-2.227692307692308E+05).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(-2.227692307692308E+05).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(-0.000569230769231E+05).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(-0.001138461538462E+05).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(-0.001707692307692E+05).epsilon(Tolerance));
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
    REQUIRE(stress(0) == Approx(-3.566322628959422E+05).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(-0.712267213640271E+05).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(-0.712267213640271E+05).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.000000).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.000000).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.000000).epsilon(Tolerance));

    // Check state variables
    REQUIRE(state_vars.empty() == false);
    REQUIRE(state_vars.at("p") == Approx(166361.9018746655).epsilon(Tolerance));
    REQUIRE(state_vars.at("q") == Approx(285405.5415319151).epsilon(Tolerance));
    REQUIRE(state_vars.at("void_ratio") ==
            Approx(jmaterial["void_ratio_initial"]).epsilon(Tolerance));
    REQUIRE(state_vars.at("e_image") ==
            Approx(0.7573915987).epsilon(Tolerance));
    REQUIRE(state_vars.at("p_image") ==
            Approx(227257.5630587415).epsilon(Tolerance));
    REQUIRE(state_vars.at("psi_image") ==
            Approx(0.092608401300000).epsilon(Tolerance));
    REQUIRE(state_vars.at("p_cohesion") ==
            Approx(9.363731120021845e+03).epsilon(Tolerance));
    REQUIRE(state_vars.at("zeta_cohesion") ==
            Approx(0.936373112002184).epsilon(Tolerance));
    REQUIRE(state_vars.at("p_dilation") ==
            Approx(1.967398016789841e+04).epsilon(Tolerance));
    REQUIRE(state_vars.at("zeta_dilation") ==
            Approx(0.983699008394921).epsilon(Tolerance));
    REQUIRE(state_vars.at("epds") == Approx(0.0032870629).epsilon(Tolerance));
    REQUIRE(state_vars.at("plastic_strain0") ==
            Approx(0.0041962007).epsilon(Tolerance));
    REQUIRE(state_vars.at("plastic_strain1") ==
            Approx(-0.0007343937).epsilon(Tolerance));
    REQUIRE(state_vars.at("plastic_strain2") ==
            Approx(-0.0007343937).epsilon(Tolerance));
    REQUIRE(state_vars.at("plastic_strain3") == Approx(0.0).epsilon(Tolerance));
    REQUIRE(state_vars.at("plastic_strain4") == Approx(0.0).epsilon(Tolerance));
    REQUIRE(state_vars.at("plastic_strain5") == Approx(0.0).epsilon(Tolerance));
  }
}
#include <limits>

#include <cmath>

#include "Eigen/Dense"
#include "catch.hpp"
#include "json.hpp"

#include "cell.h"
#include "material/material.h"
#include "node.h"
#include "particle.h"

//! Check Modified cam clay undrained condition for hardening in 3D
TEST_CASE("Modified cam clay undrained condition hardening is checked in 3D",
          "[material][modified_cam_clay][3D]") {
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
  jmaterial["p_ref"] = 100000;
  jmaterial["e_ref"] = 1.12;
  jmaterial["pc0"] = 300000;
  jmaterial["ocr"] = 1.5;
  jmaterial["m"] = 1.2;
  jmaterial["lambda"] = 0.1;
  jmaterial["kappa"] = 0.03;
  jmaterial["three_invariants"] = false;
  jmaterial["bonding"] = false;
  jmaterial["subloading"] = false;

  //! Check for id = 0
  SECTION("Modified Cam Clay id is zero") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "ModifiedCamClay3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);
  }

  //! Check for id is a positive value
  SECTION("Modified Cam Clay id is positive") {
    unsigned id = std::numeric_limits<unsigned>::max();
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "ModifiedCamClay3D", std::move(id), jmaterial);
    REQUIRE(material->id() == std::numeric_limits<unsigned>::max());
  }

  //! Check material properties
  SECTION("Modified Cam Clay check material properties") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "ModifiedCamClay3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);

    // Get material properties
    REQUIRE(material->template property<double>("density") ==
            Approx(jmaterial.at("density")).epsilon(Tolerance));
    REQUIRE(material->template property<double>("youngs_modulus") ==
            Approx(jmaterial.at("youngs_modulus")).epsilon(Tolerance));
    REQUIRE(material->template property<double>("poisson_ratio") ==
            Approx(jmaterial.at("poisson_ratio")).epsilon(Tolerance));
    REQUIRE(material->template property<double>("p_ref") ==
            Approx(jmaterial.at("p_ref")).epsilon(Tolerance));
    REQUIRE(material->template property<double>("e_ref") ==
            Approx(jmaterial.at("e_ref")).epsilon(Tolerance));
    REQUIRE(material->template property<double>("pc0") ==
            Approx(jmaterial.at("pc0")).epsilon(Tolerance));
    REQUIRE(material->template property<double>("ocr") ==
            Approx(jmaterial.at("ocr")).epsilon(Tolerance));
    REQUIRE(material->template property<double>("m") ==
            Approx(jmaterial.at("m")).epsilon(Tolerance));
    REQUIRE(material->template property<double>("lambda") ==
            Approx(jmaterial.at("lambda")).epsilon(Tolerance));
    REQUIRE(material->template property<double>("kappa") ==
            Approx(jmaterial.at("kappa")).epsilon(Tolerance));

    // Check if state variable is initialised
    SECTION("State variable is initialised") {
      REQUIRE(state_variables.at("bulk_modulus") ==
              Approx(jmaterial["bulk_modulus"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("shear_modulus") ==
              Approx(3 * jmaterial["bulk_modulus"] *
                     (1 - 2 * jmaterial["poisson_ratio"]) /
                     (2 * (1 + jmaterial["poisson_ratio"])))
                  .epsilon(Tolerance));
      REQUIRE(state_variables.at("j3") == Approx(0.0).epsilon(Tolerance));
      REQUIRE(state_variables.at("p") == Approx(0.0).epsilon(Tolerance));
      REQUIRE(state_variables.at("q") == Approx(0.0).epsilon(Tolerance));
      REQUIRE(state_variables.at("theta") == Approx(0.0).epsilon(Tolerance));
      REQUIRE(state_variables.at("pc") ==
              Approx(jmaterial["pc0"]).epsilon(Tolerance));
      REQUIRE(
          state_variables.at("void_ratio") ==
          Approx(jmaterial["e_ref"] -
                 jmaterial["lambda"] * log(jmaterial["pc0"] / jmaterial["ocr"] /
                                           jmaterial["p_ref"]) -
                 jmaterial["kappa"] * log(jmaterial["ocr"]))
              .epsilon(Tolerance));
      REQUIRE(state_variables.at("m_theta") ==
              Approx(jmaterial["m"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("f_function") ==
              Approx(0.0).epsilon(Tolerance));
      REQUIRE(state_variables.at("dpvstrain") ==
              Approx(0.0).epsilon(Tolerance));
      REQUIRE(state_variables.at("dpdstrain") ==
              Approx(0.0).epsilon(Tolerance));
      REQUIRE(state_variables.at("chi") == Approx(1.0).epsilon(Tolerance));
      REQUIRE(state_variables.at("pcd") == Approx(0.0).epsilon(Tolerance));
      REQUIRE(state_variables.at("pcc") == Approx(0.0).epsilon(Tolerance));
      REQUIRE(state_variables.at("subloading_r") ==
              Approx(1.0).epsilon(Tolerance));
    }
  }

  //! Check compute stress in elastic status
  SECTION("CamClay check stresses") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "ModifiedCamClay3D", std::move(id), jmaterial);

    auto cam_clay = std::make_shared<mpm::CamClay<Dim>>(id, jmaterial);

    REQUIRE(material->id() == 0);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    stress(0) = -200000;
    stress(1) = -200000;
    stress(2) = -200000;
    REQUIRE(stress(0) == Approx(-200000.).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(-200000.).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(-200000.).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.).epsilon(Tolerance));

    // Compute stress invariants
    mpm::dense_map state_vars = material->initialise_state_variables();
    mpm::Material<Dim>::Vector6d n;
    cam_clay->compute_stress_invariants(stress, n, &state_vars);

    // Compute elastic modulus
    cam_clay->compute_elastic_tensor(&state_vars);
    REQUIRE(state_vars.at("bulk_modulus") ==
            Approx(7370964.511104).epsilon(Tolerance));
    REQUIRE(state_vars.at("shear_modulus") ==
            Approx(3902275.329408).epsilon(Tolerance));

    // Initialise strain
    mpm::Material<Dim>::Vector6d dstrain;
    dstrain.setZero();
    dstrain(0) = 0.00050000;
    dstrain(1) = 0.00050000;
    dstrain(2) = -0.00100000;
    dstrain(3) = 0.0000000;
    dstrain(4) = 0.0000000;
    dstrain(5) = 0.0000000;

    // Compute stress
    stress =
        cam_clay->compute_stress(stress, dstrain, particle.get(), &state_vars);

    // Check stresses
    REQUIRE(stress(0) == Approx(3902.27532941).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(3902.27532941).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(-7804.55065882).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.000000).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.000000).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.000000).epsilon(Tolerance));
  }
}
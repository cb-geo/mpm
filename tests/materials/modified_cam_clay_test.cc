#include <limits>

#include <cmath>

#include "Eigen/Dense"
#include "catch.hpp"
#include "json.hpp"

#include "cell.h"
#include "material.h"
#include "node.h"
#include "particle.h"

//! Check Modified cam clay undrained condition in 3D
TEST_CASE("Modified cam clay undrained condition is checked in 3D",
          "[material][modified_cam_clay][3D]") {
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

      mpm::dense_map state_variables = material->initialise_state_variables();
      REQUIRE(state_variables.at("bulk_modulus") ==
              Approx(3846153.8460000).epsilon(Tolerance));
      REQUIRE(state_variables.at("shear_modulus") ==
              Approx(4615384.61538462).epsilon(Tolerance));
      REQUIRE(state_variables.at("p") == Approx(0.0).epsilon(Tolerance));
      REQUIRE(state_variables.at("q") == Approx(0.0).epsilon(Tolerance));
      REQUIRE(state_variables.at("theta") == Approx(0.0).epsilon(Tolerance));
      REQUIRE(state_variables.at("pc") ==
              Approx(jmaterial.at("pc0")).epsilon(Tolerance));
      REQUIRE(state_variables.at("void_ratio") ==
              Approx(1.0385213287).epsilon(Tolerance));
      REQUIRE(state_variables.at("m_theta") ==
              Approx(jmaterial.at("m")).epsilon(Tolerance));
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

      const std::vector<std::string> state_vars = {"bulk_modulus",
                                                   "shear_modulus",
                                                   "j3",
                                                   "p",
                                                   "q",
                                                   "theta",
                                                   "pc",
                                                   "void_ratio",
                                                   "delta_phi",
                                                   "m_theta",
                                                   "f_function",
                                                   "dpvstrain",
                                                   "dpdstrain",
                                                   "chi",
                                                   "pcd",
                                                   "pcc",
                                                   "subloading_r"};
      auto state_vars_test = material->state_variables();
      REQUIRE(state_vars == state_vars_test);
    }
  }

  //! Check compute stress in elastic status
  SECTION("CamClay check stresses in elastic status") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "ModifiedCamClay3D", std::move(id), jmaterial);

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
        material->compute_stress(stress, dstrain, particle.get(), &state_vars);

    // Check stresses
    REQUIRE(stress(0) == Approx(-193727.6266809207).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(-193727.6266809207).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(-212544.7466381585).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.000000).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.000000).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.000000).epsilon(Tolerance));
  }

  //! Check compute stress in plastic status
  SECTION("CamClay check stresses in plastic status") {

    jmaterial["pc0"] = 200000;
    jmaterial["ocr"] = 1.;

    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "ModifiedCamClay3D", std::move(id), jmaterial);

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
        material->compute_stress(stress, dstrain, particle.get(), &state_vars);

    // Check stresses
    REQUIRE(stress(0) == Approx(-192882.4825752268).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(-192882.4825752268).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(-211655.0108685302).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.000000).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.000000).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.000000).epsilon(Tolerance));

    // Check pc
    REQUIRE(state_vars.at("pc") ==
            Approx(200368.9146817101).epsilon(Tolerance));
  }

  //! Check compute stress in plastic status with bonded properties
  SECTION("CamClay check stresses in plastic status with bonded properties") {

    jmaterial["bonding"] = true;
    jmaterial["s_h"] = 0.5;
    jmaterial["mc_a"] = 25000;
    jmaterial["mc_b"] = 1;
    jmaterial["mc_c"] = 25000;
    jmaterial["mc_d"] = 1;
    jmaterial["m_degradation"] = 1;
    jmaterial["m_shear"] = 0;

    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "ModifiedCamClay3D", std::move(id), jmaterial);

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
        material->compute_stress(stress, dstrain, particle.get(), &state_vars);

    // Check stresses
    REQUIRE(stress(0) == Approx(-193727.6266809207).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(-193727.6266809207).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(-212544.7466381585).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.000000).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.000000).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.000000).epsilon(Tolerance));

    // Check pc
    REQUIRE(state_vars.at("pc") == Approx(300000).epsilon(Tolerance));
  }

  //! Check compute stress in plastic status with subloading properties
  SECTION("CamClay check stresses in plastic status with bonded properties") {

    jmaterial["subloading"] = true;
    jmaterial["subloading_u"] = 0.5;

    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "ModifiedCamClay3D", std::move(id), jmaterial);

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
        material->compute_stress(stress, dstrain, particle.get(), &state_vars);

    // Check stresses
    REQUIRE(stress(0) == Approx(-162537.902087049).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(-162537.902087049).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(-162537.90208594).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.000000).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.000000).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.000000).epsilon(Tolerance));

    // Check pc
    REQUIRE(state_vars.at("pc") ==
            Approx(325075.8041776047).epsilon(Tolerance));
    // Check subloading_r
    REQUIRE(state_vars.at("subloading_r") == Approx(0.5).epsilon(Tolerance));
  }
}

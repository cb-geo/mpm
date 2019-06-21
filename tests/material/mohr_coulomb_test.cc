#include <limits>

#include <cmath>

#include "Eigen/Dense"
#include "catch.hpp"
#include "json.hpp"

#include "cell.h"
#include "material/material.h"
#include "material/mohr_coulomb.h"
#include "node.h"
#include "particle.h"

//! Check MohrCoulomb class in 2D
//! Cohesion only, without softening
TEST_CASE("MohrCoulomb is checked in 2D (cohesion only, without softening)",
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
  jmaterial["tolerance"] = 0.1;

  //! Check for id = 0
  SECTION("MohrCoulomb id is zero") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "MohrCoulomb2D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);
  }

  //! Check for positive id
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
    REQUIRE(material->property("tolerance") ==
            Approx(jmaterial["tolerance"]).epsilon(Tolerance));

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

  //! Check yield correction based on trial stress
  SECTION("MohrCoulomb check yield correction based on trial stress") {
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

    // Calculate modulus values
    const double K = material->property("youngs_modulus") /
                     (3.0 * (1. - 2. * material->property("poisson_ratio")));
    const double G = material->property("youngs_modulus") /
                     (2.0 * (1. + material->property("poisson_ratio")));
    const double a1 = K + (4.0 / 3.0) * G;
    const double a2 = K - (2.0 / 3.0) * G;
    // Compute elastic tensor
    mpm::Material<Dim>::Matrix6x6 de;
    de.setZero();
    de(0, 0) = a1;
    de(0, 1) = a2;
    de(0, 2) = a2;
    de(1, 0) = a2;
    de(1, 1) = a1;
    de(1, 2) = a2;
    de(2, 0) = a2;
    de(2, 1) = a2;
    de(2, 2) = a1;
    de(3, 3) = G;
    de(4, 4) = G;
    de(5, 5) = G;

    // Initialise state variables
    mpm::dense_map state_variables = material->initialise_state_variables();
    // Check if stress invariants is computed correctly based on stress
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

    // Initialise plastic correction components
    mpm::Material<Dim>::Vector6d df_dsigma, dp_dsigma;
    df_dsigma.setZero();
    dp_dsigma.setZero();
    double softening = 0.;
    // Compute plastic correction components
    mohr_coulomb->compute_df_dp(yield_type, &state_variables, stress,
                                &df_dsigma, &dp_dsigma, &softening);
    // Check plastic correction component based on stress
    // Check dF/dSigma
    REQUIRE(df_dsigma(0) == Approx(0.3618034).epsilon(Tolerance));
    REQUIRE(df_dsigma(1) == Approx(0.1381966).epsilon(Tolerance));
    REQUIRE(df_dsigma(2) == Approx(-0.5).epsilon(Tolerance));
    REQUIRE(df_dsigma(3) == Approx(-0.2236068).epsilon(Tolerance));
    REQUIRE(df_dsigma(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(df_dsigma(5) == Approx(0.).epsilon(Tolerance));
    // Check dP/dSigma
    REQUIRE(dp_dsigma(0) == Approx(0.30618622).epsilon(Tolerance));
    REQUIRE(dp_dsigma(1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dp_dsigma(2) == Approx(-0.30618622).epsilon(Tolerance));
    REQUIRE(dp_dsigma(3) == Approx(-0.30618622).epsilon(Tolerance));
    REQUIRE(dp_dsigma(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dp_dsigma(5) == Approx(0.).epsilon(Tolerance));

    //! Check for shear failure
    SECTION("Check yield correction for shear failure") {
      // Initialise incremental of strain
      mpm::Material<Dim>::Vector6d dstrain;
      dstrain.setZero();
      dstrain(0) = -0.001;
      dstrain(1) = 0.;
      dstrain(2) = 0.;
      dstrain(3) = 0.;
      // Compute trial stress
      mpm::Material<Dim>::Vector6d trial_stress = stress + de * dstrain;
      // Check if stress invariants is computed correctly based on trial stress
      REQUIRE(mohr_coulomb->compute_stress_invariants(
                  trial_stress, &state_variables) == true);
      REQUIRE(state_variables.at("phi") ==
              Approx(jmaterial["friction"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("psi") ==
              Approx(jmaterial["dilation"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("cohesion") ==
              Approx(jmaterial["cohesion"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("j2") ==
              Approx(14031558.18540430).epsilon(Tolerance));
      REQUIRE(state_variables.at("j3") ==
              Approx(-18120349297.8641).epsilon(Tolerance));
      REQUIRE(state_variables.at("epsilon") ==
              Approx(-24826.06157515).epsilon(Tolerance));
      REQUIRE(state_variables.at("rho") ==
              Approx(5297.46320146).epsilon(Tolerance));
      REQUIRE(state_variables.at("theta") ==
              Approx(0.89359516).epsilon(Tolerance));
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
      // Initialise values of yield functions based on trial stress
      Eigen::Matrix<double, 2, 1> yield_function_trial;
      auto yield_type_trial = mohr_coulomb->compute_yield_state(
          &yield_function_trial, &state_variables);
      // Check if yield function and yield state is computed correctly
      REQUIRE(yield_function_trial(0) ==
              Approx(-11623.00067857).epsilon(Tolerance));
      REQUIRE(yield_function_trial(1) ==
              Approx(1492.38393682).epsilon(Tolerance));
      REQUIRE(yield_type_trial == mohr_coulomb->FailureState::Shear);
      // Initialise plastic correction components based on trial stress
      mpm::Material<Dim>::Vector6d df_dsigma_trial, dp_dsigma_trial;
      df_dsigma_trial.setZero();
      dp_dsigma_trial.setZero();
      double softening_trial = 0.;
      // Compute plastic correction components based on trial stress
      mohr_coulomb->compute_df_dp(yield_type_trial, &state_variables,
                                  trial_stress, &df_dsigma_trial,
                                  &dp_dsigma_trial, &softening_trial);

      // Check plastic correction component based on trial stress
      // Check dFtrial/dSigma
      REQUIRE(df_dsigma_trial(0) == Approx(-0.47906443).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(1) == Approx(0.47906443).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(2) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(3) == Approx(-0.14316868).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(5) == Approx(0.).epsilon(Tolerance));
      // Check dPtrial/dSigma
      REQUIRE(dp_dsigma_trial(0) == Approx(-0.47720936).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(1) == Approx(0.29640333).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(2) == Approx(0.18080603).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(3) == Approx(-0.11559730).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(5) == Approx(0.).epsilon(Tolerance));
      // Check compute stress
      mpm::Material<Dim>::Vector6d updated_stress =
          mohr_coulomb->compute_stress(stress, dstrain, particle.get(),
                                       &state_variables);
      // Check update stress
      REQUIRE(updated_stress(0) == Approx(-16581.86769355).epsilon(Tolerance));
      REQUIRE(updated_stress(1) == Approx(-12936.72814065).epsilon(Tolerance));
      REQUIRE(updated_stress(2) == Approx(-13481.40416580).epsilon(Tolerance));
      REQUIRE(updated_stress(3) == Approx(-772.33801257).epsilon(Tolerance));
      REQUIRE(updated_stress(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(updated_stress(5) == Approx(0.).epsilon(Tolerance));
    }

    //! Check for tensile failure
    SECTION("Check yield correction for tensile failure") {
      // Initialise incremental of strain
      mpm::Material<Dim>::Vector6d dstrain;
      dstrain.setZero();
      dstrain(0) = 0.001;
      dstrain(1) = 0.;
      dstrain(2) = 0.;
      dstrain(3) = 0.;
      // Compute trial stress
      mpm::Material<Dim>::Vector6d trial_stress = stress + de * dstrain;
      // Check if stress invariants is computed correctly based on trial stress
      REQUIRE(mohr_coulomb->compute_stress_invariants(
                  trial_stress, &state_variables) == true);
      REQUIRE(state_variables.at("phi") ==
              Approx(jmaterial["friction"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("psi") ==
              Approx(jmaterial["dilation"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("cohesion") ==
              Approx(jmaterial["cohesion"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("j2") ==
              Approx(29416173.57001970).epsilon(Tolerance));
      REQUIRE(state_variables.at("j3") ==
              Approx(59568081053.2882).epsilon(Tolerance));
      REQUIRE(state_variables.at("epsilon") ==
              Approx(4041.45188433).epsilon(Tolerance));
      REQUIRE(state_variables.at("rho") ==
              Approx(7670.22471249).epsilon(Tolerance));
      REQUIRE(state_variables.at("theta") ==
              Approx(0.08181078).epsilon(Tolerance));
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
      // Initialise values of yield functions based on trial stress
      Eigen::Matrix<double, 2, 1> yield_function_trial;
      auto yield_type_trial = mohr_coulomb->compute_yield_state(
          &yield_function_trial, &state_variables);
      // Check if yield function and yield state is computed correctly
      REQUIRE(yield_function_trial(0) ==
              Approx(8575.09909665).epsilon(Tolerance));
      REQUIRE(yield_function_trial(1) ==
              Approx(2902.93416371).epsilon(Tolerance));
      REQUIRE(yield_type_trial == mohr_coulomb->FailureState::Tensile);
      // Initialise plastic correction components based on trial stress
      mpm::Material<Dim>::Vector6d df_dsigma_trial, dp_dsigma_trial;
      df_dsigma_trial.setZero();
      dp_dsigma_trial.setZero();
      double softening_trial = 0.;
      // Compute plastic correction components based on trial stress
      mohr_coulomb->compute_df_dp(yield_type_trial, &state_variables,
                                  trial_stress, &df_dsigma_trial,
                                  &dp_dsigma_trial, &softening_trial);

      // Check plastic correction component based on trial stress
      // Check dFtrial/dSigma
      REQUIRE(df_dsigma_trial(0) == Approx(0.98726817).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(1) == Approx(0.01273183).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(2) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(3) == Approx(-0.11211480).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(5) == Approx(0.).epsilon(Tolerance));
      // Check dPtrial/dSigma
      REQUIRE(dp_dsigma_trial(0) == Approx(0.87816487).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(1) == Approx(0.06958427).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(2) == Approx(0.05225086).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(3) == Approx(-0.09302255).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(5) == Approx(0.).epsilon(Tolerance));
      // Check compute stress
      mpm::Material<Dim>::Vector6d updated_stress =
          mohr_coulomb->compute_stress(stress, dstrain, particle.get(),
                                       &state_variables);
      // Check update stress
      REQUIRE(updated_stress(0) == Approx(-1036.78511143).epsilon(Tolerance));
      REQUIRE(updated_stress(1) == Approx(-4237.89955637).epsilon(Tolerance));
      REQUIRE(updated_stress(2) == Approx(-4895.75071827).epsilon(Tolerance));
      REQUIRE(updated_stress(3) == Approx(-650.24490333).epsilon(Tolerance));
      REQUIRE(updated_stress(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(updated_stress(5) == Approx(0.).epsilon(Tolerance));
    }
  }

  //! Check yield correction based on current stress
  SECTION("MohrCoulomb check yield correction based on current stress") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "MohrCoulomb2D", std::move(id), jmaterial);

    auto mohr_coulomb = std::make_shared<mpm::MohrCoulomb<Dim>>(id, jmaterial);

    REQUIRE(material->id() == 0);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    stress(0) = -2000.;
    stress(1) = -5000.;
    stress(2) = -6000.;

    // Calculate modulus values
    const double K = material->property("youngs_modulus") /
                     (3.0 * (1. - 2. * material->property("poisson_ratio")));
    const double G = material->property("youngs_modulus") /
                     (2.0 * (1. + material->property("poisson_ratio")));
    const double a1 = K + (4.0 / 3.0) * G;
    const double a2 = K - (2.0 / 3.0) * G;
    // Compute elastic tensor
    mpm::Material<Dim>::Matrix6x6 de;
    de.setZero();
    de(0, 0) = a1;
    de(0, 1) = a2;
    de(0, 2) = a2;
    de(1, 0) = a2;
    de(1, 1) = a1;
    de(1, 2) = a2;
    de(2, 0) = a2;
    de(2, 1) = a2;
    de(2, 2) = a1;
    de(3, 3) = G;
    de(4, 4) = G;
    de(5, 5) = G;

    // Initialise state variables
    mpm::dense_map state_variables = material->initialise_state_variables();
    // Check if stress invariants is computed correctly based on stress
    REQUIRE(mohr_coulomb->compute_stress_invariants(stress, &state_variables) ==
            true);
    REQUIRE(state_variables.at("phi") ==
            Approx(jmaterial["friction"]).epsilon(Tolerance));
    REQUIRE(state_variables.at("psi") ==
            Approx(jmaterial["dilation"]).epsilon(Tolerance));
    REQUIRE(state_variables.at("cohesion") ==
            Approx(jmaterial["cohesion"]).epsilon(Tolerance));
    REQUIRE(state_variables.at("j2") ==
            Approx(4333333.33333333).epsilon(Tolerance));
    REQUIRE(state_variables.at("j3") ==
            Approx(2592592592.59259).epsilon(Tolerance));
    REQUIRE(state_variables.at("epsilon") ==
            Approx(-7505.55349947).epsilon(Tolerance));
    REQUIRE(state_variables.at("rho") ==
            Approx(2943.92028878).epsilon(Tolerance));
    REQUIRE(state_variables.at("theta") ==
            Approx(0.24256387).epsilon(Tolerance));
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
    REQUIRE(yield_function(0) == Approx(-2000.).epsilon(Tolerance));
    REQUIRE(yield_function(1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(yield_type == mohr_coulomb->FailureState::Shear);

    // Initialise plastic correction components
    mpm::Material<Dim>::Vector6d df_dsigma, dp_dsigma;
    df_dsigma.setZero();
    dp_dsigma.setZero();
    double softening = 0.;
    // Compute plastic correction components
    mohr_coulomb->compute_df_dp(yield_type, &state_variables, stress,
                                &df_dsigma, &dp_dsigma, &softening);
    // Check plastic correction component based on stress
    // Check dF/dSigma
    REQUIRE(df_dsigma(0) == Approx(0.5).epsilon(Tolerance));
    REQUIRE(df_dsigma(1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(df_dsigma(2) == Approx(-0.5).epsilon(Tolerance));
    REQUIRE(df_dsigma(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(df_dsigma(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(df_dsigma(5) == Approx(0.).epsilon(Tolerance));
    // Check dP/dSigma
    REQUIRE(dp_dsigma(0) == Approx(0.48536267).epsilon(Tolerance));
    REQUIRE(dp_dsigma(1) == Approx(-0.13867505).epsilon(Tolerance));
    REQUIRE(dp_dsigma(2) == Approx(-0.34668762).epsilon(Tolerance));
    REQUIRE(dp_dsigma(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dp_dsigma(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dp_dsigma(5) == Approx(0.).epsilon(Tolerance));

    //! Check for shear failure
    SECTION("Check yield correction for shear failure") {
      // Initialise incremental of strain
      mpm::Material<Dim>::Vector6d dstrain;
      dstrain.setZero();
      dstrain(0) = 0.001;
      dstrain(1) = 0.;
      dstrain(2) = 0.;
      dstrain(3) = 0.;
      // Compute trial stress
      mpm::Material<Dim>::Vector6d trial_stress = stress + de * dstrain;
      // Check if stress invariants is computed correctly based on trial stress
      REQUIRE(mohr_coulomb->compute_stress_invariants(
                  trial_stress, &state_variables) == true);
      REQUIRE(state_variables.at("phi") ==
              Approx(jmaterial["friction"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("psi") ==
              Approx(jmaterial["dilation"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("cohesion") ==
              Approx(jmaterial["cohesion"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("j2") ==
              Approx(42005917.1597633).epsilon(Tolerance));
      REQUIRE(state_variables.at("j3") ==
              Approx(101989076012.745).epsilon(Tolerance));
      REQUIRE(state_variables.at("epsilon") ==
              Approx(6928.20323028).epsilon(Tolerance));
      REQUIRE(state_variables.at("rho") ==
              Approx(9165.79698223).epsilon(Tolerance));
      REQUIRE(state_variables.at("theta") ==
              Approx(0.07722297).epsilon(Tolerance));
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
      // Initialise values of yield functions based on trial stress
      Eigen::Matrix<double, 2, 1> yield_function_trial;
      auto yield_type_trial = mohr_coulomb->compute_yield_state(
          &yield_function_trial, &state_variables);
      // Check if yield function and yield state is computed correctly
      REQUIRE(yield_function_trial(0) ==
              Approx(11461.53846154).epsilon(Tolerance));
      REQUIRE(yield_function_trial(1) ==
              Approx(3846.15384615).epsilon(Tolerance));
      REQUIRE(yield_type_trial == mohr_coulomb->FailureState::Tensile);
      // Initialise plastic correction components based on trial stress
      mpm::Material<Dim>::Vector6d df_dsigma_trial, dp_dsigma_trial;
      df_dsigma_trial.setZero();
      dp_dsigma_trial.setZero();
      double softening_trial = 0.;
      // Compute plastic correction components based on trial stress
      mohr_coulomb->compute_df_dp(yield_type_trial, &state_variables,
                                  trial_stress, &df_dsigma_trial,
                                  &dp_dsigma_trial, &softening_trial);

      // Check plastic correction component based on trial stress
      // Check dFtrial/dSigma
      REQUIRE(df_dsigma_trial(0) == Approx(1.).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(1) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(2) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(3) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(5) == Approx(0.).epsilon(Tolerance));
      // Check dPtrial/dSigma
      REQUIRE(dp_dsigma_trial(0) == Approx(0.88874614).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(1) == Approx(0.05882082).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(2) == Approx(0.05243303).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(3) == Approx(0.).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(5) == Approx(0.).epsilon(Tolerance));
      // Check compute stress
      mpm::Material<Dim>::Vector6d updated_stress =
          mohr_coulomb->compute_stress(stress, dstrain, particle.get(),
                                       &state_variables);
      // Check update stress
      REQUIRE(updated_stress(0) == Approx(0.).epsilon(Tolerance));
      REQUIRE(updated_stress(1) == Approx(-1348.42814081).epsilon(Tolerance));
      REQUIRE(updated_stress(2) == Approx(-488.58203346).epsilon(Tolerance));
      REQUIRE(updated_stress(3) == Approx(0.).epsilon(Tolerance));
      REQUIRE(updated_stress(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(updated_stress(5) == Approx(0.).epsilon(Tolerance));
    }
  }
}

//! Check MohrCoulomb class in 2D
//! Cohesion and friction, without softening
TEST_CASE("MohrCoulomb is checked in 2D (c & phi, without softening)",
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
  jmaterial["friction"] = 30.;
  jmaterial["dilation"] = 0.;
  jmaterial["cohesion"] = 2000.;
  jmaterial["residual_friction"] = 0.;
  jmaterial["residual_dilation"] = 0.;
  jmaterial["residual_cohesion"] = 1000.;
  jmaterial["peak_epds"] = 0.;
  jmaterial["critical_epds"] = 0.;
  jmaterial["tension_cutoff"] = 0.;
  jmaterial["tolerance"] = 0.1;

  //! Check yield correction based on trial stress
  SECTION("MohrCoulomb check yield correction based on trial stress") {
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

    // Calculate modulus values
    const double K = material->property("youngs_modulus") /
                     (3.0 * (1. - 2. * material->property("poisson_ratio")));
    const double G = material->property("youngs_modulus") /
                     (2.0 * (1. + material->property("poisson_ratio")));
    const double a1 = K + (4.0 / 3.0) * G;
    const double a2 = K - (2.0 / 3.0) * G;
    // Compute elastic tensor
    mpm::Material<Dim>::Matrix6x6 de;
    de.setZero();
    de(0, 0) = a1;
    de(0, 1) = a2;
    de(0, 2) = a2;
    de(1, 0) = a2;
    de(1, 1) = a1;
    de(1, 2) = a2;
    de(2, 0) = a2;
    de(2, 1) = a2;
    de(2, 2) = a1;
    de(3, 3) = G;
    de(4, 4) = G;
    de(5, 5) = G;

    // Initialise state variables
    mpm::dense_map state_variables = material->initialise_state_variables();
    // Check if stress invariants is computed correctly based on stress
    REQUIRE(mohr_coulomb->compute_stress_invariants(stress, &state_variables) ==
            true);
    REQUIRE(state_variables.at("phi") == Approx(0.52359878).epsilon(Tolerance));
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
    REQUIRE(yield_function(1) == Approx(-3774.1679421).epsilon(Tolerance));
    REQUIRE(yield_type == mohr_coulomb->FailureState::Elastic);

    // Initialise plastic correction components
    mpm::Material<Dim>::Vector6d df_dsigma, dp_dsigma;
    df_dsigma.setZero();
    dp_dsigma.setZero();
    double softening = 0.;
    // Compute plastic correction components
    mohr_coulomb->compute_df_dp(yield_type, &state_variables, stress,
                                &df_dsigma, &dp_dsigma, &softening);
    // Check plastic correction component based on stress
    // Check dF/dSigma
    REQUIRE(df_dsigma(0) == Approx(0.62666187).epsilon(Tolerance));
    REQUIRE(df_dsigma(1) == Approx(0.23936353).epsilon(Tolerance));
    REQUIRE(df_dsigma(2) == Approx(-0.28867513).epsilon(Tolerance));
    REQUIRE(df_dsigma(3) == Approx(-0.38729833).epsilon(Tolerance));
    REQUIRE(df_dsigma(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(df_dsigma(5) == Approx(0.).epsilon(Tolerance));
    // Check dP/dSigma
    REQUIRE(dp_dsigma(0) == Approx(0.39671222).epsilon(Tolerance));
    REQUIRE(dp_dsigma(1) == Approx(-0.05157114).epsilon(Tolerance));
    REQUIRE(dp_dsigma(2) == Approx(-0.34514108).epsilon(Tolerance));
    REQUIRE(dp_dsigma(3) == Approx(-0.44828335).epsilon(Tolerance));
    REQUIRE(dp_dsigma(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dp_dsigma(5) == Approx(0.).epsilon(Tolerance));

    //! Check for shear failure
    SECTION("Check yield correction for shear failure") {
      // Initialise incremental of strain
      mpm::Material<Dim>::Vector6d dstrain;
      dstrain.setZero();
      dstrain(0) = 0.0001;
      dstrain(1) = 0.;
      dstrain(2) = 0.;
      dstrain(3) = 0.002;
      // Compute trial stress
      mpm::Material<Dim>::Vector6d trial_stress = stress + de * dstrain;
      // Check if stress invariants is computed correctly based on trial stress
      REQUIRE(mohr_coulomb->compute_stress_invariants(
                  trial_stress, &state_variables) == true);
      REQUIRE(state_variables.at("phi") ==
              Approx(0.52359878).epsilon(Tolerance));
      REQUIRE(state_variables.at("psi") ==
              Approx(jmaterial["dilation"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("cohesion") ==
              Approx(jmaterial["cohesion"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("j2") ==
              Approx(46753451.6765286).epsilon(Tolerance));
      REQUIRE(state_variables.at("j3") ==
              Approx(56758188775.9403).epsilon(Tolerance));
      REQUIRE(state_variables.at("epsilon") ==
              Approx(-8948.92917244).epsilon(Tolerance));
      REQUIRE(state_variables.at("rho") ==
              Approx(9669.89676021).epsilon(Tolerance));
      REQUIRE(state_variables.at("theta") ==
              Approx(0.36378823).epsilon(Tolerance));
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
      // Initialise values of yield functions based on trial stress
      Eigen::Matrix<double, 2, 1> yield_function_trial;
      auto yield_type_trial = mohr_coulomb->compute_yield_state(
          &yield_function_trial, &state_variables);
      // Check if yield function and yield state is computed correctly
      REQUIRE(yield_function_trial(0) ==
              Approx(2212.05893238).epsilon(Tolerance));
      REQUIRE(yield_function_trial(1) ==
              Approx(3174.54763108).epsilon(Tolerance));
      REQUIRE(yield_type_trial == mohr_coulomb->FailureState::Shear);
      // Initialise plastic correction components based on trial stress
      mpm::Material<Dim>::Vector6d df_dsigma_trial, dp_dsigma_trial;
      df_dsigma_trial.setZero();
      dp_dsigma_trial.setZero();
      double softening_trial = 0.;
      // Compute plastic correction components based on trial stress
      mohr_coulomb->compute_df_dp(yield_type_trial, &state_variables,
                                  trial_stress, &df_dsigma_trial,
                                  &dp_dsigma_trial, &softening_trial);

      // Check plastic correction component based on trial stress
      // Check dFtrial/dSigma
      REQUIRE(df_dsigma_trial(0) == Approx(0.36433344).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(1) == Approx(0.21301683).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(2) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(3) == Approx(0.57237152).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(5) == Approx(0.).epsilon(Tolerance));
      // Check dPtrial/dSigma
      REQUIRE(dp_dsigma_trial(0) == Approx(0.21883386).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(1) == Approx(0.08262485).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(2) == Approx(-0.30145871).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(3) == Approx(0.51522540).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(5) == Approx(0.).epsilon(Tolerance));
      // Check compute stress
      mpm::Material<Dim>::Vector6d updated_stress =
          mohr_coulomb->compute_stress(stress, dstrain, particle.get(),
                                       &state_variables);
      // Check update stress
      REQUIRE(updated_stress(0) == Approx(-6491.90864415).epsilon(Tolerance));
      REQUIRE(updated_stress(1) == Approx(-6494.64090848).epsilon(Tolerance));
      REQUIRE(updated_stress(2) == Approx(-2513.45044738).epsilon(Tolerance));
      REQUIRE(updated_stress(3) == Approx(3351.32138956).epsilon(Tolerance));
      REQUIRE(updated_stress(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(updated_stress(5) == Approx(0.).epsilon(Tolerance));
    }

    //! Check for tensile failure
    SECTION("Check yield correction for tensile failure") {
      // Initialise incremental of strain
      mpm::Material<Dim>::Vector6d dstrain;
      dstrain.setZero();
      dstrain(0) = 0.001;
      dstrain(1) = 0.;
      dstrain(2) = 0.;
      dstrain(3) = 0.;
      // Compute trial stress
      mpm::Material<Dim>::Vector6d trial_stress = stress + de * dstrain;
      // Check if stress invariants is computed correctly based on trial stress
      REQUIRE(mohr_coulomb->compute_stress_invariants(
                  trial_stress, &state_variables) == true);
      REQUIRE(state_variables.at("phi") ==
              Approx(0.52359878).epsilon(Tolerance));
      REQUIRE(state_variables.at("psi") ==
              Approx(jmaterial["dilation"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("cohesion") ==
              Approx(jmaterial["cohesion"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("j2") ==
              Approx(29416173.57001970).epsilon(Tolerance));
      REQUIRE(state_variables.at("j3") ==
              Approx(59568081053.2882).epsilon(Tolerance));
      REQUIRE(state_variables.at("epsilon") ==
              Approx(4041.45188433).epsilon(Tolerance));
      REQUIRE(state_variables.at("rho") ==
              Approx(7670.22471249).epsilon(Tolerance));
      REQUIRE(state_variables.at("theta") ==
              Approx(0.08181078).epsilon(Tolerance));
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
      // Initialise values of yield functions based on trial stress
      Eigen::Matrix<double, 2, 1> yield_function_trial;
      auto yield_type_trial = mohr_coulomb->compute_yield_state(
          &yield_function_trial, &state_variables);
      // Check if yield function and yield state is computed correctly
      REQUIRE(yield_function_trial(0) ==
              Approx(8575.09909665).epsilon(Tolerance));
      REQUIRE(yield_function_trial(1) ==
              Approx(5781.54613102).epsilon(Tolerance));
      REQUIRE(yield_type_trial == mohr_coulomb->FailureState::Tensile);
      // Initialise plastic correction components based on trial stress
      mpm::Material<Dim>::Vector6d df_dsigma_trial, dp_dsigma_trial;
      df_dsigma_trial.setZero();
      dp_dsigma_trial.setZero();
      double softening_trial = 0.;
      // Compute plastic correction components based on trial stress
      mohr_coulomb->compute_df_dp(yield_type_trial, &state_variables,
                                  trial_stress, &df_dsigma_trial,
                                  &dp_dsigma_trial, &softening_trial);

      // Check plastic correction component based on trial stress
      // Check dFtrial/dSigma
      REQUIRE(df_dsigma_trial(0) == Approx(0.98726817).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(1) == Approx(0.01273183).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(2) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(3) == Approx(-0.11211480).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(5) == Approx(0.).epsilon(Tolerance));
      // Check dPtrial/dSigma
      REQUIRE(dp_dsigma_trial(0) == Approx(0.87816487).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(1) == Approx(0.06958427).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(2) == Approx(0.05225086).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(3) == Approx(-0.09302255).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(5) == Approx(0.).epsilon(Tolerance));
      // Check compute stress
      mpm::Material<Dim>::Vector6d updated_stress =
          mohr_coulomb->compute_stress(stress, dstrain, particle.get(),
                                       &state_variables);
      // Check update stress
      REQUIRE(updated_stress(0) == Approx(-140.3999116).epsilon(Tolerance));
      REQUIRE(updated_stress(1) == Approx(-4560.80574858).epsilon(Tolerance));
      REQUIRE(updated_stress(2) == Approx(-5469.2297259).epsilon(Tolerance));
      REQUIRE(updated_stress(3) == Approx(-754.27113222).epsilon(Tolerance));
      REQUIRE(updated_stress(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(updated_stress(5) == Approx(0.).epsilon(Tolerance));
    }
  }

  //! Check yield correction based on current stress
  SECTION("MohrCoulomb check yield correction based on current stress") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "MohrCoulomb2D", std::move(id), jmaterial);

    auto mohr_coulomb = std::make_shared<mpm::MohrCoulomb<Dim>>(id, jmaterial);

    REQUIRE(material->id() == 0);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    stress(0) = -1000.;
    stress(1) = -7000.;
    stress(2) = -9928.20323028;

    // Calculate modulus values
    const double K = material->property("youngs_modulus") /
                     (3.0 * (1. - 2. * material->property("poisson_ratio")));
    const double G = material->property("youngs_modulus") /
                     (2.0 * (1. + material->property("poisson_ratio")));
    const double a1 = K + (4.0 / 3.0) * G;
    const double a2 = K - (2.0 / 3.0) * G;
    // Compute elastic tensor
    mpm::Material<Dim>::Matrix6x6 de;
    de.setZero();
    de(0, 0) = a1;
    de(0, 1) = a2;
    de(0, 2) = a2;
    de(1, 0) = a2;
    de(1, 1) = a1;
    de(1, 2) = a2;
    de(2, 0) = a2;
    de(2, 1) = a2;
    de(2, 2) = a1;
    de(3, 3) = G;
    de(4, 4) = G;
    de(5, 5) = G;

    // Initialise state variables
    mpm::dense_map state_variables = material->initialise_state_variables();
    // Check if stress invariants is computed correctly based on stress
    REQUIRE(mohr_coulomb->compute_stress_invariants(stress, &state_variables) ==
            true);
    REQUIRE(state_variables.at("phi") == Approx(0.52359878).epsilon(Tolerance));
    REQUIRE(state_variables.at("psi") ==
            Approx(jmaterial["dilation"]).epsilon(Tolerance));
    REQUIRE(state_variables.at("cohesion") ==
            Approx(jmaterial["cohesion"]).epsilon(Tolerance));
    REQUIRE(state_variables.at("j2") ==
            Approx(20714531.1798341).epsilon(Tolerance));
    REQUIRE(state_variables.at("j3") ==
            Approx(20136747919.1226).epsilon(Tolerance));
    REQUIRE(state_variables.at("epsilon") ==
            Approx(-10350.85296109).epsilon(Tolerance));
    REQUIRE(state_variables.at("rho") ==
            Approx(6436.54117983).epsilon(Tolerance));
    REQUIRE(state_variables.at("theta") ==
            Approx(0.32751078).epsilon(Tolerance));
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
    REQUIRE(yield_function(0) == Approx(-1000.).epsilon(Tolerance));
    REQUIRE(yield_function(1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(yield_type == mohr_coulomb->FailureState::Shear);

    // Initialise plastic correction components
    mpm::Material<Dim>::Vector6d df_dsigma, dp_dsigma;
    df_dsigma.setZero();
    dp_dsigma.setZero();
    double softening = 0.;
    // Compute plastic correction components
    mohr_coulomb->compute_df_dp(yield_type, &state_variables, stress,
                                &df_dsigma, &dp_dsigma, &softening);
    // Check plastic correction component based on stress
    // Check dF/dSigma
    REQUIRE(df_dsigma(0) == Approx(0.86602540).epsilon(Tolerance));
    REQUIRE(df_dsigma(1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(df_dsigma(2) == Approx(-0.28867513).epsilon(Tolerance));
    REQUIRE(df_dsigma(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(df_dsigma(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(df_dsigma(5) == Approx(0.).epsilon(Tolerance));
    // Check dP/dSigma
    REQUIRE(dp_dsigma(0) == Approx(0.67187026).epsilon(Tolerance));
    REQUIRE(dp_dsigma(1) == Approx(-0.30786798).epsilon(Tolerance));
    REQUIRE(dp_dsigma(2) == Approx(-0.36400229).epsilon(Tolerance));
    REQUIRE(dp_dsigma(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dp_dsigma(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dp_dsigma(5) == Approx(0.).epsilon(Tolerance));

    //! Check for shear failure
    SECTION("Check yield correction for shear failure") {
      // Initialise incremental of strain
      mpm::Material<Dim>::Vector6d dstrain;
      dstrain.setZero();
      dstrain(0) = 0.0001;
      dstrain(1) = 0.;
      dstrain(2) = 0.;
      dstrain(3) = 0.001;
      // Compute trial stress
      mpm::Material<Dim>::Vector6d trial_stress = stress + de * dstrain;
      // Check if stress invariants is computed correctly based on trial stress
      REQUIRE(mohr_coulomb->compute_stress_invariants(
                  trial_stress, &state_variables) == true);
      REQUIRE(state_variables.at("phi") ==
              Approx(0.52359878).epsilon(Tolerance));
      REQUIRE(state_variables.at("psi") ==
              Approx(jmaterial["dilation"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("cohesion") ==
              Approx(jmaterial["cohesion"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("j2") ==
              Approx(39532413.6649157).epsilon(Tolerance));
      REQUIRE(state_variables.at("j3") ==
              Approx(91832809718.2421).epsilon(Tolerance));
      REQUIRE(state_variables.at("epsilon") ==
              Approx(-8907.47728811).epsilon(Tolerance));
      REQUIRE(state_variables.at("rho") ==
              Approx(8891.8404917).epsilon(Tolerance));
      REQUIRE(state_variables.at("theta") ==
              Approx(0.09473338).epsilon(Tolerance));
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
      // Initialise values of yield functions based on trial stress
      Eigen::Matrix<double, 2, 1> yield_function_trial;
      auto yield_type_trial = mohr_coulomb->compute_yield_state(
          &yield_function_trial, &state_variables);
      // Check if yield function and yield state is computed correctly
      REQUIRE(yield_function_trial(0) ==
              Approx(2084.86947857).epsilon(Tolerance));
      REQUIRE(yield_function_trial(1) ==
              Approx(2505.03198892).epsilon(Tolerance));
      REQUIRE(yield_type_trial == mohr_coulomb->FailureState::Shear);
      // Initialise plastic correction components based on trial stress
      mpm::Material<Dim>::Vector6d df_dsigma_trial, dp_dsigma_trial;
      df_dsigma_trial.setZero();
      dp_dsigma_trial.setZero();
      double softening_trial = 0.;
      // Compute plastic correction components based on trial stress
      mohr_coulomb->compute_df_dp(yield_type_trial, &state_variables,
                                  trial_stress, &df_dsigma_trial,
                                  &dp_dsigma_trial, &softening_trial);

      // Check plastic correction component based on trial stress
      // Check dFtrial/dSigma
      REQUIRE(df_dsigma_trial(0) == Approx(0.71907297).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(1) == Approx(0.14695243).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(2) == Approx(-0.28867513).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(3) == Approx(0.32506849).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(5) == Approx(0.).epsilon(Tolerance));
      // Check dPtrial/dSigma
      REQUIRE(dp_dsigma_trial(0) == Approx(0.50315689).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(1) == Approx(-0.1607794).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(2) == Approx(-0.34237749).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(3) == Approx(0.37723653).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(5) == Approx(0.).epsilon(Tolerance));
      // Check compute stress
      mpm::Material<Dim>::Vector6d updated_stress =
          mohr_coulomb->compute_stress(stress, dstrain, particle.get(),
                                       &state_variables);
      // Check update stress
      REQUIRE(updated_stress(0) == Approx(-2358.37307729).epsilon(Tolerance));
      REQUIRE(updated_stress(1) == Approx(-5535.40840205).epsilon(Tolerance));
      REQUIRE(updated_stress(2) == Approx(-7534.42175094).epsilon(Tolerance));
      REQUIRE(updated_stress(3) == Approx(3066.34734064).epsilon(Tolerance));
      REQUIRE(updated_stress(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(updated_stress(5) == Approx(0.).epsilon(Tolerance));
    }
  }
}

//! Check MohrCoulomb class in 2D
//! Cohesion, friction and dilation, without softening
TEST_CASE("MohrCoulomb is checked in 2D (c & phi & psi, without softening)",
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
  jmaterial["friction"] = 30.;
  jmaterial["dilation"] = 15.;
  jmaterial["cohesion"] = 2000.;
  jmaterial["residual_friction"] = 0.;
  jmaterial["residual_dilation"] = 0.;
  jmaterial["residual_cohesion"] = 1000.;
  jmaterial["peak_epds"] = 0.;
  jmaterial["critical_epds"] = 0.;
  jmaterial["tension_cutoff"] = 0.;
  jmaterial["tolerance"] = 0.1;

  //! Check yield correction based on trial stress
  SECTION("MohrCoulomb check yield correction based on trial stress") {
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

    // Calculate modulus values
    const double K = material->property("youngs_modulus") /
                     (3.0 * (1. - 2. * material->property("poisson_ratio")));
    const double G = material->property("youngs_modulus") /
                     (2.0 * (1. + material->property("poisson_ratio")));
    const double a1 = K + (4.0 / 3.0) * G;
    const double a2 = K - (2.0 / 3.0) * G;
    // Compute elastic tensor
    mpm::Material<Dim>::Matrix6x6 de;
    de.setZero();
    de(0, 0) = a1;
    de(0, 1) = a2;
    de(0, 2) = a2;
    de(1, 0) = a2;
    de(1, 1) = a1;
    de(1, 2) = a2;
    de(2, 0) = a2;
    de(2, 1) = a2;
    de(2, 2) = a1;
    de(3, 3) = G;
    de(4, 4) = G;
    de(5, 5) = G;

    // Initialise state variables
    mpm::dense_map state_variables = material->initialise_state_variables();
    // Check if stress invariants is computed correctly based on stress
    REQUIRE(mohr_coulomb->compute_stress_invariants(stress, &state_variables) ==
            true);
    REQUIRE(state_variables.at("phi") == Approx(0.52359878).epsilon(Tolerance));
    REQUIRE(state_variables.at("psi") == Approx(0.26179939).epsilon(Tolerance));
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
    REQUIRE(yield_function(1) == Approx(-3774.1679421).epsilon(Tolerance));
    REQUIRE(yield_type == mohr_coulomb->FailureState::Elastic);

    // Initialise plastic correction components
    mpm::Material<Dim>::Vector6d df_dsigma, dp_dsigma;
    df_dsigma.setZero();
    dp_dsigma.setZero();
    double softening = 0.;
    // Compute plastic correction components
    mohr_coulomb->compute_df_dp(yield_type, &state_variables, stress,
                                &df_dsigma, &dp_dsigma, &softening);
    // Check plastic correction component based on stress
    // Check dF/dSigma
    REQUIRE(df_dsigma(0) == Approx(0.62666187).epsilon(Tolerance));
    REQUIRE(df_dsigma(1) == Approx(0.23936353).epsilon(Tolerance));
    REQUIRE(df_dsigma(2) == Approx(-0.28867513).epsilon(Tolerance));
    REQUIRE(df_dsigma(3) == Approx(-0.38729833).epsilon(Tolerance));
    REQUIRE(df_dsigma(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(df_dsigma(5) == Approx(0.).epsilon(Tolerance));
    // Check dP/dSigma
    REQUIRE(dp_dsigma(0) == Approx(0.48581658).epsilon(Tolerance));
    REQUIRE(dp_dsigma(1) == Approx(0.03777283).epsilon(Tolerance));
    REQUIRE(dp_dsigma(2) == Approx(-0.25564021).epsilon(Tolerance));
    REQUIRE(dp_dsigma(3) == Approx(-0.44804375).epsilon(Tolerance));
    REQUIRE(dp_dsigma(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dp_dsigma(5) == Approx(0.).epsilon(Tolerance));

    //! Check for shear failure
    SECTION("Check yield correction for shear failure") {
      // Initialise incremental of strain
      mpm::Material<Dim>::Vector6d dstrain;
      dstrain.setZero();
      dstrain(0) = 0.0001;
      dstrain(1) = 0.;
      dstrain(2) = 0.;
      dstrain(3) = 0.002;
      // Compute trial stress
      mpm::Material<Dim>::Vector6d trial_stress = stress + de * dstrain;
      // Check if stress invariants is computed correctly based on trial stress
      REQUIRE(mohr_coulomb->compute_stress_invariants(
                  trial_stress, &state_variables) == true);
      REQUIRE(state_variables.at("phi") ==
              Approx(0.52359878).epsilon(Tolerance));
      REQUIRE(state_variables.at("psi") ==
              Approx(0.26179939).epsilon(Tolerance));
      REQUIRE(state_variables.at("cohesion") ==
              Approx(jmaterial["cohesion"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("j2") ==
              Approx(46753451.6765286).epsilon(Tolerance));
      REQUIRE(state_variables.at("j3") ==
              Approx(56758188775.9403).epsilon(Tolerance));
      REQUIRE(state_variables.at("epsilon") ==
              Approx(-8948.92917244).epsilon(Tolerance));
      REQUIRE(state_variables.at("rho") ==
              Approx(9669.89676021).epsilon(Tolerance));
      REQUIRE(state_variables.at("theta") ==
              Approx(0.36378823).epsilon(Tolerance));
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
      // Initialise values of yield functions based on trial stress
      Eigen::Matrix<double, 2, 1> yield_function_trial;
      auto yield_type_trial = mohr_coulomb->compute_yield_state(
          &yield_function_trial, &state_variables);
      // Check if yield function and yield state is computed correctly
      REQUIRE(yield_function_trial(0) ==
              Approx(2212.05893238).epsilon(Tolerance));
      REQUIRE(yield_function_trial(1) ==
              Approx(3174.54763108).epsilon(Tolerance));
      REQUIRE(yield_type_trial == mohr_coulomb->FailureState::Shear);
      // Initialise plastic correction components based on trial stress
      mpm::Material<Dim>::Vector6d df_dsigma_trial, dp_dsigma_trial;
      df_dsigma_trial.setZero();
      dp_dsigma_trial.setZero();
      double softening_trial = 0.;
      // Compute plastic correction components based on trial stress
      mohr_coulomb->compute_df_dp(yield_type_trial, &state_variables,
                                  trial_stress, &df_dsigma_trial,
                                  &dp_dsigma_trial, &softening_trial);

      // Check plastic correction component based on trial stress
      // Check dFtrial/dSigma
      REQUIRE(df_dsigma_trial(0) == Approx(0.36433344).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(1) == Approx(0.21301683).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(2) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(3) == Approx(0.57237152).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(5) == Approx(0.).epsilon(Tolerance));
      // Check dPtrial/dSigma
      REQUIRE(dp_dsigma_trial(0) == Approx(0.30814480).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(1) == Approx(0.17193919).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(2) == Approx(-0.21213479).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(3) == Approx(0.51521254).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(5) == Approx(0.).epsilon(Tolerance));
      // Check compute stress
      mpm::Material<Dim>::Vector6d updated_stress =
          mohr_coulomb->compute_stress(stress, dstrain, particle.get(),
                                       &state_variables);
      // Check update stress
      REQUIRE(updated_stress(0) == Approx(-7573.14655485).epsilon(Tolerance));
      REQUIRE(updated_stress(1) == Approx(-8293.81378005).epsilon(Tolerance));
      REQUIRE(updated_stress(2) == Approx(-6337.06362775).epsilon(Tolerance));
      REQUIRE(updated_stress(3) == Approx(4709.15490250).epsilon(Tolerance));
      REQUIRE(updated_stress(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(updated_stress(5) == Approx(0.).epsilon(Tolerance));
    }

    //! Check for tensile failure
    SECTION("Check yield correction for tensile failure") {
      // Initialise incremental of strain
      mpm::Material<Dim>::Vector6d dstrain;
      dstrain.setZero();
      dstrain(0) = 0.001;
      dstrain(1) = 0.;
      dstrain(2) = 0.;
      dstrain(3) = 0.;
      // Compute trial stress
      mpm::Material<Dim>::Vector6d trial_stress = stress + de * dstrain;
      // Check if stress invariants is computed correctly based on trial stress
      REQUIRE(mohr_coulomb->compute_stress_invariants(
                  trial_stress, &state_variables) == true);
      REQUIRE(state_variables.at("phi") ==
              Approx(0.52359878).epsilon(Tolerance));
      REQUIRE(state_variables.at("psi") ==
              Approx(0.26179939).epsilon(Tolerance));
      REQUIRE(state_variables.at("cohesion") ==
              Approx(jmaterial["cohesion"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("j2") ==
              Approx(29416173.57001970).epsilon(Tolerance));
      REQUIRE(state_variables.at("j3") ==
              Approx(59568081053.2882).epsilon(Tolerance));
      REQUIRE(state_variables.at("epsilon") ==
              Approx(4041.45188433).epsilon(Tolerance));
      REQUIRE(state_variables.at("rho") ==
              Approx(7670.22471249).epsilon(Tolerance));
      REQUIRE(state_variables.at("theta") ==
              Approx(0.08181078).epsilon(Tolerance));
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
      // Initialise values of yield functions based on trial stress
      Eigen::Matrix<double, 2, 1> yield_function_trial;
      auto yield_type_trial = mohr_coulomb->compute_yield_state(
          &yield_function_trial, &state_variables);
      // Check if yield function and yield state is computed correctly
      REQUIRE(yield_function_trial(0) ==
              Approx(8575.09909665).epsilon(Tolerance));
      REQUIRE(yield_function_trial(1) ==
              Approx(5781.54613102).epsilon(Tolerance));
      REQUIRE(yield_type_trial == mohr_coulomb->FailureState::Tensile);
      // Initialise plastic correction components based on trial stress
      mpm::Material<Dim>::Vector6d df_dsigma_trial, dp_dsigma_trial;
      df_dsigma_trial.setZero();
      dp_dsigma_trial.setZero();
      double softening_trial = 0.;
      // Compute plastic correction components based on trial stress
      mohr_coulomb->compute_df_dp(yield_type_trial, &state_variables,
                                  trial_stress, &df_dsigma_trial,
                                  &dp_dsigma_trial, &softening_trial);

      // Check plastic correction component based on trial stress
      // Check dFtrial/dSigma
      REQUIRE(df_dsigma_trial(0) == Approx(0.98726817).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(1) == Approx(0.01273183).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(2) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(3) == Approx(-0.11211480).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(5) == Approx(0.).epsilon(Tolerance));
      // Check dPtrial/dSigma
      REQUIRE(dp_dsigma_trial(0) == Approx(0.87816487).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(1) == Approx(0.06958427).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(2) == Approx(0.05225086).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(3) == Approx(-0.09302255).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(5) == Approx(0.).epsilon(Tolerance));
      // Check compute stress
      mpm::Material<Dim>::Vector6d updated_stress =
          mohr_coulomb->compute_stress(stress, dstrain, particle.get(),
                                       &state_variables);
      // Check update stress
      REQUIRE(updated_stress(0) == Approx(-140.3999116).epsilon(Tolerance));
      REQUIRE(updated_stress(1) == Approx(-4560.80574858).epsilon(Tolerance));
      REQUIRE(updated_stress(2) == Approx(-5469.2297259).epsilon(Tolerance));
      REQUIRE(updated_stress(3) == Approx(-754.27113222).epsilon(Tolerance));
      REQUIRE(updated_stress(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(updated_stress(5) == Approx(0.).epsilon(Tolerance));
    }
  }

  //! Check yield correction based on current stress
  SECTION("MohrCoulomb check yield correction based on current stress") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "MohrCoulomb2D", std::move(id), jmaterial);

    auto mohr_coulomb = std::make_shared<mpm::MohrCoulomb<Dim>>(id, jmaterial);

    REQUIRE(material->id() == 0);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    stress(0) = -1000.;
    stress(1) = -7000.;
    stress(2) = -9928.20323028;

    // Calculate modulus values
    const double K = material->property("youngs_modulus") /
                     (3.0 * (1. - 2. * material->property("poisson_ratio")));
    const double G = material->property("youngs_modulus") /
                     (2.0 * (1. + material->property("poisson_ratio")));
    const double a1 = K + (4.0 / 3.0) * G;
    const double a2 = K - (2.0 / 3.0) * G;
    // Compute elastic tensor
    mpm::Material<Dim>::Matrix6x6 de;
    de.setZero();
    de(0, 0) = a1;
    de(0, 1) = a2;
    de(0, 2) = a2;
    de(1, 0) = a2;
    de(1, 1) = a1;
    de(1, 2) = a2;
    de(2, 0) = a2;
    de(2, 1) = a2;
    de(2, 2) = a1;
    de(3, 3) = G;
    de(4, 4) = G;
    de(5, 5) = G;

    // Initialise state variables
    mpm::dense_map state_variables = material->initialise_state_variables();
    // Check if stress invariants is computed correctly based on stress
    REQUIRE(mohr_coulomb->compute_stress_invariants(stress, &state_variables) ==
            true);
    REQUIRE(state_variables.at("phi") == Approx(0.52359878).epsilon(Tolerance));
    REQUIRE(state_variables.at("psi") == Approx(0.26179939).epsilon(Tolerance));
    REQUIRE(state_variables.at("cohesion") ==
            Approx(jmaterial["cohesion"]).epsilon(Tolerance));
    REQUIRE(state_variables.at("j2") ==
            Approx(20714531.1798341).epsilon(Tolerance));
    REQUIRE(state_variables.at("j3") ==
            Approx(20136747919.1226).epsilon(Tolerance));
    REQUIRE(state_variables.at("epsilon") ==
            Approx(-10350.85296109).epsilon(Tolerance));
    REQUIRE(state_variables.at("rho") ==
            Approx(6436.54117983).epsilon(Tolerance));
    REQUIRE(state_variables.at("theta") ==
            Approx(0.32751078).epsilon(Tolerance));
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
    REQUIRE(yield_function(0) == Approx(-1000.).epsilon(Tolerance));
    REQUIRE(yield_function(1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(yield_type == mohr_coulomb->FailureState::Shear);

    // Initialise plastic correction components
    mpm::Material<Dim>::Vector6d df_dsigma, dp_dsigma;
    df_dsigma.setZero();
    dp_dsigma.setZero();
    double softening = 0.;
    // Compute plastic correction components
    mohr_coulomb->compute_df_dp(yield_type, &state_variables, stress,
                                &df_dsigma, &dp_dsigma, &softening);
    // Check plastic correction component based on stress
    // Check dF/dSigma
    REQUIRE(df_dsigma(0) == Approx(0.86602540).epsilon(Tolerance));
    REQUIRE(df_dsigma(1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(df_dsigma(2) == Approx(-0.28867513).epsilon(Tolerance));
    REQUIRE(df_dsigma(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(df_dsigma(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(df_dsigma(5) == Approx(0.).epsilon(Tolerance));
    // Check dP/dSigma
    REQUIRE(dp_dsigma(0) == Approx(0.76114953).epsilon(Tolerance));
    REQUIRE(dp_dsigma(1) == Approx(-0.21853457).epsilon(Tolerance));
    REQUIRE(dp_dsigma(2) == Approx(-0.27466577).epsilon(Tolerance));
    REQUIRE(dp_dsigma(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dp_dsigma(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dp_dsigma(5) == Approx(0.).epsilon(Tolerance));

    //! Check for shear failure
    SECTION("Check yield correction for shear failure") {
      // Initialise incremental of strain
      mpm::Material<Dim>::Vector6d dstrain;
      dstrain.setZero();
      dstrain(0) = 0.0001;
      dstrain(1) = 0.;
      dstrain(2) = 0.;
      dstrain(3) = 0.001;
      // Compute trial stress
      mpm::Material<Dim>::Vector6d trial_stress = stress + de * dstrain;
      // Check if stress invariants is computed correctly based on trial stress
      REQUIRE(mohr_coulomb->compute_stress_invariants(
                  trial_stress, &state_variables) == true);
      REQUIRE(state_variables.at("phi") ==
              Approx(0.52359878).epsilon(Tolerance));
      REQUIRE(state_variables.at("psi") ==
              Approx(0.26179939).epsilon(Tolerance));
      REQUIRE(state_variables.at("cohesion") ==
              Approx(jmaterial["cohesion"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("j2") ==
              Approx(39532413.6649157).epsilon(Tolerance));
      REQUIRE(state_variables.at("j3") ==
              Approx(91832809718.2421).epsilon(Tolerance));
      REQUIRE(state_variables.at("epsilon") ==
              Approx(-8907.47728811).epsilon(Tolerance));
      REQUIRE(state_variables.at("rho") ==
              Approx(8891.8404917).epsilon(Tolerance));
      REQUIRE(state_variables.at("theta") ==
              Approx(0.09473338).epsilon(Tolerance));
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
      // Initialise values of yield functions based on trial stress
      Eigen::Matrix<double, 2, 1> yield_function_trial;
      auto yield_type_trial = mohr_coulomb->compute_yield_state(
          &yield_function_trial, &state_variables);
      // Check if yield function and yield state is computed correctly
      REQUIRE(yield_function_trial(0) ==
              Approx(2084.86947857).epsilon(Tolerance));
      REQUIRE(yield_function_trial(1) ==
              Approx(2505.03198892).epsilon(Tolerance));
      REQUIRE(yield_type_trial == mohr_coulomb->FailureState::Shear);
      // Initialise plastic correction components based on trial stress
      mpm::Material<Dim>::Vector6d df_dsigma_trial, dp_dsigma_trial;
      df_dsigma_trial.setZero();
      dp_dsigma_trial.setZero();
      double softening_trial = 0.;
      // Compute plastic correction components based on trial stress
      mohr_coulomb->compute_df_dp(yield_type_trial, &state_variables,
                                  trial_stress, &df_dsigma_trial,
                                  &dp_dsigma_trial, &softening_trial);

      // Check plastic correction component based on trial stress
      // Check dFtrial/dSigma
      REQUIRE(df_dsigma_trial(0) == Approx(0.71907297).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(1) == Approx(0.14695243).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(2) == Approx(-0.28867513).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(3) == Approx(0.32506849).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(5) == Approx(0.).epsilon(Tolerance));
      // Check dPtrial/dSigma
      REQUIRE(dp_dsigma_trial(0) == Approx(0.59245977).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(1) == Approx(-0.07145868).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(2) == Approx(-0.25305189).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(3) == Approx(0.37722639).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(5) == Approx(0.).epsilon(Tolerance));
      // Check compute stress
      mpm::Material<Dim>::Vector6d updated_stress =
          mohr_coulomb->compute_stress(stress, dstrain, particle.get(),
                                       &state_variables);
      // Check update stress
      REQUIRE(updated_stress(0) == Approx(-2791.76726194).epsilon(Tolerance));
      REQUIRE(updated_stress(1) == Approx(-6838.46782671).epsilon(Tolerance));
      REQUIRE(updated_stress(2) == Approx(-9129.41035445).epsilon(Tolerance));
      REQUIRE(updated_stress(3) == Approx(3306.91731143).epsilon(Tolerance));
      REQUIRE(updated_stress(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(updated_stress(5) == Approx(0.).epsilon(Tolerance));
    }
  }
}

//! Check MohrCoulomb class in 2D
//! Cohesion, friction and dilation, with softening
TEST_CASE("MohrCoulomb is checked in 2D (c & phi & psi, with softening)",
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
  jmaterial["softening"] = true;
  jmaterial["friction"] = 30.;
  jmaterial["dilation"] = 15.;
  jmaterial["cohesion"] = 2000.;
  jmaterial["residual_friction"] = 0.;
  jmaterial["residual_dilation"] = 0.;
  jmaterial["residual_cohesion"] = 1000.;
  jmaterial["peak_epds"] = 1.E-16;
  jmaterial["critical_epds"] = 0.001;
  jmaterial["tension_cutoff"] = 0.;
  jmaterial["tolerance"] = 0.1;

  //! Check yield correction based on trial stress
  SECTION("MohrCoulomb check yield correction based on trial stress") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "MohrCoulomb2D", std::move(id), jmaterial);

    auto mohr_coulomb = std::make_shared<mpm::MohrCoulomb<Dim>>(id, jmaterial);

    REQUIRE(material->id() == 0);

    // Calculate modulus values
    const double K = material->property("youngs_modulus") /
                     (3.0 * (1. - 2. * material->property("poisson_ratio")));
    const double G = material->property("youngs_modulus") /
                     (2.0 * (1. + material->property("poisson_ratio")));
    const double a1 = K + (4.0 / 3.0) * G;
    const double a2 = K - (2.0 / 3.0) * G;
    // Compute elastic tensor
    mpm::Material<Dim>::Matrix6x6 de;
    de.setZero();
    de(0, 0) = a1;
    de(0, 1) = a2;
    de(0, 2) = a2;
    de(1, 0) = a2;
    de(1, 1) = a1;
    de(1, 2) = a2;
    de(2, 0) = a2;
    de(2, 1) = a2;
    de(2, 2) = a1;
    de(3, 3) = G;
    de(4, 4) = G;
    de(5, 5) = G;

    // Initialise state variables
    mpm::dense_map state_variables = material->initialise_state_variables();

    //! Check for shear failure (epde_peak < epds < epde_residual)
    SECTION("Check for shear failure (epde_peak < epds < epde_residual)") {

      // Tolerance for computation of stress
      const double Tolerance_stress = 1.E-5;

      // Initialise stress
      mpm::Material<Dim>::Vector6d stress;
      stress.setZero();
      stress(0) = -5000.;
      stress(1) = -6000.;
      stress(2) = -7000.;
      stress(3) = -1000.;

      // Check if stress invariants is computed correctly based on stress
      REQUIRE(mohr_coulomb->compute_stress_invariants(
                  stress, &state_variables) == true);
      REQUIRE(state_variables.at("phi") ==
              Approx(0.52359878).epsilon(Tolerance));
      REQUIRE(state_variables.at("psi") ==
              Approx(0.26179939).epsilon(Tolerance));
      REQUIRE(state_variables.at("cohesion") ==
              Approx(jmaterial["cohesion"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("j2") == Approx(2000000.).epsilon(Tolerance));
      REQUIRE(state_variables.at("j3") ==
              Approx(1000000000.).epsilon(Tolerance));
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

      // Define current plastic strain
      state_variables.at("plastic_strain0") = 0.00005;
      state_variables.at("plastic_strain1") = 0.00006;
      state_variables.at("plastic_strain2") = 0.00007;
      state_variables.at("plastic_strain3") = 0.00008;
      state_variables.at("epds") = 0.00004761;
      // Modified MC parameters
      state_variables.at("phi") = 0.49867048772358;
      state_variables.at("psi") = 0.24933524386179;
      state_variables.at("cohesion") = 1952.39047714;
      // Initialise values of yield functions
      Eigen::Matrix<double, 2, 1> yield_function;
      auto yield_type =
          mohr_coulomb->compute_yield_state(&yield_function, &state_variables);
      // Check if yield function and yield state is computed correctly
      REQUIRE(yield_function(0) == Approx(-4381.96601125).epsilon(Tolerance));
      REQUIRE(yield_function(1) == Approx(-3561.03580708).epsilon(Tolerance));
      REQUIRE(yield_type == mohr_coulomb->FailureState::Elastic);

      // Initialise plastic correction components
      mpm::Material<Dim>::Vector6d df_dsigma, dp_dsigma;
      df_dsigma.setZero();
      dp_dsigma.setZero();
      double softening = 0.;
      // Compute plastic correction components
      mohr_coulomb->compute_df_dp(yield_type, &state_variables, stress,
                                  &df_dsigma, &dp_dsigma, &softening);
      // Check plastic correction component based on stress
      // Check dF/dSigma
      REQUIRE(df_dsigma(0) == Approx(0.60900389).epsilon(Tolerance));
      REQUIRE(df_dsigma(1) == Approx(0.23261879).epsilon(Tolerance));
      REQUIRE(df_dsigma(2) == Approx(-0.29704523).epsilon(Tolerance));
      REQUIRE(df_dsigma(3) == Approx(-0.3763851).epsilon(Tolerance));
      REQUIRE(df_dsigma(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma(5) == Approx(0.).epsilon(Tolerance));
      // Check dP/dSigma
      REQUIRE(dp_dsigma(0) == Approx(0.47388304).epsilon(Tolerance));
      REQUIRE(dp_dsigma(1) == Approx(0.03580579).epsilon(Tolerance));
      REQUIRE(dp_dsigma(2) == Approx(-0.25505488).epsilon(Tolerance));
      REQUIRE(dp_dsigma(3) == Approx(-0.43807725).epsilon(Tolerance));
      REQUIRE(dp_dsigma(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(dp_dsigma(5) == Approx(0.).epsilon(Tolerance));

      // Initialise incremental of strain
      mpm::Material<Dim>::Vector6d dstrain;
      dstrain.setZero();
      dstrain(0) = 0.0001;
      dstrain(1) = 0.;
      dstrain(2) = 0.;
      dstrain(3) = 0.002;
      // Compute trial stress
      mpm::Material<Dim>::Vector6d trial_stress = stress + de * dstrain;
      // Check if stress invariants is computed correctly based on trial stress
      REQUIRE(mohr_coulomb->compute_stress_invariants(
                  trial_stress, &state_variables) == true);
      REQUIRE(state_variables.at("phi") ==
              Approx(0.49867048772358).epsilon(Tolerance));
      REQUIRE(state_variables.at("psi") ==
              Approx(0.24933524386179).epsilon(Tolerance));
      REQUIRE(state_variables.at("cohesion") ==
              Approx(1952.39047714).epsilon(Tolerance));
      REQUIRE(state_variables.at("j2") ==
              Approx(46753451.6765286).epsilon(Tolerance));
      REQUIRE(state_variables.at("j3") ==
              Approx(56758188775.9403).epsilon(Tolerance));
      REQUIRE(state_variables.at("epsilon") ==
              Approx(-8948.92917244).epsilon(Tolerance));
      REQUIRE(state_variables.at("rho") ==
              Approx(9669.89676021).epsilon(Tolerance));
      REQUIRE(state_variables.at("theta") ==
              Approx(0.36378823).epsilon(Tolerance));
      REQUIRE(state_variables.at("epds") ==
              Approx(0.00004761).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain0") ==
              Approx(0.00005).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain1") ==
              Approx(0.00006).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain2") ==
              Approx(0.00007).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain3") ==
              Approx(0.00008).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain4") ==
              Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain5") ==
              Approx(0.).epsilon(Tolerance));
      // Initialise values of yield functions based on trial stress
      Eigen::Matrix<double, 2, 1> yield_function_trial;
      auto yield_type_trial = mohr_coulomb->compute_yield_state(
          &yield_function_trial, &state_variables);
      // Check if yield function and yield state is computed correctly
      REQUIRE(yield_function_trial(0) ==
              Approx(2212.05893238).epsilon(Tolerance));
      REQUIRE(yield_function_trial(1) ==
              Approx(3262.66672836).epsilon(Tolerance));
      REQUIRE(yield_type_trial == mohr_coulomb->FailureState::Shear);
      // Initialise plastic correction components based on trial stress
      mpm::Material<Dim>::Vector6d df_dsigma_trial, dp_dsigma_trial;
      df_dsigma_trial.setZero();
      dp_dsigma_trial.setZero();
      double softening_trial = 0.;
      // Compute plastic correction components based on trial stress
      mohr_coulomb->compute_df_dp(yield_type_trial, &state_variables,
                                  trial_stress, &df_dsigma_trial,
                                  &dp_dsigma_trial, &softening_trial);

      // Check plastic correction component based on trial stress
      // Check dFtrial/dSigma
      REQUIRE(df_dsigma_trial(0) == Approx(0.34689653).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(1) == Approx(0.19768091).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(2) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(3) == Approx(0.56442433).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(5) == Approx(0.).epsilon(Tolerance));
      // Check dPtrial/dSigma
      REQUIRE(dp_dsigma_trial(0) == Approx(0.29664754).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(1) == Approx(0.16277439).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(2) == Approx(-0.20478799).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(3) == Approx(0.50638973).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(5) == Approx(0.).epsilon(Tolerance));
      // Check compute stress
      mpm::Material<Dim>::Vector6d updated_stress =
          mohr_coulomb->compute_stress(stress, dstrain, particle.get(),
                                       &state_variables);
      // Check update stress
      REQUIRE(updated_stress(0) ==
              Approx(-10642.32149186).epsilon(Tolerance_stress));
      REQUIRE(updated_stress(1) ==
              Approx(-10492.92006243).epsilon(Tolerance_stress));
      REQUIRE(updated_stress(2) ==
              Approx(-6225.11965155).epsilon(Tolerance_stress));
      REQUIRE(updated_stress(3) ==
              Approx(3063.5902731).epsilon(Tolerance_stress));
      REQUIRE(updated_stress(4) == Approx(0.).epsilon(Tolerance_stress));
      REQUIRE(updated_stress(5) == Approx(0.).epsilon(Tolerance_stress));

      // Check plastic strain
      REQUIRE(state_variables.at("epds") ==
              Approx(0.00080550198218).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain0") ==
              Approx(0.00060269).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain1") ==
              Approx(0.00036327).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain2") ==
              Approx(-0.00031155).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain3") ==
              Approx(0.00102347).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain4") ==
              Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain5") ==
              Approx(0.).epsilon(Tolerance));
    }

    //! Check for shear failure (epds < epde_peak)
    SECTION("Check for shear failure (epds < epde_peak)") {
      // Initialise stress
      mpm::Material<Dim>::Vector6d stress;
      stress.setZero();
      stress(0) = -5000.;
      stress(1) = -6000.;
      stress(2) = -7000.;
      stress(3) = -1000.;

      // Check if stress invariants is computed correctly based on stress
      REQUIRE(mohr_coulomb->compute_stress_invariants(
                  stress, &state_variables) == true);
      REQUIRE(state_variables.at("phi") ==
              Approx(0.52359878).epsilon(Tolerance));
      REQUIRE(state_variables.at("psi") ==
              Approx(0.26179939).epsilon(Tolerance));
      REQUIRE(state_variables.at("cohesion") ==
              Approx(jmaterial["cohesion"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("j2") == Approx(2000000.).epsilon(Tolerance));
      REQUIRE(state_variables.at("j3") ==
              Approx(1000000000.).epsilon(Tolerance));
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
      REQUIRE(yield_function(1) == Approx(-3774.1679421).epsilon(Tolerance));
      REQUIRE(yield_type == mohr_coulomb->FailureState::Elastic);

      // Initialise plastic correction components
      mpm::Material<Dim>::Vector6d df_dsigma, dp_dsigma;
      df_dsigma.setZero();
      dp_dsigma.setZero();
      double softening = 0.;
      // Compute plastic correction components
      mohr_coulomb->compute_df_dp(yield_type, &state_variables, stress,
                                  &df_dsigma, &dp_dsigma, &softening);
      // Check plastic correction component based on stress
      // Check dF/dSigma
      REQUIRE(df_dsigma(0) == Approx(0.62666187).epsilon(Tolerance));
      REQUIRE(df_dsigma(1) == Approx(0.23936353).epsilon(Tolerance));
      REQUIRE(df_dsigma(2) == Approx(-0.28867513).epsilon(Tolerance));
      REQUIRE(df_dsigma(3) == Approx(-0.38729833).epsilon(Tolerance));
      REQUIRE(df_dsigma(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma(5) == Approx(0.).epsilon(Tolerance));
      // Check dP/dSigma
      REQUIRE(dp_dsigma(0) == Approx(0.48581658).epsilon(Tolerance));
      REQUIRE(dp_dsigma(1) == Approx(0.03777283).epsilon(Tolerance));
      REQUIRE(dp_dsigma(2) == Approx(-0.25564021).epsilon(Tolerance));
      REQUIRE(dp_dsigma(3) == Approx(-0.44804375).epsilon(Tolerance));
      REQUIRE(dp_dsigma(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(dp_dsigma(5) == Approx(0.).epsilon(Tolerance));

      // Initialise incremental of strain
      mpm::Material<Dim>::Vector6d dstrain;
      dstrain.setZero();
      dstrain(0) = 0.0001;
      dstrain(1) = 0.;
      dstrain(2) = 0.;
      dstrain(3) = 0.002;
      // Compute trial stress
      mpm::Material<Dim>::Vector6d trial_stress = stress + de * dstrain;
      // Check if stress invariants is computed correctly based on trial stress
      REQUIRE(mohr_coulomb->compute_stress_invariants(
                  trial_stress, &state_variables) == true);
      REQUIRE(state_variables.at("j2") ==
              Approx(46753451.6765286).epsilon(Tolerance));
      REQUIRE(state_variables.at("j3") ==
              Approx(56758188775.9403).epsilon(Tolerance));
      REQUIRE(state_variables.at("epsilon") ==
              Approx(-8948.92917244).epsilon(Tolerance));
      REQUIRE(state_variables.at("rho") ==
              Approx(9669.89676021).epsilon(Tolerance));
      REQUIRE(state_variables.at("theta") ==
              Approx(0.36378823).epsilon(Tolerance));
      // Initialise values of yield functions based on trial stress
      Eigen::Matrix<double, 2, 1> yield_function_trial;
      auto yield_type_trial = mohr_coulomb->compute_yield_state(
          &yield_function_trial, &state_variables);
      // Check if yield function and yield state is computed correctly
      REQUIRE(yield_function_trial(0) ==
              Approx(2212.05893238).epsilon(Tolerance));
      REQUIRE(yield_function_trial(1) ==
              Approx(3174.54763108).epsilon(Tolerance));
      REQUIRE(yield_type_trial == mohr_coulomb->FailureState::Shear);
      // Initialise plastic correction components based on trial stress
      mpm::Material<Dim>::Vector6d df_dsigma_trial, dp_dsigma_trial;
      df_dsigma_trial.setZero();
      dp_dsigma_trial.setZero();
      double softening_trial = 0.;
      // Compute plastic correction components based on trial stress
      mohr_coulomb->compute_df_dp(yield_type_trial, &state_variables,
                                  trial_stress, &df_dsigma_trial,
                                  &dp_dsigma_trial, &softening_trial);

      // Check plastic correction component based on trial stress
      // Check dFtrial/dSigma
      REQUIRE(df_dsigma_trial(0) == Approx(0.36433344).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(1) == Approx(0.21301683).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(2) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(3) == Approx(0.57237152).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(5) == Approx(0.).epsilon(Tolerance));
      // Check dPtrial/dSigma
      REQUIRE(dp_dsigma_trial(0) == Approx(0.3081448).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(1) == Approx(0.17193919).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(2) == Approx(-0.21213479).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(3) == Approx(0.51521254).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(5) == Approx(0.).epsilon(Tolerance));
      // Check compute stress
      mpm::Material<Dim>::Vector6d updated_stress =
          mohr_coulomb->compute_stress(stress, dstrain, particle.get(),
                                       &state_variables);
      // Check update stress
      REQUIRE(updated_stress(0) == Approx(-7573.14655485).epsilon(Tolerance));
      REQUIRE(updated_stress(1) == Approx(-8293.81378005).epsilon(Tolerance));
      REQUIRE(updated_stress(2) == Approx(-6337.06362775).epsilon(Tolerance));
      REQUIRE(updated_stress(3) == Approx(4709.1549025).epsilon(Tolerance));
      REQUIRE(updated_stress(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(updated_stress(5) == Approx(0.).epsilon(Tolerance));

      // Check plastic strain
      REQUIRE(state_variables.at("epds") ==
              Approx(0.00043108).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain0") ==
              Approx(0.00030839).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain1") ==
              Approx(0.00017208).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain2") ==
              Approx(-0.0002123).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain3") ==
              Approx(0.00051562).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain4") ==
              Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain5") ==
              Approx(0.).epsilon(Tolerance));
    }

    //! Check for shear failure (epds > epde_residual)
    SECTION("Check for shear failure (epds > epde_residual)") {
      // Initialise stress
      mpm::Material<Dim>::Vector6d stress;
      stress.setZero();
      stress(0) = -5000.;
      stress(1) = -6000.;
      stress(2) = -6500.;
      stress(3) = 0.;

      // Check if stress invariants is computed correctly based on stress
      REQUIRE(mohr_coulomb->compute_stress_invariants(
                  stress, &state_variables) == true);
      REQUIRE(state_variables.at("phi") ==
              Approx(0.52359878).epsilon(Tolerance));
      REQUIRE(state_variables.at("psi") ==
              Approx(0.26179939).epsilon(Tolerance));
      REQUIRE(state_variables.at("cohesion") ==
              Approx(jmaterial["cohesion"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("j2") ==
              Approx(583333.33333333).epsilon(Tolerance));
      REQUIRE(state_variables.at("j3") ==
              Approx(92592592.5925928).epsilon(Tolerance));
      REQUIRE(state_variables.at("epsilon") ==
              Approx(-10103.62971082).epsilon(Tolerance));
      REQUIRE(state_variables.at("rho") ==
              Approx(1080.12344973).epsilon(Tolerance));
      REQUIRE(state_variables.at("theta") ==
              Approx(0.33347317).epsilon(Tolerance));
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

      // Define current plastic strain
      state_variables.at("plastic_strain0") = 0.001;
      state_variables.at("plastic_strain1") = 0.002;
      state_variables.at("plastic_strain2") = 0.003;
      state_variables.at("plastic_strain3") = 0.001;
      state_variables.at("epds") = 0.00129099;
      // Modified MC parameters
      state_variables.at("phi") = 0.;
      state_variables.at("psi") = 0.;
      state_variables.at("cohesion") = 1000.;

      // Initialise values of yield functions
      Eigen::Matrix<double, 2, 1> yield_function;
      auto yield_type =
          mohr_coulomb->compute_yield_state(&yield_function, &state_variables);
      // Check if yield function and yield state is computed correctly
      REQUIRE(yield_function(0) == Approx(-5000.).epsilon(Tolerance));
      REQUIRE(yield_function(1) == Approx(-250.).epsilon(Tolerance));
      REQUIRE(yield_type == mohr_coulomb->FailureState::Elastic);

      // Initialise plastic correction components
      mpm::Material<Dim>::Vector6d df_dsigma, dp_dsigma;
      df_dsigma.setZero();
      dp_dsigma.setZero();
      double softening = 0.;
      // Compute plastic correction components
      mohr_coulomb->compute_df_dp(yield_type, &state_variables, stress,
                                  &df_dsigma, &dp_dsigma, &softening);
      // Check plastic correction component based on stress
      // Check dF/dSigma
      REQUIRE(df_dsigma(0) == Approx(0.5).epsilon(Tolerance));
      REQUIRE(df_dsigma(1) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma(2) == Approx(-0.5).epsilon(Tolerance));
      REQUIRE(df_dsigma(3) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma(5) == Approx(0.).epsilon(Tolerance));
      // Check dP/dSigma
      REQUIRE(dp_dsigma(0) == Approx(0.47245559).epsilon(Tolerance));
      REQUIRE(dp_dsigma(1) == Approx(-0.09449112).epsilon(Tolerance));
      REQUIRE(dp_dsigma(2) == Approx(-0.37796447).epsilon(Tolerance));
      REQUIRE(dp_dsigma(3) == Approx(0.).epsilon(Tolerance));
      REQUIRE(dp_dsigma(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(dp_dsigma(5) == Approx(0.).epsilon(Tolerance));

      // Initialise incremental of strain
      mpm::Material<Dim>::Vector6d dstrain;
      dstrain.setZero();
      dstrain(0) = 0.0001;
      dstrain(1) = 0.;
      dstrain(2) = 0.;
      dstrain(3) = 0.002;
      // Compute trial stress
      mpm::Material<Dim>::Vector6d trial_stress = stress + de * dstrain;
      // Check if stress invariants is computed correctly based on trial stress
      REQUIRE(mohr_coulomb->compute_stress_invariants(
                  trial_stress, &state_variables) == true);
      REQUIRE(state_variables.at("phi") == Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("psi") == Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("cohesion") ==
              Approx(1000.).epsilon(Tolerance));
      REQUIRE(state_variables.at("j2") ==
              Approx(60593195.2662722).epsilon(Tolerance));
      REQUIRE(state_variables.at("j3") ==
              Approx(55145653163.4046).epsilon(Tolerance));
      REQUIRE(state_variables.at("epsilon") ==
              Approx(-8660.25403784).epsilon(Tolerance));
      REQUIRE(state_variables.at("rho") ==
              Approx(11008.46903673).epsilon(Tolerance));
      REQUIRE(state_variables.at("theta") ==
              Approx(0.42072067).epsilon(Tolerance));
      REQUIRE(state_variables.at("epds") ==
              Approx(0.00129099).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain0") ==
              Approx(0.001).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain1") ==
              Approx(0.002).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain2") ==
              Approx(0.003).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain3") ==
              Approx(0.001).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain4") ==
              Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain5") ==
              Approx(0.).epsilon(Tolerance));
      // Initialise values of yield functions based on trial stress
      Eigen::Matrix<double, 2, 1> yield_function_trial;
      auto yield_type_trial = mohr_coulomb->compute_yield_state(
          &yield_function_trial, &state_variables);
      // Check if yield function and yield state is computed correctly
      REQUIRE(yield_function_trial(0) ==
              Approx(3204.54446772).epsilon(Tolerance));
      REQUIRE(yield_function_trial(1) ==
              Approx(6743.00600619).epsilon(Tolerance));
      REQUIRE(yield_type_trial == mohr_coulomb->FailureState::Shear);
      // Initialise plastic correction components based on trial stress
      mpm::Material<Dim>::Vector6d df_dsigma_trial, dp_dsigma_trial;
      df_dsigma_trial.setZero();
      dp_dsigma_trial.setZero();
      double softening_trial = 0.;
      // Compute plastic correction components based on trial stress
      mohr_coulomb->compute_df_dp(yield_type_trial, &state_variables,
                                  trial_stress, &df_dsigma_trial,
                                  &dp_dsigma_trial, &softening_trial);

      // Check plastic correction component based on trial stress
      // Check dFtrial/dSigma
      REQUIRE(df_dsigma_trial(0) == Approx(0.05712351).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(1) == Approx(-0.05712351).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(2) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(3) == Approx(0.49672619).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(5) == Approx(0.).epsilon(Tolerance));
      // Check dPtrial/dSigma
      REQUIRE(dp_dsigma_trial(0) == Approx(0.07488303).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(1) == Approx(-0.02353467).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(2) == Approx(-0.05134836).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(3) == Approx(0.42790303).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(5) == Approx(0.).epsilon(Tolerance));
      // Check compute stress
      mpm::Material<Dim>::Vector6d updated_stress =
          mohr_coulomb->compute_stress(stress, dstrain, particle.get(),
                                       &state_variables);
      // Check update stress
      REQUIRE(updated_stress(0) == Approx(-3917.72356097).epsilon(Tolerance));
      REQUIRE(updated_stress(1) == Approx(-5340.1440237).epsilon(Tolerance));
      REQUIRE(updated_stress(2) == Approx(-5742.13241534).epsilon(Tolerance));
      REQUIRE(updated_stress(3) == Approx(529.26866296).epsilon(Tolerance));
      REQUIRE(updated_stress(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(updated_stress(5) == Approx(0.).epsilon(Tolerance));

      // Check plastic strain
      REQUIRE(state_variables.at("epds") ==
              Approx(0.00199714).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain0") ==
              Approx(0.0010343).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain1") ==
              Approx(0.00198922).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain2") ==
              Approx(0.00297648).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain3") ==
              Approx(0.00286239).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain4") ==
              Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain5") ==
              Approx(0.).epsilon(Tolerance));
    }
  }

  //! Check yield correction based on current stress
  SECTION("MohrCoulomb check yield correction based on current stress") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "MohrCoulomb2D", std::move(id), jmaterial);

    auto mohr_coulomb = std::make_shared<mpm::MohrCoulomb<Dim>>(id, jmaterial);

    REQUIRE(material->id() == 0);

    // Calculate modulus values
    const double K = material->property("youngs_modulus") /
                     (3.0 * (1. - 2. * material->property("poisson_ratio")));
    const double G = material->property("youngs_modulus") /
                     (2.0 * (1. + material->property("poisson_ratio")));
    const double a1 = K + (4.0 / 3.0) * G;
    const double a2 = K - (2.0 / 3.0) * G;
    // Compute elastic tensor
    mpm::Material<Dim>::Matrix6x6 de;
    de.setZero();
    de(0, 0) = a1;
    de(0, 1) = a2;
    de(0, 2) = a2;
    de(1, 0) = a2;
    de(1, 1) = a1;
    de(1, 2) = a2;
    de(2, 0) = a2;
    de(2, 1) = a2;
    de(2, 2) = a1;
    de(3, 3) = G;
    de(4, 4) = G;
    de(5, 5) = G;

    // Initialise state variables
    mpm::dense_map state_variables = material->initialise_state_variables();

    //! Check for shear failure (epde_peak < epds < epde_residual)
    SECTION("Check for shear failure (epde_peak < epds < epde_residual)") {

      // Tolerance for computation of stress
      const double Tolerance_stress = 1.E-5;

      // Initialise stress
      mpm::Material<Dim>::Vector6d stress;
      stress.setZero();
      stress(0) = -5000.;
      stress(1) = -6000.;
      stress(2) = -7000.;
      stress(3) = -4186.6;

      // Check if stress invariants is computed correctly based on stress
      REQUIRE(mohr_coulomb->compute_stress_invariants(
                  stress, &state_variables) == true);
      REQUIRE(state_variables.at("phi") ==
              Approx(0.52359878).epsilon(Tolerance));
      REQUIRE(state_variables.at("psi") ==
              Approx(0.26179939).epsilon(Tolerance));
      REQUIRE(state_variables.at("cohesion") ==
              Approx(jmaterial["cohesion"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("j2") ==
              Approx(18527619.56).epsilon(Tolerance));
      REQUIRE(state_variables.at("j3") ==
              Approx(17527619560.).epsilon(Tolerance));
      REQUIRE(state_variables.at("epsilon") ==
              Approx(-10392.30484541).epsilon(Tolerance));
      REQUIRE(state_variables.at("rho") ==
              Approx(6087.30146452).epsilon(Tolerance));
      REQUIRE(state_variables.at("theta") ==
              Approx(0.32101934).epsilon(Tolerance));
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

      // Define current plastic strain
      state_variables.at("plastic_strain0") = 0.0001;
      state_variables.at("plastic_strain1") = 0.0002;
      state_variables.at("plastic_strain2") = 0.00007;
      state_variables.at("plastic_strain3") = 0.00008;
      state_variables.at("epds") = 0.00009117;
      // Modified MC parameters
      state_variables.at("phi") = 0.47586473847588;
      state_variables.at("psi") = 0.23793236923794;
      state_variables.at("cohesion") = 1908.83470445882;
      // Initialise values of yield functions
      Eigen::Matrix<double, 2, 1> yield_function;
      auto yield_type =
          mohr_coulomb->compute_yield_state(&yield_function, &state_variables);
      // Check if yield function and yield state is computed correctly
      REQUIRE(yield_function(0) == Approx(-1283.64854880).epsilon(Tolerance));
      REQUIRE(yield_function(1) == Approx(0.00464538).epsilon(Tolerance));
      REQUIRE(yield_type == mohr_coulomb->FailureState::Shear);

      // Initialise plastic correction components
      mpm::Material<Dim>::Vector6d df_dsigma, dp_dsigma;
      df_dsigma.setZero();
      dp_dsigma.setZero();
      double softening = 0.;
      // Compute plastic correction components
      mohr_coulomb->compute_df_dp(yield_type, &state_variables, stress,
                                  &df_dsigma, &dp_dsigma, &softening);
      // Check plastic correction component based on stress
      // Check dF/dSigma
      REQUIRE(df_dsigma(0) == Approx(0.32438701).epsilon(Tolerance));
      REQUIRE(df_dsigma(1) == Approx(0.19097903).epsilon(Tolerance));
      REQUIRE(df_dsigma(2) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma(3) == Approx(-0.55852584).epsilon(Tolerance));
      REQUIRE(df_dsigma(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma(5) == Approx(0.).epsilon(Tolerance));
      // Check dP/dSigma
      REQUIRE(dp_dsigma(0) == Approx(0.28406168).epsilon(Tolerance));
      REQUIRE(dp_dsigma(1) == Approx(0.16523414).epsilon(Tolerance));
      REQUIRE(dp_dsigma(2) == Approx(-0.20676947).epsilon(Tolerance));
      REQUIRE(dp_dsigma(3) == Approx(-0.4974834).epsilon(Tolerance));
      REQUIRE(dp_dsigma(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(dp_dsigma(5) == Approx(0.).epsilon(Tolerance));

      // Initialise incremental of strain
      mpm::Material<Dim>::Vector6d dstrain;
      dstrain.setZero();
      dstrain(0) = 0.0005;
      dstrain(1) = -0.0005;
      dstrain(2) = 0.;
      dstrain(3) = 0.0001;
      // Compute trial stress
      mpm::Material<Dim>::Vector6d trial_stress = stress + de * dstrain;
      // Check if stress invariants is computed correctly based on trial stress
      REQUIRE(mohr_coulomb->compute_stress_invariants(
                  trial_stress, &state_variables) == true);
      REQUIRE(state_variables.at("phi") ==
              Approx(0.47586473847588).epsilon(Tolerance));
      REQUIRE(state_variables.at("psi") ==
              Approx(0.23793236923794).epsilon(Tolerance));
      REQUIRE(state_variables.at("cohesion") ==
              Approx(1908.83470445882).epsilon(Tolerance));
      REQUIRE(state_variables.at("j2") ==
              Approx(34094140.2700592).epsilon(Tolerance));
      REQUIRE(state_variables.at("j3") ==
              Approx(33094140270.0592).epsilon(Tolerance));
      REQUIRE(state_variables.at("epsilon") ==
              Approx(-10392.30484541).epsilon(Tolerance));
      REQUIRE(state_variables.at("rho") ==
              Approx(8257.6195444).epsilon(Tolerance));
      REQUIRE(state_variables.at("theta") ==
              Approx(0.3747326).epsilon(Tolerance));
      REQUIRE(state_variables.at("epds") ==
              Approx(0.00009117).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain0") ==
              Approx(0.0001).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain1") ==
              Approx(0.0002).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain2") ==
              Approx(0.00007).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain3") ==
              Approx(0.00008).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain4") ==
              Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain5") ==
              Approx(0.).epsilon(Tolerance));

      // Initialise values of yield functions based on trial stress
      Eigen::Matrix<double, 2, 1> yield_function_trial;
      auto yield_type_trial = mohr_coulomb->compute_yield_state(
          &yield_function_trial, &state_variables);
      // Check if yield function and yield state is computed correctly
      REQUIRE(yield_function_trial(0) ==
              Approx(274.43852423).epsilon(Tolerance));
      REQUIRE(yield_function_trial(1) ==
              Approx(1752.8367834).epsilon(Tolerance));
      REQUIRE(yield_type_trial == mohr_coulomb->FailureState::Shear);
      // Initialise plastic correction components based on trial stress
      mpm::Material<Dim>::Vector6d df_dsigma_trial, dp_dsigma_trial;
      df_dsigma_trial.setZero();
      dp_dsigma_trial.setZero();
      double softening_trial = 0.;
      // Compute plastic correction components based on trial stress
      mohr_coulomb->compute_df_dp(yield_type_trial, &state_variables,
                                  trial_stress, &df_dsigma_trial,
                                  &dp_dsigma_trial, &softening_trial);

      // Check plastic correction component based on trial stress
      // Check dFtrial/dSigma
      REQUIRE(df_dsigma_trial(0) == Approx(0.68104703).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(1) == Approx(-0.16568099).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(2) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(3) == Approx(-0.37035584).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(5) == Approx(0.).epsilon(Tolerance));
      // Check dPtrial/dSigma
      REQUIRE(dp_dsigma_trial(0) == Approx(0.59829333).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(1) == Approx(-0.15988784).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(2) == Approx(-0.19587913).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(3) == Approx(-0.33162576).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(5) == Approx(0.).epsilon(Tolerance));
      // Check compute stress
      mpm::Material<Dim>::Vector6d updated_stress =
          mohr_coulomb->compute_stress(stress, dstrain, particle.get(),
                                       &state_variables);
      // Check update stress
      REQUIRE(updated_stress(0) ==
              Approx(-3309.71699019).epsilon(Tolerance_stress));
      REQUIRE(updated_stress(1) ==
              Approx(-11067.2682187).epsilon(Tolerance_stress));
      REQUIRE(updated_stress(2) ==
              Approx(-6931.40841882).epsilon(Tolerance_stress));
      REQUIRE(updated_stress(3) ==
              Approx(-2820.47296506).epsilon(Tolerance_stress));
      REQUIRE(updated_stress(4) == Approx(0.).epsilon(Tolerance_stress));
      REQUIRE(updated_stress(5) == Approx(0.).epsilon(Tolerance_stress));

      // Check plastic strain
      REQUIRE(state_variables.at("epds") ==
              Approx(0.00022948).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain0") ==
              Approx(0.00028101).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain1") ==
              Approx(0.00025949).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain2") ==
              Approx(-0.00003817).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain3") ==
              Approx(-0.00017519).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain4") ==
              Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain5") ==
              Approx(0.).epsilon(Tolerance));
    }

    //! Check for shear failure (epds < epde_peak)
    SECTION("Check for shear failure (epds < epde_peak)") {
      // Initialise stress
      mpm::Material<Dim>::Vector6d stress;
      stress.setZero();
      stress(0) = -1000.;
      stress(1) = -6000.;
      stress(2) = -9350.4;
      stress(3) = -1000.;

      // Check if stress invariants is computed correctly based on stress
      REQUIRE(mohr_coulomb->compute_stress_invariants(
                  stress, &state_variables) == true);
      REQUIRE(state_variables.at("phi") ==
              Approx(0.52359878).epsilon(Tolerance));
      REQUIRE(state_variables.at("psi") ==
              Approx(0.26179939).epsilon(Tolerance));
      REQUIRE(state_variables.at("cohesion") ==
              Approx(jmaterial["cohesion"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("j2") ==
              Approx(18659060.0533333).epsilon(Tolerance));
      REQUIRE(state_variables.at("j3") ==
              Approx(13444141125.3286).epsilon(Tolerance));
      REQUIRE(state_variables.at("epsilon") ==
              Approx(-9439.90784136).epsilon(Tolerance));
      REQUIRE(state_variables.at("rho") ==
              Approx(6108.85587542).epsilon(Tolerance));
      REQUIRE(state_variables.at("theta") ==
              Approx(0.37419232).epsilon(Tolerance));
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
      REQUIRE(yield_function(0) == Approx(-807.41759643).epsilon(Tolerance));
      REQUIRE(yield_function(1) == Approx(-0.01617146).epsilon(Tolerance));
      REQUIRE(yield_type == mohr_coulomb->FailureState::Shear);

      // Initialise plastic correction components
      mpm::Material<Dim>::Vector6d df_dsigma, dp_dsigma;
      df_dsigma.setZero();
      dp_dsigma.setZero();
      double softening = 0.;
      // Compute plastic correction components
      mohr_coulomb->compute_df_dp(yield_type, &state_variables, stress,
                                  &df_dsigma, &dp_dsigma, &softening);
      // Check plastic correction component based on stress
      // Check dF/dSigma
      REQUIRE(df_dsigma(0) == Approx(0.8350549).epsilon(Tolerance));
      REQUIRE(df_dsigma(1) == Approx(0.0309705).epsilon(Tolerance));
      REQUIRE(df_dsigma(2) == Approx(-0.28867513).epsilon(Tolerance));
      REQUIRE(df_dsigma(3) == Approx(-0.16081688).epsilon(Tolerance));
      REQUIRE(df_dsigma(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma(5) == Approx(0.).epsilon(Tolerance));
      // Check dP/dSigma
      REQUIRE(dp_dsigma(0) == Approx(0.72452524).epsilon(Tolerance));
      REQUIRE(dp_dsigma(1) == Approx(-0.17541555).epsilon(Tolerance));
      REQUIRE(dp_dsigma(2) == Approx(-0.28116049).epsilon(Tolerance));
      REQUIRE(dp_dsigma(3) == Approx(-0.17998816).epsilon(Tolerance));
      REQUIRE(dp_dsigma(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(dp_dsigma(5) == Approx(0.).epsilon(Tolerance));

      // Initialise incremental of strain
      mpm::Material<Dim>::Vector6d dstrain;
      dstrain.setZero();
      dstrain(0) = 0.0001;
      dstrain(1) = 0.;
      dstrain(2) = 0.;
      dstrain(3) = 0.001;
      // Compute trial stress
      mpm::Material<Dim>::Vector6d trial_stress = stress + de * dstrain;
      // Check if stress invariants is computed correctly based on trial stress
      REQUIRE(mohr_coulomb->compute_stress_invariants(
                  trial_stress, &state_variables) == true);
      REQUIRE(state_variables.at("j2") ==
              Approx(29380069.9152663).epsilon(Tolerance));
      REQUIRE(state_variables.at("j3") ==
              Approx(50304548357.2239).epsilon(Tolerance));
      REQUIRE(state_variables.at("epsilon") ==
              Approx(-7996.53216838).epsilon(Tolerance));
      REQUIRE(state_variables.at("rho") ==
              Approx(7665.51627945).epsilon(Tolerance));
      REQUIRE(state_variables.at("theta") ==
              Approx(0.20272541).epsilon(Tolerance));
      // Initialise values of yield functions based on trial stress
      Eigen::Matrix<double, 2, 1> yield_function_trial;
      auto yield_type_trial = mohr_coulomb->compute_yield_state(
          &yield_function_trial, &state_variables);
      // Check if yield function and yield state is computed correctly
      REQUIRE(yield_function_trial(0) ==
              Approx(1513.89550814).epsilon(Tolerance));
      REQUIRE(yield_function_trial(1) ==
              Approx(1843.75660036).epsilon(Tolerance));
      REQUIRE(yield_type_trial == mohr_coulomb->FailureState::Shear);
      // Initialise plastic correction components based on trial stress
      mpm::Material<Dim>::Vector6d df_dsigma_trial, dp_dsigma_trial;
      df_dsigma_trial.setZero();
      dp_dsigma_trial.setZero();
      double softening_trial = 0.;
      // Compute plastic correction components based on trial stress
      mohr_coulomb->compute_df_dp(yield_type_trial, &state_variables,
                                  trial_stress, &df_dsigma_trial,
                                  &dp_dsigma_trial, &softening_trial);

      // Check plastic correction component based on trial stress
      // Check dFtrial/dSigma
      REQUIRE(df_dsigma_trial(0) == Approx(0.74124692).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(1) == Approx(0.12477848).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(2) == Approx(-0.28867513).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(3) == Approx(0.30412443).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(5) == Approx(0.).epsilon(Tolerance));
      // Check dPtrial/dSigma
      REQUIRE(dp_dsigma_trial(0) == Approx(0.61938609).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(1) == Approx(-0.09026886).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(2) == Approx(-0.26116803).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(3) == Approx(0.35009644).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(5) == Approx(0.).epsilon(Tolerance));
      // Check compute stress
      mpm::Material<Dim>::Vector6d updated_stress =
          mohr_coulomb->compute_stress(stress, dstrain, particle.get(),
                                       &state_variables);
      // Check update stress
      REQUIRE(updated_stress(0) == Approx(-1909.04417786).epsilon(Tolerance));
      REQUIRE(updated_stress(1) == Approx(-5722.13370553).epsilon(Tolerance));
      REQUIRE(updated_stress(2) == Approx(-8601.25290489).epsilon(Tolerance));
      REQUIRE(updated_stress(3) == Approx(2466.78909344).epsilon(Tolerance));
      REQUIRE(updated_stress(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(updated_stress(5) == Approx(0.).epsilon(Tolerance));

      // Check plastic strain
      REQUIRE(state_variables.at("epds") ==
              Approx(0.00020144).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain0") ==
              Approx(0.00022171).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain1") ==
              Approx(-0.00003258).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain2") ==
              Approx(-0.00009385).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain3") ==
              Approx(0.00009863).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain4") ==
              Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain5") ==
              Approx(0.).epsilon(Tolerance));
    }

    //! Check for shear failure (epds > epde_residual)
    SECTION("Check for shear failure (epds > epde_residual)") {
      // Initialise stress
      mpm::Material<Dim>::Vector6d stress;
      stress.setZero();
      stress(0) = -5000.;
      stress(1) = -6000.;
      stress(2) = -7000.;
      stress(3) = 0;

      // Check if stress invariants is computed correctly based on stress
      REQUIRE(mohr_coulomb->compute_stress_invariants(
                  stress, &state_variables) == true);
      REQUIRE(state_variables.at("phi") ==
              Approx(0.52359878).epsilon(Tolerance));
      REQUIRE(state_variables.at("psi") ==
              Approx(0.26179939).epsilon(Tolerance));
      REQUIRE(state_variables.at("cohesion") ==
              Approx(jmaterial["cohesion"]).epsilon(Tolerance));
      REQUIRE(state_variables.at("j2") == Approx(1000000.).epsilon(Tolerance));
      REQUIRE(state_variables.at("j3") == Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("epsilon") ==
              Approx(-10392.30484541).epsilon(Tolerance));
      REQUIRE(state_variables.at("rho") ==
              Approx(1414.21356237).epsilon(Tolerance));
      REQUIRE(state_variables.at("theta") ==
              Approx(0.52359878).epsilon(Tolerance));
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

      // Define current plastic strain
      state_variables.at("plastic_strain0") = 0.001;
      state_variables.at("plastic_strain1") = 0.002;
      state_variables.at("plastic_strain2") = 0.003;
      state_variables.at("plastic_strain3") = 0.001;
      state_variables.at("epds") = 0.00129099;
      // Modified MC parameters
      state_variables.at("phi") = 0.;
      state_variables.at("psi") = 0.;
      state_variables.at("cohesion") = 1000.;

      // Initialise values of yield functions
      Eigen::Matrix<double, 2, 1> yield_function;
      auto yield_type =
          mohr_coulomb->compute_yield_state(&yield_function, &state_variables);
      // Check if yield function and yield state is computed correctly
      REQUIRE(yield_function(0) == Approx(-5000.).epsilon(Tolerance));
      REQUIRE(yield_function(1) == Approx(0.).epsilon(Tolerance));
      REQUIRE(yield_type == mohr_coulomb->FailureState::Shear);

      // Initialise plastic correction components
      mpm::Material<Dim>::Vector6d df_dsigma, dp_dsigma;
      df_dsigma.setZero();
      dp_dsigma.setZero();
      double softening = 0.;
      // Compute plastic correction components
      mohr_coulomb->compute_df_dp(yield_type, &state_variables, stress,
                                  &df_dsigma, &dp_dsigma, &softening);
      // Check plastic correction component based on stress
      // Check dF/dSigma
      REQUIRE(df_dsigma(0) == Approx(0.5).epsilon(Tolerance));
      REQUIRE(df_dsigma(1) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma(2) == Approx(-0.5).epsilon(Tolerance));
      REQUIRE(df_dsigma(3) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma(5) == Approx(0.).epsilon(Tolerance));
      // Check dP/dSigma
      REQUIRE(dp_dsigma(0) == Approx(0.43301270).epsilon(Tolerance));
      REQUIRE(dp_dsigma(1) == Approx(0.).epsilon(Tolerance));
      REQUIRE(dp_dsigma(2) == Approx(-0.4330127).epsilon(Tolerance));
      REQUIRE(dp_dsigma(3) == Approx(0.).epsilon(Tolerance));
      REQUIRE(dp_dsigma(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(dp_dsigma(5) == Approx(0.).epsilon(Tolerance));

      // Initialise incremental of strain
      mpm::Material<Dim>::Vector6d dstrain;
      dstrain.setZero();
      dstrain(0) = 0.0001;
      dstrain(1) = 0.;
      dstrain(2) = 0.;
      dstrain(3) = 0.002;
      // Compute trial stress
      mpm::Material<Dim>::Vector6d trial_stress = stress + de * dstrain;
      // Check if stress invariants is computed correctly based on trial stress
      REQUIRE(mohr_coulomb->compute_stress_invariants(
                  trial_stress, &state_variables) == true);
      REQUIRE(state_variables.at("phi") == Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("psi") == Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("cohesion") ==
              Approx(1000.).epsilon(Tolerance));
      REQUIRE(state_variables.at("j2") ==
              Approx(61138067.061144).epsilon(Tolerance));
      REQUIRE(state_variables.at("j3") ==
              Approx(74831167079.6878).epsilon(Tolerance));
      REQUIRE(state_variables.at("epsilon") ==
              Approx(-8948.92917244).epsilon(Tolerance));
      REQUIRE(state_variables.at("rho") ==
              Approx(11057.85395646).epsilon(Tolerance));
      REQUIRE(state_variables.at("theta") ==
              Approx(0.38398831).epsilon(Tolerance));
      REQUIRE(state_variables.at("epds") ==
              Approx(0.00129099).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain0") ==
              Approx(0.001).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain1") ==
              Approx(0.002).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain2") ==
              Approx(0.003).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain3") ==
              Approx(0.001).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain4") ==
              Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain5") ==
              Approx(0.).epsilon(Tolerance));
      // Initialise values of yield functions based on trial stress
      Eigen::Matrix<double, 2, 1> yield_function_trial;
      auto yield_type_trial = mohr_coulomb->compute_yield_state(
          &yield_function_trial, &state_variables);
      // Check if yield function and yield state is computed correctly
      REQUIRE(yield_function_trial(0) ==
              Approx(3204.54446772).epsilon(Tolerance));
      REQUIRE(yield_function_trial(1) ==
              Approx(6743.00600619).epsilon(Tolerance));
      REQUIRE(yield_type_trial == mohr_coulomb->FailureState::Shear);
      // Initialise plastic correction components based on trial stress
      mpm::Material<Dim>::Vector6d df_dsigma_trial, dp_dsigma_trial;
      df_dsigma_trial.setZero();
      dp_dsigma_trial.setZero();
      double softening_trial = 0.;
      // Compute plastic correction components based on trial stress
      mohr_coulomb->compute_df_dp(yield_type_trial, &state_variables,
                                  trial_stress, &df_dsigma_trial,
                                  &dp_dsigma_trial, &softening_trial);

      // Check plastic correction component based on trial stress
      // Check dFtrial/dSigma
      REQUIRE(df_dsigma_trial(0) == Approx(0.05712351).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(1) == Approx(-0.05712351).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(2) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(3) == Approx(0.49672619).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(df_dsigma_trial(5) == Approx(0.).epsilon(Tolerance));
      // Check dPtrial/dSigma
      REQUIRE(dp_dsigma_trial(0) == Approx(0.08377842).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(1) == Approx(-0.01419973).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(2) == Approx(-0.06957869).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(3) == Approx(0.42599199).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(dp_dsigma_trial(5) == Approx(0.).epsilon(Tolerance));
      // Check compute stress
      mpm::Material<Dim>::Vector6d updated_stress =
          mohr_coulomb->compute_stress(stress, dstrain, particle.get(),
                                       &state_variables);
      // Check update stress
      REQUIRE(updated_stress(0) == Approx(-4530.74184152).epsilon(Tolerance));
      REQUIRE(updated_stress(1) == Approx(-5311.19503602).epsilon(Tolerance));
      REQUIRE(updated_stress(2) == Approx(-5658.06312247).epsilon(Tolerance));
      REQUIRE(updated_stress(3) == Approx(877.19525432).epsilon(Tolerance));
      REQUIRE(updated_stress(4) == Approx(0.).epsilon(Tolerance));
      REQUIRE(updated_stress(5) == Approx(0.).epsilon(Tolerance));

      // Check plastic strain
      REQUIRE(state_variables.at("epds") ==
              Approx(0.00190403).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain0") ==
              Approx(0.001114).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain1") ==
              Approx(0.00198546).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain2") ==
              Approx(0.00290055).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain3") ==
              Approx(0.00277193).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain4") ==
              Approx(0.).epsilon(Tolerance));
      REQUIRE(state_variables.at("plastic_strain5") ==
              Approx(0.).epsilon(Tolerance));
    }
  }
}

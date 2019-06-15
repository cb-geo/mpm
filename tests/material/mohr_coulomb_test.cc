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

  //! Check stress invariants and yield state
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

    // Initialise incremental of strain
    mpm::Material<Dim>::Vector6d dstrain;
    dstrain.setZero();
    dstrain(0) = -0.001;
    dstrain(1) = 0.;
    dstrain(2) = 0.;
    dstrain(3) = 0.;

    // Calculate modulus values
    const double K = material->property("youngs_modulus") /
                     (3.0 * (1. - 2. * material->property("poisson_ratio")));
    const double G = material->property("youngs_modulus") /
                     (2.0 * (1. + material->property("poisson_ratio")));
    const double a1 = K + (4.0 / 3.0) * G;
    const double a2 = K - (2.0 / 3.0) * G;
    // Compute elastic tensor
    mpm::Material<Dim>::Matrix6x6 de;
    de(0, 0) = a1;
    de(0, 1) = a2;
    de(0, 2) = a2;
    de(0, 3) = 0;
    de(0, 4) = 0;
    de(0, 5) = 0;
    de(1, 0) = a2;
    de(1, 1) = a1;
    de(1, 2) = a2;
    de(1, 3) = 0;
    de(1, 4) = 0;
    de(1, 5) = 0;
    de(2, 0) = a2;
    de(2, 1) = a2;
    de(2, 2) = a1;
    de(2, 3) = 0;
    de(2, 4) = 0;
    de(2, 5) = 0;
    de(3, 0) = 0;
    de(3, 1) = 0;
    de(3, 2) = 0;
    de(3, 3) = G;
    de(3, 4) = 0;
    de(3, 5) = 0;
    de(4, 0) = 0;
    de(4, 1) = 0;
    de(4, 2) = 0;
    de(4, 3) = 0;
    de(4, 4) = G;
    de(4, 5) = 0;
    de(5, 0) = 0;
    de(5, 1) = 0;
    de(5, 2) = 0;
    de(5, 3) = 0;
    de(5, 4) = 0;
    de(5, 5) = G;

    // Compute trial stress
    mpm::Material<Dim>::Vector6d trial_stress = stress + de * dstrain;

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

    // Check if stress invariants is computed correctly based on trial stress
    REQUIRE(mohr_coulomb->compute_stress_invariants(trial_stress,
                                                    &state_variables) == true);
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
    mpm::Material<Dim>::Vector6d updated_stress = mohr_coulomb->compute_stress(
        stress, dstrain, particle.get(), &state_variables);
    // Check update stress
    REQUIRE(updated_stress(0) == Approx(-16581.86769355).epsilon(Tolerance));
    REQUIRE(updated_stress(1) == Approx(-12936.72814065).epsilon(Tolerance));
    REQUIRE(updated_stress(2) == Approx(-13481.40416580).epsilon(Tolerance));
    REQUIRE(updated_stress(3) == Approx(-772.33801257).epsilon(Tolerance));
    REQUIRE(updated_stress(4) == Approx(-0.).epsilon(Tolerance));
    REQUIRE(updated_stress(5) == Approx(-0.).epsilon(Tolerance));
  }
}

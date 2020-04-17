#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

#include "material.h"

//! \brief Check material namespace functions
TEST_CASE("Material utility is checked", "[material]") {

  // Tolerance
  const double Tolerance = 1.E-6;

  SECTION("Check for zero stresses") {

    // Initialise stress
    Eigen::Matrix<double, 6, 1> stress;
    stress.setZero();

    // Check stress
    REQUIRE(stress(0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.).epsilon(Tolerance));

    // Compute mean p
    double p = mpm::material::p(stress);
    REQUIRE(p == Approx(0.).epsilon(Tolerance));

    // Compute deviatoric stress
    Eigen::Matrix<double, 6, 1> deviatoric_stress =
        mpm::material::deviatoric_stress(stress);
    REQUIRE(deviatoric_stress(0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(deviatoric_stress(1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(deviatoric_stress(2) == Approx(0.).epsilon(Tolerance));
    REQUIRE(deviatoric_stress(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(deviatoric_stress(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(deviatoric_stress(5) == Approx(0.).epsilon(Tolerance));

    // Compute J2
    double j2 = mpm::material::j2(stress);
    REQUIRE(j2 == Approx(0.).epsilon(Tolerance));

    // Compute J3
    double j3 = mpm::material::j3(stress);
    REQUIRE(j3 == Approx(0.).epsilon(Tolerance));

    // Compute deviatoric q
    double q = mpm::material::q(stress);
    REQUIRE(q == Approx(0.).epsilon(Tolerance));

    // Compute Lode angle theta
    double lode_angle = mpm::material::lode_angle(stress);
    REQUIRE(lode_angle == Approx(M_PI / 6.).epsilon(Tolerance));

    // Compute dp_dsigma
    Eigen::Matrix<double, 6, 1> dp_dsigma = mpm::material::dp_dsigma(stress);
    REQUIRE(dp_dsigma(0) == Approx(1. / 3.).epsilon(Tolerance));
    REQUIRE(dp_dsigma(1) == Approx(1. / 3.).epsilon(Tolerance));
    REQUIRE(dp_dsigma(2) == Approx(1. / 3.).epsilon(Tolerance));
    REQUIRE(dp_dsigma(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dp_dsigma(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dp_dsigma(5) == Approx(0.).epsilon(Tolerance));

    // Compute dq_disgma
    Eigen::Matrix<double, 6, 1> dq_dsigma = mpm::material::dq_dsigma(stress);
    REQUIRE(dq_dsigma(0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dq_dsigma(1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dq_dsigma(2) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dq_dsigma(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dq_dsigma(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dq_dsigma(5) == Approx(0.).epsilon(Tolerance));

    // Compute dj2_dsigma
    Eigen::Matrix<double, 6, 1> dj2_dsigma = mpm::material::dj2_dsigma(stress);
    REQUIRE(dj2_dsigma(0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dj2_dsigma(1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dj2_dsigma(2) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dj2_dsigma(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dj2_dsigma(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dj2_dsigma(5) == Approx(0.).epsilon(Tolerance));

    // Compute dj3_dsigma
    Eigen::Matrix<double, 6, 1> dj3_dsigma = mpm::material::dj3_dsigma(stress);
    REQUIRE(dj3_dsigma(0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dj3_dsigma(1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dj3_dsigma(2) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dj3_dsigma(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dj3_dsigma(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dj3_dsigma(5) == Approx(0.).epsilon(Tolerance));

    // Compute dtheta_dsigma
    Eigen::Matrix<double, 6, 1> dtheta_dsigma =
        mpm::material::dtheta_dsigma(stress);
    REQUIRE(dtheta_dsigma(0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dtheta_dsigma(1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dtheta_dsigma(2) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dtheta_dsigma(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dtheta_dsigma(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dtheta_dsigma(5) == Approx(0.).epsilon(Tolerance));
  }

  SECTION("Check for non-zero stresses") {

    // Initialise stress
    Eigen::Matrix<double, 6, 1> stress;
    stress(0) = -200.;
    stress(1) = -150.2;
    stress(2) = -150.2;
    stress(3) = 52.;
    stress(4) = -14.5;
    stress(5) = -33.;

    // Check stress
    REQUIRE(stress(0) == Approx(-200.).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(-150.2).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(-150.2).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(52.).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(-14.5).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(-33.).epsilon(Tolerance));

    // Compute mean p
    double p = mpm::material::p(stress);
    REQUIRE(p == Approx(166.8).epsilon(Tolerance));

    // Compute deviatoric stress
    Eigen::Matrix<double, 6, 1> deviatoric_stress =
        mpm::material::deviatoric_stress(stress);
    REQUIRE(deviatoric_stress(0) == Approx(-33.2).epsilon(Tolerance));
    REQUIRE(deviatoric_stress(1) == Approx(16.6).epsilon(Tolerance));
    REQUIRE(deviatoric_stress(2) == Approx(16.6).epsilon(Tolerance));
    REQUIRE(deviatoric_stress(3) == Approx(52.).epsilon(Tolerance));
    REQUIRE(deviatoric_stress(4) == Approx(-14.5).epsilon(Tolerance));
    REQUIRE(deviatoric_stress(5) == Approx(-33.).epsilon(Tolerance));

    // Compute J2
    double j2 = mpm::material::j2(stress);
    REQUIRE(j2 == Approx(4829.93).epsilon(Tolerance));

    // Compute J3
    double j3 = mpm::material::j3(stress);
    REQUIRE(j3 == Approx(-15368.092).epsilon(Tolerance));

    // Compute deviatoric q
    double q = mpm::material::q(stress);
    REQUIRE(q == Approx(120.3735436048968).epsilon(Tolerance));

    // Compute Lode angle theta
    double lode_angle = mpm::material::lode_angle(stress);
    REQUIRE(lode_angle == Approx(0.563342522771415).epsilon(Tolerance));

    // Compute dp_dsigma
    Eigen::Matrix<double, 6, 1> dp_dsigma = mpm::material::dp_dsigma(stress);
    REQUIRE(dp_dsigma(0) == Approx(1. / 3.).epsilon(Tolerance));
    REQUIRE(dp_dsigma(1) == Approx(1. / 3.).epsilon(Tolerance));
    REQUIRE(dp_dsigma(2) == Approx(1. / 3.).epsilon(Tolerance));
    REQUIRE(dp_dsigma(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dp_dsigma(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dp_dsigma(5) == Approx(0.).epsilon(Tolerance));

    // Compute dq_disgma
    Eigen::Matrix<double, 6, 1> dq_dsigma = mpm::material::dq_dsigma(stress);
    REQUIRE(dq_dsigma(0) == Approx(-0.413712170536900).epsilon(Tolerance));
    REQUIRE(dq_dsigma(1) == Approx(0.206856085268450).epsilon(Tolerance));
    REQUIRE(dq_dsigma(2) == Approx(0.206856085268450).epsilon(Tolerance));
    REQUIRE(dq_dsigma(3) == Approx(1.295965835416794).epsilon(Tolerance));
    REQUIRE(dq_dsigma(4) == Approx(-0.361375088721991).epsilon(Tolerance));
    REQUIRE(dq_dsigma(5) == Approx(-0.822439857091427).epsilon(Tolerance));

    // Compute dj2_dsigma
    Eigen::Matrix<double, 6, 1> dj2_dsigma = mpm::material::dj2_dsigma(stress);
    REQUIRE(dj2_dsigma(0) == Approx(-33.2).epsilon(Tolerance));
    REQUIRE(dj2_dsigma(1) == Approx(16.6).epsilon(Tolerance));
    REQUIRE(dj2_dsigma(2) == Approx(16.6).epsilon(Tolerance));
    REQUIRE(dj2_dsigma(3) == Approx(104.).epsilon(Tolerance));
    REQUIRE(dj2_dsigma(4) == Approx(-29.).epsilon(Tolerance));
    REQUIRE(dj2_dsigma(5) == Approx(-66.).epsilon(Tolerance));

    // Compute dj3_dsigma
    Eigen::Matrix<double, 6, 1> dj3_dsigma = mpm::material::dj3_dsigma(stress);
    REQUIRE(dj3_dsigma(0) == Approx(1675.28666666667).epsilon(Tolerance));
    REQUIRE(dj3_dsigma(1) == Approx(-30.1433333333333).epsilon(Tolerance));
    REQUIRE(dj3_dsigma(2) == Approx(-1645.14333333333).epsilon(Tolerance));
    REQUIRE(dj3_dsigma(3) == Approx(-769.4).epsilon(Tolerance));
    REQUIRE(dj3_dsigma(4) == Approx(-4394.8).epsilon(Tolerance));
    REQUIRE(dj3_dsigma(5) == Approx(-412.4).epsilon(Tolerance));

    // Compute dtheta_dsigma
    Eigen::Matrix<double, 6, 1> dtheta_dsigma =
        mpm::material::dtheta_dsigma(stress);
    REQUIRE(dtheta_dsigma(0) ==
            Approx(-0.00394140286438774).epsilon(Tolerance));
    REQUIRE(dtheta_dsigma(1) ==
            Approx(-0.00012754374859611).epsilon(Tolerance));
    REQUIRE(dtheta_dsigma(2) == Approx(0.00406894661298385).epsilon(Tolerance));
    REQUIRE(dtheta_dsigma(3) ==
            Approx(0.000709459104494452).epsilon(Tolerance));
    REQUIRE(dtheta_dsigma(4) == Approx(0.0117793023407611).epsilon(Tolerance));
    REQUIRE(dtheta_dsigma(5) == Approx(0.00189011673149042).epsilon(Tolerance));
  }
}
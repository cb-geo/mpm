#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

#include "material_utility.h"

//! \brief Check materials namespace functions
TEST_CASE("materials utility is checked", "[materials]") {

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
    double p = mpm::materials::p(stress);
    REQUIRE(p == Approx(0.).epsilon(Tolerance));

    // Compute deviatoric stress
    Eigen::Matrix<double, 6, 1> deviatoric_stress =
        mpm::materials::deviatoric_stress(stress);
    REQUIRE(deviatoric_stress(0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(deviatoric_stress(1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(deviatoric_stress(2) == Approx(0.).epsilon(Tolerance));
    REQUIRE(deviatoric_stress(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(deviatoric_stress(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(deviatoric_stress(5) == Approx(0.).epsilon(Tolerance));

    // Compute J2
    double j2 = mpm::materials::j2(stress);
    REQUIRE(j2 == Approx(0.).epsilon(Tolerance));

    // Compute J3
    double j3 = mpm::materials::j3(stress);
    REQUIRE(j3 == Approx(0.).epsilon(Tolerance));

    // Compute deviatoric q
    double q = mpm::materials::q(stress);
    REQUIRE(q == Approx(0.).epsilon(Tolerance));

    // Compute Lode angle theta
    double lode_angle = mpm::materials::lode_angle(stress);
    REQUIRE(lode_angle == Approx(M_PI / 6.).epsilon(Tolerance));

    // Compute Lode angle theta with tolerance
    double lode_angle_tolerance = mpm::materials::lode_angle(stress, Tolerance);
    REQUIRE(lode_angle_tolerance == Approx(M_PI / 6.).epsilon(Tolerance));

    // Compute dp_dsigma
    Eigen::Matrix<double, 6, 1> dp_dsigma = mpm::materials::dp_dsigma(stress);
    REQUIRE(dp_dsigma(0) == Approx(1. / 3.).epsilon(Tolerance));
    REQUIRE(dp_dsigma(1) == Approx(1. / 3.).epsilon(Tolerance));
    REQUIRE(dp_dsigma(2) == Approx(1. / 3.).epsilon(Tolerance));
    REQUIRE(dp_dsigma(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dp_dsigma(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dp_dsigma(5) == Approx(0.).epsilon(Tolerance));

    // Compute dq_disgma
    Eigen::Matrix<double, 6, 1> dq_dsigma = mpm::materials::dq_dsigma(stress);
    REQUIRE(dq_dsigma(0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dq_dsigma(1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dq_dsigma(2) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dq_dsigma(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dq_dsigma(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dq_dsigma(5) == Approx(0.).epsilon(Tolerance));

    // Compute dj2_dsigma
    Eigen::Matrix<double, 6, 1> dj2_dsigma = mpm::materials::dj2_dsigma(stress);
    REQUIRE(dj2_dsigma(0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dj2_dsigma(1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dj2_dsigma(2) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dj2_dsigma(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dj2_dsigma(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dj2_dsigma(5) == Approx(0.).epsilon(Tolerance));

    // Compute dj3_dsigma
    Eigen::Matrix<double, 6, 1> dj3_dsigma = mpm::materials::dj3_dsigma(stress);
    REQUIRE(dj3_dsigma(0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dj3_dsigma(1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dj3_dsigma(2) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dj3_dsigma(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dj3_dsigma(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dj3_dsigma(5) == Approx(0.).epsilon(Tolerance));

    // Compute dtheta_dsigma
    Eigen::Matrix<double, 6, 1> dtheta_dsigma =
        mpm::materials::dtheta_dsigma(stress);
    REQUIRE(dtheta_dsigma(0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dtheta_dsigma(1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dtheta_dsigma(2) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dtheta_dsigma(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dtheta_dsigma(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dtheta_dsigma(5) == Approx(0.).epsilon(Tolerance));

    // Compute dtheta_dsigma with tolerance
    Eigen::Matrix<double, 6, 1> dtheta_dsigma_tolerance =
        mpm::materials::dtheta_dsigma(stress, Tolerance);
    REQUIRE(dtheta_dsigma_tolerance(0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dtheta_dsigma_tolerance(1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dtheta_dsigma_tolerance(2) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dtheta_dsigma_tolerance(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dtheta_dsigma_tolerance(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dtheta_dsigma_tolerance(5) == Approx(0.).epsilon(Tolerance));

    // Initialise dstrain
    Eigen::Matrix<double, 6, 1> dstrain;
    dstrain.setZero();

    // Check dstrain
    REQUIRE(dstrain(0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dstrain(1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dstrain(2) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dstrain(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dstrain(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dstrain(5) == Approx(0.).epsilon(Tolerance));

    // Initialise de elastic stiffness matrix
    Eigen::Matrix<double, 6, 6> de;
    const double K = 1.0E+7;
    const double G = 6.0E+6;
    const double a1 = K + (4.0 / 3.0) * G;
    const double a2 = K - (2.0 / 3.0) * G;

    // clang-format off
    de(0,0)=a1;    de(0,1)=a2;    de(0,2)=a2;    de(0,3)=0;    de(0,4)=0;    de(0,5)=0;
    de(1,0)=a2;    de(1,1)=a1;    de(1,2)=a2;    de(1,3)=0;    de(1,4)=0;    de(1,5)=0;
    de(2,0)=a2;    de(2,1)=a2;    de(2,2)=a1;    de(2,3)=0;    de(2,4)=0;    de(2,5)=0;
    de(3,0)= 0;    de(3,1)= 0;    de(3,2)= 0;    de(3,3)=G;    de(3,4)=0;    de(3,5)=0;
    de(4,0)= 0;    de(4,1)= 0;    de(4,2)= 0;    de(4,3)=0;    de(4,4)=G;    de(4,5)=0;
    de(5,0)= 0;    de(5,1)= 0;    de(5,2)= 0;    de(5,3)=0;    de(5,4)=0;    de(5,5)=G;
    // clang-format on

    // Check dplastic strain
    Eigen::Matrix<double, 6, 1> dplastic_strain =
        mpm::materials::dplastic_strain(stress, dstrain, de);
    REQUIRE(dplastic_strain(0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dplastic_strain(1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dplastic_strain(2) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dplastic_strain(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dplastic_strain(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dplastic_strain(5) == Approx(0.).epsilon(Tolerance));

    // Check pdstrain
    double pdstrain = mpm::materials::pdstrain(dstrain);
    REQUIRE(pdstrain == Approx(0.).epsilon(Tolerance));
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
    double p = mpm::materials::p(stress);
    REQUIRE(p == Approx(-166.8).epsilon(Tolerance));

    // Compute deviatoric stress
    Eigen::Matrix<double, 6, 1> deviatoric_stress =
        mpm::materials::deviatoric_stress(stress);
    REQUIRE(deviatoric_stress(0) == Approx(-33.2).epsilon(Tolerance));
    REQUIRE(deviatoric_stress(1) == Approx(16.6).epsilon(Tolerance));
    REQUIRE(deviatoric_stress(2) == Approx(16.6).epsilon(Tolerance));
    REQUIRE(deviatoric_stress(3) == Approx(52.).epsilon(Tolerance));
    REQUIRE(deviatoric_stress(4) == Approx(-14.5).epsilon(Tolerance));
    REQUIRE(deviatoric_stress(5) == Approx(-33.).epsilon(Tolerance));

    // Compute J2
    double j2 = mpm::materials::j2(stress);
    REQUIRE(j2 == Approx(4829.93).epsilon(Tolerance));

    // Compute J3
    double j3 = mpm::materials::j3(stress);
    REQUIRE(j3 == Approx(-15368.092).epsilon(Tolerance));

    // Compute deviatoric q
    double q = mpm::materials::q(stress);
    REQUIRE(q == Approx(120.3735436048968).epsilon(Tolerance));

    // Compute Lode angle theta
    double lode_angle = mpm::materials::lode_angle(stress);
    REQUIRE(lode_angle == Approx(0.563342522771415).epsilon(Tolerance));

    // Compute Lode angle theta with tolerance
    double lode_angle_tolerance = mpm::materials::lode_angle(stress, Tolerance);
    REQUIRE(lode_angle_tolerance ==
            Approx(0.563342522771415).epsilon(Tolerance));

    // Compute dp_dsigma
    Eigen::Matrix<double, 6, 1> dp_dsigma = mpm::materials::dp_dsigma(stress);
    REQUIRE(dp_dsigma(0) == Approx(1. / 3.).epsilon(Tolerance));
    REQUIRE(dp_dsigma(1) == Approx(1. / 3.).epsilon(Tolerance));
    REQUIRE(dp_dsigma(2) == Approx(1. / 3.).epsilon(Tolerance));
    REQUIRE(dp_dsigma(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dp_dsigma(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(dp_dsigma(5) == Approx(0.).epsilon(Tolerance));

    // Compute dq_disgma
    Eigen::Matrix<double, 6, 1> dq_dsigma = mpm::materials::dq_dsigma(stress);
    REQUIRE(dq_dsigma(0) == Approx(-0.413712170536900).epsilon(Tolerance));
    REQUIRE(dq_dsigma(1) == Approx(0.206856085268450).epsilon(Tolerance));
    REQUIRE(dq_dsigma(2) == Approx(0.206856085268450).epsilon(Tolerance));
    REQUIRE(dq_dsigma(3) == Approx(1.295965835416794).epsilon(Tolerance));
    REQUIRE(dq_dsigma(4) == Approx(-0.361375088721991).epsilon(Tolerance));
    REQUIRE(dq_dsigma(5) == Approx(-0.822439857091427).epsilon(Tolerance));

    // Compute dj2_dsigma
    Eigen::Matrix<double, 6, 1> dj2_dsigma = mpm::materials::dj2_dsigma(stress);
    REQUIRE(dj2_dsigma(0) == Approx(-33.2).epsilon(Tolerance));
    REQUIRE(dj2_dsigma(1) == Approx(16.6).epsilon(Tolerance));
    REQUIRE(dj2_dsigma(2) == Approx(16.6).epsilon(Tolerance));
    REQUIRE(dj2_dsigma(3) == Approx(104.).epsilon(Tolerance));
    REQUIRE(dj2_dsigma(4) == Approx(-29.).epsilon(Tolerance));
    REQUIRE(dj2_dsigma(5) == Approx(-66.).epsilon(Tolerance));

    // Compute dj3_dsigma
    Eigen::Matrix<double, 6, 1> dj3_dsigma = mpm::materials::dj3_dsigma(stress);
    REQUIRE(dj3_dsigma(0) == Approx(1675.28666666667).epsilon(Tolerance));
    REQUIRE(dj3_dsigma(1) == Approx(-30.1433333333333).epsilon(Tolerance));
    REQUIRE(dj3_dsigma(2) == Approx(-1645.14333333333).epsilon(Tolerance));
    REQUIRE(dj3_dsigma(3) == Approx(-769.4).epsilon(Tolerance));
    REQUIRE(dj3_dsigma(4) == Approx(-4394.8).epsilon(Tolerance));
    REQUIRE(dj3_dsigma(5) == Approx(-412.4).epsilon(Tolerance));

    // Compute dtheta_dsigma
    Eigen::Matrix<double, 6, 1> dtheta_dsigma =
        mpm::materials::dtheta_dsigma(stress);
    REQUIRE(dtheta_dsigma(0) ==
            Approx(-0.00394140286438774).epsilon(Tolerance));
    REQUIRE(dtheta_dsigma(1) ==
            Approx(-0.00012754374859611).epsilon(Tolerance));
    REQUIRE(dtheta_dsigma(2) == Approx(0.00406894661298385).epsilon(Tolerance));
    REQUIRE(dtheta_dsigma(3) ==
            Approx(0.000709459104494452).epsilon(Tolerance));
    REQUIRE(dtheta_dsigma(4) == Approx(0.0117793023407611).epsilon(Tolerance));
    REQUIRE(dtheta_dsigma(5) == Approx(0.00189011673149042).epsilon(Tolerance));

    // Compute dtheta_dsigma with tolerance
    Eigen::Matrix<double, 6, 1> dtheta_dsigma_tolerance =
        mpm::materials::dtheta_dsigma(stress, Tolerance);
    REQUIRE(dtheta_dsigma_tolerance(0) ==
            Approx(-0.00394140286438774).epsilon(Tolerance));
    REQUIRE(dtheta_dsigma_tolerance(1) ==
            Approx(-0.00012754374859611).epsilon(Tolerance));
    REQUIRE(dtheta_dsigma_tolerance(2) ==
            Approx(0.00406894661298385).epsilon(Tolerance));
    REQUIRE(dtheta_dsigma_tolerance(3) ==
            Approx(0.000709459104494452).epsilon(Tolerance));
    REQUIRE(dtheta_dsigma_tolerance(4) ==
            Approx(0.0117793023407611).epsilon(Tolerance));
    REQUIRE(dtheta_dsigma_tolerance(5) ==
            Approx(0.00189011673149042).epsilon(Tolerance));

    // Initialise dstrain
    Eigen::Matrix<double, 6, 1> dstrain;
    dstrain(0) = 0.001;
    dstrain(1) = -0.002;
    dstrain(2) = -0.0015;
    dstrain(3) = 0.0007;
    dstrain(4) = -0.0006;
    dstrain(5) = 0.0009;

    // Check dstrain
    REQUIRE(dstrain(0) == Approx(0.001).epsilon(Tolerance));
    REQUIRE(dstrain(1) == Approx(-0.002).epsilon(Tolerance));
    REQUIRE(dstrain(2) == Approx(-0.0015).epsilon(Tolerance));
    REQUIRE(dstrain(3) == Approx(0.0007).epsilon(Tolerance));
    REQUIRE(dstrain(4) == Approx(-0.0006).epsilon(Tolerance));
    REQUIRE(dstrain(5) == Approx(0.0009).epsilon(Tolerance));

    // Initialise de elastic stiffness matrix
    Eigen::Matrix<double, 6, 6> de;
    const double K = 1.0E+7;
    const double G = 6.0E+6;
    const double a1 = K + (4.0 / 3.0) * G;
    const double a2 = K - (2.0 / 3.0) * G;

    // clang-format off
    de(0,0)=a1;    de(0,1)=a2;    de(0,2)=a2;    de(0,3)=0;    de(0,4)=0;    de(0,5)=0;
    de(1,0)=a2;    de(1,1)=a1;    de(1,2)=a2;    de(1,3)=0;    de(1,4)=0;    de(1,5)=0;
    de(2,0)=a2;    de(2,1)=a2;    de(2,2)=a1;    de(2,3)=0;    de(2,4)=0;    de(2,5)=0;
    de(3,0)= 0;    de(3,1)= 0;    de(3,2)= 0;    de(3,3)=G;    de(3,4)=0;    de(3,5)=0;
    de(4,0)= 0;    de(4,1)= 0;    de(4,2)= 0;    de(4,3)=0;    de(4,4)=G;    de(4,5)=0;
    de(5,0)= 0;    de(5,1)= 0;    de(5,2)= 0;    de(5,3)=0;    de(5,4)=0;    de(5,5)=G;
    // clang-format on

    // Check dplastic strain
    Eigen::Matrix<double, 6, 1> dplastic_strain =
        mpm::materials::dplastic_strain(stress, dstrain, de);
    REQUIRE(dplastic_strain(0) == Approx(0.001008326666667).epsilon(Tolerance));
    REQUIRE(dplastic_strain(1) ==
            Approx(-0.001995823333333).epsilon(Tolerance));
    REQUIRE(dplastic_strain(2) ==
            Approx(-0.001495823333333).epsilon(Tolerance));
    REQUIRE(dplastic_strain(3) == Approx(0.000691333333333).epsilon(Tolerance));
    REQUIRE(dplastic_strain(4) ==
            Approx(-0.000597583333333).epsilon(Tolerance));
    REQUIRE(dplastic_strain(5) == Approx(0.000905500000000).epsilon(Tolerance));

    // Check pdstrain
    double pdstrain = mpm::materials::pdstrain(dstrain);
    REQUIRE(pdstrain == Approx(0.001999444367263).epsilon(Tolerance));
  }
}
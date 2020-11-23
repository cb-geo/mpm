#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

#include "radial_basis_function.h"

TEST_CASE("Radial basis function 2D", "[RBF][2D]") {
  // Dimension
  const unsigned Dim = 2;
  // Tolerance
  const double Tolerance = 1.E-7;

  SECTION("Kernel 2D") {
    // Initialize variables
    double smoothing_length, norm_distance;
    std::string type;

    // Check error: wrong type name
    type = "test_error";
    smoothing_length = 1.0;
    norm_distance = 1.0;
    REQUIRE_THROWS(mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                         norm_distance, type));

    // Check error: wrong dimension
    REQUIRE_THROWS(mpm::RadialBasisFunction::kernel<1>(
        smoothing_length, norm_distance, "cubic_spline"));
    REQUIRE_THROWS(mpm::RadialBasisFunction::kernel<1>(
        smoothing_length, norm_distance, "quintic_spline"));
    REQUIRE_THROWS(mpm::RadialBasisFunction::kernel<1>(
        smoothing_length, norm_distance, "gaussian"));
    REQUIRE_THROWS(mpm::RadialBasisFunction::kernel<1>(
        smoothing_length, norm_distance, "super_gaussian"));

    // Cubic Spline
    type = "cubic_spline";
    norm_distance = 3.0;
    double kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                          norm_distance, type);
    REQUIRE(kernel == Approx(0.0).epsilon(Tolerance));

    norm_distance = 2.0;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.0).epsilon(Tolerance));

    norm_distance = 1.5;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.014210262776062).epsilon(Tolerance));

    norm_distance = 1.0;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.113682102208497).epsilon(Tolerance));

    norm_distance = 0.5;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.326836043849428).epsilon(Tolerance));

    norm_distance = 0.0;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.454728408833987).epsilon(Tolerance));

    // Quintic Spline
    type = "quintic_spline";
    norm_distance = 5.0;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.0).epsilon(Tolerance));

    norm_distance = 3.0;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.0).epsilon(Tolerance));

    norm_distance = 2.5;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(2.08100082494633E-05).epsilon(Tolerance));

    norm_distance = 2.0;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.000665920263983).epsilon(Tolerance));

    norm_distance = 1.5;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.004931971955123).epsilon(Tolerance));

    norm_distance = 1.0;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.017313926863554).epsilon(Tolerance));

    norm_distance = 0.5;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.035002433875597).epsilon(Tolerance));

    norm_distance = 0.0;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.043950737422867).epsilon(Tolerance));

    // Gaussian
    type = "gaussian";
    norm_distance = 3.0 + Tolerance;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.0).epsilon(Tolerance));

    norm_distance = 3.0;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(3.92825606927949E-05).epsilon(Tolerance));

    norm_distance = 1.5;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.033549615174147).epsilon(Tolerance));

    norm_distance = 0.0;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.318309886183791).epsilon(Tolerance));

    // Super Gaussian
    type = "super_gaussian";
    norm_distance = 3.0 + Tolerance;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.0).epsilon(Tolerance));

    norm_distance = 3.0;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(-0.00027497792485).epsilon(Tolerance));

    norm_distance = 1.5;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(-0.008387403793537).epsilon(Tolerance));

    norm_distance = 0.0;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.636619772367581).epsilon(Tolerance));
  }

  SECTION("Gradient 2D") {
    // Initialize variables
    double smoothing_length;
    Eigen::Matrix<double, Dim, 1> relative_distance;
    std::string type;

    // Check error: wrong type name
    type = "test_error";
    smoothing_length = 1.0;
    relative_distance.setZero();
    REQUIRE_THROWS(mpm::RadialBasisFunction::gradient<Dim>(
        smoothing_length, relative_distance, type));

    // Check error: wrong dimension
    Eigen::Matrix<double, 1, 1> error_mat;
    error_mat.setZero();
    REQUIRE_THROWS(mpm::RadialBasisFunction::gradient<1>(
        smoothing_length, error_mat, "cubic_spline"));
    REQUIRE_THROWS(mpm::RadialBasisFunction::gradient<1>(
        smoothing_length, error_mat, "quintic_spline"));
    REQUIRE_THROWS(mpm::RadialBasisFunction::gradient<1>(
        smoothing_length, error_mat, "gaussian"));
    REQUIRE_THROWS(mpm::RadialBasisFunction::gradient<1>(
        smoothing_length, error_mat, "super_gaussian"));

    // Cubic Spline
    type = "cubic_spline";
    relative_distance.setZero();
    Eigen::Matrix<double, Dim, 1> gradient =
        mpm::RadialBasisFunction::gradient<Dim>(smoothing_length,
                                                relative_distance, type);
    for (unsigned i = 0; i < gradient.size(); ++i)
      REQUIRE(gradient(i) == Approx(0.0).epsilon(Tolerance));

    relative_distance << 0.5, -0.5;
    gradient = mpm::RadialBasisFunction::gradient<Dim>(smoothing_length,
                                                       relative_distance, type);

    REQUIRE(gradient(0) == Approx(-0.320358379080714).epsilon(Tolerance));
    REQUIRE(gradient(1) == Approx(0.320358379080714).epsilon(Tolerance));

    relative_distance << 1.0, -0.5;
    gradient = mpm::RadialBasisFunction::gradient<Dim>(smoothing_length,
                                                       relative_distance, type);

    REQUIRE(gradient(0) == Approx(-0.237280496186688).epsilon(Tolerance));
    REQUIRE(gradient(1) == Approx(0.118640248093344).epsilon(Tolerance));

    relative_distance << 1.5, 1.5;
    gradient = mpm::RadialBasisFunction::gradient<Dim>(smoothing_length,
                                                       relative_distance, type);
    for (unsigned i = 0; i < gradient.size(); ++i)
      REQUIRE(gradient(i) == Approx(0.0).epsilon(Tolerance));

    // Quintic Spline
    type = "quintic_spline";
    relative_distance.setZero();
    gradient = mpm::RadialBasisFunction::gradient<Dim>(smoothing_length,
                                                       relative_distance, type);
    for (unsigned i = 0; i < gradient.size(); ++i)
      REQUIRE(gradient(i) == Approx(0.0).epsilon(Tolerance));

    relative_distance << 0.5, -0.5;
    gradient = mpm::RadialBasisFunction::gradient<Dim>(smoothing_length,
                                                       relative_distance, type);

    REQUIRE(gradient(0) == Approx(-0.025863567099382).epsilon(Tolerance));
    REQUIRE(gradient(1) == Approx(0.025863567099382).epsilon(Tolerance));

    relative_distance << 1.0, -0.5;
    gradient = mpm::RadialBasisFunction::gradient<Dim>(smoothing_length,
                                                       relative_distance, type);

    REQUIRE(gradient(0) == Approx(-0.026546314385505).epsilon(Tolerance));
    REQUIRE(gradient(1) == Approx(0.013273157192753).epsilon(Tolerance));

    relative_distance << -1.9, 1.3;
    gradient = mpm::RadialBasisFunction::gradient<Dim>(smoothing_length,
                                                       relative_distance, type);

    REQUIRE(gradient(0) == Approx(0.000651627282552).epsilon(Tolerance));
    REQUIRE(gradient(1) == Approx(-0.000445850245957).epsilon(Tolerance));

    relative_distance << 2.5, 3.5;
    gradient = mpm::RadialBasisFunction::gradient<Dim>(smoothing_length,
                                                       relative_distance, type);
    for (unsigned i = 0; i < gradient.size(); ++i)
      REQUIRE(gradient(i) == Approx(0.0).epsilon(Tolerance));

    // Gaussian
    type = "gaussian";
    relative_distance.setZero();
    gradient = mpm::RadialBasisFunction::gradient<Dim>(smoothing_length,
                                                       relative_distance, type);
    for (unsigned i = 0; i < gradient.size(); ++i)
      REQUIRE(gradient(i) == Approx(0.0).epsilon(Tolerance));

    relative_distance << -1.9, 1.3;
    gradient = mpm::RadialBasisFunction::gradient<Dim>(smoothing_length,
                                                       relative_distance, type);

    REQUIRE(gradient(0) == Approx(0.00603772001586).epsilon(Tolerance));
    REQUIRE(gradient(1) == Approx(-0.004131071589799).epsilon(Tolerance));

    relative_distance << 2.5, 3.5;
    gradient = mpm::RadialBasisFunction::gradient<Dim>(smoothing_length,
                                                       relative_distance, type);
    for (unsigned i = 0; i < gradient.size(); ++i)
      REQUIRE(gradient(i) == Approx(0.0).epsilon(Tolerance));

    // Super Gaussian
    type = "super_gaussian";
    relative_distance.setZero();
    gradient = mpm::RadialBasisFunction::gradient<Dim>(smoothing_length,
                                                       relative_distance, type);
    for (unsigned i = 0; i < gradient.size(); ++i)
      REQUIRE(gradient(i) == Approx(0.0).epsilon(Tolerance));

    relative_distance << -1.9, 1.3;
    gradient = mpm::RadialBasisFunction::gradient<Dim>(smoothing_length,
                                                       relative_distance, type);

    REQUIRE(gradient(0) == Approx(-0.013886756036479).epsilon(Tolerance));
    REQUIRE(gradient(1) == Approx(0.009501464656538).epsilon(Tolerance));

    relative_distance << 2.5, 3.5;
    gradient = mpm::RadialBasisFunction::gradient<Dim>(smoothing_length,
                                                       relative_distance, type);
    for (unsigned i = 0; i < gradient.size(); ++i)
      REQUIRE(gradient(i) == Approx(0.0).epsilon(Tolerance));
  }
}

TEST_CASE("Radial basis function 3D", "[RBF][3D]") {
  // Dimension
  const unsigned Dim = 3;
  // Tolerance
  const double Tolerance = 1.E-9;

  SECTION("Kernel 3D") {
    // Initialize variables
    double smoothing_length, norm_distance;
    std::string type;

    // Cubic Spline
    type = "cubic_spline";
    smoothing_length = 1.0;
    norm_distance = 3.0;
    double kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                          norm_distance, type);
    REQUIRE(kernel == Approx(0.0).epsilon(Tolerance));

    norm_distance = 2.0;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.0).epsilon(Tolerance));

    norm_distance = 1.5;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.009947183943243).epsilon(Tolerance));

    norm_distance = 1.0;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.079577471545948).epsilon(Tolerance));

    norm_distance = 0.5;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.2287852306946).epsilon(Tolerance));

    norm_distance = 0.0;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.318309886183791).epsilon(Tolerance));

    // Quintic Spline
    type = "quintic_spline";
    norm_distance = 5.0;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.0).epsilon(Tolerance));

    norm_distance = 3.0;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.0).epsilon(Tolerance));

    norm_distance = 2.5;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(8.31240998042629E-05).epsilon(Tolerance));

    norm_distance = 2.0;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.002659971193736).epsilon(Tolerance));

    norm_distance = 1.5;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.01970041165361).epsilon(Tolerance));

    norm_distance = 1.0;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.069159251037147).epsilon(Tolerance));

    norm_distance = 0.5;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.13981473587077).epsilon(Tolerance));

    norm_distance = 0.0;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.175558098786603).epsilon(Tolerance));

    // Gaussian
    type = "gaussian";
    norm_distance = 3.0 + Tolerance;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.0).epsilon(Tolerance));

    norm_distance = 3.0;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(2.21628115579574E-05).epsilon(Tolerance));

    norm_distance = 1.5;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.018928343413289).epsilon(Tolerance));

    norm_distance = 0.0;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.179587122125167).epsilon(Tolerance));

    // Super Gaussian
    type = "super_gaussian";
    norm_distance = 3.0 + Tolerance;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.0).epsilon(Tolerance));

    norm_distance = 3.0;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(-0.000144058275127).epsilon(Tolerance));

    norm_distance = 1.5;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.004732085853322).epsilon(Tolerance));

    norm_distance = 0.0;
    kernel = mpm::RadialBasisFunction::kernel<Dim>(smoothing_length,
                                                   norm_distance, type);
    REQUIRE(kernel == Approx(0.448967805312916).epsilon(Tolerance));
  }

  SECTION("Gradient 3D") {
    // Initialize variables
    double smoothing_length = 1.0;
    Eigen::Matrix<double, Dim, 1> relative_distance;
    std::string type;

    // Cubic Spline
    type = "cubic_spline";
    relative_distance.setZero();
    Eigen::Matrix<double, Dim, 1> gradient =
        mpm::RadialBasisFunction::gradient<Dim>(smoothing_length,
                                                relative_distance, type);
    for (unsigned i = 0; i < gradient.size(); ++i)
      REQUIRE(gradient(i) == Approx(0.0).epsilon(Tolerance));

    relative_distance << 0.5, -0.4, 0.1;
    gradient = mpm::RadialBasisFunction::gradient<Dim>(smoothing_length,
                                                       relative_distance, type);

    REQUIRE(gradient(0) == Approx(-0.245390397939789).epsilon(Tolerance));
    REQUIRE(gradient(1) == Approx(0.196312318351831).epsilon(Tolerance));
    REQUIRE(gradient(2) == Approx(-0.049078079587958).epsilon(Tolerance));

    relative_distance << 1.0, -0.6, 0.3;
    gradient = mpm::RadialBasisFunction::gradient<Dim>(smoothing_length,
                                                       relative_distance, type);

    REQUIRE(gradient(0) == Approx(-0.1255681536468).epsilon(Tolerance));
    REQUIRE(gradient(1) == Approx(0.07534089218808).epsilon(Tolerance));
    REQUIRE(gradient(2) == Approx(-0.03767044609404).epsilon(Tolerance));

    relative_distance << 1.5, 1.5, 2.0;
    gradient = mpm::RadialBasisFunction::gradient<Dim>(smoothing_length,
                                                       relative_distance, type);
    for (unsigned i = 0; i < gradient.size(); ++i)
      REQUIRE(gradient(i) == Approx(0.0).epsilon(Tolerance));

    // Quintic Spline
    type = "quintic_spline";
    relative_distance.setZero();
    gradient = mpm::RadialBasisFunction::gradient<Dim>(smoothing_length,
                                                       relative_distance, type);
    for (unsigned i = 0; i < gradient.size(); ++i)
      REQUIRE(gradient(i) == Approx(0.0).epsilon(Tolerance));

    relative_distance << 0.5, -0.5, 0.2;
    gradient = mpm::RadialBasisFunction::gradient<Dim>(smoothing_length,
                                                       relative_distance, type);

    REQUIRE(gradient(0) == Approx(-0.099803272175507).epsilon(Tolerance));
    REQUIRE(gradient(1) == Approx(0.099803272175507).epsilon(Tolerance));
    REQUIRE(gradient(2) == Approx(-0.039921308870203).epsilon(Tolerance));

    relative_distance << 0.5, -0.5, 1.3;
    gradient = mpm::RadialBasisFunction::gradient<Dim>(smoothing_length,
                                                       relative_distance, type);

    REQUIRE(gradient(0) == Approx(-0.022021780742007).epsilon(Tolerance));
    REQUIRE(gradient(1) == Approx(0.022021780742007).epsilon(Tolerance));
    REQUIRE(gradient(2) == Approx(-0.057256629929218).epsilon(Tolerance));

    relative_distance << -1.9, 1.3, 1.8;
    gradient = mpm::RadialBasisFunction::gradient<Dim>(smoothing_length,
                                                       relative_distance, type);

    REQUIRE(gradient(0) == Approx(3.14726383394464E-07).epsilon(Tolerance));
    REQUIRE(gradient(1) == Approx(-2.15339104427791E-07).epsilon(Tolerance));
    REQUIRE(gradient(2) == Approx(-2.98161836900018E-07).epsilon(Tolerance));

    relative_distance << 2.5, 3.5, 3.5;
    gradient = mpm::RadialBasisFunction::gradient<Dim>(smoothing_length,
                                                       relative_distance, type);
    for (unsigned i = 0; i < gradient.size(); ++i)
      REQUIRE(gradient(i) == Approx(0.0).epsilon(Tolerance));

    // Gaussian
    type = "gaussian";
    relative_distance.setZero();
    gradient = mpm::RadialBasisFunction::gradient<Dim>(smoothing_length,
                                                       relative_distance, type);
    for (unsigned i = 0; i < gradient.size(); ++i)
      REQUIRE(gradient(i) == Approx(0.0).epsilon(Tolerance));

    relative_distance << -1.9, 1.3, 1.0;
    gradient = mpm::RadialBasisFunction::gradient<Dim>(smoothing_length,
                                                       relative_distance, type);

    REQUIRE(gradient(0) == Approx(0.001253151422955).epsilon(Tolerance));
    REQUIRE(gradient(1) == Approx(-0.000857419394653).epsilon(Tolerance));
    REQUIRE(gradient(2) == Approx(-0.000659553380503).epsilon(Tolerance));

    relative_distance << 2.5, 3.5, 3.5;
    gradient = mpm::RadialBasisFunction::gradient<Dim>(smoothing_length,
                                                       relative_distance, type);
    for (unsigned i = 0; i < gradient.size(); ++i)
      REQUIRE(gradient(i) == Approx(0.0).epsilon(Tolerance));

    // Super Gaussian
    type = "super_gaussian";
    relative_distance.setZero();
    gradient = mpm::RadialBasisFunction::gradient<Dim>(smoothing_length,
                                                       relative_distance, type);
    for (unsigned i = 0; i < gradient.size(); ++i)
      REQUIRE(gradient(i) == Approx(0.0).epsilon(Tolerance));

    relative_distance << -1.9, 1.3, 1.0;
    gradient = mpm::RadialBasisFunction::gradient<Dim>(smoothing_length,
                                                       relative_distance, type);

    REQUIRE(gradient(0) == Approx(-0.003508823984274).epsilon(Tolerance));
    REQUIRE(gradient(1) == Approx(0.00240077430503).epsilon(Tolerance));
    REQUIRE(gradient(2) == Approx(0.001846749465407).epsilon(Tolerance));

    relative_distance << 2.5, 3.5, 3.5;
    gradient = mpm::RadialBasisFunction::gradient<Dim>(smoothing_length,
                                                       relative_distance, type);
    for (unsigned i = 0; i < gradient.size(); ++i)
      REQUIRE(gradient(i) == Approx(0.0).epsilon(Tolerance));
  }
}

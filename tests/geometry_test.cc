#include <cmath>
#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

#include "geometry.h"

//! \brief Check geometry class for 2D case
TEST_CASE("Geometry is checked for 2D case", "[geometry][2D]") {

  // Dimension
  const unsigned Dim = 2;
  // Tolerance
  const double Tolerance = 1.E-7;

  SECTION("Check inverse rotation matrix") {
    Eigen::Matrix<double, 2, 1> angles;
    // clang-format off
    angles << 10. * M_PI / 180.,   // alpha
              30. * M_PI / 180.;   // beta
    // clang-format on

    Eigen::Matrix<double, 2, 2> rotation_matrix;
    // clang-format off
    rotation_matrix <<  0.766044443118978, -0.642787609686539,
                        0.642787609686539,  0.766044443118978;
    // clang-format on
    auto check_rotation_matrix = mpm::geometry::rotation_matrix(angles);
    REQUIRE(check_rotation_matrix.cols() == 2);
    REQUIRE(check_rotation_matrix.rows() == 2);
    for (unsigned i = 0; i < check_rotation_matrix.rows(); ++i) {
      for (unsigned j = 0; j < check_rotation_matrix.cols(); ++j) {
        REQUIRE(check_rotation_matrix(i, j) ==
                Approx(rotation_matrix(i, j)).epsilon(Tolerance));
      }
    }
  }

  SECTION("Check angle between two vectors") {
    Eigen::Matrix<double, 2, 1> vector_a;
    vector_a << 3., 0.;
    Eigen::Matrix<double, 2, 1> vector_b;
    vector_b << -2., 2.;

    const double angle = 2.356194490192345;

    REQUIRE(mpm::geometry::angle_between_vectors(vector_a, vector_b) ==
            Approx(angle).epsilon(Tolerance));
  }

  SECTION("Check euler angles computations") {
    Eigen::Matrix<double, Dim, Dim> new_axes;
    // clang-format off
    new_axes << 2., -2.,
                2.,  2.;
    // clang-format on

    Eigen::Matrix<double, Dim, 1> check_euler_angles;
    check_euler_angles << 0.78539816339, 0.;

    for (unsigned i = 0; i < Dim; ++i) {
      REQUIRE(mpm::geometry::euler_angles_cartesian(new_axes)(i) ==
              Approx(check_euler_angles(i)).epsilon(Tolerance));
    }
  }
}

//! \brief Check cell class for 3D case
TEST_CASE("Geometry is checked for 3D case", "[geometry][3D]") {

  // Dimension
  const unsigned Dim = 3;
  // Tolerance
  const double Tolerance = 1.E-7;

  SECTION("Check inverse rotation matrix") {
    Eigen::Matrix<double, 3, 1> angles;
    // clang-format off
    angles << 10. * M_PI / 180.,   // alpha 
              30. * M_PI / 180.,   // beta
              60. * M_PI / 180.;   // gamma
    // clang-format on

    Eigen::Matrix<double, 3, 3> rotation_matrix;
    // clang-format off
    rotation_matrix <<  0.809456487535711,  -0.567595743096322,  0.150383733180435,
                        0.396585671433487,   0.339610177142757, -0.852868531952443,
                        0.433012701892219,   0.75,                  0.5;
    // clang-format on
    const auto check_rotation_matrix = mpm::geometry::rotation_matrix(angles);
    REQUIRE(check_rotation_matrix.cols() == 3);
    REQUIRE(check_rotation_matrix.rows() == 3);
    for (unsigned i = 0; i < check_rotation_matrix.rows(); ++i) {
      for (unsigned j = 0; j < check_rotation_matrix.cols(); ++j) {
        REQUIRE(check_rotation_matrix(i, j) ==
                Approx(rotation_matrix(i, j)).epsilon(Tolerance));
      }
    }
  }

  SECTION("Check angle between two vectors") {
    Eigen::Matrix<double, 3, 1> vector_a;
    vector_a << 3., 0., 4.;
    Eigen::Matrix<double, 3, 1> vector_b;
    vector_b << -2., 2., 4.;

    const double angle = 1.150261991510931;

    REQUIRE(mpm::geometry::angle_between_vectors(vector_a, vector_b) ==
            Approx(angle).epsilon(Tolerance));
  }

  SECTION("Check euler angles computations") {

    Eigen::Matrix<double, Dim, Dim> new_axes;
    // clang-format off
    new_axes << 2., -2., -2.,
                0.,  0., -2.,
                2.,  2.,  0.;
    // clang-format on

    Eigen::Matrix<double, Dim, 1> check_euler_angles;
    check_euler_angles << 0.78539816339, 1.0471975512, 1.57079632679;

    for (unsigned i = 0; i < Dim; ++i) {
      REQUIRE(mpm::geometry::euler_angles_cartesian(new_axes)(i) ==
              Approx(check_euler_angles(i)).epsilon(Tolerance));
    }
  }
}

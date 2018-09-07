// Hexahedron element test
#include <memory>

#include "catch.hpp"

#include "hexahedron_element.h"

//! \brief Check hexahedron element class
TEST_CASE("Hexahedron elements are checked", "[hex][element][3D]") {
  const unsigned Dim = 3;
  const double Tolerance = 1.E-7;

  //! Check for 8 noded element
  SECTION("Hexahedron element with eight nodes") {
    const unsigned nfunctions = 8;
    std::shared_ptr<mpm::Element<Dim>> hex =
        std::make_shared<mpm::HexahedronElement<Dim, nfunctions>>();

    // Check degree
    REQUIRE(hex->degree() == mpm::ElementDegree::Linear);

    // Coordinates is (0, 0, 0)
    SECTION("Eight noded hexahedron element for coordinates(0, 0, 0)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();
      auto shapefn = hex->shapefn(coords);

      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      REQUIRE(shapefn(0) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(shapefn(3) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(shapefn(4) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(shapefn(5) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(shapefn(6) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(shapefn(7) == Approx(0.125).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = hex->grad_shapefn(coords);
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      REQUIRE(gradsf(0, 0) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(3, 0) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(4, 0) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(5, 0) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(6, 0) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(7, 0) == Approx(-0.125).epsilon(Tolerance));

      REQUIRE(gradsf(0, 1) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(3, 1) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(4, 1) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(5, 1) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(6, 1) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(7, 1) == Approx(0.125).epsilon(Tolerance));

      REQUIRE(gradsf(0, 2) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(1, 2) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(2, 2) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(3, 2) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(4, 2) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(5, 2) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(6, 2) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(7, 2) == Approx(0.125).epsilon(Tolerance));
    }

    // Coordinates is (-1, -1, -1);
    SECTION("Eight noded hexahedron element for coordinates(-1, -1, -1)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << -1., -1., -1.;
      auto shapefn = hex->shapefn(coords);
      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      REQUIRE(shapefn(0) == Approx(1.0).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(3) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(4) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(5) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(6) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(7) == Approx(0.0).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = hex->grad_shapefn(coords);
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      REQUIRE(gradsf(0, 0) == Approx(-0.5).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(0.5).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(3, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(4, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(5, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(6, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(7, 0) == Approx(0.0).epsilon(Tolerance));

      REQUIRE(gradsf(0, 1) == Approx(-0.5).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(3, 1) == Approx(0.5).epsilon(Tolerance));
      REQUIRE(gradsf(4, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(5, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(6, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(7, 1) == Approx(0.0).epsilon(Tolerance));

      REQUIRE(gradsf(0, 2) == Approx(-0.5).epsilon(Tolerance));
      REQUIRE(gradsf(1, 2) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 2) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(3, 2) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(4, 2) == Approx(0.5).epsilon(Tolerance));
      REQUIRE(gradsf(5, 2) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(6, 2) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(7, 2) == Approx(0.0).epsilon(Tolerance));
    }

    // Coordinates is (1, 1, 1)
    SECTION("Eight noded hexahedron element for coordinates(1, 1, 1)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << 1., 1., 1.;
      auto shapefn = hex->shapefn(coords);

      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      REQUIRE(shapefn(0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(3) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(4) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(5) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(6) == Approx(1.0).epsilon(Tolerance));
      REQUIRE(shapefn(7) == Approx(0.0).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = hex->grad_shapefn(coords);
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      REQUIRE(gradsf(0, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(3, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(4, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(5, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(6, 0) == Approx(0.5).epsilon(Tolerance));
      REQUIRE(gradsf(7, 0) == Approx(-0.5).epsilon(Tolerance));

      REQUIRE(gradsf(0, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(3, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(4, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(5, 1) == Approx(-0.5).epsilon(Tolerance));
      REQUIRE(gradsf(6, 1) == Approx(0.5).epsilon(Tolerance));
      REQUIRE(gradsf(7, 1) == Approx(0.0).epsilon(Tolerance));

      REQUIRE(gradsf(0, 2) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(1, 2) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 2) == Approx(-0.5).epsilon(Tolerance));
      REQUIRE(gradsf(3, 2) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(4, 2) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(5, 2) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(6, 2) == Approx(0.5).epsilon(Tolerance));
      REQUIRE(gradsf(7, 2) == Approx(0.0).epsilon(Tolerance));
    }

    SECTION("Eight noded hexahedron shapefn with deformation gradient") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();
      Eigen::Matrix<double, Dim, 1> psize;
      psize.setZero();
      Eigen::Matrix<double, Dim, 1> defgrad;
      defgrad.setZero();
      auto shapefn = hex->shapefn(coords, psize, defgrad);

      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      REQUIRE(shapefn(0) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(shapefn(3) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(shapefn(4) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(shapefn(5) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(shapefn(6) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(shapefn(7) == Approx(0.125).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = hex->grad_shapefn(coords, psize, defgrad);
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      REQUIRE(gradsf(0, 0) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(3, 0) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(4, 0) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(5, 0) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(6, 0) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(7, 0) == Approx(-0.125).epsilon(Tolerance));

      REQUIRE(gradsf(0, 1) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(3, 1) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(4, 1) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(5, 1) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(6, 1) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(7, 1) == Approx(0.125).epsilon(Tolerance));

      REQUIRE(gradsf(0, 2) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(1, 2) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(2, 2) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(3, 2) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(4, 2) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(5, 2) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(6, 2) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(7, 2) == Approx(0.125).epsilon(Tolerance));
    }

    // Check Jacobian
    SECTION(
        "Eight noded hexahedron Jacobian for local coordinates(0.5,0.5,0.5)") {
      Eigen::Matrix<double, 8, Dim> coords;
      // clang-format off
      coords << 2., 1., 0.5,
                4., 2., 1.0,
                2., 4., 1.0,
                1., 3., 0.5,
                2., 1., 1.5,
                4., 2., 2.0,
                2., 4., 2.0,
                1., 3., 1.5;
      // clang-format on

      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0.5, 0.5, 0.5;

      // Jacobian result
      Eigen::Matrix<double, Dim, Dim> jacobian;
      // clang-format off
      jacobian << 0.625, 0.5, 0.25,
                 -0.875, 1.0, 0.00,
                  0.000, 0.0, 0.50;
      // clang-format on

      // Get Jacobian
      auto jac = hex->jacobian(xi, coords);

      // Check size of jacobian
      REQUIRE(jac.size() == jacobian.size());

      // Check Jacobian
      for (unsigned i = 0; i < Dim; ++i)
        for (unsigned j = 0; j < Dim; ++j)
          REQUIRE(jac(i, j) == Approx(jacobian(i, j)).epsilon(Tolerance));
    }

    // Check Jacobian
    SECTION("Eight noded hexahedron Jacobian with deformation gradient") {
      Eigen::Matrix<double, 8, Dim> coords;
      // clang-format off
      coords << 2., 1., 0.5,
                4., 2., 1.0,
                2., 4., 1.0,
                1., 3., 0.5,
                2., 1., 1.5,
                4., 2., 2.0,
                2., 4., 2.0,
                1., 3., 1.5;
      // clang-format on

      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0.5, 0.5, 0.5;
      Eigen::Matrix<double, Dim, 1> psize;
      psize.setZero();
      Eigen::Matrix<double, Dim, 1> defgrad;
      defgrad.setZero();

      // Jacobian result
      Eigen::Matrix<double, Dim, Dim> jacobian;
      // clang-format off
      jacobian << 0.625, 0.5, 0.25,
                 -0.875, 1.0, 0.00,
                  0.000, 0.0, 0.50;
      // clang-format on

      // Get Jacobian
      auto jac = hex->jacobian(xi, coords, psize, defgrad);

      // Check size of jacobian
      REQUIRE(jac.size() == jacobian.size());

      // Check Jacobian
      for (unsigned i = 0; i < Dim; ++i)
        for (unsigned j = 0; j < Dim; ++j)
          REQUIRE(jac(i, j) == Approx(jacobian(i, j)).epsilon(Tolerance));
    }

    // Coordinates is (0, 0, 0)
    SECTION("Eight noded hexahedron B-matrix for coordinates(0, 0, 0)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();
      // Get B-Matrix
      auto bmatrix = hex->bmatrix(coords);

      // Check gradient of shape functions
      auto gradsf = hex->grad_shapefn(coords);

      // Check size of B-matrix
      REQUIRE(bmatrix.size() == nfunctions);

      for (unsigned i = 0; i < nfunctions; ++i) {
        REQUIRE(bmatrix.at(i)(0, 0) == Approx(gradsf(i, 0)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(0, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(0, 2) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 1) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 2) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 2) == Approx(gradsf(i, 2)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(3, 0) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(3, 1) == Approx(gradsf(i, 0)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(3, 2) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(4, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(4, 1) == Approx(gradsf(i, 2)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(4, 2) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(5, 0) == Approx(gradsf(i, 2)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(5, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(5, 2) == Approx(gradsf(i, 0)).epsilon(Tolerance));
      }
    }

    // Coordinates is (0.5, 0.5, 0.5)
    SECTION("Eight noded hexahedron B-matrix for coordinates(0.5, 0.5, 0.5)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << 0.5, 0.5, 0.5;

      // Get B-Matrix
      auto bmatrix = hex->bmatrix(coords);

      // Check gradient of shape functions
      auto gradsf = hex->grad_shapefn(coords);

      // Check size of B-matrix
      REQUIRE(bmatrix.size() == nfunctions);

      for (unsigned i = 0; i < nfunctions; ++i) {
        REQUIRE(bmatrix.at(i)(0, 0) == Approx(gradsf(i, 0)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(0, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(0, 2) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 1) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 2) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 2) == Approx(gradsf(i, 2)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(3, 0) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(3, 1) == Approx(gradsf(i, 0)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(3, 2) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(4, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(4, 1) == Approx(gradsf(i, 2)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(4, 2) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(5, 0) == Approx(gradsf(i, 2)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(5, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(5, 2) == Approx(gradsf(i, 0)).epsilon(Tolerance));
      }
    }

    // Coordinates is (-0.5, -0.5, -0.5)
    SECTION(
        "Eight noded hexahedron B-matrix for coordinates(-0.5, -0.5, -0.5)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << -0.5, -0.5, -0.5;

      // Get B-Matrix
      auto bmatrix = hex->bmatrix(coords);

      // Check gradient of shape functions
      auto gradsf = hex->grad_shapefn(coords);

      // Check size of B-matrix
      REQUIRE(bmatrix.size() == nfunctions);

      for (unsigned i = 0; i < nfunctions; ++i) {
        REQUIRE(bmatrix.at(i)(0, 0) == Approx(gradsf(i, 0)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(0, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(0, 2) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 1) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 2) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 2) == Approx(gradsf(i, 2)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(3, 0) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(3, 1) == Approx(gradsf(i, 0)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(3, 2) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(4, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(4, 1) == Approx(gradsf(i, 2)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(4, 2) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(5, 0) == Approx(gradsf(i, 2)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(5, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(5, 2) == Approx(gradsf(i, 0)).epsilon(Tolerance));
      }
    }

    // Coordinates is (0, 0, 0)
    SECTION("Eight noded hexahedron B-matrix cell for coordinates(0, 0, 0)") {
      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0., 0., 0.;

      Eigen::Matrix<double, 8, Dim> coords;
      // clang-format off
      coords << 0., 0., 0.,
                1., 0., 0.,
                1., 1., 0.,
                0., 1., 0.,
                0., 0., 1.,
                1., 0., 1.,
                1., 1., 1.,
                0., 1., 1.;
      // clang-format on

      // Get B-Matrix
      auto bmatrix = hex->bmatrix(xi, coords);

      // Check gradient of shape functions
      auto gradsf = hex->grad_shapefn(xi);
      gradsf *= 2.;

      // Check size of B-matrix
      REQUIRE(bmatrix.size() == nfunctions);

      for (unsigned i = 0; i < nfunctions; ++i) {
        REQUIRE(bmatrix.at(i)(0, 0) == Approx(gradsf(i, 0)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(0, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(0, 2) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 1) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 2) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 2) == Approx(gradsf(i, 2)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(3, 0) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(3, 1) == Approx(gradsf(i, 0)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(3, 2) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(4, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(4, 1) == Approx(gradsf(i, 2)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(4, 2) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(5, 0) == Approx(gradsf(i, 2)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(5, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(5, 2) == Approx(gradsf(i, 0)).epsilon(Tolerance));
      }
    }

    // Coordinates is (0, 0, 0)
    SECTION("Eight noded hexahedron B-matrix cell for xi(0.5, 0.5, 0.5)") {
      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0.5, 0.5, 0.5;

      Eigen::Matrix<double, 8, Dim> coords;
      // clang-format off
      coords << 0., 0., 0.,
                1., 0., 0.,
                1., 1., 0.,
                0., 1., 0.,
                0., 0., 1.,
                1., 0., 1.,
                1., 1., 1.,
                0., 1., 1.;
      // clang-format on

      // Get B-Matrix
      auto bmatrix = hex->bmatrix(xi, coords);

      // Check gradient of shape functions
      auto gradsf = hex->grad_shapefn(xi);
      gradsf *= 2.;

      // Check size of B-matrix
      REQUIRE(bmatrix.size() == nfunctions);

      for (unsigned i = 0; i < nfunctions; ++i) {
        REQUIRE(bmatrix.at(i)(0, 0) == Approx(gradsf(i, 0)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(0, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(0, 2) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 1) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 2) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 2) == Approx(gradsf(i, 2)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(3, 0) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(3, 1) == Approx(gradsf(i, 0)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(3, 2) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(4, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(4, 1) == Approx(gradsf(i, 2)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(4, 2) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(5, 0) == Approx(gradsf(i, 2)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(5, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(5, 2) == Approx(gradsf(i, 0)).epsilon(Tolerance));
      }
    }

    // Coordinates is (-0.5, -0.5, -0.5)
    SECTION("Eight noded hexahedron B-matrix cell xi(-0.5, -0.5, -0.5)") {
      Eigen::Matrix<double, Dim, 1> xi;
      xi << -0.5, -0.5, -0.5;

      Eigen::Matrix<double, 8, Dim> coords;
      // clang-format off
      coords << 0., 0., 0.,
                1., 0., 0.,
                1., 1., 0.,
                0., 1., 0.,
                0., 0., 1.,
                1., 0., 1.,
                1., 1., 1.,
                0., 1., 1.;
      // clang-format on

      // Get B-Matrix
      auto bmatrix = hex->bmatrix(xi, coords);

      // Check gradient of shape functions
      auto gradsf = hex->grad_shapefn(xi);
      gradsf *= 2.;

      // Check size of B-matrix
      REQUIRE(bmatrix.size() == nfunctions);

      for (unsigned i = 0; i < nfunctions; ++i) {
        REQUIRE(bmatrix.at(i)(0, 0) == Approx(gradsf(i, 0)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(0, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(0, 2) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 1) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 2) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 2) == Approx(gradsf(i, 2)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(3, 0) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(3, 1) == Approx(gradsf(i, 0)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(3, 2) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(4, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(4, 1) == Approx(gradsf(i, 2)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(4, 2) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(5, 0) == Approx(gradsf(i, 2)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(5, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(5, 2) == Approx(gradsf(i, 0)).epsilon(Tolerance));
      }
    }

    SECTION("Eight noded hexahedron B-matrix with deformation gradient") {
      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0., 0., 0.;

      Eigen::Matrix<double, Dim, 1> psize;
      psize.setZero();
      Eigen::Matrix<double, Dim, 1> defgrad;
      defgrad.setZero();

      Eigen::Matrix<double, 8, Dim> coords;
      // clang-format off
      coords << 0., 0., 0.,
                1., 0., 0.,
                1., 1., 0.,
                0., 1., 0.,
                0., 0., 1.,
                1., 0., 1.,
                1., 1., 1.,
                0., 1., 1.;
      // clang-format on

      // Get B-Matrix
      auto bmatrix = hex->bmatrix(xi, coords, psize, defgrad);

      // Check gradient of shape functions
      auto gradsf = hex->grad_shapefn(xi, psize, defgrad);
      gradsf *= 2.;

      // Check size of B-matrix
      REQUIRE(bmatrix.size() == nfunctions);

      for (unsigned i = 0; i < nfunctions; ++i) {
        REQUIRE(bmatrix.at(i)(0, 0) == Approx(gradsf(i, 0)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(0, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(0, 2) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 1) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 2) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 2) == Approx(gradsf(i, 2)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(3, 0) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(3, 1) == Approx(gradsf(i, 0)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(3, 2) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(4, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(4, 1) == Approx(gradsf(i, 2)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(4, 2) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(5, 0) == Approx(gradsf(i, 2)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(5, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(5, 2) == Approx(gradsf(i, 0)).epsilon(Tolerance));
      }
    }

    SECTION("Eight noded hexahedron B-matrix and Jacobian failure") {
      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0., 0., 0.;

      Eigen::Matrix<double, 7, Dim> coords;
      // clang-format off
      coords << 0., 0., 0.,
                1., 0., 0.,
                1., 1., 0.,
                0., 1., 0.,
                0., 0., 1.,
                1., 0., 1.,
                1., 1., 1.;
      // clang-format on
      // Get B-Matrix
      auto bmatrix = hex->bmatrix(xi, coords);
      auto jacobian = hex->jacobian(xi, coords);
    }

    // Mass matrix of a cell
    SECTION("Eight noded hexahedron mass-matrix") {
      std::vector<Eigen::Matrix<double, Dim, 1>> xi_s;

      Eigen::Matrix<double, Dim, 1> xi;
      const double one_by_sqrt3 = std::fabs(1 / std::sqrt(3));
      xi << -one_by_sqrt3, -one_by_sqrt3, -one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << one_by_sqrt3, -one_by_sqrt3, -one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << one_by_sqrt3, one_by_sqrt3, -one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << -one_by_sqrt3, one_by_sqrt3, -one_by_sqrt3;
      xi_s.emplace_back(xi);

      xi << -one_by_sqrt3, -one_by_sqrt3, one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << one_by_sqrt3, -one_by_sqrt3, one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << one_by_sqrt3, one_by_sqrt3, one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << -one_by_sqrt3, one_by_sqrt3, one_by_sqrt3;
      xi_s.emplace_back(xi);

      REQUIRE(xi_s.size() == 8);

      // Get mass matrix
      const auto mass_matrix = hex->mass_matrix(xi_s);

      // Check size of mass-matrix
      REQUIRE(mass_matrix.rows() == nfunctions);
      REQUIRE(mass_matrix.cols() == nfunctions);

      // Sum should be equal to 1. * xi_s.size()
      REQUIRE(mass_matrix.sum() == Approx(1. * xi_s.size()).epsilon(Tolerance));

      Eigen::Matrix<double, 8, 8> mass;
      // clang-format off
      mass << 0.29629629629629640, 0.14814814814814810, 0.07407407407407404,
              0.14814814814814810, 0.14814814814814810, 0.07407407407407404,
              0.03703703703703701, 0.07407407407407404, 0.14814814814814810,
              0.29629629629629640, 0.14814814814814810, 0.07407407407407404,
              0.07407407407407404, 0.14814814814814810, 0.07407407407407403,
              0.03703703703703701, 0.07407407407407404, 0.14814814814814810,
              0.29629629629629640, 0.14814814814814810, 0.03703703703703701,
              0.07407407407407404, 0.14814814814814810, 0.07407407407407404,
              0.14814814814814810, 0.07407407407407404, 0.14814814814814810,
              0.29629629629629640, 0.07407407407407404, 0.03703703703703701,
              0.07407407407407404, 0.14814814814814810, 0.14814814814814810,
              0.07407407407407404, 0.03703703703703701, 0.07407407407407404,
              0.29629629629629630, 0.14814814814814810, 0.07407407407407404,
              0.14814814814814810, 0.07407407407407404, 0.14814814814814810,
              0.07407407407407404, 0.03703703703703701, 0.14814814814814810,
              0.29629629629629630, 0.14814814814814810, 0.07407407407407404,
              0.03703703703703701, 0.07407407407407403, 0.14814814814814810,
              0.07407407407407404, 0.07407407407407404, 0.14814814814814810,
              0.29629629629629630, 0.14814814814814810, 0.07407407407407404,
              0.03703703703703701, 0.07407407407407404, 0.14814814814814810,
              0.14814814814814810, 0.07407407407407404, 0.14814814814814810,
              0.2962962962962963;
      // clang-format on
      for (unsigned i = 0; i < nfunctions; ++i)
        for (unsigned j = 0; j < nfunctions; ++j)
          REQUIRE(mass_matrix(i, j) == Approx(mass(i, j)).epsilon(Tolerance));
    }

    // Laplace matrix of a cell
    SECTION("Eight noded hexahedron laplace-matrix") {
      std::vector<Eigen::Matrix<double, Dim, 1>> xi_s;

      Eigen::Matrix<double, Dim, 1> xi;
      const double one_by_sqrt3 = std::fabs(1 / std::sqrt(3));
      xi << -one_by_sqrt3, -one_by_sqrt3, -one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << one_by_sqrt3, -one_by_sqrt3, -one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << one_by_sqrt3, one_by_sqrt3, -one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << -one_by_sqrt3, one_by_sqrt3, -one_by_sqrt3;
      xi_s.emplace_back(xi);

      xi << -one_by_sqrt3, -one_by_sqrt3, one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << one_by_sqrt3, -one_by_sqrt3, one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << one_by_sqrt3, one_by_sqrt3, one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << -one_by_sqrt3, one_by_sqrt3, one_by_sqrt3;
      xi_s.emplace_back(xi);

      REQUIRE(xi_s.size() == 8);

      // Nodal coordinates
      Eigen::Matrix<double, 8, Dim> coords;
      // clang-format off
      coords << 2., 1., 0.5,
                4., 2., 1.0,
                2., 4., 1.0,
                1., 3., 0.5,
                2., 1., 1.5,
                4., 2., 2.0,
                2., 4., 2.0,
                1., 3., 1.5;
      // clang-format on

      // Get laplace matrix
      const auto laplace_matrix = hex->laplace_matrix(xi_s, coords);

      // Check size of laplace-matrix
      REQUIRE(laplace_matrix.rows() == nfunctions);
      REQUIRE(laplace_matrix.cols() == nfunctions);

      // Sum should be equal to 0.
      REQUIRE(laplace_matrix.sum() == Approx(0.).epsilon(Tolerance));

      Eigen::Matrix<double, 8, 8> laplace;
      // clang-format off
      laplace <<  0.9643677583341226,  0.261044778874057,  -0.1901695368724208,
                  0.2334579956379644, -0.6204508598313142, -0.2200815520760039,
                 -0.1638493372007791, -0.264319246865626,   0.261044778874057,
                  1.242773383095871,   0.4330321922069018,  0.1222723413095025,
                 -0.7574835440234503, -0.726073928639958,  -0.26533034028844,
                 -0.3102348825344837, -0.1901695368724208,  0.4330321922069018,
                  2.27268250189131,    0.3331516005618761, -0.6905145911387602,
                 -0.7630346561745743, -0.4913341040024276, -0.9038134064719052,
                  0.2334579956379644,  0.1222723413095025,  0.3331516005618761,
                  1.134597613242935,  -0.6628500090327518, -0.2316871676231324,
                 -0.2307894607206865, -0.6981529133757072, -0.6204508598313142,
                 -0.7574835440234503, -0.6905145911387602, -0.6628500090327518,
                  1.887162135673954,   0.4504916955937011, -0.185224986473325,
                  0.5788701592319462, -0.2200815520760039, -0.726073928639958,
                 -0.7630346561745743, -0.2316871676231324,  0.4504916955937011,
                  1.186264235677633,   0.1769044815337366,  0.1272168917085984,
                 -0.1638493372007791, -0.26533034028844, -  0.4913341040024276,
                 -0.2307894607206865, -0.185224986473325,   0.1769044815337366,
                  1.095314415432314,   0.0643093317196072, -0.264319246865626,
                 -0.3102348825344837, -0.9038134064719052, -0.6981529133757072,
                  0.5788701592319462,  0.1272168917085984,  0.0643093317196072,
                  1.40612406658757;
      // clang-format on
      for (unsigned i = 0; i < nfunctions; ++i)
        for (unsigned j = 0; j < nfunctions; ++j)
          REQUIRE(laplace_matrix(i, j) ==
                  Approx(laplace(i, j)).epsilon(Tolerance));
    }

    SECTION("Eight noded hexahedron coordinates of unit cell") {
      const unsigned nfunctions = 8;
      // Coordinates of a unit cell
      Eigen::Matrix<double, nfunctions, Dim> unit_cell;
      // clang-format off
      unit_cell << -1., -1., -1.,
                    1., -1., -1.,
                    1.,  1., -1.,
                   -1.,  1., -1.,
                   -1., -1.,  1.,
                    1., -1.,  1.,
                    1.,  1.,  1.,
                   -1.,  1.,  1.;
      // clang-format on

      auto coordinates = hex->unit_cell_coordinates();
      REQUIRE(coordinates.rows() == nfunctions);
      REQUIRE(coordinates.cols() == Dim);
      for (unsigned i = 0; i < nfunctions; ++i) {  // Iterate through nfunctions
        for (unsigned j = 0; j < Dim; ++j) {       // Dimension
          REQUIRE(coordinates(i, j) ==
                  Approx(unit_cell(i, j)).epsilon(Tolerance));
        }
      }
    }

    SECTION("Eight noded hexahedron element for sides indices") {
      // Check for sides indices
      Eigen::MatrixXi indices = hex->sides_indices();
      REQUIRE(indices.rows() == 12);
      REQUIRE(indices.cols() == 2);
      REQUIRE(indices(0, 0) == 0);
      REQUIRE(indices(0, 1) == 1);

      REQUIRE(indices(1, 0) == 1);
      REQUIRE(indices(1, 1) == 2);

      REQUIRE(indices(2, 0) == 2);
      REQUIRE(indices(2, 1) == 3);

      REQUIRE(indices(3, 0) == 3);
      REQUIRE(indices(3, 1) == 0);

      REQUIRE(indices(4, 0) == 4);
      REQUIRE(indices(4, 1) == 5);

      REQUIRE(indices(5, 0) == 5);
      REQUIRE(indices(5, 1) == 6);

      REQUIRE(indices(6, 0) == 6);
      REQUIRE(indices(6, 1) == 7);

      REQUIRE(indices(7, 0) == 7);
      REQUIRE(indices(7, 1) == 4);

      REQUIRE(indices(8, 0) == 0);
      REQUIRE(indices(8, 1) == 4);

      REQUIRE(indices(9, 0) == 1);
      REQUIRE(indices(9, 1) == 5);

      REQUIRE(indices(10, 0) == 2);
      REQUIRE(indices(10, 1) == 6);

      REQUIRE(indices(11, 0) == 3);
      REQUIRE(indices(11, 1) == 7);
    }

    SECTION("Eight noded hexahedron element for corner indices") {
      // Check for volume indices
      Eigen::VectorXi indices = hex->corner_indices();
      REQUIRE(indices.size() == 8);
      REQUIRE(indices(0) == 0);
      REQUIRE(indices(1) == 1);
      REQUIRE(indices(2) == 2);
      REQUIRE(indices(3) == 3);
      REQUIRE(indices(4) == 4);
      REQUIRE(indices(5) == 5);
      REQUIRE(indices(6) == 6);
      REQUIRE(indices(7) == 7);
    }

    SECTION("Eight noded hexahedron element for inhedron indices") {
      // Check for inhedron indices
      Eigen::MatrixXi indices = hex->inhedron_indices();
      REQUIRE(indices.rows() == 12);
      REQUIRE(indices.cols() == 3);
      REQUIRE(indices(0, 0) == 0);
      REQUIRE(indices(0, 1) == 5);
      REQUIRE(indices(0, 2) == 4);

      REQUIRE(indices(1, 0) == 0);
      REQUIRE(indices(1, 1) == 1);
      REQUIRE(indices(1, 2) == 5);

      REQUIRE(indices(2, 0) == 3);
      REQUIRE(indices(2, 1) == 6);
      REQUIRE(indices(2, 2) == 7);

      REQUIRE(indices(3, 0) == 3);
      REQUIRE(indices(3, 1) == 2);
      REQUIRE(indices(3, 2) == 6);

      REQUIRE(indices(4, 0) == 2);
      REQUIRE(indices(4, 1) == 1);
      REQUIRE(indices(4, 2) == 6);

      REQUIRE(indices(5, 0) == 6);
      REQUIRE(indices(5, 1) == 1);
      REQUIRE(indices(5, 2) == 5);

      REQUIRE(indices(6, 0) == 7);
      REQUIRE(indices(6, 1) == 6);
      REQUIRE(indices(6, 2) == 5);

      REQUIRE(indices(7, 0) == 5);
      REQUIRE(indices(7, 1) == 4);
      REQUIRE(indices(7, 2) == 7);

      REQUIRE(indices(8, 0) == 7);
      REQUIRE(indices(8, 1) == 4);
      REQUIRE(indices(8, 2) == 0);

      REQUIRE(indices(9, 0) == 7);
      REQUIRE(indices(9, 1) == 0);
      REQUIRE(indices(9, 2) == 3);

      REQUIRE(indices(10, 0) == 3);
      REQUIRE(indices(10, 1) == 0);
      REQUIRE(indices(10, 2) == 1);

      REQUIRE(indices(11, 0) == 3);
      REQUIRE(indices(11, 1) == 1);
      REQUIRE(indices(11, 2) == 2);
    }

    SECTION("Eight noded hexahedron shape function for face indices") {
      // Check for face indices
      Eigen::Matrix<int, 6, 4> indices;
      // clang-format off
      indices << 0, 1, 5, 4,
                 5, 1, 2, 0,
                 7, 6, 2, 3,
                 0, 4, 7, 3,
                 1, 0, 3, 2,
                 4, 5, 6, 7;
      // clang-format on

      // Check for all face indices
      for (unsigned i = 0; i < indices.rows(); ++i) {
        const auto check_indices = hex->face_indices(i);
        REQUIRE(check_indices.rows() == 4);
        REQUIRE(check_indices.cols() == 1);

        for (unsigned j = 0; j < indices.cols(); ++j)
          REQUIRE(check_indices(j) == indices(i, j));
      }
    }
  }

  // 20-Node (Serendipity) Hexahedron Element
  //! Check for 8 noded element
  SECTION("Hexahedron element with twenty nodes") {
    const unsigned nfunctions = 20;
    std::shared_ptr<mpm::Element<Dim>> hex =
        std::make_shared<mpm::HexahedronElement<Dim, nfunctions>>();

    // Check degree
    REQUIRE(hex->degree() == mpm::ElementDegree::Quadratic);

    // Coordinates is (0, 0, 0)
    SECTION("Twenty noded hexahedron element for coordinates(0, 0, 0)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();
      auto shapefn = hex->shapefn(coords);

      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      // Edge nodes
      // N0 = 0.125*(1 - xi(0))(1. - xi(1))(1 - xi(2))(-xi(0) - xi(1) -xi(2) -2)
      // N1 = 0.125*(1 + xi(0))(1. - xi(1))(1 - xi(2))(+xi(0) - xi(1) -xi(2) -2)
      // N2 = 0.125*(1 + xi(0))(1. + xi(1))(1 - xi(2))(+xi(0) + xi(1) -xi(2) -2)
      // N3 = 0.125*(1 - xi(0))(1. + xi(1))(1 - xi(2))(-xi(0) + xi(1) -xi(2) -2)
      // N4 = 0.125*(1 - xi(0))(1. - xi(1))(1 + xi(2))(-xi(0) - xi(1) +xi(2) -2)
      // N5 = 0.125*(1 + xi(0))(1. - xi(1))(1 + xi(2))(+xi(0) - xi(1) +xi(2) -2)
      // N6 = 0.125*(1 + xi(0))(1. + xi(1))(1 + xi(2))(+xi(0) + xi(1) +xi(2) -2)
      // N7 = 0.125*(1 - xi(0))(1. + xi(1))(1 + xi(2))(-xi(0) + xi(1) +xi(2) -2)
      REQUIRE(shapefn(0) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(shapefn(3) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(shapefn(4) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(shapefn(5) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(shapefn(6) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(shapefn(7) == Approx(-0.25).epsilon(Tolerance));

      // Midside nodes
      // N8 = 0.25*(1 - xi(0)^2)(1 - xi(1))(1 - xi(2))
      // N9 = 0.25*(1 - xi(0))(1 - xi(1)^2)(1 - xi(2))
      // N10 = 0.25*(1 - xi(0))(1 - xi(1))(1 - xi(2)^2)
      // N11 = 0.25*(1 + xi(0))(1 - xi(1)^2)(1 - xi(2))
      // N12 = 0.25*(1 + xi(0))(1 - xi(1))(1 - xi(2)^2)
      // N13 = 0.25*(1 - xi(0)^2)(1 + xi(1))(1 - xi(2))
      // N14 = 0.25*(1 + xi(0))(1 + xi(1))(1 - xi(2)^2)
      // N15 = 0.25*(1 - xi(0))(1 + xi(1))(1 - xi(2)^2)
      // N16 = 0.25*(1 - xi(0)^2)(1 - xi(1))(1 + xi(2))
      // N17 = 0.25*(1 - xi(0))(1 - xi(1)^2)(1 + xi(2))
      // N18 = 0.25*(1 + xi(0))(1 - xi(1)^2)(1 + xi(2))
      // N19 = 0.25*(1 - xi(0)^2)(1 + xi(1))(1 + xi(2))
      REQUIRE(shapefn(8) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(9) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(10) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(11) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(12) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(13) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(14) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(15) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(16) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(17) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(18) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(19) == Approx(0.25).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = hex->grad_shapefn(coords);
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      // Derivatives with respect to eta-direction
      // Edge nodes
      // G(0,0) = 0.125*(1. - xi(1))(1 - xi(2))(2*xi(0) + xi(1) + xi(2) + 1.)
      // G(1,0) = 0.125*(1. - xi(1))(1 - xi(2))(2*xi(0) - xi(1) - xi(2) - 1.)
      // G(2,0) = 0.125*(1. + xi(1))(1 - xi(2))(2*xi(0) + xi(1) - xi(2) - 1.)
      // G(3,0) = 0.125*(1. + xi(1))(1 - xi(2))(2*xi(0) - xi(1) + xi(2) + 1.)
      // G(4,0) = 0.125*(1. - xi(1))(1 + xi(2))(2*xi(0) + xi(1) - xi(2) + 1.)
      // G(5,0) = 0.125*(1. - xi(1))(1 + xi(2))(2*xi(0) - xi(1) + xi(2) - 1.)
      // G(6,0) = 0.125*(1. + xi(1))(1 + xi(2))(2*xi(0) + xi(1) + xi(2) - 1.)
      // G(7,0) = 0.125*(1. + xi(1))(1 + xi(2))(2*xi(0) - xi(1) - xi(2) + 1.)
      REQUIRE(gradsf(0, 0) == Approx(+0.125).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(3, 0) == Approx(+0.125).epsilon(Tolerance));
      REQUIRE(gradsf(4, 0) == Approx(+0.125).epsilon(Tolerance));
      REQUIRE(gradsf(5, 0) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(6, 0) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(7, 0) == Approx(+0.125).epsilon(Tolerance));

      // Midside nodes
      // G(8, 0) = -0.5*xi(0)*(1 - xi(1))(1 - xi(2));
      // G(9, 0) = -0.25*(1 - xi(1)xi(1))(1 - xi(2));
      // G(10, 0) = -0.25*(1 - xi(2)xi(2))(1 - xi(1));
      // G(11, 0) = 0.25*(1 - xi(1)xi(1))(1 - xi(2));
      // G(12, 0) = 0.25*(1 - xi(2)xi(2))(1 - xi(1));
      // G(13, 0) = -0.5*xi(0)(1 + xi(1))(1 - xi(2));
      // G(14, 0) = 0.25*(1 - xi(2)xi(2))(1 + xi(1));
      // G(15, 0) = -0.25*(1 - xi(2)xi(2))(1 + xi(1));
      // G(16, 0) = -0.5*xi(0)(1 - xi(1))(1 + xi(2));
      // G(17, 0) = -0.25*(1 - xi(1)xi(1))(1 + xi(2));
      // G(18, 0) = 0.25*(1 - xi(1)xi(1))(1 + xi(2));
      // G(19, 0) = -0.5*xi(0)(1 + xi(1))(1 + xi(2));
      REQUIRE(gradsf(8, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(9, 0) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(10, 0) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(11, 0) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(gradsf(12, 0) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(gradsf(13, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(14, 0) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(gradsf(15, 0) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(16, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(17, 0) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(18, 0) == Approx(+0.25).epsilon(Tolerance));
      REQUIRE(gradsf(19, 0) == Approx(0.0).epsilon(Tolerance));

      // Derivatives with respect to psi-direction
      // Edge nodes
      // G(0,1) = 0.125*(1. - xi(0))(1 - xi(2))(+xi(0) + 2*xi(1) + xi(2) + 1.)
      // G(1,1) = 0.125*(1. + xi(0))(1 - xi(2))(-xi(0) + 2*xi(1) + xi(2) + 1.)
      // G(2,1) = 0.125*(1. + xi(0))(1 - xi(2))(+xi(0) + 2*xi(1) - xi(2) - 1.)
      // G(3,1) = 0.125*(1. - xi(0))(1 - xi(2))(-xi(0) + 2*xi(1) - xi(2) - 1.)
      // G(4,1) = 0.125*(1. - xi(0))(1 + xi(2))(+xi(0) + 2*xi(1) - xi(2) + 1.)
      // G(5,1) = 0.125*(1. + xi(0))(1 + xi(2))(-xi(0) + 2*xi(1) - xi(2) + 1.)
      // G(6,1) = 0.125*(1. + xi(0))(1 + xi(2))(+xi(0) + 2*xi(1) + xi(2) - 1.)
      // G(7,1) = 0.125*(1. - xi(0))(1 + xi(2))(-xi(0) + 2*xi(1) + xi(2) - 1.)
      REQUIRE(gradsf(0, 1) == Approx(+0.125).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(+0.125).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(3, 1) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(4, 1) == Approx(+0.125).epsilon(Tolerance));
      REQUIRE(gradsf(5, 1) == Approx(+0.125).epsilon(Tolerance));
      REQUIRE(gradsf(6, 1) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(7, 1) == Approx(-0.125).epsilon(Tolerance));

      // Midside nodes
      // G(8, 1) = -0.25*(1 - xi(0)xi(0))(1 - xi(2));
      // G(9, 1) = -0.5*xi(1)(1 - xi(0))(1 - xi(2));
      // G(10, 1) = -0.25*(1 - xi(2)xi(2))(1 - xi(0));
      // G(11, 1) = -0.5*xi(1)(1 + xi(0))(1 - xi(2));
      // G(12, 1) = -0.25*(1 - xi(2)xi(2))(1 + xi(0));
      // G(13, 1) = 0.25*(1 - xi(0)xi(0))(1 - xi(2));
      // G(14, 1) = 0.25*(1 - xi(2)xi(2))(1 + xi(0));
      // G(15, 1) = 0.25*(1 - xi(2)xi(2))(1 - xi(0));
      // G(16, 1) = -0.25*(1 - xi(0)xi(0))(1 + xi(2));
      // G(17, 1) = -0.5*xi(1)(1 - xi(0))(1 + xi(2));
      // G(18, 1) = -0.5*xi(1)(1 + xi(0))(1 + xi(2));
      // G(19, 1) = 0.25*(1 - xi(0)xi(0))(1 + xi(2));
      REQUIRE(gradsf(8, 1) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(9, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(10, 1) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(11, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(12, 1) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(13, 1) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(gradsf(14, 1) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(gradsf(15, 1) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(gradsf(16, 1) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(17, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(18, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(19, 1) == Approx(0.25).epsilon(Tolerance));

      // Derivatives with respect to mu-direction
      // Edge nodes
      // G(0,2) = 0.125*(1. - xi(0))(1 - xi(1))(+xi(0) + xi(1) + 2*xi(2) + 1.)
      // G(1,2) = 0.125*(1. + xi(0))(1 - xi(1))(-xi(0) + xi(1) + 2*xi(2) + 1.)
      // G(2,2) = 0.125*(1. + xi(0))(1 + xi(1))(-xi(0) - xi(1) + 2*xi(2) + 1.)
      // G(3,2) = 0.125*(1. - xi(0))(1 + xi(1))(+xi(0) - xi(1) + 2*xi(2) + 1.)
      // G(4,2) = 0.125*(1. - xi(0))(1 - xi(1))(-xi(0) - xi(1) + 2*xi(2) - 1.)
      // G(5,2) = 0.125*(1. + xi(0))(1 - xi(1))(+xi(0) - xi(1) + 2*xi(2) - 1.)
      // G(6,2) = 0.125*(1. + xi(0))(1 + xi(1))(+xi(0) + xi(1) + 2*xi(2) - 1.)
      // G(7,2) = 0.125*(1. - xi(0))(1 + xi(1))(-xi(0) + xi(1) + 2*xi(2) - 1.)
      REQUIRE(gradsf(0, 2) == Approx(+0.125).epsilon(Tolerance));
      REQUIRE(gradsf(1, 2) == Approx(+0.125).epsilon(Tolerance));
      REQUIRE(gradsf(2, 2) == Approx(+0.125).epsilon(Tolerance));
      REQUIRE(gradsf(3, 2) == Approx(+0.125).epsilon(Tolerance));
      REQUIRE(gradsf(4, 2) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(5, 2) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(6, 2) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(7, 2) == Approx(-0.125).epsilon(Tolerance));

      // Midside nodes
      // G(8, 2) = -0.25*(1 - xi(0)xi(0))(1 - xi(1));
      // G(9, 2) = -0.25*(1 - xi(1)xi(1))(1 - xi(0));
      // G(10, 2) = -0.5*xi(2)(1 - xi(0))(1 - xi(1));
      // G(11, 2) = -0.25*(1 - xi(1)xi(1))(1 + xi(0));
      // G(12, 2) = -0.5*xi(2)(1 + xi(0))(1 - xi(1));
      // G(13, 2) = -0.25*(1 - xi(0)xi(0))(1 + xi(1));
      // G(14, 2) = -0.5*xi(2)(1 + xi(0))(1 + xi(1));
      // G(15, 2) = -0.5*xi(2)(1 - xi(0))(1 + xi(1));
      // G(16, 2) = 0.25*(1 - xi(0)xi(0))(1 - xi(1));
      // G(17, 2) = 0.25*(1 - xi(1)xi(1))(1 - xi(0));
      // G(18, 2) = 0.25*(1 - xi(1)xi(1))(1 + xi(0));
      // G(19, 2) = 0.25*(1 - xi(0)xi(0))(1 + xi(1));
      REQUIRE(gradsf(8, 2) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(9, 2) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(10, 2) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(11, 2) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(12, 2) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(13, 2) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(14, 2) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(15, 2) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(16, 2) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(gradsf(17, 2) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(gradsf(18, 2) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(gradsf(19, 2) == Approx(0.25).epsilon(Tolerance));
    }

    // Coordinates is (-0.5, -0.5, -0,5)
    SECTION(
        "Twenty noded hexahedron element for coordinates(-0.5, "
        "-0.5, -0.5)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << -0.5, -0.5, -0.5;
      auto shapefn = hex->shapefn(coords);

      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      // Edge nodes
      // N0 = 0.125*(1 - xi(0))(1. - xi(1))(1 - xi(2))(-xi(0) - xi(1) -xi(2) -2)
      // N1 = 0.125*(1 + xi(0))(1. - xi(1))(1 - xi(2))(+xi(0) - xi(1) -xi(2) -2)
      // N2 = 0.125*(1 + xi(0))(1. + xi(1))(1 - xi(2))(+xi(0) + xi(1) -xi(2) -2)
      // N3 = 0.125*(1 - xi(0))(1. + xi(1))(1 - xi(2))(-xi(0) + xi(1) -xi(2) -2)
      // N4 = 0.125*(1 - xi(0))(1. - xi(1))(1 + xi(2))(-xi(0) - xi(1) +xi(2) -2)
      // N5 = 0.125*(1 + xi(0))(1. - xi(1))(1 + xi(2))(+xi(0) - xi(1) +xi(2) -2)
      // N6 = 0.125*(1 + xi(0))(1. + xi(1))(1 + xi(2))(+xi(0) + xi(1) +xi(2) -2)
      // N7 = 0.125*(1 - xi(0))(1. + xi(1))(1 + xi(2))(-xi(0) + xi(1) +xi(2) -2)
      REQUIRE(shapefn(0) == Approx(-0.2109375).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(-0.2109375).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(-0.1171875).epsilon(Tolerance));
      REQUIRE(shapefn(3) == Approx(-0.2109375).epsilon(Tolerance));
      REQUIRE(shapefn(4) == Approx(-0.2109375).epsilon(Tolerance));
      REQUIRE(shapefn(5) == Approx(-0.1171875).epsilon(Tolerance));
      REQUIRE(shapefn(6) == Approx(-0.0546875).epsilon(Tolerance));
      REQUIRE(shapefn(7) == Approx(-0.1171875).epsilon(Tolerance));

      // Midside nodes
      // N8 = 0.25*(1 - xi(0)^2)(1 - xi(1))(1 - xi(2))
      // N9 = 0.25*(1 - xi(0))(1 - xi(1)^2)(1 - xi(2))
      // N10 = 0.25*(1 - xi(0))(1 - xi(1))(1 - xi(2)^2)
      // N11 = 0.25*(1 + xi(0))(1 - xi(1)^2)(1 - xi(2))
      // N12 = 0.25*(1 + xi(0))(1 - xi(1))(1 - xi(2)^2)
      // N13 = 0.25*(1 - xi(0)^2)(1 + xi(1))(1 - xi(2))
      // N14 = 0.25*(1 + xi(0))(1 + xi(1))(1 - xi(2)^2)
      // N15 = 0.25*(1 - xi(0))(1 + xi(1))(1 - xi(2)^2)
      // N16 = 0.25*(1 - xi(0)^2)(1 - xi(1))(1 + xi(2))
      // N17 = 0.25*(1 - xi(0))(1 - xi(1)^2)(1 + xi(2))
      // N18 = 0.25*(1 + xi(0))(1 - xi(1)^2)(1 + xi(2))
      // N19 = 0.25*(1 - xi(0)^2)(1 + xi(1))(1 + xi(2))
      REQUIRE(shapefn(8) == Approx(0.421875).epsilon(Tolerance));
      REQUIRE(shapefn(9) == Approx(0.421875).epsilon(Tolerance));
      REQUIRE(shapefn(10) == Approx(0.421875).epsilon(Tolerance));
      REQUIRE(shapefn(11) == Approx(0.140625).epsilon(Tolerance));
      REQUIRE(shapefn(12) == Approx(0.140625).epsilon(Tolerance));
      REQUIRE(shapefn(13) == Approx(0.140625).epsilon(Tolerance));
      REQUIRE(shapefn(14) == Approx(0.046875).epsilon(Tolerance));
      REQUIRE(shapefn(15) == Approx(0.140625).epsilon(Tolerance));
      REQUIRE(shapefn(16) == Approx(0.140625).epsilon(Tolerance));
      REQUIRE(shapefn(17) == Approx(0.140625).epsilon(Tolerance));
      REQUIRE(shapefn(18) == Approx(0.046875).epsilon(Tolerance));
      REQUIRE(shapefn(19) == Approx(0.046875).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = hex->grad_shapefn(coords);
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      // Derivatives with respect to eta-direction
      // Edge nodes
      // G(0,0) = 0.125*(1. - xi(1))(1 - xi(2))(2*xi(0) + xi(1) + xi(2) + 1.)
      // G(1,0) = 0.125*(1. - xi(1))(1 - xi(2))(2*xi(0) - xi(1) - xi(2) - 1.)
      // G(2,0) = 0.125*(1. + xi(1))(1 - xi(2))(2*xi(0) + xi(1) - xi(2) - 1.)
      // G(3,0) = 0.125*(1. + xi(1))(1 - xi(2))(2*xi(0) - xi(1) + xi(2) + 1.)
      // G(4,0) = 0.125*(1. - xi(1))(1 + xi(2))(2*xi(0) + xi(1) - xi(2) + 1.)
      // G(5,0) = 0.125*(1. - xi(1))(1 + xi(2))(2*xi(0) - xi(1) + xi(2) - 1.)
      // G(6,0) = 0.125*(1. + xi(1))(1 + xi(2))(2*xi(0) + xi(1) + xi(2) - 1.)
      // G(7,0) = 0.125*(1. + xi(1))(1 + xi(2))(2*xi(0) - xi(1) - xi(2) + 1.)
      REQUIRE(gradsf(0, 0) == Approx(-0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(-0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(-0.1875).epsilon(Tolerance));
      REQUIRE(gradsf(3, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(4, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(5, 0) == Approx(-0.1875).epsilon(Tolerance));
      REQUIRE(gradsf(6, 0) == Approx(-0.09375).epsilon(Tolerance));
      REQUIRE(gradsf(7, 0) == Approx(0.03125).epsilon(Tolerance));

      // Midside nodes
      // G(8, 0) = -0.5*xi(0)*(1 - xi(1))(1 - xi(2));
      // G(9, 0) = -0.25*(1 - xi(1)xi(1))(1 - xi(2));
      // G(10, 0) = -0.25*(1 - xi(2)xi(2))(1 - xi(1));
      // G(11, 0) = 0.25*(1 - xi(1)xi(1))(1 - xi(2));
      // G(12, 0) = 0.25*(1 - xi(2)xi(2))(1 - xi(1));
      // G(13, 0) = -0.5*xi(0)(1 + xi(1))(1 - xi(2));
      // G(14, 0) = 0.25*(1 - xi(2)xi(2))(1 + xi(1));
      // G(15, 0) = -0.25*(1 - xi(2)xi(2))(1 + xi(1));
      // G(16, 0) = -0.5*xi(0)(1 - xi(1))(1 + xi(2));
      // G(17, 0) = -0.25*(1 - xi(1)xi(1))(1 + xi(2));
      // G(18, 0) = 0.25*(1 - xi(1)xi(1))(1 + xi(2));
      // G(19, 0) = -0.5*xi(0)(1 + xi(1))(1 + xi(2));
      REQUIRE(gradsf(8, 0) == Approx(0.5625).epsilon(Tolerance));
      REQUIRE(gradsf(9, 0) == Approx(-0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(10, 0) == Approx(-0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(11, 0) == Approx(0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(12, 0) == Approx(0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(13, 0) == Approx(0.1875).epsilon(Tolerance));
      REQUIRE(gradsf(14, 0) == Approx(0.09375).epsilon(Tolerance));
      REQUIRE(gradsf(15, 0) == Approx(-0.09375).epsilon(Tolerance));
      REQUIRE(gradsf(16, 0) == Approx(0.1875).epsilon(Tolerance));
      REQUIRE(gradsf(17, 0) == Approx(-0.09375).epsilon(Tolerance));
      REQUIRE(gradsf(18, 0) == Approx(0.09375).epsilon(Tolerance));
      REQUIRE(gradsf(19, 0) == Approx(0.0625).epsilon(Tolerance));

      // Derivatives with respect to psi-direction
      // Edge nodes
      // G(0,1) = 0.125*(1. - xi(0))(1 - xi(2))(+xi(0) + 2*xi(1) + xi(2) + 1.)
      // G(1,1) = 0.125*(1. + xi(0))(1 - xi(2))(-xi(0) + 2*xi(1) + xi(2) + 1.)
      // G(2,1) = 0.125*(1. + xi(0))(1 - xi(2))(+xi(0) + 2*xi(1) - xi(2) - 1.)
      // G(3,1) = 0.125*(1. - xi(0))(1 - xi(2))(-xi(0) + 2*xi(1) - xi(2) - 1.)
      // G(4,1) = 0.125*(1. - xi(0))(1 + xi(2))(+xi(0) + 2*xi(1) - xi(2) + 1.)
      // G(5,1) = 0.125*(1. + xi(0))(1 + xi(2))(-xi(0) + 2*xi(1) - xi(2) + 1.)
      // G(6,1) = 0.125*(1. + xi(0))(1 + xi(2))(+xi(0) + 2*xi(1) + xi(2) - 1.)
      // G(7,1) = 0.125*(1. - xi(0))(1 + xi(2))(-xi(0) + 2*xi(1) + xi(2) - 1.)
      REQUIRE(gradsf(0, 1) == Approx(-0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(0.).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(-0.1875).epsilon(Tolerance));
      REQUIRE(gradsf(3, 1) == Approx(-0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(4, 1) == Approx(0.).epsilon(Tolerance));
      REQUIRE(gradsf(5, 1) == Approx(0.03125).epsilon(Tolerance));
      REQUIRE(gradsf(6, 1) == Approx(-0.09375).epsilon(Tolerance));
      REQUIRE(gradsf(7, 1) == Approx(-0.1875).epsilon(Tolerance));

      // Midside nodes
      // G(8, 1) = -0.25*(1 - xi(0)xi(0))(1 - xi(2));
      // G(9, 1) = -0.5*xi(1)(1 - xi(0))(1 - xi(2));
      // G(10, 1) = -0.25*(1 - xi(2)xi(2))(1 - xi(0));
      // G(11, 1) = -0.5*xi(1)(1 + xi(0))(1 - xi(2));
      // G(12, 1) = -0.25*(1 - xi(2)xi(2))(1 + xi(0));
      // G(13, 1) = 0.25*(1 - xi(0)xi(0))(1 - xi(2));
      // G(14, 1) = 0.25*(1 - xi(2)xi(2))(1 + xi(0));
      // G(15, 1) = 0.25*(1 - xi(2)xi(2))(1 - xi(0));
      // G(16, 1) = -0.25*(1 - xi(0)xi(0))(1 + xi(2));
      // G(17, 1) = -0.5*xi(1)(1 - xi(0))(1 + xi(2));
      // G(18, 1) = -0.5*xi(1)(1 + xi(0))(1 + xi(2));
      // G(19, 1) = 0.25*(1 - xi(0)xi(0))(1 + xi(2));
      REQUIRE(gradsf(8, 1) == Approx(-0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(9, 1) == Approx(0.5625).epsilon(Tolerance));
      REQUIRE(gradsf(10, 1) == Approx(-0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(11, 1) == Approx(0.1875).epsilon(Tolerance));
      REQUIRE(gradsf(12, 1) == Approx(-0.09375).epsilon(Tolerance));
      REQUIRE(gradsf(13, 1) == Approx(0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(14, 1) == Approx(0.09375).epsilon(Tolerance));
      REQUIRE(gradsf(15, 1) == Approx(0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(16, 1) == Approx(-0.09375).epsilon(Tolerance));
      REQUIRE(gradsf(17, 1) == Approx(0.1875).epsilon(Tolerance));
      REQUIRE(gradsf(18, 1) == Approx(0.0625).epsilon(Tolerance));
      REQUIRE(gradsf(19, 1) == Approx(0.09375).epsilon(Tolerance));

      // Derivatives with respect to mu-direction
      // Edge nodes
      // G(0,2) = 0.125*(1. - xi(0))(1 - xi(1))(+xi(0) + xi(1) + 2*xi(2) + 1.)
      // G(1,2) = 0.125*(1. + xi(0))(1 - xi(1))(-xi(0) + xi(1) + 2*xi(2) + 1.)
      // G(2,2) = 0.125*(1. + xi(0))(1 + xi(1))(-xi(0) - xi(1) + 2*xi(2) + 1.)
      // G(3,2) = 0.125*(1. - xi(0))(1 + xi(1))(+xi(0) - xi(1) + 2*xi(2) + 1.)
      // G(4,2) = 0.125*(1. - xi(0))(1 - xi(1))(-xi(0) - xi(1) + 2*xi(2) - 1.)
      // G(5,2) = 0.125*(1. + xi(0))(1 - xi(1))(+xi(0) - xi(1) + 2*xi(2) - 1.)
      // G(6,2) = 0.125*(1. + xi(0))(1 + xi(1))(+xi(0) + xi(1) + 2*xi(2) - 1.)
      // G(7,2) = 0.125*(1. - xi(0))(1 + xi(1))(-xi(0) + xi(1) + 2*xi(2) - 1.)
      REQUIRE(gradsf(0, 2) == Approx(-0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(1, 2) == Approx(0.).epsilon(Tolerance));
      REQUIRE(gradsf(2, 2) == Approx(0.03125).epsilon(Tolerance));
      REQUIRE(gradsf(3, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(4, 2) == Approx(-0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(5, 2) == Approx(-0.1875).epsilon(Tolerance));
      REQUIRE(gradsf(6, 2) == Approx(-0.09375).epsilon(Tolerance));
      REQUIRE(gradsf(7, 2) == Approx(-0.1875).epsilon(Tolerance));

      // Midside nodes
      // G(8, 2) = -0.25*(1 - xi(0)xi(0))(1 - xi(1));
      // G(9, 2) = -0.25*(1 - xi(1)xi(1))(1 - xi(0));
      // G(10, 2) = -0.5*xi(2)(1 - xi(0))(1 - xi(1));
      // G(11, 2) = -0.25*(1 - xi(1)xi(1))(1 + xi(0));
      // G(12, 2) = -0.5*xi(2)(1 + xi(0))(1 - xi(1));
      // G(13, 2) = -0.25*(1 - xi(0)xi(0))(1 + xi(1));
      // G(14, 2) = -0.5*xi(2)(1 + xi(0))(1 + xi(1));
      // G(15, 2) = -0.5*xi(2)(1 - xi(0))(1 + xi(1));
      // G(16, 2) = 0.25*(1 - xi(0)xi(0))(1 - xi(1));
      // G(17, 2) = 0.25*(1 - xi(1)xi(1))(1 - xi(0));
      // G(18, 2) = 0.25*(1 - xi(1)xi(1))(1 + xi(0));
      // G(19, 2) = 0.25*(1 - xi(0)xi(0))(1 + xi(1));
      REQUIRE(gradsf(8, 2) == Approx(-0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(9, 2) == Approx(-0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(10, 2) == Approx(0.5625).epsilon(Tolerance));
      REQUIRE(gradsf(11, 2) == Approx(-0.09375).epsilon(Tolerance));
      REQUIRE(gradsf(12, 2) == Approx(0.1875).epsilon(Tolerance));
      REQUIRE(gradsf(13, 2) == Approx(-0.09375).epsilon(Tolerance));
      REQUIRE(gradsf(14, 2) == Approx(0.0625).epsilon(Tolerance));
      REQUIRE(gradsf(15, 2) == Approx(0.1875).epsilon(Tolerance));
      REQUIRE(gradsf(16, 2) == Approx(0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(17, 2) == Approx(0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(18, 2) == Approx(0.09375).epsilon(Tolerance));
      REQUIRE(gradsf(19, 2) == Approx(0.09375).epsilon(Tolerance));
    }
    // Coordinates is (0.5, 0.5, 0,5)
    SECTION(
        "Twenty noded hexahedron element for coordinates(0.5, "
        "0.5, 0.5)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << 0.5, 0.5, 0.5;
      auto shapefn = hex->shapefn(coords);

      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      // Edge nodes
      // N0 = 0.125*(1 - xi(0))(1. - xi(1))(1 - xi(2))(-xi(0) - xi(1) -xi(2) -2)
      // N1 = 0.125*(1 + xi(0))(1. - xi(1))(1 - xi(2))(+xi(0) - xi(1) -xi(2) -2)
      // N2 = 0.125*(1 + xi(0))(1. + xi(1))(1 - xi(2))(+xi(0) + xi(1) -xi(2) -2)
      // N3 = 0.125*(1 - xi(0))(1. + xi(1))(1 - xi(2))(-xi(0) + xi(1) -xi(2) -2)
      // N4 = 0.125*(1 - xi(0))(1. - xi(1))(1 + xi(2))(-xi(0) - xi(1) +xi(2) -2)
      // N5 = 0.125*(1 + xi(0))(1. - xi(1))(1 + xi(2))(+xi(0) - xi(1) +xi(2) -2)
      // N6 = 0.125*(1 + xi(0))(1. + xi(1))(1 + xi(2))(+xi(0) + xi(1) +xi(2) -2)
      // N7 = 0.125*(1 - xi(0))(1. + xi(1))(1 + xi(2))(-xi(0) + xi(1) +xi(2) -2)
      REQUIRE(shapefn(0) == Approx(-0.0546875).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(-0.1171875).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(-0.2109375).epsilon(Tolerance));
      REQUIRE(shapefn(3) == Approx(-0.1171875).epsilon(Tolerance));
      REQUIRE(shapefn(4) == Approx(-0.1171875).epsilon(Tolerance));
      REQUIRE(shapefn(5) == Approx(-0.2109375).epsilon(Tolerance));
      REQUIRE(shapefn(6) == Approx(-0.2109375).epsilon(Tolerance));
      REQUIRE(shapefn(7) == Approx(-0.2109375).epsilon(Tolerance));

      // Midside nodes
      // N8 = 0.25*(1 - xi(0)^2)(1 - xi(1))(1 - xi(2))
      // N9 = 0.25*(1 - xi(0))(1 - xi(1)^2)(1 - xi(2))
      // N10 = 0.25*(1 - xi(0))(1 - xi(1))(1 - xi(2)^2)
      // N11 = 0.25*(1 + xi(0))(1 - xi(1)^2)(1 - xi(2))
      // N12 = 0.25*(1 + xi(0))(1 - xi(1))(1 - xi(2)^2)
      // N13 = 0.25*(1 - xi(0)^2)(1 + xi(1))(1 - xi(2))
      // N14 = 0.25*(1 + xi(0))(1 + xi(1))(1 - xi(2)^2)
      // N15 = 0.25*(1 - xi(0))(1 + xi(1))(1 - xi(2)^2)
      // N16 = 0.25*(1 - xi(0)^2)(1 - xi(1))(1 + xi(2))
      // N17 = 0.25*(1 - xi(0))(1 - xi(1)^2)(1 + xi(2))
      // N18 = 0.25*(1 + xi(0))(1 - xi(1)^2)(1 + xi(2))
      // N19 = 0.25*(1 - xi(0)^2)(1 + xi(1))(1 + xi(2))
      REQUIRE(shapefn(8) == Approx(0.046875).epsilon(Tolerance));
      REQUIRE(shapefn(9) == Approx(0.046875).epsilon(Tolerance));
      REQUIRE(shapefn(10) == Approx(0.046875).epsilon(Tolerance));
      REQUIRE(shapefn(11) == Approx(0.140625).epsilon(Tolerance));
      REQUIRE(shapefn(12) == Approx(0.140625).epsilon(Tolerance));
      REQUIRE(shapefn(13) == Approx(0.140625).epsilon(Tolerance));
      REQUIRE(shapefn(14) == Approx(0.421875).epsilon(Tolerance));
      REQUIRE(shapefn(15) == Approx(0.140625).epsilon(Tolerance));
      REQUIRE(shapefn(16) == Approx(0.140625).epsilon(Tolerance));
      REQUIRE(shapefn(17) == Approx(0.140625).epsilon(Tolerance));
      REQUIRE(shapefn(18) == Approx(0.421875).epsilon(Tolerance));
      REQUIRE(shapefn(19) == Approx(0.421875).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = hex->grad_shapefn(coords);
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      // Derivatives with respect to eta-direction
      // Edge nodes
      // G(0,0) = 0.125*(1. - xi(1))(1 - xi(2))(2*xi(0) + xi(1) + xi(2) + 1.)
      // G(1,0) = 0.125*(1. - xi(1))(1 - xi(2))(2*xi(0) - xi(1) - xi(2) - 1.)
      // G(2,0) = 0.125*(1. + xi(1))(1 - xi(2))(2*xi(0) + xi(1) - xi(2) - 1.)
      // G(3,0) = 0.125*(1. + xi(1))(1 - xi(2))(2*xi(0) - xi(1) + xi(2) + 1.)
      // G(4,0) = 0.125*(1. - xi(1))(1 + xi(2))(2*xi(0) + xi(1) - xi(2) + 1.)
      // G(5,0) = 0.125*(1. - xi(1))(1 + xi(2))(2*xi(0) - xi(1) + xi(2) - 1.)
      // G(6,0) = 0.125*(1. + xi(1))(1 + xi(2))(2*xi(0) + xi(1) + xi(2) - 1.)
      // G(7,0) = 0.125*(1. + xi(1))(1 + xi(2))(2*xi(0) - xi(1) - xi(2) + 1.)
      REQUIRE(gradsf(0, 0) == Approx(0.09375).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(-0.03125).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(3, 0) == Approx(0.1875).epsilon(Tolerance));
      REQUIRE(gradsf(4, 0) == Approx(0.1875).epsilon(Tolerance));
      REQUIRE(gradsf(5, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(6, 0) == Approx(0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(7, 0) == Approx(0.28125).epsilon(Tolerance));

      // Midside nodes
      // G(8, 0) = -0.5*xi(0)*(1 - xi(1))(1 - xi(2));
      // G(9, 0) = -0.25*(1 - xi(1)xi(1))(1 - xi(2));
      // G(10, 0) = -0.25*(1 - xi(2)xi(2))(1 - xi(1));
      // G(11, 0) = 0.25*(1 - xi(1)xi(1))(1 - xi(2));
      // G(12, 0) = 0.25*(1 - xi(2)xi(2))(1 - xi(1));
      // G(13, 0) = -0.5*xi(0)(1 + xi(1))(1 - xi(2));
      // G(14, 0) = 0.25*(1 - xi(2)xi(2))(1 + xi(1));
      // G(15, 0) = -0.25*(1 - xi(2)xi(2))(1 + xi(1));
      // G(16, 0) = -0.5*xi(0)(1 - xi(1))(1 + xi(2));
      // G(17, 0) = -0.25*(1 - xi(1)xi(1))(1 + xi(2));
      // G(18, 0) = 0.25*(1 - xi(1)xi(1))(1 + xi(2));
      // G(19, 0) = -0.5*xi(0)(1 + xi(1))(1 + xi(2));
      REQUIRE(gradsf(8, 0) == Approx(-0.0625).epsilon(Tolerance));
      REQUIRE(gradsf(9, 0) == Approx(-0.09375).epsilon(Tolerance));
      REQUIRE(gradsf(10, 0) == Approx(-0.09375).epsilon(Tolerance));
      REQUIRE(gradsf(11, 0) == Approx(0.09375).epsilon(Tolerance));
      REQUIRE(gradsf(12, 0) == Approx(0.09375).epsilon(Tolerance));
      REQUIRE(gradsf(13, 0) == Approx(-0.1875).epsilon(Tolerance));
      REQUIRE(gradsf(14, 0) == Approx(0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(15, 0) == Approx(-0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(16, 0) == Approx(-0.1875).epsilon(Tolerance));
      REQUIRE(gradsf(17, 0) == Approx(-0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(18, 0) == Approx(0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(19, 0) == Approx(-0.5625).epsilon(Tolerance));

      // Derivatives with respect to psi-direction
      // Edge nodes
      // G(0,1) = 0.125*(1. - xi(0))(1 - xi(2))(+xi(0) + 2*xi(1) + xi(2) + 1.)
      // G(1,1) = 0.125*(1. + xi(0))(1 - xi(2))(-xi(0) + 2*xi(1) + xi(2) + 1.)
      // G(2,1) = 0.125*(1. + xi(0))(1 - xi(2))(+xi(0) + 2*xi(1) - xi(2) - 1.)
      // G(3,1) = 0.125*(1. - xi(0))(1 - xi(2))(-xi(0) + 2*xi(1) - xi(2) - 1.)
      // G(4,1) = 0.125*(1. - xi(0))(1 + xi(2))(+xi(0) + 2*xi(1) - xi(2) + 1.)
      // G(5,1) = 0.125*(1. + xi(0))(1 + xi(2))(-xi(0) + 2*xi(1) - xi(2) + 1.)
      // G(6,1) = 0.125*(1. + xi(0))(1 + xi(2))(+xi(0) + 2*xi(1) + xi(2) - 1.)
      // G(7,1) = 0.125*(1. - xi(0))(1 + xi(2))(-xi(0) + 2*xi(1) + xi(2) - 1.)
      REQUIRE(gradsf(0, 1) == Approx(0.09375).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(0.1875).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(3, 1) == Approx(-0.03125).epsilon(Tolerance));
      REQUIRE(gradsf(4, 1) == Approx(0.1875).epsilon(Tolerance));
      REQUIRE(gradsf(5, 1) == Approx(0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(6, 1) == Approx(0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(7, 1) == Approx(0.0).epsilon(Tolerance));

      // Midside nodes
      // G(8, 1) = -0.25*(1 - xi(0)xi(0))(1 - xi(2));
      // G(9, 1) = -0.5*xi(1)(1 - xi(0))(1 - xi(2));
      // G(10, 1) = -0.25*(1 - xi(2)xi(2))(1 - xi(0));
      // G(11, 1) = -0.5*xi(1)(1 + xi(0))(1 - xi(2));
      // G(12, 1) = -0.25*(1 - xi(2)xi(2))(1 + xi(0));
      // G(13, 1) = 0.25*(1 - xi(0)xi(0))(1 - xi(2));
      // G(14, 1) = 0.25*(1 - xi(2)xi(2))(1 + xi(0));
      // G(15, 1) = 0.25*(1 - xi(2)xi(2))(1 - xi(0));
      // G(16, 1) = -0.25*(1 - xi(0)xi(0))(1 + xi(2));
      // G(17, 1) = -0.5*xi(1)(1 - xi(0))(1 + xi(2));
      // G(18, 1) = -0.5*xi(1)(1 + xi(0))(1 + xi(2));
      // G(19, 1) = 0.25*(1 - xi(0)xi(0))(1 + xi(2));
      REQUIRE(gradsf(8, 1) == Approx(-0.09375).epsilon(Tolerance));
      REQUIRE(gradsf(9, 1) == Approx(-0.0625).epsilon(Tolerance));
      REQUIRE(gradsf(10, 1) == Approx(-0.09375).epsilon(Tolerance));
      REQUIRE(gradsf(11, 1) == Approx(-0.1875).epsilon(Tolerance));
      REQUIRE(gradsf(12, 1) == Approx(-0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(13, 1) == Approx(0.09375).epsilon(Tolerance));
      REQUIRE(gradsf(14, 1) == Approx(0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(15, 1) == Approx(0.09375).epsilon(Tolerance));
      REQUIRE(gradsf(16, 1) == Approx(-0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(17, 1) == Approx(-0.1875).epsilon(Tolerance));
      REQUIRE(gradsf(18, 1) == Approx(-0.5625).epsilon(Tolerance));
      REQUIRE(gradsf(19, 1) == Approx(0.28125).epsilon(Tolerance));

      // Derivatives with respect to mu-direction
      // Edge nodes
      // G(0,2) = 0.125*(1. - xi(0))(1 - xi(1))(+xi(0) + xi(1) + 2*xi(2) + 1.)
      // G(1,2) = 0.125*(1. + xi(0))(1 - xi(1))(-xi(0) + xi(1) + 2*xi(2) + 1.)
      // G(2,2) = 0.125*(1. + xi(0))(1 + xi(1))(-xi(0) - xi(1) + 2*xi(2) + 1.)
      // G(3,2) = 0.125*(1. - xi(0))(1 + xi(1))(+xi(0) - xi(1) + 2*xi(2) + 1.)
      // G(4,2) = 0.125*(1. - xi(0))(1 - xi(1))(-xi(0) - xi(1) + 2*xi(2) - 1.)
      // G(5,2) = 0.125*(1. + xi(0))(1 - xi(1))(+xi(0) - xi(1) + 2*xi(2) - 1.)
      // G(6,2) = 0.125*(1. + xi(0))(1 + xi(1))(+xi(0) + xi(1) + 2*xi(2) - 1.)
      // G(7,2) = 0.125*(1. - xi(0))(1 + xi(1))(-xi(0) + xi(1) + 2*xi(2) - 1.)
      REQUIRE(gradsf(0, 2) == Approx(0.09375).epsilon(Tolerance));
      REQUIRE(gradsf(1, 2) == Approx(0.1875).epsilon(Tolerance));
      REQUIRE(gradsf(2, 2) == Approx(0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(3, 2) == Approx(0.1875).epsilon(Tolerance));
      REQUIRE(gradsf(4, 2) == Approx(-0.03125).epsilon(Tolerance));
      REQUIRE(gradsf(5, 2) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(6, 2) == Approx(0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(7, 2) == Approx(0.0).epsilon(Tolerance));

      // Midside nodes
      // G(8, 2) = -0.25*(1 - xi(0)xi(0))(1 - xi(1));
      // G(9, 2) = -0.25*(1 - xi(1)xi(1))(1 - xi(0));
      // G(10, 2) = -0.5*xi(2)(1 - xi(0))(1 - xi(1));
      // G(11, 2) = -0.25*(1 - xi(1)xi(1))(1 + xi(0));
      // G(12, 2) = -0.5*xi(2)(1 + xi(0))(1 - xi(1));
      // G(13, 2) = -0.25*(1 - xi(0)xi(0))(1 + xi(1));
      // G(14, 2) = -0.5*xi(2)(1 + xi(0))(1 + xi(1));
      // G(15, 2) = -0.5*xi(2)(1 - xi(0))(1 + xi(1));
      // G(16, 2) = 0.25*(1 - xi(0)xi(0))(1 - xi(1));
      // G(17, 2) = 0.25*(1 - xi(1)xi(1))(1 - xi(0));
      // G(18, 2) = 0.25*(1 - xi(1)xi(1))(1 + xi(0));
      // G(19, 2) = 0.25*(1 - xi(0)xi(0))(1 + xi(1));
      REQUIRE(gradsf(8, 2) == Approx(-0.09375).epsilon(Tolerance));
      REQUIRE(gradsf(9, 2) == Approx(-0.09375).epsilon(Tolerance));
      REQUIRE(gradsf(10, 2) == Approx(-0.0625).epsilon(Tolerance));
      REQUIRE(gradsf(11, 2) == Approx(-0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(12, 2) == Approx(-0.1875).epsilon(Tolerance));
      REQUIRE(gradsf(13, 2) == Approx(-0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(14, 2) == Approx(-0.5625).epsilon(Tolerance));
      REQUIRE(gradsf(15, 2) == Approx(-0.1875).epsilon(Tolerance));
      REQUIRE(gradsf(16, 2) == Approx(0.09375).epsilon(Tolerance));
      REQUIRE(gradsf(17, 2) == Approx(0.09375).epsilon(Tolerance));
      REQUIRE(gradsf(18, 2) == Approx(0.28125).epsilon(Tolerance));
      REQUIRE(gradsf(19, 2) == Approx(0.28125).epsilon(Tolerance));
    }

    SECTION("Twenty noded hexahedron element with grad deformation") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();
      Eigen::Matrix<double, Dim, 1> psize;
      psize.setZero();
      Eigen::Matrix<double, Dim, 1> defgrad;
      defgrad.setZero();

      auto shapefn = hex->shapefn(coords, psize, defgrad);

      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      // Edge nodes
      REQUIRE(shapefn(0) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(shapefn(3) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(shapefn(4) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(shapefn(5) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(shapefn(6) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(shapefn(7) == Approx(-0.25).epsilon(Tolerance));

      // Midside nodes
      REQUIRE(shapefn(8) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(9) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(10) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(11) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(12) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(13) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(14) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(15) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(16) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(17) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(18) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(19) == Approx(0.25).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = hex->grad_shapefn(coords, psize, defgrad);
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      // Derivatives with respect to eta-direction
      // Edge nodes
      REQUIRE(gradsf(0, 0) == Approx(+0.125).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(3, 0) == Approx(+0.125).epsilon(Tolerance));
      REQUIRE(gradsf(4, 0) == Approx(+0.125).epsilon(Tolerance));
      REQUIRE(gradsf(5, 0) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(6, 0) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(7, 0) == Approx(+0.125).epsilon(Tolerance));

      // Midside nodes
      REQUIRE(gradsf(8, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(9, 0) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(10, 0) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(11, 0) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(gradsf(12, 0) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(gradsf(13, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(14, 0) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(gradsf(15, 0) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(16, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(17, 0) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(18, 0) == Approx(+0.25).epsilon(Tolerance));
      REQUIRE(gradsf(19, 0) == Approx(0.0).epsilon(Tolerance));

      // Derivatives with respect to psi-direction
      // Edge nodes
      REQUIRE(gradsf(0, 1) == Approx(+0.125).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(+0.125).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(3, 1) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(4, 1) == Approx(+0.125).epsilon(Tolerance));
      REQUIRE(gradsf(5, 1) == Approx(+0.125).epsilon(Tolerance));
      REQUIRE(gradsf(6, 1) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(7, 1) == Approx(-0.125).epsilon(Tolerance));

      // Midside nodes
      REQUIRE(gradsf(8, 1) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(9, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(10, 1) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(11, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(12, 1) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(13, 1) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(gradsf(14, 1) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(gradsf(15, 1) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(gradsf(16, 1) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(17, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(18, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(19, 1) == Approx(0.25).epsilon(Tolerance));

      // Derivatives with respect to mu-direction
      // Edge nodes
      REQUIRE(gradsf(0, 2) == Approx(+0.125).epsilon(Tolerance));
      REQUIRE(gradsf(1, 2) == Approx(+0.125).epsilon(Tolerance));
      REQUIRE(gradsf(2, 2) == Approx(+0.125).epsilon(Tolerance));
      REQUIRE(gradsf(3, 2) == Approx(+0.125).epsilon(Tolerance));
      REQUIRE(gradsf(4, 2) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(5, 2) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(6, 2) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(7, 2) == Approx(-0.125).epsilon(Tolerance));

      // Midside nodes
      REQUIRE(gradsf(8, 2) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(9, 2) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(10, 2) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(11, 2) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(12, 2) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(13, 2) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(14, 2) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(15, 2) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(16, 2) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(gradsf(17, 2) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(gradsf(18, 2) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(gradsf(19, 2) == Approx(0.25).epsilon(Tolerance));
    }

    // Check Jacobian
    SECTION("20 noded hexahedron Jacobian for local coordinates(0.5,0.5,0.5)") {
      Eigen::Matrix<double, 20, Dim> coords;
      // clang-format off
      coords << 2.0, 1.0, 0.50,
                4.0, 2.0, 1.00,
                2.0, 4.0, 1.00,
                1.0, 3.0, 0.50,
                2.0, 1.0, 1.50,
                4.0, 2.0, 2.00,
                2.0, 4.0, 2.00,
                1.0, 3.0, 1.50,
                3.0, 1.5, 0.75,
                1.5, 2.0, 0.50,
                2.0, 1.0, 1.00,
                3.0, 3.0, 1.00,
                4.0, 2.0, 1.50,
                1.5, 3.5, 0.75,
                2.0, 4.0, 1.50,
                1.0, 3.0, 1.00,
                3.0, 1.5, 1.75,
                1.5, 2.0, 1.50,
                4.0, 2.0, 1.50,
                1.5, 3.5, 1.75;
      // clang-format on

      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0.5, 0.5, 0.5;

      // Jacobian result
      Eigen::Matrix<double, Dim, Dim> jacobian;
      // clang-format off
      jacobian << 0.90625, 0.21875, 0.109375,
                 -1.43750, 1.56250, 0.281250,
                  0.28125,-0.28125, 0.359375;
      // clang-format on

      // Get Jacobian
      auto jac = hex->jacobian(xi, coords);

      // Check size of jacobian
      REQUIRE(jac.size() == jacobian.size());

      // Check Jacobian
      for (unsigned i = 0; i < Dim; ++i)
        for (unsigned j = 0; j < Dim; ++j)
          REQUIRE(jac(i, j) == Approx(jacobian(i, j)).epsilon(Tolerance));
    }

    // Check Jacobian with deformation gradient
    SECTION("20 noded hexahedron Jacobian with deformation gradient") {
      Eigen::Matrix<double, 20, Dim> coords;
      // clang-format off
      coords << 2.0, 1.0, 0.50,
                4.0, 2.0, 1.00,
                2.0, 4.0, 1.00,
                1.0, 3.0, 0.50,
                2.0, 1.0, 1.50,
                4.0, 2.0, 2.00,
                2.0, 4.0, 2.00,
                1.0, 3.0, 1.50,
                3.0, 1.5, 0.75,
                1.5, 2.0, 0.50,
                2.0, 1.0, 1.00,
                3.0, 3.0, 1.00,
                4.0, 2.0, 1.50,
                1.5, 3.5, 0.75,
                2.0, 4.0, 1.50,
                1.0, 3.0, 1.00,
                3.0, 1.5, 1.75,
                1.5, 2.0, 1.50,
                4.0, 2.0, 1.50,
                1.5, 3.5, 1.75;
      // clang-format on

      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0.5, 0.5, 0.5;

      Eigen::Matrix<double, Dim, 1> psize;
      psize << 0.25, 0.5, 0.75;
      Eigen::Matrix<double, Dim, 1> defgrad;
      defgrad.setZero();

      // Jacobian result
      Eigen::Matrix<double, Dim, Dim> jacobian;
      // clang-format off
      jacobian << 0.90625, 0.21875, 0.109375,
                 -1.43750, 1.56250, 0.281250,
                  0.28125,-0.28125, 0.359375;
      // clang-format on

      // Get Jacobian
      auto jac = hex->jacobian(xi, coords, psize, defgrad);

      // Check size of jacobian
      REQUIRE(jac.size() == jacobian.size());

      // Check Jacobian
      for (unsigned i = 0; i < Dim; ++i)
        for (unsigned j = 0; j < Dim; ++j)
          REQUIRE(jac(i, j) == Approx(jacobian(i, j)).epsilon(Tolerance));
    }

    // Coordinates is (0, 0, 0)
    SECTION("20-noded hexahedron B-matrix for coordinates(0, 0, 0)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();
      // Get B-Matrix
      auto bmatrix = hex->bmatrix(coords);

      // Check gradient of shape functions
      auto gradsf = hex->grad_shapefn(coords);

      // Check size of B-matrix
      REQUIRE(bmatrix.size() == nfunctions);

      for (unsigned i = 0; i < nfunctions; ++i) {
        REQUIRE(bmatrix.at(i)(0, 0) == Approx(gradsf(i, 0)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(0, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(0, 2) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 1) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 2) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 2) == Approx(gradsf(i, 2)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(3, 0) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(3, 1) == Approx(gradsf(i, 0)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(3, 2) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(4, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(4, 1) == Approx(gradsf(i, 2)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(4, 2) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(5, 0) == Approx(gradsf(i, 2)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(5, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(5, 2) == Approx(gradsf(i, 0)).epsilon(Tolerance));
      }
    }

    // Coordinates is (0.5, 0.5, 0.5)
    SECTION("20-noded hexahedron B-matrix for coordinates(0.5, 0.5, 0.5)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << 0.5, 0.5, 0.5;

      // Get B-Matrix
      auto bmatrix = hex->bmatrix(coords);

      // Check gradient of shape functions
      auto gradsf = hex->grad_shapefn(coords);

      // Check size of B-matrix
      REQUIRE(bmatrix.size() == nfunctions);

      for (unsigned i = 0; i < nfunctions; ++i) {
        REQUIRE(bmatrix.at(i)(0, 0) == Approx(gradsf(i, 0)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(0, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(0, 2) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 1) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 2) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 2) == Approx(gradsf(i, 2)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(3, 0) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(3, 1) == Approx(gradsf(i, 0)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(3, 2) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(4, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(4, 1) == Approx(gradsf(i, 2)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(4, 2) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(5, 0) == Approx(gradsf(i, 2)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(5, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(5, 2) == Approx(gradsf(i, 0)).epsilon(Tolerance));
      }
    }

    // Coordinates is (-0.5, -0.5, -0.5)
    SECTION("20-noded hexahedron B-matrix for coordinates(-0.5, -0.5, -0.5)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << -0.5, -0.5, -0.5;

      // Get B-Matrix
      auto bmatrix = hex->bmatrix(coords);

      // Check gradient of shape functions
      auto gradsf = hex->grad_shapefn(coords);

      // Check size of B-matrix
      REQUIRE(bmatrix.size() == nfunctions);

      for (unsigned i = 0; i < nfunctions; ++i) {
        REQUIRE(bmatrix.at(i)(0, 0) == Approx(gradsf(i, 0)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(0, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(0, 2) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 1) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 2) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 2) == Approx(gradsf(i, 2)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(3, 0) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(3, 1) == Approx(gradsf(i, 0)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(3, 2) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(4, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(4, 1) == Approx(gradsf(i, 2)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(4, 2) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(5, 0) == Approx(gradsf(i, 2)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(5, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(5, 2) == Approx(gradsf(i, 0)).epsilon(Tolerance));
      }
    }

    SECTION("20-noded hexahedron B-matrix and Jacobian failure") {
      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0., 0., 0.;

      Eigen::Matrix<double, 7, Dim> coords;
      // clang-format off
      coords << 0., 0., 0.,
                1., 0., 0.,
                1., 1., 0.,
                0., 1., 0.,
                0., 0., 1.,
                1., 0., 1.,
                1., 1., 1.;
      // clang-format on
      // Get B-Matrix
      auto bmatrix = hex->bmatrix(xi, coords);
      auto jacobian = hex->jacobian(xi, coords);
    }

    // Mass matrix of a cell
    SECTION("20-noded hexahedron mass-matrix") {
      std::vector<Eigen::Matrix<double, Dim, 1>> xi_s;

      Eigen::Matrix<double, Dim, 1> xi;
      const double one_by_sqrt3 = std::fabs(1 / std::sqrt(3));
      xi << -one_by_sqrt3, -one_by_sqrt3, -one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << one_by_sqrt3, -one_by_sqrt3, -one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << one_by_sqrt3, one_by_sqrt3, -one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << -one_by_sqrt3, one_by_sqrt3, -one_by_sqrt3;
      xi_s.emplace_back(xi);

      xi << -one_by_sqrt3, -one_by_sqrt3, one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << one_by_sqrt3, -one_by_sqrt3, one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << one_by_sqrt3, one_by_sqrt3, one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << -one_by_sqrt3, one_by_sqrt3, one_by_sqrt3;
      xi_s.emplace_back(xi);

      REQUIRE(xi_s.size() == 8);

      // Get mass matrix
      const auto mass_matrix = hex->mass_matrix(xi_s);

      // Check size of mass-matrix
      REQUIRE(mass_matrix.rows() == nfunctions);
      REQUIRE(mass_matrix.cols() == nfunctions);

      // Sum should be equal to 1. * xi_s.size()
      REQUIRE(mass_matrix.sum() == Approx(1. * xi_s.size()).epsilon(Tolerance));

      Eigen::Matrix<double, 20, 20> mass;
      //! clang-format off
      mass << 0.1481481481481481, 0.1234567901234567, 0.1234567901234567,
          0.1234567901234567, 0.1234567901234567, 0.1234567901234567,
          0.111111111111111, 0.1234567901234567, -0.1975308641975307,
          -0.1975308641975307, -0.1975308641975307, -0.1728395061728394,
          -0.1728395061728394, -0.1728395061728394, -0.1234567901234567,
          -0.1728395061728394, -0.1728395061728394, -0.1728395061728394,
          -0.1234567901234567, -0.1234567901234567, 0.1234567901234567,
          0.1481481481481481, 0.1234567901234567, 0.1234567901234567,
          0.1234567901234567, 0.1234567901234567, 0.1234567901234567,
          0.111111111111111, -0.1975308641975307, -0.1728395061728394,
          -0.1728395061728394, -0.1975308641975307, -0.1975308641975308,
          -0.1728395061728394, -0.1728395061728394, -0.1234567901234567,
          -0.1728395061728394, -0.1234567901234567, -0.1728395061728394,
          -0.1234567901234567, 0.1234567901234567, 0.1234567901234567,
          0.1481481481481481, 0.1234567901234567, 0.111111111111111,
          0.1234567901234567, 0.1234567901234567, 0.1234567901234567,
          -0.1728395061728394, -0.1728395061728394, -0.1234567901234567,
          -0.1975308641975307, -0.1728395061728394, -0.1975308641975307,
          -0.1975308641975307, -0.1728395061728394, -0.1234567901234567,
          -0.1234567901234567, -0.1728395061728394, -0.1728395061728394,
          0.1234567901234567, 0.1234567901234567, 0.1234567901234567,
          0.1481481481481481, 0.1234567901234567, 0.111111111111111,
          0.1234567901234567, 0.1234567901234567, -0.1728395061728394,
          -0.1975308641975307, -0.1728395061728394, -0.1728395061728394,
          -0.1234567901234567, -0.1975308641975307, -0.1728395061728394,
          -0.1975308641975308, -0.1234567901234567, -0.1728395061728394,
          -0.1234567901234567, -0.1728395061728394, 0.1234567901234567,
          0.1234567901234567, 0.111111111111111, 0.1234567901234567,
          0.1481481481481481, 0.1234567901234567, 0.1234567901234567,
          0.1234567901234567, -0.1728395061728394, -0.1728395061728394,
          -0.1975308641975308, -0.1234567901234567, -0.1728395061728394,
          -0.1234567901234567, -0.1234567901234567, -0.1728395061728394,
          -0.1975308641975307, -0.1975308641975307, -0.1728395061728394,
          -0.1728395061728394, 0.1234567901234567, 0.1234567901234567,
          0.1234567901234567, 0.111111111111111, 0.1234567901234567,
          0.1481481481481481, 0.1234567901234567, 0.1234567901234567,
          -0.1728395061728394, -0.1234567901234567, -0.1728395061728394,
          -0.1728395061728394, -0.1975308641975308, -0.1234567901234567,
          -0.1728395061728394, -0.1234567901234567, -0.1975308641975307,
          -0.1728395061728394, -0.1975308641975307, -0.1728395061728394,
          0.111111111111111, 0.1234567901234567, 0.1234567901234567,
          0.1234567901234567, 0.1234567901234567, 0.1234567901234567,
          0.1481481481481481, 0.1234567901234567, -0.1234567901234567,
          -0.1234567901234567, -0.1234567901234567, -0.1728395061728394,
          -0.1728395061728394, -0.1728395061728394, -0.1975308641975307,
          -0.1728395061728394, -0.1728395061728394, -0.1728395061728394,
          -0.1975308641975307, -0.1975308641975307, 0.1234567901234567,
          0.111111111111111, 0.1234567901234567, 0.1234567901234567,
          0.1234567901234567, 0.1234567901234567, 0.1234567901234567,
          0.1481481481481481, -0.1234567901234567, -0.1728395061728394,
          -0.1728395061728394, -0.1234567901234567, -0.1234567901234567,
          -0.1728395061728394, -0.1728395061728394, -0.1975308641975307,
          -0.1728395061728394, -0.1975308641975307, -0.1728395061728394,
          -0.1975308641975307, -0.1975308641975307, -0.1975308641975307,
          -0.1728395061728394, -0.1728395061728394, -0.1728395061728394,
          -0.1728395061728394, -0.1234567901234567, -0.1234567901234567,
          0.3950617283950616, 0.2962962962962961, 0.2962962962962961,
          0.2962962962962961, 0.2962962962962961, 0.1975308641975307,
          0.148148148148148, 0.148148148148148, 0.1975308641975307,
          0.148148148148148, 0.148148148148148, 0.09876543209876533,
          -0.1975308641975307, -0.1728395061728394, -0.1728395061728394,
          -0.1975308641975307, -0.1728395061728394, -0.1234567901234567,
          -0.1234567901234567, -0.1728395061728394, 0.2962962962962961,
          0.3950617283950616, 0.2962962962962961, 0.1975308641975307,
          0.148148148148148, 0.2962962962962961, 0.148148148148148,
          0.2962962962962961, 0.148148148148148, 0.1975308641975307,
          0.09876543209876534, 0.148148148148148, -0.1975308641975307,
          -0.1728395061728394, -0.1234567901234567, -0.1728395061728394,
          -0.1975308641975308, -0.1728395061728394, -0.1234567901234567,
          -0.1728395061728394, 0.2962962962962961, 0.2962962962962961,
          0.3950617283950616, 0.148148148148148, 0.1975308641975307,
          0.148148148148148, 0.09876543209876534, 0.1975308641975307,
          0.2962962962962961, 0.2962962962962961, 0.148148148148148,
          0.148148148148148, -0.1728395061728394, -0.1975308641975307,
          -0.1975308641975307, -0.1728395061728394, -0.1234567901234567,
          -0.1728395061728394, -0.1728395061728394, -0.1234567901234567,
          0.2962962962962961, 0.1975308641975307, 0.148148148148148,
          0.3950617283950616, 0.2962962962962961, 0.2962962962962961,
          0.2962962962962961, 0.148148148148148, 0.148148148148148,
          0.09876543209876534, 0.1975308641975307, 0.148148148148148,
          -0.1728395061728394, -0.1975308641975308, -0.1728395061728394,
          -0.1234567901234567, -0.1728395061728394, -0.1975308641975308,
          -0.1728395061728394, -0.1234567901234567, 0.2962962962962961,
          0.148148148148148, 0.1975308641975307, 0.2962962962962961,
          0.3950617283950615, 0.148148148148148, 0.1975308641975307,
          0.09876543209876534, 0.2962962962962961, 0.148148148148148,
          0.2962962962962961, 0.148148148148148, -0.1728395061728394,
          -0.1728395061728394, -0.1975308641975307, -0.1975308641975307,
          -0.1234567901234567, -0.1234567901234567, -0.1728395061728394,
          -0.1728395061728394, 0.1975308641975307, 0.2962962962962961,
          0.148148148148148, 0.2962962962962961, 0.148148148148148,
          0.3950617283950616, 0.2962962962962961, 0.2962962962962961,
          0.09876543209876534, 0.148148148148148, 0.148148148148148,
          0.1975308641975307, -0.1234567901234567, -0.1728395061728394,
          -0.1975308641975307, -0.1728395061728394, -0.1234567901234567,
          -0.1728395061728394, -0.1975308641975307, -0.1728395061728394,
          0.148148148148148, 0.148148148148148, 0.09876543209876534,
          0.2962962962962961, 0.1975308641975307, 0.2962962962962961,
          0.3950617283950615, 0.1975308641975307, 0.148148148148148,
          0.148148148148148, 0.2962962962962961, 0.2962962962962961,
          -0.1728395061728394, -0.1234567901234567, -0.1728395061728394,
          -0.1975308641975308, -0.1728395061728394, -0.1234567901234567,
          -0.1728395061728394, -0.1975308641975307, 0.148148148148148,
          0.2962962962962961, 0.1975308641975307, 0.148148148148148,
          0.09876543209876534, 0.2962962962962961, 0.1975308641975307,
          0.3950617283950615, 0.148148148148148, 0.2962962962962961,
          0.148148148148148, 0.2962962962962961, -0.1728395061728394,
          -0.1728395061728394, -0.1234567901234567, -0.1234567901234567,
          -0.1975308641975307, -0.1975308641975307, -0.1728395061728394,
          -0.1728395061728394, 0.1975308641975307, 0.148148148148148,
          0.2962962962962961, 0.148148148148148, 0.2962962962962961,
          0.09876543209876534, 0.148148148148148, 0.148148148148148,
          0.3950617283950615, 0.2962962962962961, 0.2962962962962961,
          0.1975308641975307, -0.1728395061728394, -0.1234567901234567,
          -0.1234567901234567, -0.1728395061728394, -0.1975308641975307,
          -0.1728395061728394, -0.1728395061728394, -0.1975308641975307,
          0.148148148148148, 0.1975308641975307, 0.2962962962962961,
          0.09876543209876534, 0.148148148148148, 0.148148148148148,
          0.148148148148148, 0.2962962962962961, 0.2962962962962961,
          0.3950617283950615, 0.1975308641975307, 0.2962962962962961,
          -0.1234567901234567, -0.1728395061728394, -0.1728395061728394,
          -0.1234567901234567, -0.1728395061728394, -0.1975308641975307,
          -0.1975308641975307, -0.1728395061728394, 0.148148148148148,
          0.09876543209876534, 0.148148148148148, 0.1975308641975307,
          0.2962962962962961, 0.148148148148148, 0.2962962962962961,
          0.148148148148148, 0.2962962962962961, 0.1975308641975307,
          0.3950617283950615, 0.2962962962962961, -0.1234567901234567,
          -0.1234567901234567, -0.1728395061728394, -0.1728395061728394,
          -0.1728395061728394, -0.1728395061728394, -0.1975308641975307,
          -0.1975308641975307, 0.09876543209876533, 0.148148148148148,
          0.148148148148148, 0.148148148148148, 0.148148148148148,
          0.1975308641975307, 0.2962962962962961, 0.2962962962962961,
          0.1975308641975307, 0.2962962962962961, 0.2962962962962961,
          0.3950617283950615;
      // clang-format on
      for (unsigned i = 0; i < nfunctions; ++i)
        for (unsigned j = 0; j < nfunctions; ++j)
          REQUIRE(mass_matrix(i, j) == Approx(mass(i, j)).epsilon(Tolerance));
    }

    // Laplace matrix of a cell
    SECTION("20-noded hexahedron laplace-matrix") {
      std::vector<Eigen::Matrix<double, Dim, 1>> xi_s;

      Eigen::Matrix<double, Dim, 1> xi;
      const double one_by_sqrt3 = std::fabs(1 / std::sqrt(3));
      xi << -one_by_sqrt3, -one_by_sqrt3, -one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << one_by_sqrt3, -one_by_sqrt3, -one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << one_by_sqrt3, one_by_sqrt3, -one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << -one_by_sqrt3, one_by_sqrt3, -one_by_sqrt3;
      xi_s.emplace_back(xi);

      xi << -one_by_sqrt3, -one_by_sqrt3, one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << one_by_sqrt3, -one_by_sqrt3, one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << one_by_sqrt3, one_by_sqrt3, one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << -one_by_sqrt3, one_by_sqrt3, one_by_sqrt3;
      xi_s.emplace_back(xi);

      REQUIRE(xi_s.size() == 8);

      // Nodal coordinates
      Eigen::Matrix<double, 20, Dim> coords;
      // clang-format off
      coords << 2.0, 1.0, 0.50,
                4.0, 2.0, 1.00,
                2.0, 4.0, 1.00,
                1.0, 3.0, 0.50,
                2.0, 1.0, 1.50,
                4.0, 2.0, 2.00,
                2.0, 4.0, 2.00,
                1.0, 3.0, 1.50,
                3.0, 1.5, 0.75,
                1.5, 2.0, 0.50,
                2.0, 1.0, 1.00,
                3.0, 3.0, 1.00,
                4.0, 2.0, 1.50,
                1.5, 3.5, 0.75,
                2.0, 4.0, 1.50,
                1.0, 3.0, 1.00,
                3.0, 1.5, 1.75,
                1.5, 2.0, 1.50,
                4.0, 2.0, 1.50,
                1.5, 3.5, 1.75;
      // clang-format on

      // Get laplace matrix
      const auto laplace_matrix = hex->laplace_matrix(xi_s, coords);

      // Check size of laplace-matrix
      REQUIRE(laplace_matrix.rows() == nfunctions);
      REQUIRE(laplace_matrix.cols() == nfunctions);

      Eigen::Matrix<double, 20, 20> laplace;
      laplace << 6.15933929330602, 4.938676360208379, 4.058798842761456,
          3.776996153514026, 4.976923615044938, 4.348568809939851,
          3.907765270436546, 5.794922003041152, -2.378685014427916,
          -5.536048018007822, -8.84093622113221, -1.100857511792741,
          -5.522691796058115, 0.01790980108445428, -6.030343075250605,
          -14.68581000074875, -1.224340104438079, -0.005114463255402324,
          -2.234642701444805, -1.356471975117627, 4.938676360208379,
          11.21484802784489, 4.081716008429723, 3.42897599003046,
          5.025879299946409, 8.170918200545696, 4.95877517199694,
          7.782823613402732, -2.858753754386522, -3.612270433156081,
          -7.501753263993352, 2.264122966861815, -12.42024109318314,
          -2.014924579198793, -6.892202977180818, -7.611220041618338,
          -3.049877388674338, -0.06047605860149874, -4.571120776660473,
          -1.65973684528307, 4.058798842761456, 4.081716008429723,
          8.528392031420996, 4.523672525376581, 0.6731108109138284,
          4.001455589844431, 4.158805328105587, 5.922297826135225,
          -3.128651652280993, -5.056094482013514, -3.383579335130428,
          -3.449854529765026, -4.509018906055925, -6.611884017194679,
          -9.161348600196341, -11.07427442942752, 1.994885866560606,
          5.019500172357599, -1.507880835141334, 0.2231785103095052,
          3.776996153514026, 3.42897599003046, 4.523672525376581,
          5.659178156732192, 2.059164452549481, 4.760061473171758,
          3.323047914967523, 4.366168456560783, -2.87529736009937,
          -4.150936493636122, -4.089638778195734, -4.496272990998748,
          -3.34471447300587, -1.105338197190417, -5.574589174004736,
          -8.273216409082542, 2.850808750123276, -0.4187973830786254,
          -1.62047095951107, -1.289711081566152, 4.976923615044938,
          5.025879299946409, 0.6731108109138284, 2.059164452549481,
          7.596193475939065, 4.5620913238822, 3.835574427742611,
          4.393331161268728, -1.092981827654348, -3.343689344221354,
          -9.630640335753997, 2.429841121276548, -6.215116692494317,
          6.644839781301256, -3.223741811124502, -10.90764512285358,
          -4.114786139361408, -5.749284093800984, -3.564920648690664,
          -3.832477684274747, 4.348568809939851, 8.170918200545696,
          4.001455589844431, 4.760061473171758, 4.5620913238822,
          10.27495321916873, 5.094648956715952, 7.432212216606898,
          -2.486249474018963, -2.437380936216606, -6.718037799849757,
          -0.2648377020835075, -11.32726336676532, 1.262993685950383,
          -6.298130156830595, -6.311463206932629, -0.1224764752776564,
          -3.137369670859187, -5.428028023658495, -3.368706027054427,
          3.907765270436546, 4.95877517199694, 4.158805328105587,
          3.323047914967523, 3.835574427742611, 5.094648956715952,
          4.868117032930769, 5.376372541761476, -1.727190081128829,
          -3.45683439863675, -5.605169446019096, -0.2591007294201401,
          -6.196235967511092, 0.8700445595441595, -6.693580082984033,
          -9.973298940902511, -1.228546058258308, -1.778473873127997,
          -2.863588535286143, -2.306783229826141, 5.794922003041152,
          7.782823613402732, 5.922297826135225, 4.366168456560783,
          4.393331161268728, 7.432212216606898, 5.376372541761476,
          8.921045899194265, -2.387907021436512, -5.390523498842208,
          -7.820615208588086, -0.1830238819466619, -9.026098123982962,
          -3.124740398165035, -8.267048789751669, -13.99740974932118,
          -1.311995151763001, 1.03924667861088, -4.282897589581626,
          -1.802131403274742, -2.378685014427916, -2.858753754386522,
          -3.128651652280993, -2.87529736009937, -1.092981827654348,
          -2.486249474018963, -1.727190081128829, -2.387907021436512,
          3.478537433168269, 3.162486654381813, 1.929362069412624,
          2.554816111258024, 2.450802325183254, 2.005097580329037,
          3.077233484452764, 5.289675192086496, -1.306252274523052,
          -1.924540455797193, 0.8089307276172231, 0.1745311932139862,
          -5.536048018007822, -3.612270433156081, -5.056094482013514,
          -4.150936493636122, -3.343689344221354, -2.437380936216606,
          -3.45683439863675, -5.390523498842208, 3.162486654381813,
          7.694270889376506, 6.325116353879784, 2.795579123304787,
          2.805784366472655, 2.094313819397243, 6.266075584923926,
          16.85227124419685, -0.3218247904602536, -2.325448469162336,
          1.216345925982085, 0.3342907581563925, -8.84093622113221,
          -7.501753263993352, -3.383579335130428, -4.089638778195734,
          -9.630640335753997, -6.718037799849757, -5.605169446019096,
          -7.820615208588086, 1.929362069412624, 6.325116353879784,
          16.65481772392231, -1.52751034833257, 9.56845082209581,
          -4.97377468321884, 7.715613883285052, 20.62390415224664,
          3.85201132461631, 3.008046123097777, 3.771832424523221,
          3.072975093350886, -1.100857511792741, 2.264122966861815,
          -3.449854529765026, -4.496272990998748, 2.429841121276548,
          -0.2648377020835075, -0.2591007294201401, -0.1830238819466619,
          2.554816111258024, 2.795579123304787, -1.52751034833257,
          8.309579819087396, -3.434090554679442, 4.008645820931794,
          2.207691804916842, 3.532193269393903, -6.355249485048906,
          -3.10295257534748, -1.764661695165005, -1.48764845435588,
          -5.522691796058115, -12.42024109318314, -4.509018906055925,
          -3.34471447300587, -6.215116692494317, -11.32726336676532,
          -6.196235967511092, -9.026098123982962, 2.450802325183254,
          2.805784366472655, 9.56845082209581, -3.434090554679442,
          17.31995670162867, -0.08715670847230617, 8.583409530767243,
          8.191938770150042, 2.639743911640359, 1.425972706741691,
          6.068607056478453, 2.618649377999702, 0.01790980108445428,
          -2.014924579198793, -6.611884017194679, -1.105338197190417,
          6.644839781301256, 1.262993685950383, 0.8700445595441595,
          -3.124740398165035, 2.005097580329037, 2.094313819397243,
          -4.97377468321884, 4.008645820931794, -0.08715670847230617,
          21.60969173425412, 3.78272085197206, -1.306125286215504,
          -3.983060122830408, -16.57598866067371, -2.703454821116731,
          -6.994707556124125, -6.030343075250605, -6.892202977180818,
          -9.161348600196341, -5.574589174004736, -3.223741811124502,
          -6.298130156830595, -6.693580082984033, -8.267048789751669,
          3.077233484452764, 6.266075584923926, 7.715613883285052,
          2.207691804916842, 8.583409530767243, 3.78272085197206,
          14.76681367477093, 17.55516221554868, -1.044520834401446,
          -3.276239931573811, 1.824319371886537, 0.4318785484492591,
          -14.68581000074875, -7.611220041618338, -11.07427442942752,
          -8.273216409082542, -10.90764512285358, -6.311463206932629,
          -9.973298940902511, -13.99740974932118, 5.289675192086496,
          16.85227124419685, 20.62390415224664, 3.532193269393903,
          8.191938770150042, -1.306125286215504, 17.55516221554868,
          48.97962793203256, 2.464278846470994, -1.995148786476195,
          3.552023642690095, 2.608844340193837, -1.224340104438079,
          -3.049877388674338, 1.994885866560606, 2.850808750123276,
          -4.114786139361408, -0.1224764752776564, -1.228546058258308,
          -1.311995151763001, -1.306252274523052, -0.3218247904602536,
          3.85201132461631, -6.355249485048906, 2.639743911640359,
          -3.983060122830408, -1.044520834401446, 2.464278846470994,
          7.986302123593876, 3.502372446512387, 2.127198741799675,
          2.142680789218174, -0.005114463255402324, -0.06047605860149874,
          5.019500172357599, -0.4187973830786254, -5.749284093800984,
          -3.137369670859187, -1.778473873127997, 1.03924667861088,
          -1.924540455797193, -2.325448469162336, 3.008046123097777,
          -3.10295257534748, 1.425972706741691, -16.57598866067371,
          -3.276239931573811, -1.995148786476195, 3.502372446512387,
          16.50062670924684, 3.479348267192873, 6.813822313911404,
          -2.234642701444805, -4.571120776660473, -1.507880835141334,
          -1.62047095951107, -3.564920648690664, -5.428028023658495,
          -2.863588535286143, -4.282897589581626, 0.8089307276172231,
          1.216345925982085, 3.771832424523221, -1.764661695165005,
          6.068607056478453, -2.703454821116731, 1.824319371886537,
          3.552023642690095, 2.127198741799675, 3.479348267192873,
          4.73094360865391, 3.277690027890099, -1.356471975117627,
          -1.65973684528307, 0.2231785103095052, -1.289711081566152,
          -3.832477684274747, -3.368706027054427, -2.306783229826141,
          -1.802131403274742, 0.1745311932139862, 0.3342907581563925,
          3.072975093350886, -1.48764845435588, 2.618649377999702,
          -6.994707556124125, 0.4318785484492591, 2.608844340193837,
          2.142680789218174, 6.813822313911404, 3.277690027890099,
          4.484446465585116;
      // clang-format on
      for (unsigned i = 0; i < nfunctions; ++i)
        for (unsigned j = 0; j < nfunctions; ++j)
          REQUIRE(laplace_matrix(i, j) ==
                  Approx(laplace(i, j)).epsilon(Tolerance));
    }

    SECTION("20-noded hexahedron coordinates of unit cell") {
      const unsigned nfunctions = 20;

      // Coordinates of a unit cell
      Eigen::Matrix<double, nfunctions, Dim> unit_cell;
      // clang-format off
      unit_cell << -1., -1., -1.,
                    1., -1., -1.,
                    1.,  1., -1.,
                   -1.,  1., -1.,
                   -1., -1.,  1.,
                    1., -1.,  1.,
                    1.,  1.,  1.,
                   -1.,  1.,  1.,
                    0., -1., -1.,
                   -1.,  0., -1.,
                   -1., -1.,  0.,
                    1.,  0., -1.,
                    1., -1.,  0.,
                    0.,  1., -1.,
                    1.,  1.,  0.,
                   -1.,  1.,  0.,
                    0., -1.,  1.,
                   -1.,  0.,  1.,
                    1.,  0.,  1.,
                    0.,  1.,  1.;
      // clang-format on

      auto coordinates = hex->unit_cell_coordinates();
      REQUIRE(coordinates.rows() == nfunctions);
      REQUIRE(coordinates.cols() == Dim);
      for (unsigned i = 0; i < nfunctions; ++i) {  // Iterate through nfunctions
        for (unsigned j = 0; j < Dim; ++j) {       // Dimension
          REQUIRE(coordinates(i, j) ==
                  Approx(unit_cell(i, j)).epsilon(Tolerance));
        }
      }
    }

    SECTION("20-noded hexahedron element for sides indices") {
      // Check for sides indices
      Eigen::MatrixXi indices = hex->sides_indices();
      REQUIRE(indices.rows() == 12);
      REQUIRE(indices.cols() == 2);
      REQUIRE(indices(0, 0) == 0);
      REQUIRE(indices(0, 1) == 1);

      REQUIRE(indices(1, 0) == 1);
      REQUIRE(indices(1, 1) == 2);

      REQUIRE(indices(2, 0) == 2);
      REQUIRE(indices(2, 1) == 3);

      REQUIRE(indices(3, 0) == 3);
      REQUIRE(indices(3, 1) == 0);

      REQUIRE(indices(4, 0) == 4);
      REQUIRE(indices(4, 1) == 5);

      REQUIRE(indices(5, 0) == 5);
      REQUIRE(indices(5, 1) == 6);

      REQUIRE(indices(6, 0) == 6);
      REQUIRE(indices(6, 1) == 7);

      REQUIRE(indices(7, 0) == 7);
      REQUIRE(indices(7, 1) == 4);

      REQUIRE(indices(8, 0) == 0);
      REQUIRE(indices(8, 1) == 4);

      REQUIRE(indices(9, 0) == 1);
      REQUIRE(indices(9, 1) == 5);

      REQUIRE(indices(10, 0) == 2);
      REQUIRE(indices(10, 1) == 6);

      REQUIRE(indices(11, 0) == 3);
      REQUIRE(indices(11, 1) == 7);
    }

    SECTION("20-noded hexahedron element for corner indices") {
      // Check for volume indices
      Eigen::VectorXi indices = hex->corner_indices();
      REQUIRE(indices.size() == 8);
      REQUIRE(indices(0) == 0);
      REQUIRE(indices(1) == 1);
      REQUIRE(indices(2) == 2);
      REQUIRE(indices(3) == 3);
      REQUIRE(indices(4) == 4);
      REQUIRE(indices(5) == 5);
      REQUIRE(indices(6) == 6);
      REQUIRE(indices(7) == 7);
    }

    SECTION("20-noded hexahedron element for inhedron indices") {
      // Check for inhedron indices
      Eigen::MatrixXi indices = hex->inhedron_indices();
      REQUIRE(indices.rows() == 12);
      REQUIRE(indices.cols() == 3);
      REQUIRE(indices(0, 0) == 0);
      REQUIRE(indices(0, 1) == 5);
      REQUIRE(indices(0, 2) == 4);

      REQUIRE(indices(1, 0) == 0);
      REQUIRE(indices(1, 1) == 1);
      REQUIRE(indices(1, 2) == 5);

      REQUIRE(indices(2, 0) == 3);
      REQUIRE(indices(2, 1) == 6);
      REQUIRE(indices(2, 2) == 7);

      REQUIRE(indices(3, 0) == 3);
      REQUIRE(indices(3, 1) == 2);
      REQUIRE(indices(3, 2) == 6);

      REQUIRE(indices(4, 0) == 2);
      REQUIRE(indices(4, 1) == 1);
      REQUIRE(indices(4, 2) == 6);

      REQUIRE(indices(5, 0) == 6);
      REQUIRE(indices(5, 1) == 1);
      REQUIRE(indices(5, 2) == 5);

      REQUIRE(indices(6, 0) == 7);
      REQUIRE(indices(6, 1) == 6);
      REQUIRE(indices(6, 2) == 5);

      REQUIRE(indices(7, 0) == 5);
      REQUIRE(indices(7, 1) == 4);
      REQUIRE(indices(7, 2) == 7);

      REQUIRE(indices(8, 0) == 7);
      REQUIRE(indices(8, 1) == 4);
      REQUIRE(indices(8, 2) == 0);

      REQUIRE(indices(9, 0) == 7);
      REQUIRE(indices(9, 1) == 0);
      REQUIRE(indices(9, 2) == 3);

      REQUIRE(indices(10, 0) == 3);
      REQUIRE(indices(10, 1) == 0);
      REQUIRE(indices(10, 2) == 1);

      REQUIRE(indices(11, 0) == 3);
      REQUIRE(indices(11, 1) == 1);
      REQUIRE(indices(11, 2) == 2);
    }

    SECTION("20-noded noded hexahedron shape function for face indices") {
      // Check for face indices
      Eigen::Matrix<int, 6, 8> indices;
      // clang-format off
      indices << 0, 1, 5, 4,  8, 12, 16, 10,
                 5, 1, 2, 0, 12, 11, 14, 18,
                 7, 6, 2, 3, 19, 14, 13, 15,
                 0, 4, 7, 3, 10, 17, 15,  9,
                 1, 0, 3, 2,  8,  9, 13, 11,
                 4, 5, 6, 7, 16, 18, 19, 17;
      // clang-format on

      // Check for all face indices
      for (unsigned i = 0; i < indices.rows(); ++i) {
        const auto check_indices = hex->face_indices(i);
        REQUIRE(check_indices.rows() == 8);
        REQUIRE(check_indices.cols() == 1);

        for (unsigned j = 0; j < indices.cols(); ++j)
          REQUIRE(check_indices(j) == indices(i, j));
      }
    }
  }
}

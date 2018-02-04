// Hexahedron shape function test
#include <memory>

#include "catch.hpp"

#include "hex_shapefn.h"

//! \brief Check hexahedron shapefn class
TEST_CASE("Hexahedron shape functions are checked",
          "[hexsf][hex][shapefn][3D]") {
  const unsigned Dim = 3;
  const double Tolerance = 1.E-7;

  //! Check for 8 noded shape function
  SECTION("Hexahedron shape function with eight nodes") {
    const unsigned nfunctions = 8;
    std::shared_ptr<mpm::ShapeFn<Dim>> hexsf =
        std::make_shared<mpm::HexahedronShapeFn<Dim, nfunctions>>();

    // Coordinates is (0, 0, 0)
    SECTION("Eight noded hexahedron shape function for coordinates(0, 0, 0)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();
      auto shapefn = hexsf->shapefn(coords);

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
      auto gradsf = hexsf->grad_shapefn(coords);
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
    SECTION(
        "Eight noded hexahedron shape function for coordinates(-1, -1, -1)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << -1., -1., -1.;
      auto shapefn = hexsf->shapefn(coords);
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
      auto gradsf = hexsf->grad_shapefn(coords);
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
    SECTION("Eight noded hexahedron shape function for coordinates(1, 1, 1)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << 1., 1., 1.;
      auto shapefn = hexsf->shapefn(coords);

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
      auto gradsf = hexsf->grad_shapefn(coords);
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

    SECTION("Eight noded hexahedron shape function for volume indices") {
      // Check for volume indices
      Eigen::VectorXi indices = hexsf->volume_indices();
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

    SECTION("Eight noded hexahedron shape function for inhedron indices") {
      // Check for inhedron indices
      Eigen::MatrixXi indices = hexsf->inhedron_indices();
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
  }

  // 20-Node (Serendipity) Hexahedron Element
  //! Check for 8 noded shape function
  SECTION("Hexahedron shape function with twenty nodes") {
    const unsigned nfunctions = 20;
    std::shared_ptr<mpm::ShapeFn<Dim>> hexsf =
        std::make_shared<mpm::HexahedronShapeFn<Dim, nfunctions>>();

    // Coordinates is (0, 0, 0)
    SECTION("Twenty noded hexahedron shape function for coordinates(0, 0, 0)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();
      auto shapefn = hexsf->shapefn(coords);

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
      auto gradsf = hexsf->grad_shapefn(coords);
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
        "Twenty noded hexahedron shape function for coordinates(-0.5, "
        "-0.5, -0.5)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << -0.5, -0.5, -0.5;
      auto shapefn = hexsf->shapefn(coords);

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
      auto gradsf = hexsf->grad_shapefn(coords);
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
        "Twenty noded hexahedron shape function for coordinates(0.5, "
        "0.5, 0.5)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << 0.5, 0.5, 0.5;
      auto shapefn = hexsf->shapefn(coords);

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
      auto gradsf = hexsf->grad_shapefn(coords);
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

    SECTION("20-noded hexahedron shape function for volume indices") {
      // Check for volume indices
      Eigen::VectorXi indices = hexsf->volume_indices();
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

    SECTION("20-noded hexahedron shape function for inhedron indices") {
      // Check for inhedron indices
      Eigen::MatrixXi indices = hexsf->inhedron_indices();
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
  }
}

// Hexahedron quadrature test
#include <cmath>

#include <limits>
#include <memory>

#include "catch.hpp"
#include "hexahedron_quadrature.h"

//! \brief Check HexahedronQuadratures class
TEST_CASE("Hexahedron quadratures are checked",
          "[hexquadrature][hex][quadrature][3D]") {
  const unsigned Dim = 3;
  const double Tolerance = 1.E-7;

  //! Check for a single point quadrature function
  SECTION("Hexahedron with a single quadrature") {
    const unsigned Nquadratures = 1;

    auto quad =
        std::make_shared<mpm::HexahedronQuadrature<Dim, Nquadratures>>();

    // Check quadratures
    auto points = quad->quadratures();

    // Check size
    REQUIRE(points.rows() == 3);
    REQUIRE(points.cols() == 1);

    // Check quadrature points
    REQUIRE(points(0, 0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(points(1, 0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(points(2, 0) == Approx(0.).epsilon(Tolerance));

    // Check weights
    auto weights = quad->weights();

    // Check size
    REQUIRE(weights.size() == 1);

    // Check weights
    REQUIRE(weights(0) == Approx(8.0).epsilon(Tolerance));
  }

  //! Check for eight quadrature points
  SECTION("Hexahedron with eight quadratures") {
    const unsigned Nquadratures = 8;

    auto quad =
        std::make_shared<mpm::HexahedronQuadrature<Dim, Nquadratures>>();

    // Check quadratures
    auto points = quad->quadratures();

    // Check size
    REQUIRE(points.rows() == 3);
    REQUIRE(points.cols() == 8);

    // Check quadrature points
    REQUIRE(points(0, 0) == Approx(-1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(1, 0) == Approx(-1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(2, 0) == Approx(-1. / std::sqrt(3.)).epsilon(Tolerance));

    REQUIRE(points(0, 1) == Approx(+1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(1, 1) == Approx(-1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(2, 1) == Approx(-1. / std::sqrt(3.)).epsilon(Tolerance));

    REQUIRE(points(0, 2) == Approx(+1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(1, 2) == Approx(+1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(2, 2) == Approx(-1. / std::sqrt(3.)).epsilon(Tolerance));

    REQUIRE(points(0, 3) == Approx(-1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(1, 3) == Approx(+1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(2, 3) == Approx(-1. / std::sqrt(3.)).epsilon(Tolerance));

    REQUIRE(points(0, 4) == Approx(-1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(1, 4) == Approx(-1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(2, 4) == Approx(+1. / std::sqrt(3.)).epsilon(Tolerance));

    REQUIRE(points(0, 5) == Approx(+1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(1, 5) == Approx(-1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(2, 5) == Approx(+1. / std::sqrt(3.)).epsilon(Tolerance));

    REQUIRE(points(0, 6) == Approx(+1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(1, 6) == Approx(+1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(2, 6) == Approx(+1. / std::sqrt(3.)).epsilon(Tolerance));

    REQUIRE(points(0, 7) == Approx(-1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(1, 7) == Approx(+1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(2, 7) == Approx(+1. / std::sqrt(3.)).epsilon(Tolerance));

    // Check weights
    auto weights = quad->weights();

    // Check size
    REQUIRE(weights.size() == 8);

    // Check weights
    REQUIRE(weights(0) == Approx(1.0).epsilon(Tolerance));
    REQUIRE(weights(1) == Approx(1.0).epsilon(Tolerance));
    REQUIRE(weights(2) == Approx(1.0).epsilon(Tolerance));
    REQUIRE(weights(3) == Approx(1.0).epsilon(Tolerance));
    REQUIRE(weights(4) == Approx(1.0).epsilon(Tolerance));
    REQUIRE(weights(5) == Approx(1.0).epsilon(Tolerance));
    REQUIRE(weights(6) == Approx(1.0).epsilon(Tolerance));
    REQUIRE(weights(7) == Approx(1.0).epsilon(Tolerance));
  }

  //! Check for twenty seven quadrature points
  SECTION("Hexahedron with twenty seven quadratures") {
    const unsigned Nquadratures = 27;

    auto quad =
        std::make_shared<mpm::HexahedronQuadrature<Dim, Nquadratures>>();

    // Check quadratures
    auto points = quad->quadratures();

    // Check size
    REQUIRE(points.rows() == 3);
    REQUIRE(points.cols() == 27);

    // Check quadrature points
    REQUIRE(points(0, 0) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(1, 0) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(2, 0) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(0, 1) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(1, 1) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(2, 1) == Approx(+0.).epsilon(Tolerance));

    REQUIRE(points(0, 2) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(1, 2) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(2, 2) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(0, 3) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(1, 3) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(2, 3) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(0, 4) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(1, 4) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(2, 4) == Approx(+0.).epsilon(Tolerance));

    REQUIRE(points(0, 5) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(1, 5) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(2, 5) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(0, 6) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(1, 6) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(2, 6) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(0, 7) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(1, 7) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(2, 7) == Approx(0.).epsilon(Tolerance));

    REQUIRE(points(0, 8) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(1, 8) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(2, 8) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(0, 9) == Approx(0.).epsilon(Tolerance));
    REQUIRE(points(1, 9) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(2, 9) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(0, 10) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(1, 10) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(2, 10) == Approx(+0.).epsilon(Tolerance));

    REQUIRE(points(0, 11) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(1, 11) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(2, 11) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(0, 12) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(1, 12) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(2, 12) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(0, 13) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(1, 13) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(2, 13) == Approx(+0.).epsilon(Tolerance));

    REQUIRE(points(0, 14) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(1, 14) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(2, 14) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(0, 15) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(1, 15) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(2, 15) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(0, 16) == Approx(0.).epsilon(Tolerance));
    REQUIRE(points(1, 16) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(2, 16) == Approx(+0.).epsilon(Tolerance));

    REQUIRE(points(0, 17) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(1, 17) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(2, 17) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(0, 18) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(1, 18) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(2, 18) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(0, 19) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(1, 19) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(2, 19) == Approx(0.).epsilon(Tolerance));

    REQUIRE(points(0, 20) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(1, 20) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(2, 20) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(0, 21) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(1, 21) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(2, 21) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(0, 22) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(1, 22) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(2, 22) == Approx(+0.).epsilon(Tolerance));

    REQUIRE(points(0, 23) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(1, 23) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(2, 23) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(0, 24) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(1, 24) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(2, 24) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(0, 25) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(1, 25) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(2, 25) == Approx(+0.).epsilon(Tolerance));

    REQUIRE(points(0, 26) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(1, 26) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(2, 26) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));

    // Check weights
    auto weights = quad->weights();

    // Check size
    REQUIRE(weights.size() == 27);

    // Check weights
    REQUIRE(weights(0) ==
            Approx((5. / 9.) * (5. / 9.) * (5. / 9.)).epsilon(Tolerance));
    REQUIRE(weights(1) ==
            Approx((5. / 9.) * (5. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights(2) ==
            Approx((5. / 9.) * (5. / 9.) * (5. / 9.)).epsilon(Tolerance));
    REQUIRE(weights(3) ==
            Approx((5. / 9.) * (5. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights(4) ==
            Approx((5. / 9.) * (8. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights(5) ==
            Approx((5. / 9.) * (5. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights(6) ==
            Approx((5. / 9.) * (5. / 9.) * (5. / 9.)).epsilon(Tolerance));
    REQUIRE(weights(7) ==
            Approx((5. / 9.) * (5. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights(8) ==
            Approx((5. / 9.) * (5. / 9.) * (5. / 9.)).epsilon(Tolerance));

    REQUIRE(weights(9) ==
            Approx((5. / 9.) * (5. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights(10) ==
            Approx((5. / 9.) * (8. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights(11) ==
            Approx((5. / 9.) * (5. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights(12) ==
            Approx((5. / 9.) * (8. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights(13) ==
            Approx((8. / 9.) * (8. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights(14) ==
            Approx((5. / 9.) * (8. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights(15) ==
            Approx((5. / 9.) * (5. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights(16) ==
            Approx((5. / 9.) * (8. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights(17) ==
            Approx((5. / 9.) * (5. / 9.) * (8. / 9.)).epsilon(Tolerance));

    REQUIRE(weights(18) ==
            Approx((5. / 9.) * (5. / 9.) * (5. / 9.)).epsilon(Tolerance));
    REQUIRE(weights(19) ==
            Approx((5. / 9.) * (5. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights(20) ==
            Approx((5. / 9.) * (5. / 9.) * (5. / 9.)).epsilon(Tolerance));
    REQUIRE(weights(21) ==
            Approx((5. / 9.) * (5. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights(22) ==
            Approx((5. / 9.) * (8. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights(23) ==
            Approx((5. / 9.) * (5. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights(24) ==
            Approx((5. / 9.) * (5. / 9.) * (5. / 9.)).epsilon(Tolerance));
    REQUIRE(weights(25) ==
            Approx((5. / 9.) * (5. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights(26) ==
            Approx((5. / 9.) * (5. / 9.) * (5. / 9.)).epsilon(Tolerance));
  }

  //! Check for sixty four quadrature points
  SECTION("Hexahedron with sixty four quadratures") {
    const unsigned Nquadratures = 64;

    auto quad =
        std::make_shared<mpm::HexahedronQuadrature<Dim, Nquadratures>>();

    // Check quadratures
    auto points = quad->quadratures();

    // Check size
    REQUIRE(points.rows() == 3);
    REQUIRE(points.cols() == 64);

    // Check quadrature points
    REQUIRE(points(0, 0) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 0) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 0) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 1) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 1) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 1) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 2) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 2) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 2) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 3) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 3) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 3) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 4) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 4) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 4) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 5) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 5) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 5) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 6) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 6) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 6) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 7) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 7) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 7) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 8) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 8) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 8) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 9) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 9) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 9) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 10) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 10) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 10) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 11) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 11) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 11) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 12) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 12) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 12) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 13) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 13) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 13) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 14) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 14) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 14) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 15) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 15) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 15) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 16) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 16) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 16) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 17) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 17) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 17) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 18) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 18) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 18) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 19) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 19) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 19) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 20) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 20) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 20) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 21) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 21) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 21) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 22) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 22) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 22) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 23) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 23) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 23) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 24) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 24) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 24) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 25) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 25) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 25) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 26) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 26) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 26) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 27) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 27) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 27) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 28) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 28) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 28) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 29) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 29) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 29) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 30) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 30) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 30) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 31) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 31) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 31) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 32) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 32) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 32) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 33) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 33) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 33) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 34) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 34) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 34) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 35) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 35) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 35) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 36) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 36) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 36) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 37) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 37) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 37) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 38) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 38) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 38) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 39) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 39) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 39) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 40) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 40) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 40) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 41) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 41) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 41) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 42) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 42) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 42) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 43) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 43) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 43) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 44) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 44) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 44) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 45) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 45) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 45) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 46) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 46) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 46) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 47) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 47) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 47) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 48) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 48) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 48) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 49) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 49) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 49) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 50) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 50) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 50) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 51) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 51) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 51) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 52) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 52) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 52) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 53) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 53) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 53) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 54) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 54) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 54) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 55) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 55) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 55) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 56) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 56) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 56) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 57) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 57) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 57) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 58) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 58) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 58) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 59) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 59) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 59) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 60) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 60) ==
            Approx(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 60) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 61) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 61) ==
            Approx(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 61) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 62) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 62) ==
            Approx(+std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 62) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    REQUIRE(points(0, 63) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(1, 63) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));
    REQUIRE(points(2, 63) ==
            Approx(+std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))
                .epsilon(Tolerance));

    // Check weights
    auto weights = quad->weights();

    // Check size
    REQUIRE(weights.size() == 64);

    // Check weights
    REQUIRE(weights(0) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. - std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(1) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. - std::sqrt(30.)) / 36. *
                   (18. + std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(2) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. - std::sqrt(30.)) / 36. *
                   (18. + std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(3) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. - std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(4) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(5) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. + std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(6) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. + std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(7) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(8) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(9) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. + std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(10) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. + std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(11) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(12) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. - std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(13) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(14) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(15) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. - std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(16) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. - std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(17) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. - std::sqrt(30.)) / 36. *
                   (18. + std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(18) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. - std::sqrt(30.)) / 36. *
                   (18. + std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(19) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. - std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(20) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(21) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. + std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(22) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. + std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(23) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(24) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(25) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. + std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(26) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. + std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(27) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(28) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. - std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(29) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(30) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(31) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. - std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(32) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. - std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(33) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. - std::sqrt(30.)) / 36. *
                   (18. + std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(34) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. - std::sqrt(30.)) / 36. *
                   (18. + std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(35) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. - std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(36) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(37) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. + std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(38) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. + std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(39) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(40) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(41) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. + std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(42) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. + std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(43) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(44) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. - std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(45) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(46) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(47) ==
            Approx((18. + std::sqrt(30.)) / 36. * (18. - std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(48) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. - std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(49) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. - std::sqrt(30.)) / 36. *
                   (18. + std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(50) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. - std::sqrt(30.)) / 36. *
                   (18. + std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(51) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. - std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(52) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(53) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. + std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(54) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. + std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(55) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(56) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(57) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. + std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(58) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. + std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(59) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(60) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. - std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(61) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(62) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. + std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
    REQUIRE(weights(63) ==
            Approx((18. - std::sqrt(30.)) / 36. * (18. - std::sqrt(30.)) / 36. *
                   (18. - std::sqrt(30.)) / 36.)
                .epsilon(Tolerance));
  }
}

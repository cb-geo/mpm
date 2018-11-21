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
  SECTION("Hexahedron with four quadratures") {
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
  SECTION("Hexahedron with nine quadratures") {
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
}

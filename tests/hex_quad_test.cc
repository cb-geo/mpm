// Hexahedron quadrature test
#include <cmath>

#include <limits>
#include <memory>

#include "catch.hpp"
#include "hex_quad.h"

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
    REQUIRE(points.rows() == 1);
    REQUIRE(points.cols() == 3);

    // Check quadrature points
    REQUIRE(points(0, 0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(points(0, 1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(points(0, 2) == Approx(0.).epsilon(Tolerance));

    // Check weights
    auto weights = quad->weights();

    // Check size
    REQUIRE(weights.size() == 1);

    // Check weights
    REQUIRE(weights.at(0) == Approx(8.0).epsilon(Tolerance));
  }

  //! Check for eight quadrature points
  SECTION("Hexahedron with four quadratures") {
    const unsigned Nquadratures = 8;

    auto quad =
        std::make_shared<mpm::HexahedronQuadrature<Dim, Nquadratures>>();

    // Check quadratures
    auto points = quad->quadratures();

    // Check size
    REQUIRE(points.rows() == 8);
    REQUIRE(points.cols() == 3);

    // Check quadrature points
    REQUIRE(points(0, 0) == Approx(-1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(0, 1) == Approx(-1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(0, 2) == Approx(-1. / std::sqrt(3.)).epsilon(Tolerance));

    REQUIRE(points(1, 0) == Approx(+1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(1, 1) == Approx(-1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(1, 2) == Approx(-1. / std::sqrt(3.)).epsilon(Tolerance));

    REQUIRE(points(2, 0) == Approx(+1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(2, 1) == Approx(+1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(2, 2) == Approx(-1. / std::sqrt(3.)).epsilon(Tolerance));

    REQUIRE(points(3, 0) == Approx(-1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(3, 1) == Approx(+1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(3, 2) == Approx(-1. / std::sqrt(3.)).epsilon(Tolerance));

    REQUIRE(points(4, 0) == Approx(-1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(4, 1) == Approx(-1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(4, 2) == Approx(+1. / std::sqrt(3.)).epsilon(Tolerance));

    REQUIRE(points(5, 0) == Approx(+1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(5, 1) == Approx(-1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(5, 2) == Approx(+1. / std::sqrt(3.)).epsilon(Tolerance));

    REQUIRE(points(6, 0) == Approx(+1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(6, 1) == Approx(+1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(6, 2) == Approx(+1. / std::sqrt(3.)).epsilon(Tolerance));

    REQUIRE(points(7, 0) == Approx(-1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(7, 1) == Approx(+1. / std::sqrt(3.)).epsilon(Tolerance));
    REQUIRE(points(7, 2) == Approx(+1. / std::sqrt(3.)).epsilon(Tolerance));

    // Check weights
    auto weights = quad->weights();

    // Check size
    REQUIRE(weights.size() == 8);

    // Check weights
    REQUIRE(weights.at(0) == Approx(1.0).epsilon(Tolerance));
    REQUIRE(weights.at(1) == Approx(1.0).epsilon(Tolerance));
    REQUIRE(weights.at(2) == Approx(1.0).epsilon(Tolerance));
    REQUIRE(weights.at(3) == Approx(1.0).epsilon(Tolerance));
    REQUIRE(weights.at(4) == Approx(1.0).epsilon(Tolerance));
    REQUIRE(weights.at(5) == Approx(1.0).epsilon(Tolerance));
    REQUIRE(weights.at(6) == Approx(1.0).epsilon(Tolerance));
    REQUIRE(weights.at(7) == Approx(1.0).epsilon(Tolerance));
  }

  //! Check for twenty seven quadrature points
  SECTION("Hexahedron with nine quadratures") {
    const unsigned Nquadratures = 27;

    auto quad =
        std::make_shared<mpm::HexahedronQuadrature<Dim, Nquadratures>>();

    // Check quadratures
    auto points = quad->quadratures();

    // Check size
    REQUIRE(points.rows() == 27);
    REQUIRE(points.cols() == 3);

    // Check quadrature points
    REQUIRE(points(0, 0) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(0, 1) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(0, 1) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(1, 0) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(1, 1) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(1, 2) == Approx(+0.).epsilon(Tolerance));

    REQUIRE(points(2, 0) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(2, 1) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(2, 2) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(3, 0) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(3, 1) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(3, 2) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(4, 0) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(4, 1) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(4, 2) == Approx(+0.).epsilon(Tolerance));

    REQUIRE(points(5, 0) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(5, 1) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(5, 2) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(6, 0) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(6, 1) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(6, 2) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(7, 0) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(7, 1) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(7, 2) == Approx(0.).epsilon(Tolerance));

    REQUIRE(points(8, 0) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(8, 1) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(8, 2) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(9, 0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(points(9, 1) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(9, 2) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(10, 0) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(10, 1) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(10, 2) == Approx(+0.).epsilon(Tolerance));

    REQUIRE(points(11, 0) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(11, 1) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(11, 2) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(12, 0) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(12, 1) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(12, 2) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(13, 0) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(13, 1) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(13, 2) == Approx(+0.).epsilon(Tolerance));

    REQUIRE(points(14, 0) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(14, 1) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(14, 2) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(15, 0) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(15, 1) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(15, 2) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(16, 0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(points(16, 1) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(16, 2) == Approx(+0.).epsilon(Tolerance));

    REQUIRE(points(17, 0) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(17, 1) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(17, 2) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(18, 0) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(18, 1) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(18, 2) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(19, 0) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(19, 1) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(19, 2) == Approx(0.).epsilon(Tolerance));

    REQUIRE(points(20, 0) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(20, 1) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(20, 2) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(21, 0) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(21, 1) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(21, 2) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(22, 0) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(22, 1) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(22, 2) == Approx(+0.).epsilon(Tolerance));

    REQUIRE(points(23, 0) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(23, 1) == Approx(+0.).epsilon(Tolerance));
    REQUIRE(points(23, 2) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(24, 0) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(24, 1) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(24, 2) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(25, 0) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(25, 1) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(25, 2) == Approx(+0.).epsilon(Tolerance));

    REQUIRE(points(26, 0) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(26, 1) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(26, 2) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));

    // Check weights
    auto weights = quad->weights();

    // Check size
    REQUIRE(weights.size() == 27);

    // Check weights
    REQUIRE(weights.at(0) ==
            Approx((5. / 9.) * (5. / 9.) * (5. / 9.)).epsilon(Tolerance));
    REQUIRE(weights.at(1) ==
            Approx((5. / 9.) * (5. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights.at(2) ==
            Approx((5. / 9.) * (5. / 9.) * (5. / 9.)).epsilon(Tolerance));
    REQUIRE(weights.at(3) ==
            Approx((5. / 9.) * (5. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights.at(4) ==
            Approx((5. / 9.) * (8. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights.at(5) ==
            Approx((5. / 9.) * (5. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights.at(6) ==
            Approx((5. / 9.) * (5. / 9.) * (5. / 9.)).epsilon(Tolerance));
    REQUIRE(weights.at(7) ==
            Approx((5. / 9.) * (5. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights.at(8) ==
            Approx((5. / 9.) * (5. / 9.) * (5. / 9.)).epsilon(Tolerance));

    REQUIRE(weights.at(9) ==
            Approx((5. / 9.) * (5. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights.at(10) ==
            Approx((5. / 9.) * (8. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights.at(11) ==
            Approx((5. / 9.) * (5. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights.at(12) ==
            Approx((5. / 9.) * (8. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights.at(13) ==
            Approx((8. / 9.) * (8. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights.at(14) ==
            Approx((5. / 9.) * (8. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights.at(15) ==
            Approx((5. / 9.) * (5. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights.at(16) ==
            Approx((5. / 9.) * (8. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights.at(17) ==
            Approx((5. / 9.) * (5. / 9.) * (8. / 9.)).epsilon(Tolerance));

    REQUIRE(weights.at(18) ==
            Approx((5. / 9.) * (5. / 9.) * (5. / 9.)).epsilon(Tolerance));
    REQUIRE(weights.at(19) ==
            Approx((5. / 9.) * (5. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights.at(20) ==
            Approx((5. / 9.) * (5. / 9.) * (5. / 9.)).epsilon(Tolerance));
    REQUIRE(weights.at(21) ==
            Approx((5. / 9.) * (5. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights.at(22) ==
            Approx((5. / 9.) * (8. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights.at(23) ==
            Approx((5. / 9.) * (5. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights.at(24) ==
            Approx((5. / 9.) * (5. / 9.) * (5. / 9.)).epsilon(Tolerance));
    REQUIRE(weights.at(25) ==
            Approx((5. / 9.) * (5. / 9.) * (8. / 9.)).epsilon(Tolerance));
    REQUIRE(weights.at(26) ==
            Approx((5. / 9.) * (5. / 9.) * (5. / 9.)).epsilon(Tolerance));
  }
}
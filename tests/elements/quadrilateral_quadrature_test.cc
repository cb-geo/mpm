// Quadrilateral quadrature test
#include <cmath>

#include <limits>
#include <memory>

#include "catch.hpp"
#include "quadrilateral_quadrature.h"

//! \brief Check QuadrilateralQuadratures class
TEST_CASE("Quadrilateral quadratures are checked",
          "[quadquadrature][quad][quadrature][2D]") {
  const unsigned Dim = 2;
  const double Tolerance = 1.E-7;

  //! Check for a single point quadrature function
  SECTION("Quadrilateral with a single quadrature") {
    const unsigned Nquadratures = 1;

    auto quad =
        std::make_shared<mpm::QuadrilateralQuadrature<Dim, Nquadratures>>();

    // Che quadratures
    auto points = quad->quadratures();

    // Check size
    REQUIRE(points.rows() == 2);
    REQUIRE(points.cols() == 1);

    // Check quadrature points
    REQUIRE(points(0, 0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(points(1, 0) == Approx(0.).epsilon(Tolerance));

    // Check weights
    auto weights = quad->weights();

    // Check size
    REQUIRE(weights.size() == 1);

    // Check weights
    REQUIRE(weights(0) == Approx(4.0).epsilon(Tolerance));
  }

  //! Check for four quadrature points
  SECTION("Quadrilateral with four quadratures") {
    const unsigned Nquadratures = 4;

    auto quad =
        std::make_shared<mpm::QuadrilateralQuadrature<Dim, Nquadratures>>();

    // Check quadratures
    auto points = quad->quadratures();

    // Check size
    REQUIRE(points.rows() == 2);
    REQUIRE(points.cols() == 4);

    // Check quadrature points
    REQUIRE(points(0, 0) == Approx(-std::sqrt(3.) / 3.).epsilon(Tolerance));
    REQUIRE(points(1, 0) == Approx(-std::sqrt(3.) / 3.).epsilon(Tolerance));
    REQUIRE(points(0, 1) == Approx(+std::sqrt(3.) / 3.).epsilon(Tolerance));
    REQUIRE(points(1, 1) == Approx(-std::sqrt(3.) / 3.).epsilon(Tolerance));
    REQUIRE(points(0, 2) == Approx(+std::sqrt(3.) / 3.).epsilon(Tolerance));
    REQUIRE(points(1, 2) == Approx(+std::sqrt(3.) / 3.).epsilon(Tolerance));
    REQUIRE(points(0, 3) == Approx(-std::sqrt(3.) / 3.).epsilon(Tolerance));
    REQUIRE(points(1, 3) == Approx(+std::sqrt(3.) / 3.).epsilon(Tolerance));

    // Check weights
    auto weights = quad->weights();

    // Check size
    REQUIRE(weights.size() == 4);

    // Check weights
    REQUIRE(weights(0) == Approx(1.0).epsilon(Tolerance));
    REQUIRE(weights(1) == Approx(1.0).epsilon(Tolerance));
    REQUIRE(weights(2) == Approx(1.0).epsilon(Tolerance));
    REQUIRE(weights(3) == Approx(1.0).epsilon(Tolerance));
  }

  //! Check for four quadrature points
  SECTION("Quadrilateral with nine quadratures") {
    const unsigned Nquadratures = 9;

    auto quad =
        std::make_shared<mpm::QuadrilateralQuadrature<Dim, Nquadratures>>();

    // Check quadratures
    auto points = quad->quadratures();

    // Check size
    REQUIRE(points.rows() == 2);
    REQUIRE(points.cols() == 9);

    // Check quadrature points
    REQUIRE(points(0, 0) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(1, 0) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(0, 1) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(1, 1) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(0, 2) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(1, 2) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(0, 3) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(1, 3) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(0, 4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(points(1, 4) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(0, 5) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(1, 5) == Approx(0.).epsilon(Tolerance));

    REQUIRE(points(0, 6) == Approx(0.).epsilon(Tolerance));
    REQUIRE(points(1, 6) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(0, 7) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(1, 7) == Approx(0.).epsilon(Tolerance));

    REQUIRE(points(0, 8) == Approx(0.).epsilon(Tolerance));
    REQUIRE(points(1, 8) == Approx(0.).epsilon(Tolerance));

    // Check weights
    auto weights = quad->weights();

    // Check size
    REQUIRE(weights.size() == 9);

    // Check weights
    REQUIRE(weights(0) == Approx(25. / 81.).epsilon(Tolerance));
    REQUIRE(weights(1) == Approx(25. / 81.).epsilon(Tolerance));
    REQUIRE(weights(2) == Approx(25. / 81.).epsilon(Tolerance));
    REQUIRE(weights(3) == Approx(25. / 81.).epsilon(Tolerance));
    REQUIRE(weights(4) == Approx(40. / 81.).epsilon(Tolerance));
    REQUIRE(weights(5) == Approx(40. / 81.).epsilon(Tolerance));
    REQUIRE(weights(6) == Approx(40. / 81.).epsilon(Tolerance));
    REQUIRE(weights(7) == Approx(40. / 81.).epsilon(Tolerance));
    REQUIRE(weights(8) == Approx(64. / 81.).epsilon(Tolerance));
  }
}

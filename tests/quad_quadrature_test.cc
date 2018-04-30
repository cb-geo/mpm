// Quadrilateral quadrature test
#include <cmath>

#include <limits>
#include <memory>

#include "catch.hpp"
#include "quad_quadrature.h"

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
    REQUIRE(points.rows() == 1);
    REQUIRE(points.cols() == 2);

    // Check quadrature points
    REQUIRE(points(0, 0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(points(0, 1) == Approx(0.).epsilon(Tolerance));

    // Check weights
    auto weights = quad->weights();

    // Check size
    REQUIRE(weights.size() == 1);

    // Check weights
    REQUIRE(weights.at(0) == Approx(8.0).epsilon(Tolerance));
  }

  //! Check for four quadrature points
  SECTION("Quadrilateral with four quadratures") {
    const unsigned Nquadratures = 4;

    auto quad =
        std::make_shared<mpm::QuadrilateralQuadrature<Dim, Nquadratures>>();

    // Check quadratures
    auto points = quad->quadratures();

    // Check size
    REQUIRE(points.rows() == 4);
    REQUIRE(points.cols() == 2);

    // Check quadrature points
    REQUIRE(points(0, 0) == Approx(-std::sqrt(3.) / 3.).epsilon(Tolerance));
    REQUIRE(points(0, 1) == Approx(-std::sqrt(3.) / 3.).epsilon(Tolerance));
    REQUIRE(points(1, 0) == Approx(+std::sqrt(3.) / 3.).epsilon(Tolerance));
    REQUIRE(points(1, 1) == Approx(-std::sqrt(3.) / 3.).epsilon(Tolerance));
    REQUIRE(points(2, 0) == Approx(+std::sqrt(3.) / 3.).epsilon(Tolerance));
    REQUIRE(points(2, 1) == Approx(+std::sqrt(3.) / 3.).epsilon(Tolerance));
    REQUIRE(points(3, 0) == Approx(-std::sqrt(3.) / 3.).epsilon(Tolerance));
    REQUIRE(points(3, 1) == Approx(+std::sqrt(3.) / 3.).epsilon(Tolerance));

    // Check weights
    auto weights = quad->weights();

    // Check size
    REQUIRE(weights.size() == 4);

    // Check weights
    REQUIRE(weights.at(0) == Approx(1.0).epsilon(Tolerance));
    REQUIRE(weights.at(1) == Approx(1.0).epsilon(Tolerance));
    REQUIRE(weights.at(2) == Approx(1.0).epsilon(Tolerance));
    REQUIRE(weights.at(3) == Approx(1.0).epsilon(Tolerance));
  }

  //! Check for four quadrature points
  SECTION("Quadrilateral with nine quadratures") {
    const unsigned Nquadratures = 9;

    auto quad =
        std::make_shared<mpm::QuadrilateralQuadrature<Dim, Nquadratures>>();

    // Check quadratures
    auto points = quad->quadratures();

    // Check size
    REQUIRE(points.rows() == 9);
    REQUIRE(points.cols() == 2);

    // Check quadrature points
    REQUIRE(points(0, 0) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(0, 1) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(1, 0) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(1, 1) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(2, 0) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(2, 1) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(3, 0) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(3, 1) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(4, 0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(points(4, 1) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(5, 0) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(5, 1) == Approx(0.).epsilon(Tolerance));

    REQUIRE(points(6, 0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(points(6, 1) == Approx(+std::sqrt(3. / 5.)).epsilon(Tolerance));

    REQUIRE(points(7, 0) == Approx(-std::sqrt(3. / 5.)).epsilon(Tolerance));
    REQUIRE(points(7, 1) == Approx(0.).epsilon(Tolerance));

    REQUIRE(points(8, 0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(points(8, 1) == Approx(0.).epsilon(Tolerance));

    // Check weights
    auto weights = quad->weights();

    // Check size
    REQUIRE(weights.size() == 9);

    // Check weights
    REQUIRE(weights.at(0) == Approx(25. / 81.).epsilon(Tolerance));
    REQUIRE(weights.at(1) == Approx(25. / 81.).epsilon(Tolerance));
    REQUIRE(weights.at(2) == Approx(25. / 81.).epsilon(Tolerance));
    REQUIRE(weights.at(3) == Approx(25. / 81.).epsilon(Tolerance));
    REQUIRE(weights.at(4) == Approx(40. / 81.).epsilon(Tolerance));
    REQUIRE(weights.at(5) == Approx(40. / 81.).epsilon(Tolerance));
    REQUIRE(weights.at(6) == Approx(40. / 81.).epsilon(Tolerance));
    REQUIRE(weights.at(7) == Approx(40. / 81.).epsilon(Tolerance));
    REQUIRE(weights.at(8) == Approx(64. / 81.).epsilon(Tolerance));
  }
}
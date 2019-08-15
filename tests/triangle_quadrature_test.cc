// Quadrilateral quadrature test
#include <cmath>

#include <limits>
#include <memory>

#include "catch.hpp"
#include "triangle_quadrature.h"

//! \brief Check TriangleQuadrature class
TEST_CASE("Triangle quadratures are checked",
          "[triquadrature][tri][quadrature][2D]") {
  const unsigned Dim = 2;
  const double Tolerance = 1.E-7;

  //! Check for a single point quadrature function
  SECTION("Triangle with a single quadrature") {
    const unsigned Nquadratures = 1;

    auto tri =
        std::make_shared<mpm::TriangleQuadrature<Dim, Nquadratures>>();

    // Check quadratures
    auto points = tri->quadratures();

    // Check size
    REQUIRE(points.rows() == 2);
    REQUIRE(points.cols() == 1);

    // Check quadrature points
    REQUIRE(points(0, 0) == Approx(1. / 3.).epsilon(Tolerance));
    REQUIRE(points(1, 0) == Approx(1. / 3.).epsilon(Tolerance));

    // Check weights
    auto weights = tri->weights();

    // Check size
    REQUIRE(weights.size() == 1);

    // Check weights
    REQUIRE(weights(0) == Approx(0.5).epsilon(Tolerance));
  }

  //! Check for three quadrature points
  SECTION("Triangle with three quadratures") {
    const unsigned Nquadratures = 3;

    auto tri =
        std::make_shared<mpm::TriangleQuadrature<Dim, Nquadratures>>();

    // Check quadratures
    auto points = tri->quadratures();

    // Check size
    REQUIRE(points.rows() == 2);
    REQUIRE(points.cols() == 3);

    // Check quadrature points
    REQUIRE(points(0, 0) == Approx(1. / 6.).epsilon(Tolerance));
    REQUIRE(points(1, 0) == Approx(1. / 6.).epsilon(Tolerance));
    REQUIRE(points(0, 1) == Approx(2. / 3.).epsilon(Tolerance));
    REQUIRE(points(1, 1) == Approx(1. / 6.).epsilon(Tolerance));
    REQUIRE(points(0, 2) == Approx(1. / 6.).epsilon(Tolerance));
    REQUIRE(points(1, 2) == Approx(2. / 3.).epsilon(Tolerance));

    // Check weights
    auto weights = tri->weights();

    // Check size
    REQUIRE(weights.size() == 3);

    // Check weights
    REQUIRE(weights(0) == Approx(1. / 3.).epsilon(Tolerance));
    REQUIRE(weights(1) == Approx(1. / 3.).epsilon(Tolerance));
    REQUIRE(weights(2) == Approx(1. / 3.).epsilon(Tolerance));
  }
}

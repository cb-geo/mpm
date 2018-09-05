// Quadrilateral element test
#include <memory>

#include "catch.hpp"

#include "quadrilateral_gimp_element.h"

//! \brief Check quadrilateral element class
TEST_CASE("Quadrilateral gimp elements are checked", "[gimp]") {
  const unsigned Dim = 2;
  const double Tolerance = 1.E-7;

  //! Check for center element nodes
  SECTION("16 Node Quadrilateral GIMP Element") {
    const unsigned nfunctions = 16;
    std::shared_ptr<mpm::Element<Dim>> quad =
        std::make_shared<mpm::QuadrilateralGIMPElement<Dim, nfunctions>>();

    // Check degree
    REQUIRE(quad->degree() == mpm::ElementDegree::Linear);

    // Coordinates is (0,0) Size is (0,0)
    SECTION(
        "16 Node quadrilateral element matrix for coordinate(0,0), Size "
        "(0,0)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();
      Eigen::Matrix<double, Dim, 1> psize;
      psize.setZero();
      Eigen::Matrix<double, Dim, 1> defgrad;
      defgrad.setZero();

      auto shapefn = quad->shapefn(coords, psize, defgrad);

      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      REQUIRE(shapefn(0) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(3) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(4) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(5) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(6) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(7) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(8) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(9) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(10) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(11) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(12) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(13) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(14) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(15) == Approx(0.0).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(coords, psize, defgrad);
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      REQUIRE(gradsf(0, 0) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(gradsf(3, 0) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(4, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(5, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(6, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(7, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(8, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(9, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(10, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(11, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(12, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(13, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(14, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(15, 0) == Approx(0.0).epsilon(Tolerance));

      REQUIRE(gradsf(0, 1) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(gradsf(3, 1) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(gradsf(4, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(5, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(6, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(7, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(8, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(9, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(10, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(11, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(12, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(13, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(14, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(15, 1) == Approx(0.0).epsilon(Tolerance));
    }

    // Coordinates is (-1,-1) Size is (0,0)
    SECTION(
        "16 Node quadrilateral element matrix for coordinate(-1,-1), Size "
        "(0,0)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << -1, -1;
      Eigen::Matrix<double, Dim, 1> psize;
      psize.setZero();
      Eigen::Matrix<double, Dim, 1> defgrad;
      defgrad.setZero();

      auto shapefn = quad->shapefn(coords, psize, defgrad);

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
      REQUIRE(shapefn(8) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(9) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(10) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(11) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(12) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(13) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(14) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(15) == Approx(0.0).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(coords, psize, defgrad);
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      REQUIRE(gradsf(0, 0) == Approx(0.5).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(3, 0) == Approx(0.0).epsilon(Tolerance));

      REQUIRE(gradsf(4, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(5, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(6, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(7, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(8, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(9, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(10, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(11, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(12, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(13, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(14, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(15, 0) == Approx(-0.5).epsilon(Tolerance));

      REQUIRE(gradsf(0, 1) == Approx(0.5).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(3, 1) == Approx(0.0).epsilon(Tolerance));

      REQUIRE(gradsf(4, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(5, 1) == Approx(-0.5).epsilon(Tolerance));
      REQUIRE(gradsf(6, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(7, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(8, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(9, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(10, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(11, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(12, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(13, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(14, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(15, 1) == Approx(0.0).epsilon(Tolerance));
    }

    // Coordinates is (1,1) Size is (0,0)
    SECTION(
        "16 Node quadrilateral element matrix for coordinate(1,1), Size "
        "(1,1)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << 1, 1;
      Eigen::Matrix<double, Dim, 1> psize;
      psize.setZero();
      Eigen::Matrix<double, Dim, 1> defgrad;
      defgrad.setZero();

      auto shapefn = quad->shapefn(coords, psize, defgrad);

      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      REQUIRE(shapefn(0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(1.0).epsilon(Tolerance));
      REQUIRE(shapefn(3) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(4) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(5) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(6) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(7) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(8) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(9) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(10) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(11) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(12) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(13) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(14) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(15) == Approx(0.0).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(coords, psize, defgrad);
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      REQUIRE(gradsf(0, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(0.5).epsilon(Tolerance));
      REQUIRE(gradsf(3, 0) == Approx(-0.5).epsilon(Tolerance));

      REQUIRE(gradsf(4, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(5, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(6, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(7, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(8, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(9, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(10, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(11, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(12, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(13, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(14, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(15, 0) == Approx(0.0).epsilon(Tolerance));

      REQUIRE(gradsf(0, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(-0.5).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(0.5).epsilon(Tolerance));
      REQUIRE(gradsf(3, 1) == Approx(0.0).epsilon(Tolerance));

      REQUIRE(gradsf(4, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(5, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(6, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(7, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(8, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(9, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(10, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(11, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(12, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(13, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(14, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(15, 1) == Approx(0.0).epsilon(Tolerance));
    }
  }
}

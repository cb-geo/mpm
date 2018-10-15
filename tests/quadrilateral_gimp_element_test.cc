// Quadrilateral element test
#include <memory>

#include "catch.hpp"

#include "quadrilateral_gimp_element.h"

//! \brief Check quadrilateral element class
TEST_CASE("Quadrilateral gimp elements are checked",
          "[quad][element][2D][gimp]") {
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

      // Coordinate location of point (x,y)
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();
      // Particle size (x,y)
      Eigen::Matrix<double, Dim, 1> psize;
      psize.setZero();
      // Deformation gradient
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
      // Coordinate location of point (x,y)
      Eigen::Matrix<double, Dim, 1> coords;
      coords << -1, -1;
      // Particle size (x,y)
      Eigen::Matrix<double, Dim, 1> psize;
      psize.setZero();
      // Deformation gradient
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
      // Coordinate location of point (x,y)
      Eigen::Matrix<double, Dim, 1> coords;
      coords << 1, 1;
      // Particle size (x,y)
      Eigen::Matrix<double, Dim, 1> psize;
      psize.setZero();
      // Deformation gradient
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

    // Coordinates is (-0.8,-0.8) Size is (0.25,0.25)
    SECTION(
        "16 Node quadrilateral element matrix for coordinate(-0.8,-0.8), Size "
        "(0.25,0.25)") {
      // Location of point (x,y)
      Eigen::Matrix<double, Dim, 1> coords;
      coords << -0.8, -0.8;
      // Size of particle (x,y)
      Eigen::Matrix<double, Dim, 1> psize;
      psize << 0.25, 0.25;
      // Deformarion gradient
      Eigen::Matrix<double, Dim, 1> defgrad;
      defgrad.setZero();

      auto shapefn = quad->shapefn(coords, psize, defgrad);

      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      REQUIRE(shapefn(0) == Approx(0.80550625).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(0.090871875).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(0.0102515625).epsilon(Tolerance));
      REQUIRE(shapefn(3) == Approx(0.090871875).epsilon(Tolerance));
      REQUIRE(shapefn(4) == Approx(1.5625e-06).epsilon(Tolerance));
      REQUIRE(shapefn(5) == Approx(0.001121875).epsilon(Tolerance));
      REQUIRE(shapefn(6) == Approx(0.0001265625).epsilon(Tolerance));
      REQUIRE(shapefn(7) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(8) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(9) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(10) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(11) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(12) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(13) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(14) == Approx(0.0001265625).epsilon(Tolerance));
      REQUIRE(shapefn(15) == Approx(0.001121875).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(coords, psize, defgrad);
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      REQUIRE(gradsf(0, 0) == Approx(-0.359).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(0.403875).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(0.0455625).epsilon(Tolerance));
      REQUIRE(gradsf(3, 0) == Approx(-0.0405).epsilon(Tolerance));

      REQUIRE(gradsf(4, 0) == Approx(-6.25e-05).epsilon(Tolerance));
      REQUIRE(gradsf(5, 0) == Approx(-0.0005).epsilon(Tolerance));
      REQUIRE(gradsf(6, 0) == Approx(0.0005625).epsilon(Tolerance));
      REQUIRE(gradsf(7, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(8, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(9, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(10, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(11, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(12, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(13, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(14, 0) == Approx(-0.0050625).epsilon(Tolerance));
      REQUIRE(gradsf(15, 0) == Approx(-0.044875).epsilon(Tolerance));

      REQUIRE(gradsf(0, 1) == Approx(-0.359).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(-0.0405).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(0.0455625).epsilon(Tolerance));
      REQUIRE(gradsf(3, 1) == Approx(0.403875).epsilon(Tolerance));

      REQUIRE(gradsf(4, 1) == Approx(-6.25e-05).epsilon(Tolerance));
      REQUIRE(gradsf(5, 1) == Approx(-0.044875).epsilon(Tolerance));
      REQUIRE(gradsf(6, 1) == Approx(-0.0050625).epsilon(Tolerance));
      REQUIRE(gradsf(7, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(8, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(9, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(10, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(11, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(12, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(13, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(14, 1) == Approx(0.0005625).epsilon(Tolerance));
      REQUIRE(gradsf(15, 1) == Approx(-0.0005).epsilon(Tolerance));
    }
    // Coordinates is (0.8,0.8) Size is (0.25,0.25)
    SECTION(
        "16 Node quadrilateral element matrix for coordinate(0.8,0.8), Size "
        "(0.25,0.25)") {
      // Location of point (x,y)
      Eigen::Matrix<double, Dim, 1> coords;
      coords << 0.8, 0.8;
      // Size of particle (x,y)
      Eigen::Matrix<double, Dim, 1> psize;
      psize << 0.25, 0.25;
      // Deformarion gradient
      Eigen::Matrix<double, Dim, 1> defgrad;
      defgrad.setZero();

      auto shapefn = quad->shapefn(coords, psize, defgrad);

      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      REQUIRE(shapefn(0) == Approx(0.0102515625).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(0.090871875).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(0.80550625).epsilon(Tolerance));
      REQUIRE(shapefn(3) == Approx(0.090871875).epsilon(Tolerance));
      REQUIRE(shapefn(4) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(5) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(6) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(7) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(8) == Approx(0.0001265625).epsilon(Tolerance));
      REQUIRE(shapefn(9) == Approx(0.001121875).epsilon(Tolerance));
      REQUIRE(shapefn(10) == Approx(1.5625e-06).epsilon(Tolerance));
      REQUIRE(shapefn(11) == Approx(0.001121875).epsilon(Tolerance));
      REQUIRE(shapefn(12) == Approx(0.0001265625).epsilon(Tolerance));
      REQUIRE(shapefn(13) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(14) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(15) == Approx(0.0).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(coords, psize, defgrad);
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      REQUIRE(gradsf(0, 0) == Approx(-0.0455625).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(0.0405).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(0.359).epsilon(Tolerance));
      REQUIRE(gradsf(3, 0) == Approx(-0.403875).epsilon(Tolerance));

      REQUIRE(gradsf(4, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(5, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(6, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(7, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(8, 0) == Approx(0.0050625).epsilon(Tolerance));
      REQUIRE(gradsf(9, 0) == Approx(0.044875).epsilon(Tolerance));
      REQUIRE(gradsf(10, 0) == Approx(6.25e-05).epsilon(Tolerance));
      REQUIRE(gradsf(11, 0) == Approx(0.0005).epsilon(Tolerance));
      REQUIRE(gradsf(12, 0) == Approx(-0.0005625).epsilon(Tolerance));
      REQUIRE(gradsf(13, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(14, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(15, 0) == Approx(0.0).epsilon(Tolerance));

      REQUIRE(gradsf(0, 1) == Approx(-0.0455625).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(-0.403875).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(0.359).epsilon(Tolerance));
      REQUIRE(gradsf(3, 1) == Approx(0.0405).epsilon(Tolerance));

      REQUIRE(gradsf(4, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(5, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(6, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(7, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(8, 1) == Approx(-0.0005625).epsilon(Tolerance));
      REQUIRE(gradsf(9, 1) == Approx(0.0005).epsilon(Tolerance));
      REQUIRE(gradsf(10, 1) == Approx(6.25e-05).epsilon(Tolerance));
      REQUIRE(gradsf(11, 1) == Approx(0.044875).epsilon(Tolerance));
      REQUIRE(gradsf(12, 1) == Approx(0.0050625).epsilon(Tolerance));
      REQUIRE(gradsf(13, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(14, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(15, 1) == Approx(0.0).epsilon(Tolerance));
    }

    // Coordinates is (0,0)
    SECTION("Four noded local sf quadrilateral element for coordinates(0,0)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();
      auto shapefn = quad->shapefn_local(coords, Eigen::Vector2d::Zero(),
                                         Eigen::Vector2d::Zero());

      // Check shape function
      REQUIRE(shapefn.size() == 4);

      REQUIRE(shapefn(0) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(3) == Approx(0.25).epsilon(Tolerance));
    }

    // Check Jacobian
    SECTION("16-noded quadrilateral Jacobian with deformation gradient") {
      Eigen::Matrix<double, 16, Dim> coords;
      // clang-format off
      // clang-format off
      coords << -1., -1.,
                 1., -1.,
                 1.,  1.,
                -1.,  1.,
                -3., -3.,
                -1., -3.,
                 1., -3.,
                 3., -3.,
                 3., -1.,
                 3.,  1.,
                 3.,  3.,
                 1.,  3.,
                -1.,  3.,
                -3.,  3.,
                -3.,  1.,
                -3., -1.;
      // clang-format on

      Eigen::Matrix<double, Dim, 1> psize;
      psize.setZero();
      Eigen::Matrix<double, Dim, 1> defgrad;
      defgrad.setZero();

      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0., 0.;

      Eigen::Matrix<double, Dim, Dim> jacobian;
      // clang-format off
      jacobian << 1.0, 0,
                 0, 1.0;
      // clang-format on

      // Get Jacobian
      auto jac = quad->jacobian(xi, coords, psize, defgrad);

      // Check size of jacobian
      REQUIRE(jac.size() == jacobian.size());

      // Check Jacobian
      for (unsigned i = 0; i < Dim; ++i)
        for (unsigned j = 0; j < Dim; ++j)
          REQUIRE(jac(i, j) == Approx(jacobian(i, j)).epsilon(Tolerance));
    }

    // Check local Jacobian
    SECTION(
        "Four noded quadrilateral local Jacobian for local "
        "coordinates(0.5,0.5)") {
      Eigen::Matrix<double, 4, Dim> coords;
      // clang-format off
      coords << 2., 1.,
                4., 2.,
                2., 4.,
                1., 3.;
      // clang-format on

      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0.5, 0.5;

      Eigen::Matrix<double, Dim, Dim> jacobian;
      // clang-format off
      jacobian << 0.625, 0.5,
                 -0.875, 1.0;
      // clang-format on

      // Get Jacobian
      auto jac = quad->jacobian_local(xi, coords, Eigen::Vector2d::Zero(),
                                      Eigen::Vector2d::Zero());

      // Check size of jacobian
      REQUIRE(jac.size() == jacobian.size());

      // Check Jacobian
      for (unsigned i = 0; i < Dim; ++i)
        for (unsigned j = 0; j < Dim; ++j)
          REQUIRE(jac(i, j) == Approx(jacobian(i, j)).epsilon(Tolerance));
    }

    // Check BMatrix with deformation gradient
    SECTION(
        "Four noded quadrilateral B-matrix cell with deformation gradient") {
      // Reference coordinates
      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0.5, 0.5;

      // Nodal coordinates
      Eigen::Matrix<double, 16, Dim> coords;
      // clang-format off
      coords <<  1.,1.,
                 2.,1.,
                 2.,2.,
                 1.,2.,
                 0,0,
                 1.,0,
                 2.,0,
                 3.,0,
                 3.,1.,
                 3.,2.,
                 3.,3.,
                 2.,3.,
                 1.,3.,
                 0,3.,
                 0,2.,
                 0,1.;
      // clang-format on

      Eigen::Matrix<double, Dim, 1> psize;
      psize.setZero();
      Eigen::Matrix<double, Dim, 1> defgrad;
      defgrad.setZero();

      // Get B-Matrix
      auto bmatrix = quad->bmatrix(xi, coords, psize, defgrad);

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(xi, psize, defgrad);
      gradsf *= 2;

      // Check size of B-matrix
      REQUIRE(bmatrix.size() == nfunctions);

      for (unsigned i = 0; i < nfunctions; ++i) {
        REQUIRE(bmatrix.at(i)(0, 0) == Approx(gradsf(i, 0)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(0, 1) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 0) == Approx(0.).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(1, 1) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 0) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(bmatrix.at(i)(2, 1) == Approx(gradsf(i, 0)).epsilon(Tolerance));
      }
    }
    SECTION("Center cell gimp element volume") {
      // Check element volume
      REQUIRE(quad->unit_element_volume() == Approx(4).epsilon(Tolerance));
    }
  }
}

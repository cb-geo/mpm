// hexrilateral element test
#include <memory>

#include "catch.hpp"

#include "hexahedron_gimp_element.h"

//! \brief Check hexrilateral element class
TEST_CASE("Hexahedron gimp elements are checked", "[hex][element][3D][gimp]") {
  const unsigned Dim = 3;
  const double Tolerance = 1.E-7;

  //! Check for center element nodes
  SECTION("64 Node hexrilateral GIMP Element") {
    const unsigned nfunctions = 64;
    std::shared_ptr<mpm::Element<Dim>> hex =
        std::make_shared<mpm::HexahedronGIMPElement<Dim, nfunctions>>();

    // Check degree
    REQUIRE(hex->degree() == mpm::ElementDegree::Linear);

    // Coordinates is (0,0,0) Size is (0,0,0)
    SECTION(
        "64 Node hexrilateral element matrix for coordinate(0.5,0.5,0.5), Size "
        "(0,0,0)") {

      // Coordinate location of point (x,y)
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();
      // Particle size (x,y)
      Eigen::Matrix<double, Dim, 1> psize;
      psize.setZero();
      // Deformation gradient
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
      REQUIRE(shapefn(8) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(9) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(10) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(11) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(12) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(13) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(14) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(15) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(16) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(17) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(18) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(19) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(20) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(21) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(22) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(23) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(24) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(25) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(26) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(27) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(28) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(29) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(30) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(31) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(32) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(33) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(34) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(35) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(36) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(37) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(38) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(39) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(40) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(41) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(42) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(43) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(44) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(45) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(46) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(47) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(48) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(49) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(50) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(51) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(52) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(53) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(54) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(55) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(56) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(57) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(58) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(59) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(60) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(61) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(62) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(63) == Approx(0).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = hex->grad_shapefn(coords, psize, defgrad);
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      REQUIRE(gradsf(0, 0) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(0, 1) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(0, 2) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(1, 2) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(2, 2) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(3, 0) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(3, 1) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(3, 2) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(4, 0) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(4, 1) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(4, 2) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(5, 0) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(5, 1) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(5, 2) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(6, 0) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(6, 1) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(6, 2) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(7, 0) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(7, 1) == Approx(-0.125).epsilon(Tolerance));
      REQUIRE(gradsf(7, 2) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(gradsf(8, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(8, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(8, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(9, 0) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(9, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(9, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(10, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(10, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(10, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(11, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(11, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(11, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(12, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(12, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(12, 2) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(13, 0) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(13, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(13, 2) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(14, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(14, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(14, 2) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(15, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(15, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(15, 2) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(16, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(16, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(16, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(17, 0) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(17, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(17, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(18, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(18, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(18, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(19, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(19, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(19, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(20, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(20, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(20, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(21, 0) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(21, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(21, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(22, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(22, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(22, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(23, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(23, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(23, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(24, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(24, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(24, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(25, 0) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(25, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(25, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(26, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(26, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(26, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(27, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(27, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(27, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(28, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(28, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(28, 2) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(29, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(29, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(29, 2) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(30, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(30, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(30, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(31, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(31, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(31, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(32, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(32, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(32, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(33, 0) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(33, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(33, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(34, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(34, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(34, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(35, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(35, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(35, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(36, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(36, 1) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(36, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(37, 0) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(37, 1) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(37, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(38, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(38, 1) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(38, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(39, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(39, 1) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(39, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(40, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(40, 1) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(40, 2) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(41, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(41, 1) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(41, 2) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(42, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(42, 1) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(42, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(43, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(43, 1) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(43, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(44, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(44, 1) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(44, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(45, 0) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(45, 1) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(45, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(46, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(46, 1) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(46, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(47, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(47, 1) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(47, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(48, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(48, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(48, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(49, 0) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(49, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(49, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(50, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(50, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(50, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(51, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(51, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(51, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(52, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(52, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(52, 2) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(53, 0) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(53, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(53, 2) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(54, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(54, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(54, 2) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(55, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(55, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(55, 2) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(56, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(56, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(56, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(57, 0) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(57, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(57, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(58, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(58, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(58, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(59, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(59, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(59, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(60, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(60, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(60, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(61, 0) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(61, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(61, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(62, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(62, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(62, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(63, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(63, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(63, 2) == Approx(0).epsilon(Tolerance));
    }

    // Coordinates is (0.8,0.8,0.8 Size is (0.5,0.5,0.5)
    SECTION(
        "64 Node hexrilateral element matrix for coordinate(0.8,0.8,0.8), Size "
        "(0.5,0.5,0.5)") {

      // Coordinate location of point (x,y)
      Eigen::Matrix<double, Dim, 1> coords;
      coords << 0.8, 0.8, 0.8;
      // Size of particle (x,y)
      Eigen::Matrix<double, Dim, 1> psize;
      psize << 0.5, 0.5, 0.5;
      // Deformation gradient
      Eigen::Matrix<double, Dim, 1> defgrad;
      defgrad.setZero();

      auto shapefn = hex->shapefn(coords, psize, defgrad);

      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      REQUIRE(shapefn(0) == Approx(0.00920077734).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(0.0815575078).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(0.722941859).epsilon(Tolerance));
      REQUIRE(shapefn(3) == Approx(0.0815575078).epsilon(Tolerance));
      REQUIRE(shapefn(4) == Approx(0.0010379707).epsilon(Tolerance));
      REQUIRE(shapefn(5) == Approx(0.00920077734).epsilon(Tolerance));
      REQUIRE(shapefn(6) == Approx(0.0815575078).epsilon(Tolerance));
      REQUIRE(shapefn(7) == Approx(0.00920077734).epsilon(Tolerance));
      REQUIRE(shapefn(8) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(9) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(10) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(11) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(12) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(13) == Approx(1.28144531e-05).epsilon(Tolerance));
      REQUIRE(shapefn(14) == Approx(0.000113589844).epsilon(Tolerance));
      REQUIRE(shapefn(15) == Approx(1.58203125e-07).epsilon(Tolerance));
      REQUIRE(shapefn(16) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(17) == Approx(0.000113589844).epsilon(Tolerance));
      REQUIRE(shapefn(18) == Approx(0.00100688281).epsilon(Tolerance));
      REQUIRE(shapefn(19) == Approx(1.40234375e-06).epsilon(Tolerance));
      REQUIRE(shapefn(20) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(21) == Approx(1.58203125e-07).epsilon(Tolerance));
      REQUIRE(shapefn(22) == Approx(1.40234375e-06).epsilon(Tolerance));
      REQUIRE(shapefn(23) == Approx(1.953125e-09).epsilon(Tolerance));
      REQUIRE(shapefn(24) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(25) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(26) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(27) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(28) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(29) == Approx(0.000113589844).epsilon(Tolerance));
      REQUIRE(shapefn(30) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(31) == Approx(0.00100688281).epsilon(Tolerance));
      REQUIRE(shapefn(32) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(33) == Approx(0.000113589844).epsilon(Tolerance));
      REQUIRE(shapefn(34) == Approx(0.00100688281).epsilon(Tolerance));
      REQUIRE(shapefn(35) == Approx(1.40234375e-06).epsilon(Tolerance));
      REQUIRE(shapefn(36) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(37) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(38) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(39) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(40) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(41) == Approx(1.28144531e-05).epsilon(Tolerance));
      REQUIRE(shapefn(42) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(43) == Approx(0.000113589844).epsilon(Tolerance));
      REQUIRE(shapefn(44) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(45) == Approx(1.28144531e-05).epsilon(Tolerance));
      REQUIRE(shapefn(46) == Approx(0.000113589844).epsilon(Tolerance));
      REQUIRE(shapefn(47) == Approx(1.58203125e-07).epsilon(Tolerance));
      REQUIRE(shapefn(48) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(49) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(50) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(51) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(52) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(53) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(54) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(55) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(56) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(57) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(58) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(59) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(60) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(61) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(62) == Approx(0).epsilon(Tolerance));
      REQUIRE(shapefn(63) == Approx(0).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = hex->grad_shapefn(coords, psize, defgrad);
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      REQUIRE(gradsf(0, 0) == Approx(-0.0408923437).epsilon(Tolerance));
      REQUIRE(gradsf(0, 1) == Approx(0.004100625).epsilon(Tolerance));
      REQUIRE(gradsf(0, 2) == Approx(-0.0408923437).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(0.03634875).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(0.03634875).epsilon(Tolerance));
      REQUIRE(gradsf(1, 2) == Approx(-0.362477812).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(0.3222025).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(0.3222025).epsilon(Tolerance));
      REQUIRE(gradsf(2, 2) == Approx(0.3222025).epsilon(Tolerance));
      REQUIRE(gradsf(3, 0) == Approx(-0.362477812).epsilon(Tolerance));
      REQUIRE(gradsf(3, 1) == Approx(0.03634875).epsilon(Tolerance));
      REQUIRE(gradsf(3, 2) == Approx(0.03634875).epsilon(Tolerance));
      REQUIRE(gradsf(4, 0) == Approx(-0.00461320312).epsilon(Tolerance));
      REQUIRE(gradsf(4, 1) == Approx(-0.00461320312).epsilon(Tolerance));
      REQUIRE(gradsf(4, 2) == Approx(-0.00461320312).epsilon(Tolerance));
      REQUIRE(gradsf(5, 0) == Approx(0.004100625).epsilon(Tolerance));
      REQUIRE(gradsf(5, 1) == Approx(-0.0408923437).epsilon(Tolerance));
      REQUIRE(gradsf(5, 2) == Approx(-0.0408923437).epsilon(Tolerance));
      REQUIRE(gradsf(6, 0) == Approx(0.03634875).epsilon(Tolerance));
      REQUIRE(gradsf(6, 1) == Approx(-0.362477812).epsilon(Tolerance));
      REQUIRE(gradsf(6, 2) == Approx(0.03634875).epsilon(Tolerance));
      REQUIRE(gradsf(7, 0) == Approx(-0.0408923437).epsilon(Tolerance));
      REQUIRE(gradsf(7, 1) == Approx(-0.0408923437).epsilon(Tolerance));
      REQUIRE(gradsf(7, 2) == Approx(0.004100625).epsilon(Tolerance));
      REQUIRE(gradsf(8, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(8, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(8, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(9, 0) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(9, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(9, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(10, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(10, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(10, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(11, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(11, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(11, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(12, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(12, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(12, 2) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(13, 0) == Approx(-5.6953125e-05).epsilon(Tolerance));
      REQUIRE(gradsf(13, 1) == Approx(0.000512578125).epsilon(Tolerance));
      REQUIRE(gradsf(13, 2) == Approx(-5.6953125e-05).epsilon(Tolerance));
      REQUIRE(gradsf(14, 0) == Approx(5.0625e-05).epsilon(Tolerance));
      REQUIRE(gradsf(14, 1) == Approx(0.00454359375).epsilon(Tolerance));
      REQUIRE(gradsf(14, 2) == Approx(-0.00050484375).epsilon(Tolerance));
      REQUIRE(gradsf(15, 0) == Approx(6.328125e-06).epsilon(Tolerance));
      REQUIRE(gradsf(15, 1) == Approx(6.328125e-06).epsilon(Tolerance));
      REQUIRE(gradsf(15, 2) == Approx(-7.03125e-07).epsilon(Tolerance));
      REQUIRE(gradsf(16, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(16, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(16, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(17, 0) == Approx(-0.00050484375).epsilon(Tolerance));
      REQUIRE(gradsf(17, 1) == Approx(0.00454359375).epsilon(Tolerance));
      REQUIRE(gradsf(17, 2) == Approx(5.0625e-05).epsilon(Tolerance));
      REQUIRE(gradsf(18, 0) == Approx(0.00044875).epsilon(Tolerance));
      REQUIRE(gradsf(18, 1) == Approx(0.0402753125).epsilon(Tolerance));
      REQUIRE(gradsf(18, 2) == Approx(0.00044875).epsilon(Tolerance));
      REQUIRE(gradsf(19, 0) == Approx(5.609375e-05).epsilon(Tolerance));
      REQUIRE(gradsf(19, 1) == Approx(5.609375e-05).epsilon(Tolerance));
      REQUIRE(gradsf(19, 2) == Approx(6.25e-07).epsilon(Tolerance));
      REQUIRE(gradsf(20, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(20, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(20, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(21, 0) == Approx(-7.03125e-07).epsilon(Tolerance));
      REQUIRE(gradsf(21, 1) == Approx(6.328125e-06).epsilon(Tolerance));
      REQUIRE(gradsf(21, 2) == Approx(6.328125e-06).epsilon(Tolerance));
      REQUIRE(gradsf(22, 0) == Approx(6.25e-07).epsilon(Tolerance));
      REQUIRE(gradsf(22, 1) == Approx(5.609375e-05).epsilon(Tolerance));
      REQUIRE(gradsf(22, 2) == Approx(5.609375e-05).epsilon(Tolerance));
      REQUIRE(gradsf(23, 0) == Approx(7.8125e-08).epsilon(Tolerance));
      REQUIRE(gradsf(23, 1) == Approx(7.8125e-08).epsilon(Tolerance));
      REQUIRE(gradsf(23, 2) == Approx(7.8125e-08).epsilon(Tolerance));
      REQUIRE(gradsf(24, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(24, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(24, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(25, 0) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(25, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(25, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(26, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(26, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(26, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(27, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(27, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(27, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(28, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(28, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(28, 2) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(29, 0) == Approx(0.00454359375).epsilon(Tolerance));
      REQUIRE(gradsf(29, 1) == Approx(5.0625e-05).epsilon(Tolerance));
      REQUIRE(gradsf(29, 2) == Approx(-0.00050484375).epsilon(Tolerance));
      REQUIRE(gradsf(30, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(30, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(30, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(31, 0) == Approx(0.0402753125).epsilon(Tolerance));
      REQUIRE(gradsf(31, 1) == Approx(0.00044875).epsilon(Tolerance));
      REQUIRE(gradsf(31, 2) == Approx(0.00044875).epsilon(Tolerance));
      REQUIRE(gradsf(32, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(32, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(32, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(33, 0) == Approx(-0.00050484375).epsilon(Tolerance));
      REQUIRE(gradsf(33, 1) == Approx(5.0625e-05).epsilon(Tolerance));
      REQUIRE(gradsf(33, 2) == Approx(0.00454359375).epsilon(Tolerance));
      REQUIRE(gradsf(34, 0) == Approx(0.00044875).epsilon(Tolerance));
      REQUIRE(gradsf(34, 1) == Approx(0.00044875).epsilon(Tolerance));
      REQUIRE(gradsf(34, 2) == Approx(0.0402753125).epsilon(Tolerance));
      REQUIRE(gradsf(35, 0) == Approx(5.609375e-05).epsilon(Tolerance));
      REQUIRE(gradsf(35, 1) == Approx(6.25e-07).epsilon(Tolerance));
      REQUIRE(gradsf(35, 2) == Approx(5.609375e-05).epsilon(Tolerance));
      REQUIRE(gradsf(36, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(36, 1) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(36, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(37, 0) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(37, 1) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(37, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(38, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(38, 1) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(38, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(39, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(39, 1) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(39, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(40, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(40, 1) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(40, 2) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(41, 0) == Approx(0.000512578125).epsilon(Tolerance));
      REQUIRE(gradsf(41, 1) == Approx(-5.6953125e-05).epsilon(Tolerance));
      REQUIRE(gradsf(41, 2) == Approx(-5.6953125e-05).epsilon(Tolerance));
      REQUIRE(gradsf(42, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(42, 1) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(42, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(43, 0) == Approx(0.00454359375).epsilon(Tolerance));
      REQUIRE(gradsf(43, 1) == Approx(-0.00050484375).epsilon(Tolerance));
      REQUIRE(gradsf(43, 2) == Approx(5.0625e-05).epsilon(Tolerance));
      REQUIRE(gradsf(44, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(44, 1) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(44, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(45, 0) == Approx(-5.6953125e-05).epsilon(Tolerance));
      REQUIRE(gradsf(45, 1) == Approx(-5.6953125e-05).epsilon(Tolerance));
      REQUIRE(gradsf(45, 2) == Approx(0.000512578125).epsilon(Tolerance));
      REQUIRE(gradsf(46, 0) == Approx(5.0625e-05).epsilon(Tolerance));
      REQUIRE(gradsf(46, 1) == Approx(-0.00050484375).epsilon(Tolerance));
      REQUIRE(gradsf(46, 2) == Approx(0.00454359375).epsilon(Tolerance));
      REQUIRE(gradsf(47, 0) == Approx(6.328125e-06).epsilon(Tolerance));
      REQUIRE(gradsf(47, 1) == Approx(-7.03125e-07).epsilon(Tolerance));
      REQUIRE(gradsf(47, 2) == Approx(6.328125e-06).epsilon(Tolerance));
      REQUIRE(gradsf(48, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(48, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(48, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(49, 0) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(49, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(49, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(50, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(50, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(50, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(51, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(51, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(51, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(52, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(52, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(52, 2) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(53, 0) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(53, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(53, 2) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(54, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(54, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(54, 2) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(55, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(55, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(55, 2) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(56, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(56, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(56, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(57, 0) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(57, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(57, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(58, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(58, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(58, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(59, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(59, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(59, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(60, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(60, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(60, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(61, 0) == Approx(-0).epsilon(Tolerance));
      REQUIRE(gradsf(61, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(61, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(62, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(62, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(62, 2) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(63, 0) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(63, 1) == Approx(0).epsilon(Tolerance));
      REQUIRE(gradsf(63, 2) == Approx(0).epsilon(Tolerance));
    }

    // Coordinates is (0,0,0)
    SECTION("Eight noded local sf hexahedron element for coordinates(0,0,0)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();
      auto shapefn = hex->shapefn_local(coords, Eigen::Vector3d::Zero(),
                                        Eigen::Vector3d::Zero());

      // Check shape function
      REQUIRE(shapefn.size() == 8);

      REQUIRE(shapefn(0) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(shapefn(3) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(shapefn(4) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(shapefn(5) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(shapefn(6) == Approx(0.125).epsilon(Tolerance));
      REQUIRE(shapefn(7) == Approx(0.125).epsilon(Tolerance));
    }

    // Check Jacobian
    SECTION("64-noded hexrilateral Jacobian with deformation gradient") {

      Eigen::Matrix<double, 64, Dim> coords;
      coords << -1., 1., -1, 1., 1., -1, 1., 1., 1, -1., 1., 1, -1., -1., -1,
          1., -1., -1, 1., -1., 1, -1., -1., 1, -3., 3, -3, -1., 3, -3, 1., 3,
          -3, 3., 3, -3, -3., 3, -1, -1., 3, -1, 1., 3, -1, 3., 3, -1, -3., 3,
          1, -1., 3, 1, 1., 3, 1, 3., 3, 1, -3., 3, 3, -1., 3, 3, 1., 3, 3, 3.,
          3, 3, -3., 1., -3, -1., 1., -3, 1., 1., -3, 3., 1., -3, -3., 1., -1,
          3., 1., -1, -3., 1., 1, 3., 1., 1, -3., 1., 3, -1., 1., 3, 1., 1., 3,
          3., 1., 3, -3., -1., -3, -1., -1., -3, 1., -1., -3, 3., -1., -3, -3.,
          -1., -1, 3., -1., -1, -3., -1., 1, 3., -1., 1, -3., -1., 3, -1., -1.,
          3, 1., -1., 3, 3., -1., 3, -3., -3., -3, -1., -3., -3, 1., -3., -3,
          3., -3., -3, -3., -3., -1, -1., -3., -1, 1., -3., -1, 3., -3., -1,
          -3., -3., 1, -1., -3., 1, 1., -3., 1, 3., -3., 1, -3., -3., 3, -1.,
          -3., 3, 1., -3., 3, 3., -3., 3;

      Eigen::Matrix<double, Dim, 1> psize;
      psize.setZero();
      Eigen::Matrix<double, Dim, 1> defgrad;
      defgrad.setZero();

      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0., 0., 0.;

      Eigen::Matrix<double, Dim, Dim> jacobian;
      // clang-format off
      jacobian << 1., 0., 0.,
                  0., 1., 0.,
                  0., 0., 1;
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
    SECTION("64 noded hexahedron B-matrix cell for coordinates(0, 0, 0)") {
      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0.0, 0.0, 0.0;

      Eigen::Matrix<double, 64, Dim> coords;
      coords << -1., 1., -1, 1., 1., -1, 1., 1., 1, -1., 1., 1, -1., -1., -1,
          1., -1., -1, 1., -1., 1, -1., -1., 1, -3., 3, -3, -1., 3, -3, 1., 3,
          -3, 3., 3, -3, -3., 3, -1, -1., 3, -1, 1., 3, -1, 3., 3, -1, -3., 3,
          1, -1., 3, 1, 1., 3, 1, 3., 3, 1, -3., 3, 3, -1., 3, 3, 1., 3, 3, 3.,
          3, 3, -3., 1., -3, -1., 1., -3, 1., 1., -3, 3., 1., -3, -3., 1., -1,
          3., 1., -1, -3., 1., 1, 3., 1., 1, -3., 1., 3, -1., 1., 3, 1., 1., 3,
          3., 1., 3, -3., -1., -3, -1., -1., -3, 1., -1., -3, 3., -1., -3, -3.,
          -1., -1, 3., -1., -1, -3., -1., 1, 3., -1., 1, -3., -1., 3, -1., -1.,
          3, 1., -1., 3, 3., -1., 3, -3., -3., -3, -1., -3., -3, 1., -3., -3,
          3., -3., -3, -3., -3., -1, -1., -3., -1, 1., -3., -1, 3., -3., -1,
          -3., -3., 1, -1., -3., 1, 1., -3., 1, 3., -3., 1, -3., -3., 3, -1.,
          -3., 3, 1., -3., 3, 3., -3., 3;

      // Get B-Matrix
      auto bmatrix = hex->bmatrix(xi, coords, Eigen::Vector3d::Zero(),
                                  Eigen::Vector3d::Zero());

      // Check gradient of shape functions
      auto gradsf = hex->grad_shapefn(xi, Eigen::Vector3d::Zero(),
                                      Eigen::Vector3d::Zero());
      //  gradsf *= 2.;

      // Check dN/dx
      auto dn_dx = hex->dn_dx(xi, coords, Eigen::Vector3d::Zero(),
                              Eigen::Vector3d::Zero());
      REQUIRE(dn_dx.rows() == nfunctions);
      REQUIRE(dn_dx.cols() == Dim);
      for (unsigned i = 0; i < nfunctions; ++i) {
        REQUIRE(dn_dx(i, 0) == Approx(gradsf(i, 0)).epsilon(Tolerance));
        REQUIRE(dn_dx(i, 1) == Approx(gradsf(i, 1)).epsilon(Tolerance));
        REQUIRE(dn_dx(i, 2) == Approx(gradsf(i, 2)).epsilon(Tolerance));
      }

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

    SECTION("Center cell gimp element length") {
      // Check element length
      REQUIRE(hex->unit_element_length() == Approx(2).epsilon(Tolerance));
    }
  }
}

// Quadrilateral shape function test
#include <memory>

#include "catch.hpp"

#include "quad_shapefn.h"

//! \brief Check quadrilateral shapefn class
TEST_CASE("Quadrilateral shape functions are checked",
          "[quadsf][quad][shapefn][2D]") {
  const unsigned Dim = 2;
  const double Tolerance = 1.E-7;

  //! Check for 4 noded shape function
  SECTION("Quadrilateral shape function with four nodes") {
    const unsigned nfunctions = 4;
    std::shared_ptr<mpm::ShapeFn<Dim>> quadsf =
        std::make_shared<mpm::QuadrilateralShapeFn<Dim, nfunctions>>();

    // Check degree
    REQUIRE(quadsf->degree() == mpm::ShapeFnDegree::Linear);

    // Coordinates is (0,0)
    SECTION("Four noded quadrilateral shape function for coordinates(0,0)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();
      auto shapefn = quadsf->shapefn(coords);

      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      REQUIRE(shapefn(0) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(3) == Approx(0.25).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = quadsf->grad_shapefn(coords);
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      REQUIRE(gradsf(0, 0) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(gradsf(3, 0) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(0, 1) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(gradsf(3, 1) == Approx(0.25).epsilon(Tolerance));
    }

    // Coordinates is (-1, -1);
    SECTION("Four noded quadrilateral shape function for coordinates(-1, -1)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << -1., -1.;
      auto shapefn = quadsf->shapefn(coords);

      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      REQUIRE(shapefn(0) == Approx(1.0).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(3) == Approx(0.0).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = quadsf->grad_shapefn(coords);
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      REQUIRE(gradsf(0, 0) == Approx(-0.5).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(0.5).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(3, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(0, 1) == Approx(-0.5).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(3, 1) == Approx(0.5).epsilon(Tolerance));
    }

    // Coordinates is (1,1)
    SECTION("Four noded quadrilateral shape function for coordinates(1,1)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << 1., 1.;
      auto shapefn = quadsf->shapefn(coords);

      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      REQUIRE(shapefn(0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(1.0).epsilon(Tolerance));
      REQUIRE(shapefn(3) == Approx(0.0).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = quadsf->grad_shapefn(coords);
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      REQUIRE(gradsf(0, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(0.5).epsilon(Tolerance));
      REQUIRE(gradsf(3, 0) == Approx(-0.5).epsilon(Tolerance));

      REQUIRE(gradsf(0, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(-0.5).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(0.5).epsilon(Tolerance));
      REQUIRE(gradsf(3, 1) == Approx(0.0).epsilon(Tolerance));
    }

    // Coordinates is (0,0)
    SECTION("Four noded quadrilateral B-matrix for coordinates(0,0)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();

      // Get B-Matrix
      auto bmatrix = quadsf->bmatrix(coords);

      // Check gradient of shape functions
      auto gradsf = quadsf->grad_shapefn(coords);

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

    // Coordinates is (0.5,0.5)
    SECTION("Four noded quadrilateral B-matrix for coordinates(0.5,0.5)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << 0.5, 0.5;

      // Get B-Matrix
      auto bmatrix = quadsf->bmatrix(coords);

      // Check gradient of shape functions
      auto gradsf = quadsf->grad_shapefn(coords);

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

    // Coordinates is (-0.5,-0.5)
    SECTION("Four noded quadrilateral B-matrix for coordinates(-0.5,-0.5)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << -0.5, -0.5;

      // Get B-Matrix
      auto bmatrix = quadsf->bmatrix(coords);

      // Check gradient of shape functions
      auto gradsf = quadsf->grad_shapefn(coords);

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

    SECTION("Four noded quadrilateral coordinates of unit cell") {
      const unsigned nfunctions = 4;

      // Coordinates of a unit cell
      Eigen::Matrix<double, nfunctions, Dim> unit_cell;
      // clang-format off
      unit_cell << -1., -1.,
                    1., -1.,
                    1.,  1.,
                   -1.,  1.;
      // clang-format on

      auto coordinates = quadsf->unit_cell_coordinates();
      REQUIRE(coordinates.rows() == nfunctions);
      REQUIRE(coordinates.cols() == Dim);
      for (unsigned i = 0; i < nfunctions; ++i) {  // Iterate through nfunctions
        for (unsigned j = 0; j < Dim; ++j) {       // Dimension
          REQUIRE(coordinates(i, j) ==
                  Approx(unit_cell(i, j)).epsilon(Tolerance));
        }
      }
    }

    SECTION("Four noded quadrilateral shape function for side indices") {
      // Check for sides indices
      Eigen::MatrixXi indices = quadsf->sides_indices();
      REQUIRE(indices.rows() == 4);
      REQUIRE(indices.cols() == 2);

      REQUIRE(indices(0, 0) == 0);
      REQUIRE(indices(0, 1) == 1);

      REQUIRE(indices(1, 0) == 1);
      REQUIRE(indices(1, 1) == 2);

      REQUIRE(indices(2, 0) == 2);
      REQUIRE(indices(2, 1) == 3);

      REQUIRE(indices(3, 0) == 3);
      REQUIRE(indices(3, 1) == 0);
    }

    SECTION("Four noded quadrilateral shape function for corner indices") {
      // Check for corner indices
      Eigen::VectorXi indices = quadsf->corner_indices();
      REQUIRE(indices.size() == 4);
      REQUIRE(indices(0) == 0);
      REQUIRE(indices(1) == 1);
      REQUIRE(indices(2) == 2);
      REQUIRE(indices(3) == 3);
    }

    SECTION("Four noded quadrilateral shape function for inhedron indices") {
      // Check for inhedron indices
      Eigen::MatrixXi indices = quadsf->inhedron_indices();
      REQUIRE(indices.rows() == 4);
      REQUIRE(indices.cols() == 2);

      REQUIRE(indices(0, 0) == 0);
      REQUIRE(indices(0, 1) == 1);

      REQUIRE(indices(1, 0) == 1);
      REQUIRE(indices(1, 1) == 2);

      REQUIRE(indices(2, 0) == 2);
      REQUIRE(indices(2, 1) == 3);

      REQUIRE(indices(3, 0) == 3);
      REQUIRE(indices(3, 1) == 0);
    }
  }

  //! Check for 8 noded shape function
  SECTION("Quadrilateral shape function with eight nodes") {
    const unsigned nfunctions = 8;
    std::shared_ptr<mpm::ShapeFn<Dim>> quadsf =
        std::make_shared<mpm::QuadrilateralShapeFn<Dim, nfunctions>>();

    // Check degree
    REQUIRE(quadsf->degree() == mpm::ShapeFnDegree::Quadratic);

    // Coordinates is (0,0)
    SECTION("Eight noded quadrilateral shape function for coordinates(0,0)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();
      auto shapefn = quadsf->shapefn(coords);

      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      REQUIRE(shapefn(0) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(shapefn(3) == Approx(-0.25).epsilon(Tolerance));
      REQUIRE(shapefn(4) == Approx(0.5).epsilon(Tolerance));
      REQUIRE(shapefn(5) == Approx(0.5).epsilon(Tolerance));
      REQUIRE(shapefn(6) == Approx(0.5).epsilon(Tolerance));
      REQUIRE(shapefn(7) == Approx(0.5).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = quadsf->grad_shapefn(coords);
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      REQUIRE(gradsf(0, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(3, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(4, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(5, 0) == Approx(0.5).epsilon(Tolerance));
      REQUIRE(gradsf(6, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(7, 0) == Approx(-0.5).epsilon(Tolerance));

      REQUIRE(gradsf(0, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(3, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(4, 1) == Approx(-0.5).epsilon(Tolerance));
      REQUIRE(gradsf(5, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(6, 1) == Approx(0.5).epsilon(Tolerance));
      REQUIRE(gradsf(7, 1) == Approx(0.0).epsilon(Tolerance));
    }

    // Coordinates is (-1,-1)
    SECTION("Eight noded quadrilateral shape function for coordinates(-1,-1)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << -1., -1.;
      auto shapefn = quadsf->shapefn(coords);

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
      auto gradsf = quadsf->grad_shapefn(coords);
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      REQUIRE(gradsf(0, 0) == Approx(-1.5).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(-0.5).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(3, 0) == Approx(0.0).epsilon(Tolerance));

      REQUIRE(gradsf(4, 0) == Approx(2.0).epsilon(Tolerance));
      REQUIRE(gradsf(5, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(6, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(7, 0) == Approx(0.0).epsilon(Tolerance));

      REQUIRE(gradsf(0, 1) == Approx(-1.5).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(3, 1) == Approx(-0.5).epsilon(Tolerance));

      REQUIRE(gradsf(4, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(5, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(6, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(7, 1) == Approx(2.0).epsilon(Tolerance));
    }

    // Coordinates is (1,1)
    SECTION("Eight noded quadrilateral shape function for coordinates(1, 1)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << 1., 1.;
      auto shapefn = quadsf->shapefn(coords);

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

      // Check gradient of shape functions
      auto gradsf = quadsf->grad_shapefn(coords);
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      REQUIRE(gradsf(0, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(1.5).epsilon(Tolerance));
      REQUIRE(gradsf(3, 0) == Approx(0.5).epsilon(Tolerance));
      REQUIRE(gradsf(4, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(5, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(6, 0) == Approx(-2.0).epsilon(Tolerance));
      REQUIRE(gradsf(7, 0) == Approx(0.0).epsilon(Tolerance));

      REQUIRE(gradsf(0, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(0.5).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(1.5).epsilon(Tolerance));
      REQUIRE(gradsf(3, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(4, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(5, 1) == Approx(-2.0).epsilon(Tolerance));
      REQUIRE(gradsf(6, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(7, 1) == Approx(0.0).epsilon(Tolerance));
    }

    // Coordinates is (0,0)
    SECTION("Eight noded quadrilateral B-matrix for coordinates(0,0)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();

      // Get B-Matrix
      auto bmatrix = quadsf->bmatrix(coords);

      // Check gradient of shape functions
      auto gradsf = quadsf->grad_shapefn(coords);

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

    // Coordinates is (0.5,0.5)
    SECTION("Eight noded quadrilateral B-matrix for coordinates(0.5,0.5)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << 0.5, 0.5;

      // Get B-Matrix
      auto bmatrix = quadsf->bmatrix(coords);

      // Check gradient of shape functions
      auto gradsf = quadsf->grad_shapefn(coords);

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

    // Coordinates is (-0.5,-0.5)
    SECTION("Eight noded quadrilateral B-matrix for coordinates(-0.5,-0.5)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << -0.5, -0.5;

      // Get B-Matrix
      auto bmatrix = quadsf->bmatrix(coords);

      // Check gradient of shape functions
      auto gradsf = quadsf->grad_shapefn(coords);

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

    SECTION("Eight noded quadrilateral coordinates of unit cell") {
      const unsigned nfunctions = 8;

      // Coordinates of a unit cell
      Eigen::Matrix<double, nfunctions, Dim> unit_cell;
      // clang-format off
      unit_cell << -1., -1.,
                    1., -1.,
                    1.,  1.,
                   -1.,  1.,
                    0., -1.,
                    1.,  0.,
                    0.,  1.,
                   -1.,  0.;
      // clang-format on

      auto coordinates = quadsf->unit_cell_coordinates();
      REQUIRE(coordinates.rows() == nfunctions);
      REQUIRE(coordinates.cols() == Dim);
      for (unsigned i = 0; i < nfunctions; ++i) {  // Iterate through nfunctions
        for (unsigned j = 0; j < Dim; ++j) {       // Dimension
          REQUIRE(coordinates(i, j) ==
                  Approx(unit_cell(i, j)).epsilon(Tolerance));
        }
      }
    }

    SECTION("Eight noded quadrilateral shape function for side indices") {
      // Check for sides indices
      Eigen::MatrixXi indices = quadsf->sides_indices();
      REQUIRE(indices.rows() == 4);
      REQUIRE(indices.cols() == 2);

      REQUIRE(indices(0, 0) == 0);
      REQUIRE(indices(0, 1) == 1);

      REQUIRE(indices(1, 0) == 1);
      REQUIRE(indices(1, 1) == 2);

      REQUIRE(indices(2, 0) == 2);
      REQUIRE(indices(2, 1) == 3);

      REQUIRE(indices(3, 0) == 3);
      REQUIRE(indices(3, 1) == 0);
    }

    SECTION("Eight noded quadrilateral shape function for corner indices") {
      // Check for corner indices
      Eigen::VectorXi indices = quadsf->corner_indices();
      REQUIRE(indices.size() == 4);
      REQUIRE(indices(0) == 0);
      REQUIRE(indices(1) == 1);
      REQUIRE(indices(2) == 2);
      REQUIRE(indices(3) == 3);
    }

    SECTION("Eight noded quadrilateral shape function for inhedron indices") {
      // Check for inhedron indices
      Eigen::MatrixXi indices = quadsf->inhedron_indices();
      REQUIRE(indices.rows() == 4);
      REQUIRE(indices.cols() == 2);

      REQUIRE(indices(0, 0) == 0);
      REQUIRE(indices(0, 1) == 1);

      REQUIRE(indices(1, 0) == 1);
      REQUIRE(indices(1, 1) == 2);

      REQUIRE(indices(2, 0) == 2);
      REQUIRE(indices(2, 1) == 3);

      REQUIRE(indices(3, 0) == 3);
      REQUIRE(indices(3, 1) == 0);
    }
  }

  //! Check for 9 noded shape function
  SECTION("Quadrilateral shape function with nine nodes") {
    const unsigned nfunctions = 9;
    std::shared_ptr<mpm::ShapeFn<Dim>> quadsf =
        std::make_shared<mpm::QuadrilateralShapeFn<Dim, nfunctions>>();

    // Check degree
    REQUIRE(quadsf->degree() == mpm::ShapeFnDegree::Quadratic);

    // Coordinates is (0,0)
    SECTION("Nine noded quadrilateral shape function for coordinates(0,0)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();
      auto shapefn = quadsf->shapefn(coords);

      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      REQUIRE(shapefn(0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(3) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(4) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(5) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(6) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(7) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(8) == Approx(1.0).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = quadsf->grad_shapefn(coords);
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      REQUIRE(gradsf(0, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(3, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(4, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(5, 0) == Approx(0.5).epsilon(Tolerance));
      REQUIRE(gradsf(6, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(7, 0) == Approx(-0.5).epsilon(Tolerance));
      REQUIRE(gradsf(8, 0) == Approx(0.0).epsilon(Tolerance));

      REQUIRE(gradsf(0, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(3, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(4, 1) == Approx(-0.5).epsilon(Tolerance));
      REQUIRE(gradsf(5, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(6, 1) == Approx(0.5).epsilon(Tolerance));
      REQUIRE(gradsf(7, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(8, 1) == Approx(0.0).epsilon(Tolerance));
    }

    // Coordinates is (-1,-1)
    SECTION("Nine noded quadrilateral shape function for coordinates(-1,-1)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << -1., -1.;
      auto shapefn = quadsf->shapefn(coords);

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

      // Check gradient of shape functions
      auto gradsf = quadsf->grad_shapefn(coords);
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      REQUIRE(gradsf(0, 0) == Approx(-1.5).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(-0.5).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(3, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(4, 0) == Approx(2.0).epsilon(Tolerance));
      REQUIRE(gradsf(5, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(6, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(7, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(8, 0) == Approx(0.0).epsilon(Tolerance));

      REQUIRE(gradsf(0, 1) == Approx(-1.5).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(3, 1) == Approx(-0.5).epsilon(Tolerance));
      REQUIRE(gradsf(4, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(5, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(6, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(7, 1) == Approx(2.0).epsilon(Tolerance));
      REQUIRE(gradsf(8, 1) == Approx(0.0).epsilon(Tolerance));
    }

    // Coordinates is (1,1)
    SECTION("Nine noded quadrilateral shape function for coordinates(1, 1)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << 1., 1.;
      auto shapefn = quadsf->shapefn(coords);

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

      // Check gradient of shape functions
      auto gradsf = quadsf->grad_shapefn(coords);
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      REQUIRE(gradsf(0, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(1.5).epsilon(Tolerance));
      REQUIRE(gradsf(3, 0) == Approx(0.5).epsilon(Tolerance));
      REQUIRE(gradsf(4, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(5, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(6, 0) == Approx(-2.0).epsilon(Tolerance));
      REQUIRE(gradsf(7, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(8, 0) == Approx(0.0).epsilon(Tolerance));

      REQUIRE(gradsf(0, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(0.5).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(1.5).epsilon(Tolerance));
      REQUIRE(gradsf(3, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(4, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(5, 1) == Approx(-2.0).epsilon(Tolerance));
      REQUIRE(gradsf(6, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(7, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(8, 1) == Approx(0.0).epsilon(Tolerance));
    }

    // Coordinates is (0,0)
    SECTION("Nine noded quadrilateral B-matrix for coordinates(0,0)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();

      // Get B-Matrix
      auto bmatrix = quadsf->bmatrix(coords);

      // Check gradient of shape functions
      auto gradsf = quadsf->grad_shapefn(coords);

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

    // Coordinates is (0.5,0.5)
    SECTION("Nine noded quadrilateral B-matrix for coordinates(0.5,0.5)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << 0.5, 0.5;

      // Get B-Matrix
      auto bmatrix = quadsf->bmatrix(coords);

      // Check gradient of shape functions
      auto gradsf = quadsf->grad_shapefn(coords);

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

    // Coordinates is (-0.5,-0.5)
    SECTION("Nine noded quadrilateral B-matrix for coordinates(-0.5,-0.5)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << -0.5, -0.5;

      // Get B-Matrix
      auto bmatrix = quadsf->bmatrix(coords);

      // Check gradient of shape functions
      auto gradsf = quadsf->grad_shapefn(coords);

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

    SECTION("Nine noded quadrilateral coordinates of unit cell") {
      const unsigned nfunctions = 9;

      // Coordinates of a unit cell
      Eigen::Matrix<double, nfunctions, Dim> unit_cell;
      // clang-format off
      unit_cell << -1., -1.,
                    1., -1.,
                    1.,  1.,
                   -1.,  1.,
                    0., -1.,
                    1.,  0.,
                    0.,  1.,
                   -1.,  0.,
                    0.,  0.;
      // clang-format on

      auto coordinates = quadsf->unit_cell_coordinates();
      REQUIRE(coordinates.rows() == nfunctions);
      REQUIRE(coordinates.cols() == Dim);
      for (unsigned i = 0; i < nfunctions; ++i) {  // Iterate through nfunctions
        for (unsigned j = 0; j < Dim; ++j) {       // Dimension
          REQUIRE(coordinates(i, j) ==
                  Approx(unit_cell(i, j)).epsilon(Tolerance));
        }
      }
    }

    SECTION("Nine noded quadrilateral shape function for side indices") {
      // Check for sides indices
      Eigen::MatrixXi indices = quadsf->sides_indices();
      REQUIRE(indices.rows() == 4);
      REQUIRE(indices.cols() == 2);

      REQUIRE(indices(0, 0) == 0);
      REQUIRE(indices(0, 1) == 1);

      REQUIRE(indices(1, 0) == 1);
      REQUIRE(indices(1, 1) == 2);

      REQUIRE(indices(2, 0) == 2);
      REQUIRE(indices(2, 1) == 3);

      REQUIRE(indices(3, 0) == 3);
      REQUIRE(indices(3, 1) == 0);
    }

    SECTION("Nine noded quadrilateral shape function for corner indices") {
      // Check for corner indices
      Eigen::VectorXi indices = quadsf->corner_indices();
      REQUIRE(indices.size() == 4);
      REQUIRE(indices(0) == 0);
      REQUIRE(indices(1) == 1);
      REQUIRE(indices(2) == 2);
      REQUIRE(indices(3) == 3);
    }

    SECTION("Nine noded quadrilateral shape function for inhedron indices") {
      // Check for inhedron indices
      Eigen::MatrixXi indices = quadsf->inhedron_indices();
      REQUIRE(indices.rows() == 4);
      REQUIRE(indices.cols() == 2);

      REQUIRE(indices(0, 0) == 0);
      REQUIRE(indices(0, 1) == 1);

      REQUIRE(indices(1, 0) == 1);
      REQUIRE(indices(1, 1) == 2);

      REQUIRE(indices(2, 0) == 2);
      REQUIRE(indices(2, 1) == 3);

      REQUIRE(indices(3, 0) == 3);
      REQUIRE(indices(3, 1) == 0);
    }
  }
}

// Triangle element test
#include <memory>

#include "catch.hpp"

#include "triangle_element.h"

//! \brief Check triangle element class
TEST_CASE("Triangle elements are checked", "[tri][element][2D]") {
  const unsigned Dim = 2;
  const double Tolerance = 1.E-7;

  //! Check for 3 noded element
  SECTION("Triangle element with three nodes") {
    const unsigned nfunctions = 3;
    std::shared_ptr<mpm::Element<Dim>> tri =
        std::make_shared<mpm::TriangleElement<Dim, nfunctions>>();

    // Check degree
    REQUIRE(tri->degree() == mpm::ElementDegree::Linear);

    // Coordinates is (0,0)
    SECTION("Three noded triangle element for coordinates(0,0)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();
      auto shapefn = tri->shapefn(coords, Eigen::Vector2d::Zero(),
                                  Eigen::Vector2d::Zero());

      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      REQUIRE(shapefn(0) == Approx(1.0).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(0.0).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = tri->grad_shapefn(coords, Eigen::Vector2d::Zero(),
                                      Eigen::Vector2d::Zero());
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      REQUIRE(gradsf(0, 0) == Approx(-1.0).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(1.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(0, 1) == Approx(-1.0).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(1.0).epsilon(Tolerance));
    }

    // Coordinates is (0.333, 0.333);
    SECTION("Three noded triangle element for coordinates(0.333, 0.333)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << 0.333, 0.333;
      auto shapefn = tri->shapefn(coords, Eigen::Vector2d::Zero(),
                                  Eigen::Vector2d::Zero());

      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      REQUIRE(shapefn(0) == Approx(0.334).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(0.333).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(0.333).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = tri->grad_shapefn(coords, Eigen::Vector2d::Zero(),
                                      Eigen::Vector2d::Zero());
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      REQUIRE(gradsf(0, 0) == Approx(-1.0).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(1.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(0, 1) == Approx(-1.0).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(1.0).epsilon(Tolerance));
    }

    // Coordinates is (0.5,0.4)
    SECTION("Three noded triangle element for coordinates(0.5,0.4)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << 0.5, 0.4;
      auto shapefn = tri->shapefn(coords, Eigen::Vector2d::Zero(),
                                  Eigen::Vector2d::Zero());

      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      REQUIRE(shapefn(0) == Approx(0.1).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(0.5).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(0.4).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = tri->grad_shapefn(coords, Eigen::Vector2d::Zero(),
                                      Eigen::Vector2d::Zero());
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      REQUIRE(gradsf(0, 0) == Approx(-1.0).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(1.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(0.0).epsilon(Tolerance));

      REQUIRE(gradsf(0, 1) == Approx(-1.0).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(1.0).epsilon(Tolerance));
    }

    // Coordinates is (0,0)
    SECTION("Three noded local sf triangle element for coordinates(0,0)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();
      auto shapefn = tri->shapefn_local(coords, Eigen::Vector2d::Zero(),
                                        Eigen::Vector2d::Zero());

      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      REQUIRE(shapefn(0) == Approx(1.0).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(0.0).epsilon(Tolerance));
    }

    // Check shapefn with deformation gradient
    SECTION("Three noded triangle element shapefn with deformation gradient") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();
      Eigen::Matrix<double, Dim, 1> psize;
      psize.setZero();
      Eigen::Matrix<double, Dim, 1> defgrad;
      defgrad.setZero();
      auto shapefn = tri->shapefn(coords, psize, defgrad);

      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      REQUIRE(shapefn(0) == Approx(1.0).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(0.0).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = tri->grad_shapefn(coords, psize, defgrad);
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      REQUIRE(gradsf(0, 0) == Approx(-1.0).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(1.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(0, 1) == Approx(-1.0).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(1.0).epsilon(Tolerance));
    }

    // Check Jacobian
    SECTION(
        "Three noded triangle Jacobian for local coordinates(0.333,0.333)") {
      Eigen::Matrix<double, 3, Dim> coords;
      // clang-format off
      coords << 2., 1.,
                4., 2.,
                2., 4.;
      // clang-format on

      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0.333, 0.333;

      Eigen::Matrix<double, Dim, Dim> jacobian;
      // clang-format off
      jacobian << 2.0, 1.0,
                  0.0, 3.0;
      // clang-format on

      // Get Jacobian
      auto jac = tri->jacobian(xi, coords, Eigen::Vector2d::Zero(),
                               Eigen::Vector2d::Zero());

      // Check size of jacobian
      REQUIRE(jac.size() == jacobian.size());

      // Check Jacobian
      for (unsigned i = 0; i < Dim; ++i)
        for (unsigned j = 0; j < Dim; ++j)
          REQUIRE(jac(i, j) == Approx(jacobian(i, j)).epsilon(Tolerance));
    }

    // Check local Jacobian
    SECTION(
        "Three noded Triangle local Jacobian for local "
        "coordinates(0.333,0.333)") {
      Eigen::Matrix<double, 3, Dim> coords;
      // clang-format off
      coords << 2., 1.,
                4., 2.,
                2., 4.;
      // clang-format on

      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0.333, 0.333;

      Eigen::Matrix<double, Dim, Dim> jacobian;
      // clang-format off
      jacobian << 2.0, 1.0,
                  0.0, 3.0;
      // clang-format on

      // Get Jacobian
      auto jac = tri->jacobian_local(xi, coords, Eigen::Vector2d::Zero(),
                                     Eigen::Vector2d::Zero());

      // Check size of jacobian
      REQUIRE(jac.size() == jacobian.size());

      // Check Jacobian
      for (unsigned i = 0; i < Dim; ++i)
        for (unsigned j = 0; j < Dim; ++j)
          REQUIRE(jac(i, j) == Approx(jacobian(i, j)).epsilon(Tolerance));
    }

    // Check Jacobian
    SECTION("Three noded triangle Jacobian with deformation gradient") {
      Eigen::Matrix<double, 3, Dim> coords;
      // clang-format off
      coords << 2., 1.,
                4., 2.,
                2., 4.;
      // clang-format on

      Eigen::Matrix<double, Dim, 1> psize;
      psize.setZero();
      Eigen::Matrix<double, Dim, 1> defgrad;
      defgrad.setZero();

      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0.333, 0.333;

      Eigen::Matrix<double, Dim, Dim> jacobian;
      // clang-format off
      jacobian << 2.0, 1.0,
                  0.0, 3.0;
      // clang-format on

      // Get Jacobian
      auto jac = tri->jacobian(xi, coords, psize, defgrad);

      // Check size of jacobian
      REQUIRE(jac.size() == jacobian.size());

      // Check Jacobian
      for (unsigned i = 0; i < Dim; ++i)
        for (unsigned j = 0; j < Dim; ++j)
          REQUIRE(jac(i, j) == Approx(jacobian(i, j)).epsilon(Tolerance));
    }

    // Coordinates is (0,0)
    SECTION("Three noded triangle B-matrix cell for coordinates(0,0)") {
      // Reference coordinates
      Eigen::Matrix<double, Dim, 1> xi;
      xi.setZero();

      // Nodal coordinates
      Eigen::Matrix<double, 3, Dim> coords;
      // clang-format off
      coords << 2., 1.,
                4., 2.,
                2., 4.;
      // clang-format on

      Eigen::Matrix<double, Dim, Dim> jacobian;
      // clang-format off
      jacobian << 2.0, 1.0,
                  0.0, 3.0;
      // clang-format on

      // Get B-Matrix
      auto bmatrix = tri->bmatrix(xi, coords, Eigen::Vector2d::Zero(),
                                  Eigen::Vector2d::Zero());

      // Check gradient of shape functions
      auto gradsf = tri->grad_shapefn(xi, Eigen::Vector2d::Zero(),
                                      Eigen::Vector2d::Zero());
      gradsf *= ((jacobian.inverse()).transpose());

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

    // Coordinates is (0.333,0.333)
    SECTION("Three noded triangle B-matrix cell for coordinates(0.333,0.333)") {
      // Reference coordinates
      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0.333, 0.333;

      // Nodal coordinates
      Eigen::Matrix<double, 3, Dim> coords;
      // clang-format off
      coords << 2., 1.,
                4., 2.,
                2., 4.;
      // clang-format on

      Eigen::Matrix<double, Dim, Dim> jacobian;
      // clang-format off
      jacobian << 2.0, 1.0,
                  0.0, 3.0;
      // clang-format on

      // Get B-Matrix
      auto bmatrix = tri->bmatrix(xi, coords, Eigen::Vector2d::Zero(),
                                  Eigen::Vector2d::Zero());

      // Check gradient of shape functions
      auto gradsf = tri->grad_shapefn(xi, Eigen::Vector2d::Zero(),
                                      Eigen::Vector2d::Zero());
      gradsf *= ((jacobian.inverse()).transpose());

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
    SECTION("Three noded triangle B-matrix cell for coordinates(0.5,0.5)") {
      // Reference coordinates
      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0.5, 0.5;

      // Nodal coordinates
      Eigen::Matrix<double, 3, Dim> coords;
      // clang-format off
      coords << 2., 1.,
                4., 2.,
                2., 4.;
      // clang-format on

      Eigen::Matrix<double, Dim, Dim> jacobian;
      // clang-format off
      jacobian << 2.0, 1.0,
                  0.0, 3.0;
      // clang-format on

      // Get B-Matrix
      auto bmatrix = tri->bmatrix(xi, coords, Eigen::Vector2d::Zero(),
                                  Eigen::Vector2d::Zero());

      // Check gradient of shape functions
      auto gradsf = tri->grad_shapefn(xi, Eigen::Vector2d::Zero(),
                                      Eigen::Vector2d::Zero());
      gradsf *= ((jacobian.inverse()).transpose());

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

    // Check BMatrix with deformation gradient
    SECTION("Three noded triangle B-matrix cell with deformation gradient") {
      // Reference coordinates
      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0.333, 0.333;

      // Nodal coordinates
      Eigen::Matrix<double, 3, Dim> coords;
      // clang-format off
      coords << 2., 1.,
                4., 2.,
                2., 4.;
      // clang-format on

      Eigen::Matrix<double, Dim, Dim> jacobian;
      // clang-format off
      jacobian << 2.0, 1.0,
                  0.0, 3.0;
      // clang-format on

      Eigen::Matrix<double, Dim, 1> psize;
      psize.setZero();
      Eigen::Matrix<double, Dim, 1> defgrad;
      defgrad.setZero();

      // Get B-Matrix
      auto bmatrix = tri->bmatrix(xi, coords, psize, defgrad);

      // Check gradient of shape functions
      auto gradsf = tri->grad_shapefn(xi, psize, defgrad);
      gradsf *= ((jacobian.inverse()).transpose());

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

    SECTION("Three noded triangle B-matrix and Jacobian failure") {
      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0., 0.;

      Eigen::Matrix<double, 3, Dim> coords;
      // clang-format off
      coords << 2., 1.,
                4., 2.,
                2., 4.;
      // clang-format on
      // Get B-Matrix
      auto bmatrix = tri->bmatrix(xi, coords, Eigen::Vector2d::Zero(),
                                  Eigen::Vector2d::Zero());
      auto jacobian = tri->jacobian(xi, coords, Eigen::Vector2d::Zero(),
                                    Eigen::Vector2d::Zero());
    }

    // Ni Nj matrix of a cell
    SECTION("Three noded triangle ni-nj-matrix") {
      std::vector<Eigen::Matrix<double, Dim, 1>> xi_s;

      Eigen::Matrix<double, Dim, 1> xi;
      xi << 1.0 / 6, 1.0 / 6;
      xi_s.emplace_back(xi);
      xi << 2.0 / 3, 1.0 / 6;
      xi_s.emplace_back(xi);
      xi << 1.0 / 6, 2.0 / 3;
      xi_s.emplace_back(xi);

      REQUIRE(xi_s.size() == 3);

      // Get Ni Nj matrix
      const auto ni_nj_matrix = tri->ni_nj_matrix(xi_s);

      // Check size of ni_nj_matrix
      REQUIRE(ni_nj_matrix.rows() == nfunctions);
      REQUIRE(ni_nj_matrix.cols() == nfunctions);

      // Sum should be equal to 1. * xi_s.size()
      REQUIRE(ni_nj_matrix.sum() ==
              Approx(1. * xi_s.size()).epsilon(Tolerance));

      Eigen::Matrix<double, 3, 3> mass;
      // clang-format off
      mass <<  0.50, 0.25, 0.25,
               0.25, 0.50, 0.25,
               0.25, 0.25, 0.50;
      // clang-format on

      // auxiliary matrices for checking its multiplication by scalar
      auto ni_nj_matrix_unit = 1.0 * ni_nj_matrix;
      auto ni_nj_matrix_zero = 0.0 * ni_nj_matrix;
      auto ni_nj_matrix_negative = -2.0 * ni_nj_matrix;
      double scalar = 21.65489;
      auto ni_nj_matrix_scalar = scalar * ni_nj_matrix;

      for (unsigned i = 0; i < nfunctions; ++i) {
        for (unsigned j = 0; j < nfunctions; ++j) {
          REQUIRE(ni_nj_matrix(i, j) == Approx(mass(i, j)).epsilon(Tolerance));
          // check multiplication by unity;
          REQUIRE(ni_nj_matrix_unit(i, j) ==
                  Approx(1.0 * mass(i, j)).epsilon(Tolerance));
          // check multiplication by zero;
          REQUIRE(ni_nj_matrix_zero(i, j) ==
                  Approx(0.0 * mass(i, j)).epsilon(Tolerance));
          // check multiplication by negative number;
          REQUIRE(ni_nj_matrix_negative(i, j) ==
                  Approx(-2.0 * mass(i, j)).epsilon(Tolerance));
          // check multiplication by an arbitrary scalar;
          REQUIRE(ni_nj_matrix_scalar(i, j) ==
                  Approx(scalar * mass(i, j)).epsilon(Tolerance));
        }
      }
    }

    // Laplace matrix of a cell
    SECTION("Three noded triangle laplace-matrix") {
      std::vector<Eigen::Matrix<double, Dim, 1>> xi_s;

      Eigen::Matrix<double, Dim, 1> xi;
      xi << 1.0 / 6, 1.0 / 6;
      xi_s.emplace_back(xi);
      xi << 2.0 / 3, 1.0 / 6;
      xi_s.emplace_back(xi);
      xi << 1.0 / 6, 2.0 / 3;
      xi_s.emplace_back(xi);

      REQUIRE(xi_s.size() == 3);

      Eigen::Matrix<double, 3, Dim> coords;
      // clang-format off
      coords << 2., 1.,
                4., 2.,
                2., 4.;
      // clang-format on

      // Get laplace matrix
      const auto laplace_matrix = tri->laplace_matrix(xi_s, coords);

      // Check size of laplace-matrix
      REQUIRE(laplace_matrix.rows() == nfunctions);
      REQUIRE(laplace_matrix.cols() == nfunctions);

      // Sum should be equal to 0.
      REQUIRE(laplace_matrix.sum() == Approx(0.).epsilon(Tolerance));

      Eigen::Matrix<double, 3, 3> laplace;
      // clang-format off
      laplace <<  0.8333333333333333, -0.6666666666666667, -0.1666666666666667,
                 -0.6666666666666667,  0.8333333333333333, -0.1666666666666667,
                 -0.1666666666666667, -0.1666666666666667,  0.3333333333333333;
      // clang-format on
      for (unsigned i = 0; i < nfunctions; ++i)
        for (unsigned j = 0; j < nfunctions; ++j)
          REQUIRE(laplace_matrix(i, j) ==
                  Approx(laplace(i, j)).epsilon(Tolerance));
    }

    SECTION("Three noded triangle coordinates of unit cell") {
      const unsigned nfunctions = 3;

      // Coordinates of a unit cell
      Eigen::Matrix<double, nfunctions, Dim> unit_cell;
      // clang-format off
      unit_cell <<  0.,  0.,
                    1.,  0.,
                    0.,  1.;
      // clang-format on

      auto coordinates = tri->unit_cell_coordinates();
      REQUIRE(coordinates.rows() == nfunctions);
      REQUIRE(coordinates.cols() == Dim);
      for (unsigned i = 0; i < nfunctions; ++i) {  // Iterate through nfunctions
        for (unsigned j = 0; j < Dim; ++j) {       // Dimension
          REQUIRE(coordinates(i, j) ==
                  Approx(unit_cell(i, j)).epsilon(Tolerance));
        }
      }
    }

    SECTION("Three noded triangle element for side indices") {
      // Check for sides indices
      Eigen::MatrixXi indices = tri->sides_indices();
      REQUIRE(indices.rows() == 3);
      REQUIRE(indices.cols() == 2);

      REQUIRE(indices(0, 0) == 0);
      REQUIRE(indices(0, 1) == 1);

      REQUIRE(indices(1, 0) == 1);
      REQUIRE(indices(1, 1) == 2);

      REQUIRE(indices(2, 0) == 2);
      REQUIRE(indices(2, 1) == 0);
    }

    SECTION("Three noded triangle element for corner indices") {
      // Check for corner indices
      Eigen::VectorXi indices = tri->corner_indices();
      REQUIRE(indices.size() == 3);
      REQUIRE(indices(0) == 0);
      REQUIRE(indices(1) == 1);
      REQUIRE(indices(2) == 2);
    }

    SECTION("Three noded triangle element for inhedron indices") {
      // Check for inhedron indices
      Eigen::MatrixXi indices = tri->inhedron_indices();
      REQUIRE(indices.rows() == 3);
      REQUIRE(indices.cols() == 2);

      REQUIRE(indices(0, 0) == 0);
      REQUIRE(indices(0, 1) == 1);

      REQUIRE(indices(1, 0) == 1);
      REQUIRE(indices(1, 1) == 2);

      REQUIRE(indices(2, 0) == 2);
      REQUIRE(indices(2, 1) == 0);
    }

    SECTION("Four noded triangle shape function for face indices") {
      // Check for face indices
      Eigen::Matrix<int, 3, 2> indices;
      // clang-format off
      indices << 0, 1,
                 1, 2,
                 2, 0;
      // clang-format on
      // Check for all face indices
      for (unsigned i = 0; i < indices.rows(); ++i) {
        const auto check_indices = tri->face_indices(i);
        REQUIRE(check_indices.rows() == 2);
        REQUIRE(check_indices.cols() == 1);

        for (unsigned j = 0; j < indices.cols(); ++j)
          REQUIRE(check_indices(j) == indices(i, j));
      }

      // Check number of faces
      REQUIRE(tri->nfaces() == 3);
    }
  }

  //! Check for 6 noded element
  SECTION("Triangle element with six nodes") {
    const unsigned nfunctions = 6;
    std::shared_ptr<mpm::Element<Dim>> tri =
        std::make_shared<mpm::TriangleElement<Dim, nfunctions>>();

    // Check degree
    REQUIRE(tri->degree() == mpm::ElementDegree::Quadratic);

    // Coordinates is (0,0)
    SECTION("Six noded triangle element for coordinates(0,0)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();
      auto shapefn = tri->shapefn(coords, Eigen::Vector2d::Zero(),
                                  Eigen::Vector2d::Zero());

      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      REQUIRE(shapefn(0) == Approx(1.0).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(3) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(4) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(5) == Approx(0.0).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = tri->grad_shapefn(coords, Eigen::Vector2d::Zero(),
                                      Eigen::Vector2d::Zero());
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      REQUIRE(gradsf(0, 0) == Approx(-3.0).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(-1.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(3, 0) == Approx(4.0).epsilon(Tolerance));
      REQUIRE(gradsf(4, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(5, 0) == Approx(0.0).epsilon(Tolerance));

      REQUIRE(gradsf(0, 1) == Approx(-3.0).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(-1.0).epsilon(Tolerance));
      REQUIRE(gradsf(3, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(4, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(5, 1) == Approx(4.0).epsilon(Tolerance));
    }

    // Coordinates is (0.5,0.5)
    SECTION("Six noded triangle element for coordinates(0.5,0.5)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << 0.5, 0.5;
      auto shapefn = tri->shapefn(coords, Eigen::Vector2d::Zero(),
                                  Eigen::Vector2d::Zero());

      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      REQUIRE(shapefn(0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(3) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(4) == Approx(1.0).epsilon(Tolerance));
      REQUIRE(shapefn(5) == Approx(0.0).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = tri->grad_shapefn(coords, Eigen::Vector2d::Zero(),
                                      Eigen::Vector2d::Zero());
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      REQUIRE(gradsf(0, 0) == Approx(1.0).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(1.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(3, 0) == Approx(-2.0).epsilon(Tolerance));
      REQUIRE(gradsf(4, 0) == Approx(2.0).epsilon(Tolerance));
      REQUIRE(gradsf(5, 0) == Approx(-2.0).epsilon(Tolerance));

      REQUIRE(gradsf(0, 1) == Approx(1.0).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(1.0).epsilon(Tolerance));
      REQUIRE(gradsf(3, 1) == Approx(-2.0).epsilon(Tolerance));
      REQUIRE(gradsf(4, 1) == Approx(2.0).epsilon(Tolerance));
      REQUIRE(gradsf(5, 1) == Approx(-2.0).epsilon(Tolerance));
    }

    // Coordinates is (0.2,0.4)
    SECTION("Six noded triangle element for coordinates(0.2, 0.4)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << 0.2, 0.4;
      auto shapefn = tri->shapefn(coords, Eigen::Vector2d::Zero(),
                                  Eigen::Vector2d::Zero());

      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      REQUIRE(shapefn(0) == Approx(-0.08).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(-0.12).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(-0.08).epsilon(Tolerance));
      REQUIRE(shapefn(3) == Approx(0.32).epsilon(Tolerance));
      REQUIRE(shapefn(4) == Approx(0.32).epsilon(Tolerance));
      REQUIRE(shapefn(5) == Approx(0.64).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = tri->grad_shapefn(coords, Eigen::Vector2d::Zero(),
                                      Eigen::Vector2d::Zero());
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      REQUIRE(gradsf(0, 0) == Approx(-0.6).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(-0.2).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(3, 0) == Approx(0.8).epsilon(Tolerance));
      REQUIRE(gradsf(4, 0) == Approx(1.6).epsilon(Tolerance));
      REQUIRE(gradsf(5, 0) == Approx(-1.6).epsilon(Tolerance));

      REQUIRE(gradsf(0, 1) == Approx(-0.6).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(0.6).epsilon(Tolerance));
      REQUIRE(gradsf(3, 1) == Approx(-0.8).epsilon(Tolerance));
      REQUIRE(gradsf(4, 1) == Approx(0.8).epsilon(Tolerance));
      REQUIRE(gradsf(5, 1) == Approx(0.0).epsilon(Tolerance));
    }

    // Coordinates is (0,0)
    SECTION("Six noded local sf triangle element for coordinates(0,0)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();
      auto shapefn = tri->shapefn_local(coords, Eigen::Vector2d::Zero(),
                                        Eigen::Vector2d::Zero());

      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      REQUIRE(shapefn(0) == Approx(1.0).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(3) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(4) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(5) == Approx(0.0).epsilon(Tolerance));
    }

    // Coordinates is (0,0)
    SECTION("Six noded triangle Jacobian for local coordinates(0.5,0.5)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();
      Eigen::Matrix<double, Dim, 1> psize;
      psize.setZero();
      Eigen::Matrix<double, Dim, 1> defgrad;
      defgrad.setZero();
      auto shapefn = tri->shapefn(coords, psize, defgrad);

      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      REQUIRE(shapefn(0) == Approx(1.0).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(3) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(4) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(5) == Approx(0.0).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = tri->grad_shapefn(coords, psize, defgrad);
      REQUIRE(gradsf.rows() == nfunctions);
      REQUIRE(gradsf.cols() == Dim);

      REQUIRE(gradsf(0, 0) == Approx(-3.0).epsilon(Tolerance));
      REQUIRE(gradsf(1, 0) == Approx(-1.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(3, 0) == Approx(4.0).epsilon(Tolerance));
      REQUIRE(gradsf(4, 0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(5, 0) == Approx(0.0).epsilon(Tolerance));

      REQUIRE(gradsf(0, 1) == Approx(-3.0).epsilon(Tolerance));
      REQUIRE(gradsf(1, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(2, 1) == Approx(-1.0).epsilon(Tolerance));
      REQUIRE(gradsf(3, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(4, 1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(gradsf(5, 1) == Approx(4.0).epsilon(Tolerance));
    }

    // Check Jacobian
    SECTION(
        "Six noded triangle Jacobian for local coordinates(0.5,0.5)") {
      Eigen::Matrix<double, 6, Dim> coords;
      // clang-format off
      coords << 2.0, 1.0,
                4.0, 2.0,
                2.0, 4.0,
                3.0, 1.5,
                3.0, 3.0,
                2.0, 2.5;
      // clang-format on

      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0.5, 0.5;

      // Jacobian result
      Eigen::Matrix<double, Dim, Dim> jacobian;
      // clang-format off
      jacobian <<  2.0, 1.0,
                   0.0, 3.0;
      // clang-format on

      // Get Jacobian
      auto jac = tri->jacobian(xi, coords, Eigen::Vector2d::Zero(),
                               Eigen::Vector2d::Zero());

      // Check size of jacobian
      REQUIRE(jac.size() == jacobian.size());

      // Check Jacobian
      for (unsigned i = 0; i < Dim; ++i)
        for (unsigned j = 0; j < Dim; ++j)
          REQUIRE(jac(i, j) == Approx(jacobian(i, j)).epsilon(Tolerance));
    }

    // Check local Jacobian
    SECTION(
        "Six noded triangle local Jacobian for local "
        "coordinates(0.5,0.5)") {
      Eigen::Matrix<double, 6, Dim> coords;
      // clang-format off
      coords << 2.0, 1.0,
                4.0, 2.0,
                2.0, 4.0,
                3.0, 1.5,
                3.0, 3.0,
                2.0, 2.5;
      // clang-format on

      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0.5, 0.5;

      // Jacobian result
      Eigen::Matrix<double, Dim, Dim> jacobian;
      // clang-format off
      jacobian <<  2.0, 1.0,
                   0.0, 3.0;
      // clang-format on

      // Get Jacobian
      auto jac = tri->jacobian_local(xi, coords, Eigen::Vector2d::Zero(),
                                     Eigen::Vector2d::Zero());

      // Check size of jacobian
      REQUIRE(jac.size() == jacobian.size());

      // Check Jacobian
      for (unsigned i = 0; i < Dim; ++i)
        for (unsigned j = 0; j < Dim; ++j)
          REQUIRE(jac(i, j) == Approx(jacobian(i, j)).epsilon(Tolerance));
    }

    // Check Jacobian
    SECTION("Six noded triangle Jacobian with deformation gradient") {
      Eigen::Matrix<double, 6, Dim> coords;
      // clang-format off
      coords << 2.0, 1.0,
                4.0, 2.0,
                2.0, 4.0,
                3.0, 1.5,
                3.0, 3.0,
                2.0, 2.5;
      // clang-format on

      Eigen::Matrix<double, Dim, 1> psize;
      psize.setZero();
      Eigen::Matrix<double, Dim, 1> defgrad;
      defgrad.setZero();

      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0.5, 0.5;

      // Jacobian result
      Eigen::Matrix<double, Dim, Dim> jacobian;
      // clang-format off
      jacobian <<  2.0, 1.0,
                   0.0, 3.0;
      // clang-format on

      // Get Jacobian
      auto jac = tri->jacobian(xi, coords, psize, defgrad);

      // Check size of jacobian
      REQUIRE(jac.size() == jacobian.size());

      // Check Jacobian
      for (unsigned i = 0; i < Dim; ++i)
        for (unsigned j = 0; j < Dim; ++j)
          REQUIRE(jac(i, j) == Approx(jacobian(i, j)).epsilon(Tolerance));
    }

    // Coordinates is (0,0)
    SECTION("Six noded triangle B-matrix cell for coordinates(0,0)") {
      // Reference coordinates
      Eigen::Matrix<double, Dim, 1> xi;
      xi.setZero();

      // Nodal coordinates
      Eigen::Matrix<double, 6, Dim> coords;
      // clang-format off
      coords << 2.0, 1.0,
                4.0, 2.0,
                2.0, 4.0,
                3.0, 1.5,
                3.0, 3.0,
                2.0, 2.5;
      // clang-format on

      Eigen::Matrix<double, Dim, Dim> jacobian;
      // clang-format off
      jacobian <<  2.0, 1.0,
                   0.0, 3.0;
      // clang-format on

      // Get B-Matrix
      auto bmatrix = tri->bmatrix(xi, coords, Eigen::Vector2d::Zero(),
                                  Eigen::Vector2d::Zero());

      // Check gradient of shape functions
      auto gradsf = tri->grad_shapefn(xi, Eigen::Vector2d::Zero(),
                                      Eigen::Vector2d::Zero());
      gradsf *= ((jacobian.inverse()).transpose());

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
    SECTION("Six noded triangle B-matrix cell for coordinates(0.5,0.5)") {
      // Reference coordinates
      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0.5, 0.5;

      // Nodal coordinates
      Eigen::Matrix<double, 6, Dim> coords;
      // clang-format off
      coords << 2.0, 1.0,
                4.0, 2.0,
                2.0, 4.0,
                3.0, 1.5,
                3.0, 3.0,
                2.0, 2.5;
      // clang-format on

      Eigen::Matrix<double, Dim, Dim> jacobian;
      // clang-format off
      jacobian <<  2.0, 1.0,
                   0.0, 3.0;
      // clang-format on

      // Get B-Matrix
      auto bmatrix = tri->bmatrix(xi, coords, Eigen::Vector2d::Zero(),
                                  Eigen::Vector2d::Zero());

      // Check gradient of shape functions
      auto gradsf = tri->grad_shapefn(xi, Eigen::Vector2d::Zero(),
                                      Eigen::Vector2d::Zero());
      gradsf *= ((jacobian.inverse()).transpose());

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

    // Coordinates is (0.2,0.4)
    SECTION("Six noded triangle B-matrix cell for coordinates(0.2,0.4)") {
      // Reference coordinates
      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0.2, 0.4;

      // Nodal coordinates
      Eigen::Matrix<double, 6, Dim> coords;
      // clang-format off
      coords << 2.0, 1.0,
                4.0, 2.0,
                2.0, 4.0,
                3.0, 1.5,
                3.0, 3.0,
                2.0, 2.5;
      // clang-format on

      Eigen::Matrix<double, Dim, Dim> jacobian;
      // clang-format off
      jacobian <<  2.0, 1.0,
                   0.0, 3.0;
      // clang-format on

      // Get B-Matrix
      auto bmatrix = tri->bmatrix(xi, coords, Eigen::Vector2d::Zero(),
                                  Eigen::Vector2d::Zero());

      // Check gradient of shape functions
      auto gradsf = tri->grad_shapefn(xi, Eigen::Vector2d::Zero(),
                                      Eigen::Vector2d::Zero());
      gradsf *= ((jacobian.inverse()).transpose());

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

    // Check Bmatrix with deformation gradient
    SECTION("Six noded triangle B-matrix with deformation gradient") {
      // Reference coordinates
      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0.5, 0.5;

      // Nodal coordinates
      Eigen::Matrix<double, 6, Dim> coords;
      // clang-format off
      coords << 2.0, 1.0,
                4.0, 2.0,
                2.0, 4.0,
                3.0, 1.5,
                3.0, 3.0,
                2.0, 2.5;
      // clang-format on

      Eigen::Matrix<double, Dim, Dim> jacobian;
      // clang-format off
      jacobian <<  2.0, 1.0,
                   0.0, 3.0;
      // clang-format on

      Eigen::Matrix<double, Dim, 1> psize;
      psize << 0.5, 0.5;
      Eigen::Matrix<double, Dim, 1> defgrad;
      defgrad << 0.5, 0.5;

      // Get B-Matrix
      auto bmatrix = tri->bmatrix(xi, coords, psize, defgrad);

      // Check gradient of shape functions
      auto gradsf = tri->grad_shapefn(xi, psize, defgrad);
      gradsf *= ((jacobian.inverse()).transpose());

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

    SECTION("Six noded triangle B-matrix and Jacobian failure") {
      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0., 0.;

      Eigen::Matrix<double, 6, Dim> coords;
      // clang-format off
      coords << 2.0, 1.0,
                4.0, 2.0,
                2.0, 4.0,
                3.0, 1.5,
                3.0, 3.0,
                2.0, 2.5;
      // clang-format on
      // Get B-Matrix
      auto bmatrix = tri->bmatrix(xi, coords, Eigen::Vector2d::Zero(),
                                  Eigen::Vector2d::Zero());
      auto jacobian = tri->jacobian(xi, coords, Eigen::Vector2d::Zero(),
                                    Eigen::Vector2d::Zero());
    }

    // Ni Nj matrix of a cell
    SECTION("Six noded triangle ni-nj-matrix") {
      std::vector<Eigen::Matrix<double, Dim, 1>> xi_s;

      Eigen::Matrix<double, Dim, 1> xi;
      xi << 1.0 / 6, 1.0 / 6;
      xi_s.emplace_back(xi);
      xi << 2.0 / 3, 1.0 / 6;
      xi_s.emplace_back(xi);
      xi << 1.0 / 6, 2.0 / 3;
      xi_s.emplace_back(xi);

      REQUIRE(xi_s.size() == 3);

      // Get Ni Nj matrix
      const auto ni_nj_matrix = tri->ni_nj_matrix(xi_s);

      // Check size of ni_nj_matrix
      REQUIRE(ni_nj_matrix.rows() == nfunctions);
      REQUIRE(ni_nj_matrix.cols() == nfunctions);

      // Sum should be equal to 1. * xi_s.size()
      REQUIRE(ni_nj_matrix.sum() ==
              Approx(1. * xi_s.size()).epsilon(Tolerance));

      Eigen::Matrix<double, 6, 6> mass;
      // clang-format off
      mass <<  0.07407407407407, -0.03703703703704, -0.03703703703704,
               0.03703703703704, -0.07407407407407,  0.03703703703704,
              -0.03703703703704,  0.07407407407407, -0.03703703703704,
	       0.03703703703704,  0.03703703703704, -0.07407407407407,
              -0.03703703703704, -0.03703703703704,  0.07407407407407,
              -0.07407407407407,  0.03703703703704,  0.03703703703704,
               0.03703703703704,  0.03703703703704, -0.07407407407407,
               0.40740740740741,  0.29629629629630,  0.29629629629630,
              -0.07407407407407,  0.03703703703704,  0.03703703703704,
               0.29629629629630,  0.40740740740741,  0.29629629629630,
               0.03703703703704, -0.07407407407407,  0.03703703703704,
               0.29629629629630,  0.29629629629630,  0.40740740740741;

      // clang-format on

      // auxiliary matrices for checking its multiplication by scalar
      auto ni_nj_matrix_unit = 1.0 * ni_nj_matrix;
      auto ni_nj_matrix_zero = 0.0 * ni_nj_matrix;
      auto ni_nj_matrix_negative = -2.0 * ni_nj_matrix;
      double scalar = 21.65489;
      auto ni_nj_matrix_scalar = scalar * ni_nj_matrix;

      for (unsigned i = 0; i < nfunctions; ++i) {
        for (unsigned j = 0; j < nfunctions; ++j) {
          REQUIRE(ni_nj_matrix(i, j) == Approx(mass(i, j)).epsilon(Tolerance));
          // check multiplication by unity;
          REQUIRE(ni_nj_matrix_unit(i, j) ==
                  Approx(1.0 * mass(i, j)).epsilon(Tolerance));
          // check multiplication by zero;
          REQUIRE(ni_nj_matrix_zero(i, j) ==
                  Approx(0.0 * mass(i, j)).epsilon(Tolerance));
          // check multiplication by negative number;
          REQUIRE(ni_nj_matrix_negative(i, j) ==
                  Approx(-2.0 * mass(i, j)).epsilon(Tolerance));
          // check multiplication by an arbitrary scalar;
          REQUIRE(ni_nj_matrix_scalar(i, j) ==
                  Approx(scalar * mass(i, j)).epsilon(Tolerance));
        }
      }
    }

    // Laplace matrix of a cell
    SECTION("Six noded triangle laplace-matrix") {
      std::vector<Eigen::Matrix<double, Dim, 1>> xi_s;

      Eigen::Matrix<double, Dim, 1> xi;
      const double one_by_sqrt3 = std::fabs(1 / std::sqrt(3));
      xi << 1.0 / 6, 1.0 / 6;
      xi_s.emplace_back(xi);
      xi << 2.0 / 3, 1.0 / 6;
      xi_s.emplace_back(xi);
      xi << 1.0 / 6, 2.0 / 3;
      xi_s.emplace_back(xi);

      REQUIRE(xi_s.size() == 3);

      Eigen::Matrix<double, 6, Dim> coords;
      // clang-format off
      coords << 2.0, 1.0,
                4.0, 2.0,
                2.0, 4.0,
                3.0, 1.5,
                3.0, 3.0,
                2.0, 2.5;
      // clang-format on

      // Get laplace matrix
      const auto laplace_matrix = tri->laplace_matrix(xi_s, coords);

      // Check size of laplace-matrix
      REQUIRE(laplace_matrix.rows() == nfunctions);
      REQUIRE(laplace_matrix.cols() == nfunctions);

      // Sum should be equal to 0.
      REQUIRE(laplace_matrix.sum() == Approx(0.).epsilon(Tolerance));

      Eigen::Matrix<double, 6, 6> laplace;
      // clang-format off
      laplace <<  0.83333333333333,  0.22222222222222,  0.05555555555556,
                 -0.88888888888889,  0.00000000000000, -0.22222222222222,
                  0.22222222222222,  0.83333333333333,  0.05555555555556,
                 -0.88888888888889, -0.22222222222222,  0.00000000000000,
                  0.05555555555556,  0.05555555555556,  0.33333333333333,
                  0.00000000000000, -0.22222222222222, -0.22222222222222,
                 -0.88888888888889, -0.88888888888889,  0.00000000000000,
                  2.66666666666667, -0.44444444444444, -0.44444444444444,
                  0.00000000000000, -0.22222222222222, -0.22222222222222,
                 -0.44444444444444,  2.66666666666667, -1.77777777777778,
                 -0.22222222222222,  0.00000000000000, -0.22222222222222,
                 -0.44444444444444, -1.77777777777778,  2.66666666666667;
      // clang-format on
      for (unsigned i = 0; i < nfunctions; ++i)
        for (unsigned j = 0; j < nfunctions; ++j)
          REQUIRE(laplace_matrix(i, j) ==
                  Approx(laplace(i, j)).epsilon(Tolerance));
    }

    SECTION("Six noded triangle coordinates of unit cell") {
      const unsigned nfunctions = 6;

      // Coordinates of a unit cell
      Eigen::Matrix<double, nfunctions, Dim> unit_cell;
      // clang-format off
      unit_cell <<  0.0,  0.0,
                    1.0,  0.0,
                    0.0,  1.0,
                    0.5,  0.0,
                    0.5,  0.5,
                    0.0,  0.5;
      // clang-format on

      auto coordinates = tri->unit_cell_coordinates();
      REQUIRE(coordinates.rows() == nfunctions);
      REQUIRE(coordinates.cols() == Dim);
      for (unsigned i = 0; i < nfunctions; ++i) {  // Iterate through nfunctions
        for (unsigned j = 0; j < Dim; ++j) {       // Dimension
          REQUIRE(coordinates(i, j) ==
                  Approx(unit_cell(i, j)).epsilon(Tolerance));
        }
      }
    }

    SECTION("Six noded triangle element for side indices") {
      // Check for sides indices
      Eigen::MatrixXi indices = tri->sides_indices();
      REQUIRE(indices.rows() == 3);
      REQUIRE(indices.cols() == 2);

      REQUIRE(indices(0, 0) == 0);
      REQUIRE(indices(0, 1) == 1);

      REQUIRE(indices(1, 0) == 1);
      REQUIRE(indices(1, 1) == 2);

      REQUIRE(indices(2, 0) == 2);
      REQUIRE(indices(2, 1) == 0);
    }

    SECTION("Six noded triangle element for corner indices") {
      // Check for corner indices
      Eigen::VectorXi indices = tri->corner_indices();
      REQUIRE(indices.size() == 3);
      REQUIRE(indices(0) == 0);
      REQUIRE(indices(1) == 1);
      REQUIRE(indices(2) == 2);
    }

    SECTION("Six noded triangle element for inhedron indices") {
      // Check for inhedron indices
      Eigen::MatrixXi indices = tri->inhedron_indices();
      REQUIRE(indices.rows() == 3);
      REQUIRE(indices.cols() == 2);

      REQUIRE(indices(0, 0) == 0);
      REQUIRE(indices(0, 1) == 1);

      REQUIRE(indices(1, 0) == 1);
      REQUIRE(indices(1, 1) == 2);

      REQUIRE(indices(2, 0) == 2);
      REQUIRE(indices(2, 1) == 0);
    }

    SECTION("Six noded triangle shape function for face indices") {
      // Check for face indices
      Eigen::Matrix<int, 3, 3> indices;
      // clang-format off
      indices << 0, 1, 3,
                 1, 2, 4,
                 2, 0, 5;
      // clang-format on
      // Check for all face indices
      for (unsigned i = 0; i < indices.rows(); ++i) {
        const auto check_indices = tri->face_indices(i);
        REQUIRE(check_indices.rows() == 3);
        REQUIRE(check_indices.cols() == 1);

        for (unsigned j = 0; j < indices.cols(); ++j)
          REQUIRE(check_indices(j) == indices(i, j));
      }

      // Check number of faces
      REQUIRE(tri->nfaces() == 3);
    }
  }
}

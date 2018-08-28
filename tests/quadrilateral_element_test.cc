// Quadrilateral element test
#include <memory>

#include "catch.hpp"

#include "quadrilateral_element.h"

//! \brief Check quadrilateral element class
TEST_CASE("Quadrilateral elements are checked", "[quad][element][2D]") {
  const unsigned Dim = 2;
  const double Tolerance = 1.E-7;

  //! Check for 4 noded element
  SECTION("Quadrilateral element with four nodes") {
    const unsigned nfunctions = 4;
    std::shared_ptr<mpm::Element<Dim>> quad =
        std::make_shared<mpm::QuadrilateralElement<Dim, nfunctions>>();

    // Check degree
    REQUIRE(quad->degree() == mpm::ElementDegree::Linear);

    // Coordinates is (0,0)
    SECTION("Four noded quadrilateral element for coordinates(0,0)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();
      auto shapefn = quad->shapefn(coords);

      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      REQUIRE(shapefn(0) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(0.25).epsilon(Tolerance));
      REQUIRE(shapefn(3) == Approx(0.25).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(coords);
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
    SECTION("Four noded quadrilateral element for coordinates(-1, -1)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << -1., -1.;
      auto shapefn = quad->shapefn(coords);

      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      REQUIRE(shapefn(0) == Approx(1.0).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(3) == Approx(0.0).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(coords);
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
    SECTION("Four noded quadrilateral element for coordinates(1,1)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << 1., 1.;
      auto shapefn = quad->shapefn(coords);

      // Check shape function
      REQUIRE(shapefn.size() == nfunctions);

      REQUIRE(shapefn(0) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(1) == Approx(0.0).epsilon(Tolerance));
      REQUIRE(shapefn(2) == Approx(1.0).epsilon(Tolerance));
      REQUIRE(shapefn(3) == Approx(0.0).epsilon(Tolerance));

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(coords);
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

    // Check shapefn with deformation gradient
    SECTION(
        "Four noded quadrilateral element shapefn with deformation gradient") {
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

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(coords, psize, defgrad);
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

    // Check Jacobian
    SECTION(
        "Four noded quadrilateral Jacobian for local coordinates(0.5,0.5)") {
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
      auto jac = quad->jacobian(xi, coords);

      // Check size of jacobian
      REQUIRE(jac.size() == jacobian.size());

      // Check Jacobian
      for (unsigned i = 0; i < Dim; ++i)
        for (unsigned j = 0; j < Dim; ++j)
          REQUIRE(jac(i, j) == Approx(jacobian(i, j)).epsilon(Tolerance));
    }

    // Check Jacobian
    SECTION("Four noded quadrilateral Jacobian with deformation gradient") {
      Eigen::Matrix<double, 4, Dim> coords;
      // clang-format off
      coords << 2., 1.,
                4., 2.,
                2., 4.,
                1., 3.;
      // clang-format on

      Eigen::Matrix<double, Dim, 1> psize;
      psize.setZero();
      Eigen::Matrix<double, Dim, 1> defgrad;
      defgrad.setZero();

      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0.5, 0.5;

      Eigen::Matrix<double, Dim, Dim> jacobian;
      // clang-format off
      jacobian << 0.625, 0.5,
                 -0.875, 1.0;
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

    // Coordinates is (0,0)
    SECTION("Four noded quadrilateral B-matrix for coordinates(0,0)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();

      // Get B-Matrix
      auto bmatrix = quad->bmatrix(coords);

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(coords);

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
      auto bmatrix = quad->bmatrix(coords);

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(coords);

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
      auto bmatrix = quad->bmatrix(coords);

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(coords);

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

    // Coordinates is (0,0)
    SECTION("Four noded quadrilateral B-matrix cell for coordinates(0,0)") {
      // Reference coordinates
      Eigen::Matrix<double, Dim, 1> xi;
      xi.setZero();

      // Nodal coordinates
      Eigen::Matrix<double, 4, Dim> coords;
      // clang-format off
      coords << 0., 0.,
                1., 0.,
                1., 1.,
                0., 1.;
      // clang-format on

      // Get B-Matrix
      auto bmatrix = quad->bmatrix(xi, coords);

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(xi);
      gradsf *= 2.;

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
    SECTION("Four noded quadrilateral B-matrix cell for coordinates(0.5,0.5)") {
      // Reference coordinates
      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0.5, 0.5;

      // Nodal coordinates
      Eigen::Matrix<double, 4, Dim> coords;
      // clang-format off
      coords << 0., 0.,
                1., 0.,
                1., 1.,
                0., 1.;
      // clang-format on

      // Get B-Matrix
      auto bmatrix = quad->bmatrix(xi, coords);

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(xi);
      gradsf *= 2.;

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
    SECTION(
        "Four noded quadrilateral B-matrix cell for coordinates(-0.5,-0.5)") {
      // Reference coordinates
      Eigen::Matrix<double, Dim, 1> xi;
      xi << -0.5, -0.5;

      // Nodal coordinates
      Eigen::Matrix<double, 4, Dim> coords;
      // clang-format off
      coords << 0., 0.,
                1., 0.,
                1., 1.,
                0., 1.;
      // clang-format on

      // Get B-Matrix
      auto bmatrix = quad->bmatrix(xi, coords);

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(xi);
      gradsf *= 2.;

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
    SECTION(
        "Four noded quadrilateral B-matrix cell with deformation gradient") {
      // Reference coordinates
      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0.5, 0.5;

      // Nodal coordinates
      Eigen::Matrix<double, 4, Dim> coords;
      // clang-format off
      coords << 0., 0.,
                1., 0.,
                1., 1.,
                0., 1.;
      // clang-format on

      Eigen::Matrix<double, Dim, 1> psize;
      psize.setZero();
      Eigen::Matrix<double, Dim, 1> defgrad;
      defgrad.setZero();

      // Get B-Matrix
      auto bmatrix = quad->bmatrix(xi, coords, psize, defgrad);

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(xi, psize, defgrad);
      gradsf *= 2.;

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

    SECTION("Four noded quadrilateral B-matrix and Jacobian failure") {
      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0., 0.;

      Eigen::Matrix<double, 3, Dim> coords;
      // clang-format off
      coords << 0., 0.,
                1., 0., 
                1., 1.;
      // clang-format on
      // Get B-Matrix
      auto bmatrix = quad->bmatrix(xi, coords);
      auto jacobian = quad->jacobian(xi, coords);
    }

    // Mass matrix of a cell
    SECTION("Four noded quadrilateral mass-matrix") {
      std::vector<Eigen::Matrix<double, Dim, 1>> xi_s;

      Eigen::Matrix<double, Dim, 1> xi;
      const double one_by_sqrt3 = std::fabs(1 / std::sqrt(3));
      xi << -one_by_sqrt3, -one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << one_by_sqrt3, -one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << one_by_sqrt3, one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << -one_by_sqrt3, one_by_sqrt3;
      xi_s.emplace_back(xi);

      REQUIRE(xi_s.size() == 4);

      // Get mass matrix
      const auto mass_matrix = quad->mass_matrix(xi_s);

      // Check size of mass-matrix
      REQUIRE(mass_matrix.rows() == nfunctions);
      REQUIRE(mass_matrix.cols() == nfunctions);

      // Sum should be equal to 1. * xi_s.size()
      REQUIRE(mass_matrix.sum() == Approx(1. * xi_s.size()).epsilon(Tolerance));

      Eigen::Matrix<double, 4, 4> mass;
      // clang-format off
      mass <<  0.4444444444444445, 0.2222222222222222, 0.1111111111111111, 0.2222222222222222,
               0.2222222222222222, 0.4444444444444445, 0.2222222222222222, 0.1111111111111111,
               0.1111111111111111, 0.2222222222222222, 0.4444444444444445, 0.2222222222222222,
               0.2222222222222222, 0.1111111111111111, 0.2222222222222222, 0.4444444444444445;
      // clang-format on

      for (unsigned i = 0; i < nfunctions; ++i)
        for (unsigned j = 0; j < nfunctions; ++j)
          REQUIRE(mass_matrix(i, j) == Approx(mass(i, j)).epsilon(Tolerance));
    }

    // Laplace matrix of a cell
    SECTION("Four noded quadrilateral laplace-matrix") {
      std::vector<Eigen::Matrix<double, Dim, 1>> xi_s;

      Eigen::Matrix<double, Dim, 1> xi;
      const double one_by_sqrt3 = std::fabs(1 / std::sqrt(3));
      xi << -one_by_sqrt3, -one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << one_by_sqrt3, -one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << one_by_sqrt3, one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << -one_by_sqrt3, one_by_sqrt3;
      xi_s.emplace_back(xi);

      REQUIRE(xi_s.size() == 4);

      Eigen::Matrix<double, 4, Dim> coords;
      // clang-format off
      coords << 2., 1.,
                4., 2.,
                2., 4.,
                1., 3.;
      // clang-format on

      // Get laplace matrix
      const auto laplace_matrix = quad->laplace_matrix(xi_s, coords);

      // Check size of laplace-matrix
      REQUIRE(laplace_matrix.rows() == nfunctions);
      REQUIRE(laplace_matrix.cols() == nfunctions);

      // Sum should be equal to 0.
      REQUIRE(laplace_matrix.sum() == Approx(0.).epsilon(Tolerance));

      Eigen::Matrix<double, 4, 4> laplace;
      // clang-format off
      laplace <<  0.6572032405702867,-0.1005029077617645,-0.4834678486646353,-0.0732324841438869,
                 -0.1005029077617645, 0.4548818121031032,-0.1858117822979408,-0.1685671220433979,
                 -0.4834678486646353,-0.1858117822979408, 0.9689758463959470,-0.2996962154333709,
                 -0.0732324841438869,-0.1685671220433979,-0.2996962154333709, 0.5414958216206558;
      // clang-format on
      for (unsigned i = 0; i < nfunctions; ++i)
        for (unsigned j = 0; j < nfunctions; ++j)
          REQUIRE(laplace_matrix(i, j) ==
                  Approx(laplace(i, j)).epsilon(Tolerance));
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

      auto coordinates = quad->unit_cell_coordinates();
      REQUIRE(coordinates.rows() == nfunctions);
      REQUIRE(coordinates.cols() == Dim);
      for (unsigned i = 0; i < nfunctions; ++i) {  // Iterate through nfunctions
        for (unsigned j = 0; j < Dim; ++j) {       // Dimension
          REQUIRE(coordinates(i, j) ==
                  Approx(unit_cell(i, j)).epsilon(Tolerance));
        }
      }
    }

    SECTION("Four noded quadrilateral element for side indices") {
      // Check for sides indices
      Eigen::MatrixXi indices = quad->sides_indices();
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

    SECTION("Four noded quadrilateral element for corner indices") {
      // Check for corner indices
      Eigen::VectorXi indices = quad->corner_indices();
      REQUIRE(indices.size() == 4);
      REQUIRE(indices(0) == 0);
      REQUIRE(indices(1) == 1);
      REQUIRE(indices(2) == 2);
      REQUIRE(indices(3) == 3);
    }

    SECTION("Four noded quadrilateral element for inhedron indices") {
      // Check for inhedron indices
      Eigen::MatrixXi indices = quad->inhedron_indices();
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

    SECTION("Four noded quadrilateral shape function for face indices") {
      // Check for face indices
      Eigen::Matrix<int, 4, 2> indices;
      // clang-format off
      indices << 0, 1, 
                 1, 2, 
                 2, 3, 
                 3, 0;
      // clang-format on
      // Check for all face indices
      for (unsigned i = 0; i < indices.rows(); ++i) {
        const auto check_indices = quad->face_indices(i);
        REQUIRE(check_indices.rows() == 2);
        REQUIRE(check_indices.cols() == 1);

        for (unsigned j = 0; j < indices.cols(); ++j)
          REQUIRE(check_indices(j) == indices(i, j));
      }
    }
  }

  //! Check for 8 noded element
  SECTION("Quadrilateral element with eight nodes") {
    const unsigned nfunctions = 8;
    std::shared_ptr<mpm::Element<Dim>> quad =
        std::make_shared<mpm::QuadrilateralElement<Dim, nfunctions>>();

    // Check degree
    REQUIRE(quad->degree() == mpm::ElementDegree::Quadratic);

    // Coordinates is (0,0)
    SECTION("Eight noded quadrilateral element for coordinates(0,0)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();
      auto shapefn = quad->shapefn(coords);

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
      auto gradsf = quad->grad_shapefn(coords);
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
    SECTION("Eight noded quadrilateral element for coordinates(-1,-1)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << -1., -1.;
      auto shapefn = quad->shapefn(coords);

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
      auto gradsf = quad->grad_shapefn(coords);
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
    SECTION("Eight noded quadrilateral element for coordinates(1, 1)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << 1., 1.;
      auto shapefn = quad->shapefn(coords);

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
      auto gradsf = quad->grad_shapefn(coords);
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
    SECTION(
        "Eight noded quadrilateral element shapefn with deformation gradient") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();
      Eigen::Matrix<double, Dim, 1> psize;
      psize.setZero();
      Eigen::Matrix<double, Dim, 1> defgrad;
      defgrad.setZero();
      auto shapefn = quad->shapefn(coords, psize, defgrad);

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
      auto gradsf = quad->grad_shapefn(coords, psize, defgrad);
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

    // Check Jacobian
    SECTION(
        "Eight noded quadrilateral Jacobian for local coordinates(0.5,0.5)") {
      Eigen::Matrix<double, 8, Dim> coords;
      // clang-format off
      coords << 2.0, 1.0,
                4.0, 2.0,
                2.0, 4.0,
                1.0, 3.0,
                3.0, 1.5,
                3.0, 3.0,
                1.5, 3.5,
                1.5, 2.0;
      // clang-format on

      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0.5, 0.5;

      // Jacobian result
      Eigen::Matrix<double, Dim, Dim> jacobian;
      // clang-format off
      jacobian << 0.625, 0.5,
                 -0.875, 1.0;
      // clang-format on

      // Get Jacobian
      auto jac = quad->jacobian(xi, coords);

      // Check size of jacobian
      REQUIRE(jac.size() == jacobian.size());

      // Check Jacobian
      for (unsigned i = 0; i < Dim; ++i)
        for (unsigned j = 0; j < Dim; ++j)
          REQUIRE(jac(i, j) == Approx(jacobian(i, j)).epsilon(Tolerance));
    }

    // Check Jacobian
    SECTION("Eight noded quadrilateral Jacobian with deformation gradient") {
      Eigen::Matrix<double, 8, Dim> coords;
      // clang-format off
      coords << 2.0, 1.0,
                4.0, 2.0,
                2.0, 4.0,
                1.0, 3.0,
                3.0, 1.5,
                3.0, 3.0,
                1.5, 3.5,
                1.5, 2.0;
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
      jacobian << 0.625, 0.5,
                 -0.875, 1.0;
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

    // Coordinates is (0,0)
    SECTION("Eight noded quadrilateral B-matrix for coordinates(0,0)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();

      // Get B-Matrix
      auto bmatrix = quad->bmatrix(coords);

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(coords);

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
      auto bmatrix = quad->bmatrix(coords);

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(coords);

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
      auto bmatrix = quad->bmatrix(coords);

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(coords);

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

    // Coordinates is (0,0)
    SECTION("Eight noded quadrilateral B-matrix cell for coordinates(0,0)") {
      // Reference coordinates
      Eigen::Matrix<double, Dim, 1> xi;
      xi.setZero();

      // Nodal coordinates
      Eigen::Matrix<double, 8, Dim> coords;
      // clang-format off
      coords << 0., 0.,
                1., 0.,
                1., 1.,
                0., 1.,
                0.5, 0.,
                1.0, 0.5,
                0.5, 1.,
                0., 0.5;
      // clang-format on

      // Get B-Matrix
      auto bmatrix = quad->bmatrix(xi, coords);

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(xi);
      gradsf *= 2.;

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
    SECTION(
        "Eight noded quadrilateral B-matrix cell for coordinates(0.5,0.5)") {
      // Reference coordinates
      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0.5, 0.5;

      // Nodal coordinates
      Eigen::Matrix<double, 8, Dim> coords;
      // clang-format off
      coords << 0., 0.,
                1., 0.,
                1., 1.,
                0., 1.,
                0.5, 0.,
                1.0, 0.5,
                0.5, 1.,
                0., 0.5;
      // clang-format on

      // Get B-Matrix
      auto bmatrix = quad->bmatrix(xi, coords);

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(xi);
      gradsf *= 2.;

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
    SECTION(
        "Eight noded quadrilateral B-matrix cell for coordinates(-0.5,-0.5)") {
      // Reference coordinates
      Eigen::Matrix<double, Dim, 1> xi;
      xi << -0.5, -0.5;

      // Nodal coordinates
      Eigen::Matrix<double, 8, Dim> coords;
      // clang-format off
      coords << 0., 0.,
                1., 0.,
                1., 1.,
                0., 1.,
                0.5, 0.,
                1.0, 0.5,
                0.5, 1.,
                0., 0.5;
      // clang-format on

      // Get B-Matrix
      auto bmatrix = quad->bmatrix(xi, coords);

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(xi);
      gradsf *= 2.;

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
    SECTION("Eight noded quadrilateral B-matrix with deformation gradient") {
      // Reference coordinates
      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0.5, 0.5;

      // Nodal coordinates
      Eigen::Matrix<double, 8, Dim> coords;
      // clang-format off
      coords << 0., 0.,
                1., 0.,
                1., 1.,
                0., 1.,
                0.5, 0.,
                1.0, 0.5,
                0.5, 1.,
                0., 0.5;
      // clang-format on

      Eigen::Matrix<double, Dim, 1> psize;
      psize << 0.5, 0.5;
      Eigen::Matrix<double, Dim, 1> defgrad;
      defgrad << 0.5, 0.5;

      // Get B-Matrix
      auto bmatrix = quad->bmatrix(xi, coords, psize, defgrad);

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(xi, psize, defgrad);
      gradsf *= 2.;

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

    SECTION("Eight noded quadrilateral B-matrix and Jacobian failure") {
      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0., 0.;

      Eigen::Matrix<double, 3, Dim> coords;
      // clang-format off
      coords << 0., 0.,
                1., 0., 
                1., 1.;
      // clang-format on
      // Get B-Matrix
      auto bmatrix = quad->bmatrix(xi, coords);
      auto jacobian = quad->jacobian(xi, coords);
    }

    // Mass matrix of a cell
    SECTION("Eight noded quadrilateral mass-matrix") {
      std::vector<Eigen::Matrix<double, Dim, 1>> xi_s;

      Eigen::Matrix<double, Dim, 1> xi;
      const double one_by_sqrt3 = std::fabs(1 / std::sqrt(3));
      xi << -one_by_sqrt3, -one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << one_by_sqrt3, -one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << one_by_sqrt3, one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << -one_by_sqrt3, one_by_sqrt3;
      xi_s.emplace_back(xi);

      REQUIRE(xi_s.size() == 4);

      // Get mass matrix
      const auto mass_matrix = quad->mass_matrix(xi_s);

      // Check size of mass-matrix
      REQUIRE(mass_matrix.rows() == nfunctions);
      REQUIRE(mass_matrix.cols() == nfunctions);

      // Sum should be equal to 1. * xi_s.size()
      REQUIRE(mass_matrix.sum() == Approx(1. * xi_s.size()).epsilon(Tolerance));

      Eigen::Matrix<double, 8, 8> mass;
      mass << 0.07407407407407406, 0.00000000000000, 0.037037037037037,
          -0.00000000000000000, -0.07407407407407397, -0.1481481481481481,
          -0.1481481481481481, -0.07407407407407397, 0.0, 0.07407407407407407,
          0.0, 0.037037037037037, -0.07407407407407397, -0.07407407407407397,
          -0.1481481481481481, -0.1481481481481481, 0.037037037037037,
          0.000000000000000, 0.07407407407407407, 0.000000000000000,
          -0.1481481481481481, -0.07407407407407396, -0.07407407407407397,
          -0.1481481481481481, 0.0000000000000, 0.037037037037037,
          0.0000000000000, 0.07407407407407406, -0.1481481481481481,
          -0.1481481481481481, -0.07407407407407399, -0.07407407407407396,
          -0.07407407407407397, -0.07407407407407397, -0.1481481481481481,
          -0.1481481481481481, 0.5925925925925922, 0.4444444444444442,
          0.2962962962962961, 0.4444444444444442, -0.1481481481481481,
          -0.07407407407407397, -0.07407407407407396, -0.1481481481481481,
          0.4444444444444442, 0.5925925925925923, 0.4444444444444442,
          0.2962962962962961, -0.1481481481481481, -0.1481481481481481,
          -0.07407407407407397, -0.07407407407407399, 0.2962962962962961,
          0.4444444444444442, 0.5925925925925923, 0.4444444444444442,
          -0.07407407407407397, -0.1481481481481481, -0.1481481481481481,
          -0.07407407407407396, 0.4444444444444442, 0.2962962962962961,
          0.4444444444444442, 0.5925925925925923;

      for (unsigned i = 0; i < nfunctions; ++i)
        for (unsigned j = 0; j < nfunctions; ++j)
          REQUIRE(mass_matrix(i, j) == Approx(mass(i, j)).epsilon(Tolerance));
    }

    // Laplace matrix of a cell
    SECTION("Four noded quadrilateral laplace-matrix") {
      std::vector<Eigen::Matrix<double, Dim, 1>> xi_s;

      Eigen::Matrix<double, Dim, 1> xi;
      const double one_by_sqrt3 = std::fabs(1 / std::sqrt(3));
      xi << -one_by_sqrt3, -one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << one_by_sqrt3, -one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << one_by_sqrt3, one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << -one_by_sqrt3, one_by_sqrt3;
      xi_s.emplace_back(xi);

      REQUIRE(xi_s.size() == 4);

      Eigen::Matrix<double, 8, Dim> coords;
      // clang-format off
      coords << 2.0, 1.0,
                4.0, 2.0,
                2.0, 4.0,
                1.0, 3.0,
                3.0, 1.5,
                3.0, 3.0,
                1.5, 3.5,
                1.5, 2.0;
      // clang-format on

      // Get laplace matrix
      const auto laplace_matrix = quad->laplace_matrix(xi_s, coords);

      // Check size of laplace-matrix
      REQUIRE(laplace_matrix.rows() == nfunctions);
      REQUIRE(laplace_matrix.cols() == nfunctions);

      // Sum should be equal to 0.
      REQUIRE(laplace_matrix.sum() == Approx(0.).epsilon(Tolerance));

      Eigen::Matrix<double, 8, 8> laplace;
      laplace << 1.072493643107602, 0.4065987931213615, 0.6655397867711725,
          0.4263799884212531, -0.7648755792622943, -0.5379335339843152,
          -0.5366002581899634, -0.7316028399848165, 0.4065987931213615,
          0.6876728555050236, 0.5389592865753069, 0.4340001313936471,
          -0.5139612336052258, -0.6612353999512298, -0.5371839904579891,
          -0.3548504425808939, 0.6655397867711725, 0.5389592865753069,
          1.69222960654119, 0.6057921007929679, -0.7389668317060483,
          -0.9851551573410983, -1.134868473557333, -0.643530318076158,
          0.4263799884212531, 0.4340001313936471, 0.6057921007929679,
          0.9409579548861762, -0.4852376064304632, -0.57346791647202,
          -0.9386660952285415, -0.4097585573630192, -0.7648755792622943,
          -0.5139612336052258, -0.7389668317060483, -0.4852376064304632,
          1.972465581627769, -0.08358715573264081, 0.5368025726730385,
          0.07736025243586532, -0.5379335339843152, -0.6612353999512298,
          -0.9851551573410983, -0.57346791647202, -0.08358715573264081,
          1.964331090455655, 0.5284519760236881, 0.3485960970019614,
          -0.5366002581899634, -0.5371839904579891, -1.134868473557333,
          -0.9386660952285415, 0.5368025726730385, 0.5284519760236881,
          2.56584444050972, -0.4837801717726197, -0.7316028399848165,
          -0.3548504425808939, -0.643530318076158, -0.4097585573630192,
          0.07736025243586532, 0.3485960970019614, -0.4837801717726197,
          2.197565980339681;
      for (unsigned i = 0; i < nfunctions; ++i)
        for (unsigned j = 0; j < nfunctions; ++j)
          REQUIRE(laplace_matrix(i, j) ==
                  Approx(laplace(i, j)).epsilon(Tolerance));
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

      auto coordinates = quad->unit_cell_coordinates();
      REQUIRE(coordinates.rows() == nfunctions);
      REQUIRE(coordinates.cols() == Dim);
      for (unsigned i = 0; i < nfunctions; ++i) {  // Iterate through nfunctions
        for (unsigned j = 0; j < Dim; ++j) {       // Dimension
          REQUIRE(coordinates(i, j) ==
                  Approx(unit_cell(i, j)).epsilon(Tolerance));
        }
      }
    }

    SECTION("Eight noded quadrilateral element for side indices") {
      // Check for sides indices
      Eigen::MatrixXi indices = quad->sides_indices();
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

    SECTION("Eight noded quadrilateral element for corner indices") {
      // Check for corner indices
      Eigen::VectorXi indices = quad->corner_indices();
      REQUIRE(indices.size() == 4);
      REQUIRE(indices(0) == 0);
      REQUIRE(indices(1) == 1);
      REQUIRE(indices(2) == 2);
      REQUIRE(indices(3) == 3);
    }

    SECTION("Eight noded quadrilateral element for inhedron indices") {
      // Check for inhedron indices
      Eigen::MatrixXi indices = quad->inhedron_indices();
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

    SECTION("Eight noded quadrilateral shape function for face indices") {
      // Check for face indices
      Eigen::Matrix<int, 4, 3> indices;
      // clang-format off
      indices << 0, 1, 4, 
                 1, 2, 5, 
                 2, 3, 6, 
                 3, 0, 7;
      // clang-format on
      // Check for all face indices
      for (unsigned i = 0; i < indices.rows(); ++i) {
        const auto check_indices = quad->face_indices(i);
        REQUIRE(check_indices.rows() == 3);
        REQUIRE(check_indices.cols() == 1);

        for (unsigned j = 0; j < indices.cols(); ++j)
          REQUIRE(check_indices(j) == indices(i, j));
      }
    }
  }

  //! Check for 9 noded element
  SECTION("Quadrilateral element with nine nodes") {
    const unsigned nfunctions = 9;
    std::shared_ptr<mpm::Element<Dim>> quad =
        std::make_shared<mpm::QuadrilateralElement<Dim, nfunctions>>();

    // Check degree
    REQUIRE(quad->degree() == mpm::ElementDegree::Quadratic);

    // Coordinates is (0,0)
    SECTION("Nine noded quadrilateral element for coordinates(0,0)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();
      auto shapefn = quad->shapefn(coords);

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
      auto gradsf = quad->grad_shapefn(coords);
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
    SECTION("Nine noded quadrilateral element for coordinates(-1,-1)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << -1., -1.;
      auto shapefn = quad->shapefn(coords);

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
      auto gradsf = quad->grad_shapefn(coords);
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
    SECTION("Nine noded quadrilateral element for coordinates(1, 1)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << 1., 1.;
      auto shapefn = quad->shapefn(coords);

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
      auto gradsf = quad->grad_shapefn(coords);
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

    SECTION("Nine noded quadrilateral element with deformation gradient") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();
      Eigen::Matrix<double, Dim, 1> psize;
      psize << 0.5, 0.1;
      Eigen::Matrix<double, Dim, 1> defgrad;
      defgrad << -0.25, 0.1;

      auto shapefn = quad->shapefn(coords, psize, defgrad);

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
      auto gradsf = quad->grad_shapefn(coords, psize, defgrad);
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

    // Coordinates is (0,0)
    SECTION("Nine noded quadrilateral B-matrix for coordinates(0,0)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords.setZero();

      // Get B-Matrix
      auto bmatrix = quad->bmatrix(coords);

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(coords);

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
      auto bmatrix = quad->bmatrix(coords);

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(coords);

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

    // Mass matrix of a cell
    SECTION("Nine noded quadrilateral mass-matrix") {
      std::vector<Eigen::Matrix<double, Dim, 1>> xi_s;

      Eigen::Matrix<double, Dim, 1> xi;
      const double one_by_sqrt3 = std::fabs(1 / std::sqrt(3));
      xi << -one_by_sqrt3, -one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << one_by_sqrt3, -one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << one_by_sqrt3, one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << -one_by_sqrt3, one_by_sqrt3;
      xi_s.emplace_back(xi);

      REQUIRE(xi_s.size() == 4);

      // Get mass matrix
      const auto mass_matrix = quad->mass_matrix(xi_s);

      // Check size of mass-matrix
      REQUIRE(mass_matrix.rows() == nfunctions);
      REQUIRE(mass_matrix.cols() == nfunctions);

      // Sum should be equal to 1. * xi_s.size()
      REQUIRE(mass_matrix.sum() == Approx(1. * xi_s.size()).epsilon(Tolerance));

      Eigen::Matrix<double, 9, 9> mass;
      mass << 0.04938271604938273, -0.02469135802469136, 0.01234567901234568,
          -0.02469135802469136, 0.04938271604938271, -0.02469135802469136,
          -0.02469135802469136, 0.04938271604938272, 0.04938271604938271,
          -0.02469135802469136, 0.04938271604938273, -0.02469135802469136,
          0.01234567901234568, 0.04938271604938271, 0.04938271604938272,
          -0.02469135802469136, -0.02469135802469136, 0.04938271604938271,
          0.01234567901234568, -0.02469135802469136, 0.04938271604938273,
          -0.02469135802469136, -0.02469135802469136, 0.04938271604938271,
          0.04938271604938272, -0.02469135802469135, 0.0493827160493827,
          -0.02469135802469136, 0.01234567901234568, -0.02469135802469136,
          0.04938271604938273, -0.02469135802469135, -0.02469135802469135,
          0.04938271604938271, 0.04938271604938271, 0.0493827160493827,
          0.04938271604938271, 0.04938271604938271, -0.02469135802469136,
          -0.02469135802469135, 0.1975308641975308, 0.04938271604938271,
          -0.09876543209876538, 0.04938271604938271, 0.1975308641975307,
          -0.02469135802469136, 0.04938271604938272, 0.04938271604938271,
          -0.02469135802469135, 0.04938271604938271, 0.1975308641975308,
          0.0493827160493827, -0.09876543209876538, 0.1975308641975307,
          -0.02469135802469136, -0.02469135802469136, 0.04938271604938272,
          0.04938271604938271, -0.09876543209876538, 0.0493827160493827,
          0.1975308641975308, 0.0493827160493827, 0.1975308641975307,
          0.04938271604938272, -0.02469135802469136, -0.02469135802469135,
          0.04938271604938271, 0.04938271604938271, -0.09876543209876538,
          0.0493827160493827, 0.1975308641975308, 0.1975308641975307,
          0.04938271604938271, 0.04938271604938271, 0.0493827160493827,
          0.0493827160493827, 0.1975308641975307, 0.1975308641975307,
          0.1975308641975307, 0.1975308641975307, 0.7901234567901227;

      for (unsigned i = 0; i < nfunctions; ++i)
        for (unsigned j = 0; j < nfunctions; ++j)
          REQUIRE(mass_matrix(i, j) == Approx(mass(i, j)).epsilon(Tolerance));
    }

    // Laplace matrix of a cell
    SECTION("Nine noded quadrilateral laplace-matrix") {
      std::vector<Eigen::Matrix<double, Dim, 1>> xi_s;

      Eigen::Matrix<double, Dim, 1> xi;
      const double one_by_sqrt3 = std::fabs(1 / std::sqrt(3));
      xi << -one_by_sqrt3, -one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << one_by_sqrt3, -one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << one_by_sqrt3, one_by_sqrt3;
      xi_s.emplace_back(xi);
      xi << -one_by_sqrt3, one_by_sqrt3;
      xi_s.emplace_back(xi);

      REQUIRE(xi_s.size() == 4);

      Eigen::Matrix<double, 9, Dim> coords;
      // clang-format off
      coords << 2.0, 1.0,
                4.0, 2.0,
                2.0, 4.0,
                1.0, 3.0,
                3.0, 1.5,
                3.0, 3.0,
                1.5, 3.5,
                1.5, 2.0,
                2.25, 2.5;
      // clang-format on

      // Get laplace matrix
      const auto laplace_matrix = quad->laplace_matrix(xi_s, coords);

      // Check size of laplace-matrix
      REQUIRE(laplace_matrix.rows() == nfunctions);
      REQUIRE(laplace_matrix.cols() == nfunctions);

      // Sum should be equal to 0.
      REQUIRE(laplace_matrix.sum() == Approx(0.).epsilon(Tolerance));

      Eigen::Matrix<double, 9, 9> laplace;
      laplace << 0.5084866346946783, -0.0734446911538873, -0.05371864985162622,
          -0.1103133473371132, -0.07669209107623806, 0.1927084136591803,
          0.2589628244010591, -0.10396920063195, -0.5420198927041027,
          -0.0734446911538873, 0.2915928953674497, -0.09633562590981673,
          -0.01872968022704423, 0.006295206305480323, -0.09852050058308423,
          0.09045204385768324, 0.104856148496622, -0.2061657961534028,
          -0.05371864985162622, -0.09633562590981673, 0.8177197417085185,
          -0.08615266317527349, 0.2597195128997584, 0.05598964672214593,
          -0.02880253454656023, 0.2946061776964577, -1.163025605543604,
          -0.1103133473371132, -0.01872968022704423, -0.08615266317527349,
          0.4315782917823672, 0.1483185364464784, 0.1025466858623602,
          -0.1977303579466335, 0.163247736680731, -0.4327652020858722,
          -0.07669209107623806, 0.006295206305480323, 0.2597195128997584,
          0.1483185364464784, 1.47575966253524, -0.6652099937400512,
          -0.1746625352294245, -0.298245968990285, -0.675282329150958,
          0.1927084136591803, -0.09852050058308423, 0.05598964672214593,
          0.1025466858623602, -0.6652099937400512, 1.297791333533371,
          -0.2679300507936517, -0.1119270433390657, -0.5054484913212045,
          0.2589628244010591, 0.09045204385768324, -0.02880253454656023,
          -0.1977303579466335, -0.1746625352294245, -0.2679300507936517,
          1.639620143797324, -1.074145582008699, -0.2457639515310967,
          -0.10396920063195, 0.104856148496622, 0.2946061776964577,
          0.163247736680731, -0.298245968990285, -0.1119270433390657,
          -1.074145582008699, 1.943059456579911, -0.917481724483721,
          -0.5420198927041027, -0.2061657961534028, -1.163025605543604,
          -0.4327652020858722, -0.675282329150958, -0.5054484913212045,
          -0.2457639515310967, -0.917481724483721, 4.687952992973963;

      for (unsigned i = 0; i < nfunctions; ++i)
        for (unsigned j = 0; j < nfunctions; ++j)
          REQUIRE(laplace_matrix(i, j) ==
                  Approx(laplace(i, j)).epsilon(Tolerance));
    }

    // Check Jacobian
    SECTION(
        "Nine noded quadrilateral Jacobian for local coordinates(0.5,0.5)") {
      Eigen::Matrix<double, 9, Dim> coords;
      // clang-format off
      coords << 2.0, 1.0,
                4.0, 2.0,
                2.0, 4.0,
                1.0, 3.0,
                3.0, 1.5,
                3.0, 3.0,
                1.5, 3.5,
                1.5, 2.0,
                2.25, 2.5;
      // clang-format on

      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0.5, 0.5;

      // Jacobian result
      Eigen::Matrix<double, Dim, Dim> jacobian;
      // clang-format off
      jacobian << 0.625, 0.5,
                 -0.875, 1.0;
      // clang-format on

      // Get Jacobian
      auto jac = quad->jacobian(xi, coords);

      // Check size of jacobian
      REQUIRE(jac.size() == jacobian.size());

      // Check Jacobian
      for (unsigned i = 0; i < Dim; ++i)
        for (unsigned j = 0; j < Dim; ++j)
          REQUIRE(jac(i, j) == Approx(jacobian(i, j)).epsilon(Tolerance));
    }

    // Check Jacobian with deformation gradient
    SECTION("Nine noded quadrilateral Jacobian with deformation gradient") {
      Eigen::Matrix<double, 9, Dim> coords;
      // clang-format off
      coords << 2.0, 1.0,
                4.0, 2.0,
                2.0, 4.0,
                1.0, 3.0,
                3.0, 1.5,
                3.0, 3.0,
                1.5, 3.5,
                1.5, 2.0,
                2.25, 2.5;
      // clang-format on

      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0.5, 0.5;

      Eigen::Matrix<double, Dim, 1> psize;
      psize << -0.5, -0.5;
      Eigen::Matrix<double, Dim, 1> defgrad;
      defgrad << -0.5, -0.5;

      // Jacobian result
      Eigen::Matrix<double, Dim, Dim> jacobian;
      // clang-format off
      jacobian << 0.625, 0.5,
                 -0.875, 1.0;
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

    // Coordinates is (-0.5,-0.5)
    SECTION("Nine noded quadrilateral B-matrix for coordinates(-0.5,-0.5)") {
      Eigen::Matrix<double, Dim, 1> coords;
      coords << -0.5, -0.5;

      // Get B-Matrix
      auto bmatrix = quad->bmatrix(coords);

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(coords);

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

    // Coordinates is (0,0)
    SECTION("Nine noded quadrilateral B-matrix cell for coordinates(0,0)") {
      // Reference coordinates
      Eigen::Matrix<double, Dim, 1> xi;
      xi.setZero();

      // Nodal coordinates
      Eigen::Matrix<double, 9, Dim> coords;
      // clang-format off
      coords << 0., 0.,
                1., 0.,
                1., 1.,
                0., 1.,
                0.5, 0.,
                1.0, 0.5,
                0.5, 1.,
                0., 0.5,
                0.5, 0.5;
      // clang-format on

      // Get B-Matrix
      auto bmatrix = quad->bmatrix(xi, coords);

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(xi);
      gradsf *= 2.;

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
    SECTION("Nine noded quadrilateral B-matrix cell for coordinates(0.5,0.5)") {
      // Reference coordinates
      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0.5, 0.5;

      // Nodal coordinates
      Eigen::Matrix<double, 9, Dim> coords;
      // clang-format off
      coords << 0., 0.,
                1., 0.,
                1., 1.,
                0., 1.,
                0.5, 0.,
                1.0, 0.5,
                0.5, 1.,
                0., 0.5,
                0.5, 0.5;
      // clang-format on

      // Get B-Matrix
      auto bmatrix = quad->bmatrix(xi, coords);

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(xi);
      gradsf *= 2.;

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
    SECTION(
        "Nine noded quadrilateral B-matrix cell for coordinates(-0.5,-0.5)") {
      // Reference coordinates
      Eigen::Matrix<double, Dim, 1> xi;
      xi << -0.5, -0.5;

      // Nodal coordinates
      Eigen::Matrix<double, 9, Dim> coords;
      // clang-format off
      coords << 0., 0.,
                1., 0.,
                1., 1.,
                0., 1.,
                0.5, 0.,
                1.0, 0.5,
                0.5, 1.,
                0., 0.5,
                0.5, 0.5;
      // clang-format on

      // Get B-Matrix
      auto bmatrix = quad->bmatrix(xi, coords);

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(xi);
      gradsf *= 2.;

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

    SECTION("Nine noded quadrilateral B-matrix with deformation gradient)") {
      // Reference coordinates
      Eigen::Matrix<double, Dim, 1> xi;
      xi << -0.5, -0.5;

      Eigen::Matrix<double, Dim, 1> psize;
      psize << -0.5, -0.5;
      Eigen::Matrix<double, Dim, 1> defgrad;
      defgrad << -0.5, -0.5;

      // Nodal coordinates
      Eigen::Matrix<double, 9, Dim> coords;
      // clang-format off
      coords << 0., 0.,
                1., 0.,
                1., 1.,
                0., 1.,
                0.5, 0.,
                1.0, 0.5,
                0.5, 1.,
                0., 0.5,
                0.5, 0.5;
      // clang-format on

      // Get B-Matrix
      auto bmatrix = quad->bmatrix(xi, coords, psize, defgrad);

      // Check gradient of shape functions
      auto gradsf = quad->grad_shapefn(xi, psize, defgrad);
      gradsf *= 2.;

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

    SECTION("Nine noded quadrilateral B-matrix and Jacobian failure") {
      Eigen::Matrix<double, Dim, 1> xi;
      xi << 0., 0.;

      Eigen::Matrix<double, 3, Dim> coords;
      // clang-format off
      coords << 0., 0.,
                1., 0., 
                1., 1.;
      // clang-format on
      // Get B-Matrix
      auto bmatrix = quad->bmatrix(xi, coords);
      auto jacobian = quad->jacobian(xi, coords);
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

      auto coordinates = quad->unit_cell_coordinates();
      REQUIRE(coordinates.rows() == nfunctions);
      REQUIRE(coordinates.cols() == Dim);
      for (unsigned i = 0; i < nfunctions; ++i) {  // Iterate through nfunctions
        for (unsigned j = 0; j < Dim; ++j) {       // Dimension
          REQUIRE(coordinates(i, j) ==
                  Approx(unit_cell(i, j)).epsilon(Tolerance));
        }
      }
    }

    SECTION("Nine noded quadrilateral element for side indices") {
      // Check for sides indices
      Eigen::MatrixXi indices = quad->sides_indices();
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

    SECTION("Nine noded quadrilateral element for corner indices") {
      // Check for corner indices
      Eigen::VectorXi indices = quad->corner_indices();
      REQUIRE(indices.size() == 4);
      REQUIRE(indices(0) == 0);
      REQUIRE(indices(1) == 1);
      REQUIRE(indices(2) == 2);
      REQUIRE(indices(3) == 3);
    }

    SECTION("Nine noded quadrilateral element for inhedron indices") {
      // Check for inhedron indices
      Eigen::MatrixXi indices = quad->inhedron_indices();
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

    SECTION("Nine noded quadrilateral shape function for face indices") {
      // Check for face indices
      Eigen::Matrix<int, 4, 3> indices;
      // clang-format off
      indices << 0, 1, 4, 
                 1, 2, 5, 
                 2, 3, 6, 
                 3, 0, 7;
      // clang-format on

      // Check for all face indices
      for (unsigned i = 0; i < indices.rows(); ++i) {
        const auto check_indices = quad->face_indices(i);
        REQUIRE(check_indices.rows() == 3);
        REQUIRE(check_indices.cols() == 1);

        for (unsigned j = 0; j < indices.cols(); ++j)
          REQUIRE(check_indices(j) == indices(i, j));
      }
    }
  }
}

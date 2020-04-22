#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"
#include <Eigen/Sparse>

#include "cg_eigen.h"
#include "solver_base.h"

// Generate 3x3 test matrix
Eigen::SparseMatrix<double> CreateSymmetricTestMatrix3x3() {
  Eigen::SparseMatrix<double> matrix(3, 3);
  matrix.setZero();
  matrix.coeffRef(0, 0) = 4.0;
  matrix.coeffRef(1, 1) = 4.0;
  matrix.coeffRef(2, 2) = 4.0;
  matrix.coeffRef(0, 1) = 1.0;
  matrix.coeffRef(0, 2) = 1.0;
  matrix.coeffRef(1, 2) = 1.0;
  matrix.coeffRef(1, 0) = 1.0;
  matrix.coeffRef(2, 0) = 1.0;
  matrix.coeffRef(2, 1) = 1.0;
  return matrix;
}

TEST_CASE("Linear solver conjugate gradient", "[linear_solver]") {

  // Maximum iteration possible
  unsigned max_iter = 100;
  // Allowed solving tolerance
  double solve_tolerance = 10e-5;
  // Tolerance
  const double Tolerance = 1.E-7;

  SECTION("Eigen solver") {

    // Construct Eigen solver
    std::shared_ptr<mpm::SolverBase<Eigen::SparseMatrix<double>>>
        eigen_matrix_solver =
            Factory<mpm::SolverBase<Eigen::SparseMatrix<double>>, unsigned,
                    double>::instance()
                ->create("EigenCG", std::move(max_iter),
                         std::move(solve_tolerance));

    SECTION("Eigen 3x3 solver") {
      // Initiate solver_type
      std::string solver_type = "cg";

      // Construct matrix and vector
      const auto& A = CreateSymmetricTestMatrix3x3();
      Eigen::VectorXd b(3);
      b << 6.0, 3.0, -9.0;

      // Solve
      const auto& x_eigen = eigen_matrix_solver->solve(A, b, "cg");

      // Check
      REQUIRE(x_eigen(0) == Approx(2.).epsilon(Tolerance));
      REQUIRE(x_eigen(1) == Approx(1.).epsilon(Tolerance));
      REQUIRE(x_eigen(2) == Approx(-3.).epsilon(Tolerance));
    }
  }

#if USE_PETSC
  SECTION("PETSC solver") {
    std::shared_ptr<mpm::SolverBase<Eigen::SparseMatrix<double>>>
        petsc_matrix_solver =
            Factory<mpm::SolverBase<Eigen::SparseMatrix<double>>, unsigned,
                    double>::instance()
                ->create("KrylovPETSC", std::move(&&max_iter),
                         std::move(&&tolerance));

    SECTION("PETSC 3x3 solver") {
      // Initiate solver_type
      std::string solver_type = "cg";

      // Construct matrix and vector
      const auto& A = CreateSymmetricTestMatrix3x3();
      Eigen::VectorXd b(3);
      b << 6.0, 3.0, -9.0;

      // Solve
      const auto& x_eigen = petsc_matrix_solver->solve(A, b, "cg");

      // Check
      REQUIRE(x_eigen(0) == Approx(2.).epsilon(Tolerance));
      REQUIRE(x_eigen(1) == Approx(1.).epsilon(Tolerance));
      REQUIRE(x_eigen(2) == Approx(-3.).epsilon(Tolerance));
    }
  }
#endif
}
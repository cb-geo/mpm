#include <limits>
#include <memory>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>

#include "catch.hpp"
#include "cg_eigen.h"
#include "mpi_datatypes.h"
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

// Generate random SPD coefficient matrix with specified dimension and seed
Eigen::SparseMatrix<double> CreateRandomSymmetricTestMatrix(unsigned dim,
                                                            unsigned seed = 0) {
  std::srand(seed);
  Eigen::MatrixXd m = Eigen::MatrixXd::Random(dim, dim);
  m = (m + Eigen::MatrixXd::Constant(dim, dim, 1.2)) * 50;
  Eigen::SparseMatrix<double> A = m.sparseView();
  Eigen::SparseMatrix<double> AT = A.transpose();
  Eigen::SparseMatrix<double> matrix = A + AT;
  return matrix;
}

// Generate random RHS vector with specified dimension and seed
Eigen::VectorXd CreateRandomRHSVector(unsigned dim, unsigned seed = 0) {
  std::srand(seed);
  Eigen::VectorXd vector = Eigen::VectorXd::Random(dim);
  return vector;
}

TEST_CASE("Linear solver test", "[linear_solver]") {
  // Maximum iteration possible
  unsigned max_iter = 100;
  // Allowed solving tolerance
  double solve_tolerance = 1.E-5;
  // Tolerance
  const double Tolerance = 1.E-6;

  SECTION("Eigen solver") {

    // Construct Eigen solver
    std::shared_ptr<mpm::SolverBase<Eigen::SparseMatrix<double>>>
        eigen_matrix_solver =
            Factory<mpm::SolverBase<Eigen::SparseMatrix<double>>, unsigned,
                    double>::instance()
                ->create("CGEigen", std::move(max_iter),
                         std::move(solve_tolerance));

    SECTION("Eigen CG 3x3 solver") {
      // Initiate solver_type
      std::string solver_type = "cg";

      // Construct matrix and vector
      const auto& A = CreateSymmetricTestMatrix3x3();
      Eigen::VectorXd b(3);
      b << 6.0, 3.0, -9.0;

      // Solve
      const auto& x_eigen = eigen_matrix_solver->solve(A, b, solver_type);

      // Check
      REQUIRE(x_eigen(0) == Approx(2.).epsilon(Tolerance));
      REQUIRE(x_eigen(1) == Approx(1.).epsilon(Tolerance));
      REQUIRE(x_eigen(2) == Approx(-3.).epsilon(Tolerance));
    }

    SECTION("Eigen CG 5x5 solver") {
      // Initiate solver_type
      std::string solver_type = "cg";

      // Construct matrix and vector
      const auto& A = CreateRandomSymmetricTestMatrix(5);
      const auto& b = CreateRandomRHSVector(5);

      // Solve
      const auto& x_eigen = eigen_matrix_solver->solve(A, b, solver_type);

      // Check
      REQUIRE(x_eigen(0) == Approx(0.0120451).epsilon(Tolerance));
      REQUIRE(x_eigen(1) == Approx(0.0059998).epsilon(Tolerance));
      REQUIRE(x_eigen(2) == Approx(-0.0002562).epsilon(Tolerance));
      REQUIRE(x_eigen(3) == Approx(-0.00243321).epsilon(Tolerance));
      REQUIRE(x_eigen(4) == Approx(-0.0137986).epsilon(Tolerance));
    }

    SECTION("Eigen CG 50x50 solver") {
      // Initiate solver_type
      std::string solver_type = "cg";

      unsigned dim = 50;

      // Construct matrix and vector
      const auto& A = CreateRandomSymmetricTestMatrix(dim, 1);
      const auto& b = CreateRandomRHSVector(dim, 1);

      // Solve
      const auto& x_eigen = eigen_matrix_solver->solve(A, b, solver_type);

      // Check vector
      Eigen::VectorXd check(dim);
      check << 0.00953086, 0.0026714, 0.000905192, -0.00516796, -0.0118203,
          -0.0052406, 0.00351175, -0.0022648, 0.00881958, 0.00402084,
          -0.00374926, 0.00389573, 0.0278947, 0.0008647, 0.0108372, 0.00994053,
          -0.0152397, -0.0131816, 0.00237902, 0.0107677, 0.00205941, 0.0167911,
          -0.00424342, 0.000822485, 0.00284694, -0.000252545, 0.00401197,
          -0.0193391, -0.00497734, 0.00414726, -0.0191149, -0.00433217,
          -0.00845597, 0.00625816, 0.00608805, 0.00172621, -0.0319379,
          -0.0120693, -0.0199266, -0.0136789, 0.011165, 0.000247912,
          -0.00396261, -0.00169119, 0.00979176, 0.0168016, 0.0106939,
          -0.00551589, 0.0140779, 0.00212647;

      // Check
      for (unsigned i = 0; i < dim; i++)
        REQUIRE(x_eigen(i) == Approx(check(i)).epsilon(Tolerance));
    }
  }

#if USE_PETSC
  SECTION("PETSC solver") {
    std::shared_ptr<mpm::SolverBase<Eigen::SparseMatrix<double>>>
        petsc_matrix_solver =
            Factory<mpm::SolverBase<Eigen::SparseMatrix<double>>, unsigned,
                    double>::instance()
                ->create("KrylovPETSC", std::move(max_iter),
                         std::move(solve_tolerance));

    SECTION("PETSC CG 3x3 solver") {
      // Initiate solver_type and necessary mapping arrays
      std::string solver_type = "cg";
      std::vector<int> mapper{0, 1, 2};
      petsc_matrix_solver->assign_global_active_dof(3);
      petsc_matrix_solver->assign_rank_global_mapper(mapper);

      // Construct matrix and vector
      const auto& A = CreateSymmetricTestMatrix3x3();
      Eigen::VectorXd b(3);
      b << 6.0, 3.0, -9.0;

      // Solve
      const auto& x_eigen = petsc_matrix_solver->solve(A, b, solver_type);

      // Check
      REQUIRE(x_eigen(0) == Approx(2.).epsilon(Tolerance));
      REQUIRE(x_eigen(1) == Approx(1.).epsilon(Tolerance));
      REQUIRE(x_eigen(2) == Approx(-3.).epsilon(Tolerance));
    }

    SECTION("PETSC CG 5x5 solver") {
      // Initiate solver_type and necessary mapping arrays
      std::string solver_type = "cg";
      std::vector<int> mapper{0, 1, 2, 3, 4};
      petsc_matrix_solver->assign_global_active_dof(5);
      petsc_matrix_solver->assign_rank_global_mapper(mapper);

      // Construct matrix and vector
      const auto& A = CreateRandomSymmetricTestMatrix(5);
      const auto& b = CreateRandomRHSVector(5);

      // Solve
      const auto& x_eigen = petsc_matrix_solver->solve(A, b, solver_type);

      // Check
      REQUIRE(x_eigen(0) == Approx(0.0120451).epsilon(Tolerance));
      REQUIRE(x_eigen(1) == Approx(0.0059998).epsilon(Tolerance));
      REQUIRE(x_eigen(2) == Approx(-0.0002562).epsilon(Tolerance));
      REQUIRE(x_eigen(3) == Approx(-0.00243321).epsilon(Tolerance));
      REQUIRE(x_eigen(4) == Approx(-0.0137986).epsilon(Tolerance));
    }

    SECTION("PETSC CG 50x50 solver") {
      // Initiate solver_type and necessary mapping arrays
      std::string solver_type = "cg";
      unsigned dim = 50;

      std::vector<int> mapper;
      for (unsigned i = 0; i < dim; i++) mapper.push_back(i);
      petsc_matrix_solver->assign_global_active_dof(dim);
      petsc_matrix_solver->assign_rank_global_mapper(mapper);

      // Construct matrix and vector
      const auto& A = CreateRandomSymmetricTestMatrix(dim, 1);
      const auto& b = CreateRandomRHSVector(dim, 1);

      // Solve
      const auto& x_eigen = petsc_matrix_solver->solve(A, b, solver_type);

      // Check vector
      Eigen::VectorXd check(dim);
      check << 0.00953086, 0.0026714, 0.000905192, -0.00516796, -0.0118203,
          -0.0052406, 0.00351175, -0.0022648, 0.00881958, 0.00402084,
          -0.00374926, 0.00389573, 0.0278947, 0.0008647, 0.0108372, 0.00994053,
          -0.0152397, -0.0131816, 0.00237902, 0.0107677, 0.00205941, 0.0167911,
          -0.00424342, 0.000822485, 0.00284694, -0.000252545, 0.00401197,
          -0.0193391, -0.00497734, 0.00414726, -0.0191149, -0.00433217,
          -0.00845597, 0.00625816, 0.00608805, 0.00172621, -0.0319379,
          -0.0120693, -0.0199266, -0.0136789, 0.011165, 0.000247912,
          -0.00396261, -0.00169119, 0.00979176, 0.0168016, 0.0106939,
          -0.00551589, 0.0140779, 0.00212647;

      // Check
      for (unsigned i = 0; i < dim; i++)
        REQUIRE(x_eigen(i) == Approx(check(i)).epsilon(Tolerance));
    }
  }
#endif
}
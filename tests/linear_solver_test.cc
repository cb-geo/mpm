#include <limits>
#include <memory>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>

#include "catch.hpp"
#include "factory.h"
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
  // Allowed relative tolerance
  double rel_tolerance = 1.E-5;
  // Tolerance
  const double Tolerance = 1.E-6;

  SECTION("Eigen solver") {
    // Construct Eigen solver
    std::shared_ptr<mpm::SolverBase<Eigen::SparseMatrix<double>>>
        eigen_matrix_solver =
            Factory<mpm::SolverBase<Eigen::SparseMatrix<double>>, unsigned,
                    double>::instance()
                ->create("IterativeEigen", std::move(max_iter),
                         std::move(rel_tolerance));

    // Initiate sub_solver_type
    std::string sub_solver_type = "cg";
    REQUIRE_NOTHROW(eigen_matrix_solver->set_sub_solver_type(sub_solver_type));

    SECTION("Eigen CG 3x3 solver") {
      // Construct matrix and vector
      const auto& A = CreateSymmetricTestMatrix3x3();
      Eigen::VectorXd b(3);
      b << 6.0, 3.0, -9.0;

      // Solve
      auto x_eigen = eigen_matrix_solver->solve(A, b);

      // Check
      REQUIRE(x_eigen(0) == Approx(2.).epsilon(Tolerance));
      REQUIRE(x_eigen(1) == Approx(1.).epsilon(Tolerance));
      REQUIRE(x_eigen(2) == Approx(-3.).epsilon(Tolerance));

      SECTION("Eigen CG verbosity") {
        // Verbosity 1
        REQUIRE_NOTHROW(eigen_matrix_solver->set_verbosity(1));

        // Solve lscg
        REQUIRE_NOTHROW(eigen_matrix_solver->set_sub_solver_type("lscg"));
        x_eigen.setZero();
        x_eigen = eigen_matrix_solver->solve(A, b);

        // Solve bicgstab
        REQUIRE_NOTHROW(eigen_matrix_solver->set_sub_solver_type("bicgstab"));
        x_eigen = eigen_matrix_solver->solve(A, b);

        // Solve cg
        REQUIRE_NOTHROW(eigen_matrix_solver->set_sub_solver_type("cg"));
        x_eigen = eigen_matrix_solver->solve(A, b);

        // Verbosity 3
        REQUIRE_NOTHROW(eigen_matrix_solver->set_verbosity(3));
        x_eigen = eigen_matrix_solver->solve(A, b);
      }
    }

    SECTION("Eigen CG 5x5 solver") {
      // Construct matrix and vector
      const auto& A = CreateRandomSymmetricTestMatrix(5);
      const auto& b = CreateRandomRHSVector(5);

      // Solve
      const auto& x_eigen = eigen_matrix_solver->solve(A, b);

      // Check
      REQUIRE(x_eigen(0) == Approx(0.0120451).epsilon(Tolerance));
      REQUIRE(x_eigen(1) == Approx(0.0059998).epsilon(Tolerance));
      REQUIRE(x_eigen(2) == Approx(-0.0002562).epsilon(Tolerance));
      REQUIRE(x_eigen(3) == Approx(-0.00243321).epsilon(Tolerance));
      REQUIRE(x_eigen(4) == Approx(-0.0137986).epsilon(Tolerance));
    }

    SECTION("Eigen CG 50x50 solver") {
      unsigned dim = 50;

      // Construct matrix and vector
      const auto& A = CreateRandomSymmetricTestMatrix(dim, 1);
      const auto& b = CreateRandomRHSVector(dim, 1);

      // Solve
      auto x_eigen = eigen_matrix_solver->solve(A, b);

      // Check vector
      Eigen::VectorXd check(dim);
      check << 0.00953087, 0.00267145, 0.00090515, -0.0051681, -0.0118201,
          -0.00524066, 0.00351178, -0.00226483, 0.00881954, 0.00402078,
          -0.00374931, 0.00389581, 0.0278946, 0.00086471, 0.0108373, 0.00994061,
          -0.0152397, -0.0131817, 0.00237895, 0.0107676, 0.0020594, 0.0167911,
          -0.00424347, 0.000822635, 0.00284699, -0.000252711, 0.004012,
          -0.0193392, -0.00497733, 0.00414727, -0.0191148, -0.00433217,
          -0.00845599, 0.00625809, 0.00608814, 0.00172608, -0.0319379,
          -0.0120693, -0.0199266, -0.0136787, 0.011165, 0.000247912,
          -0.00396258, -0.00169115, 0.00979188, 0.0168016, 0.0106938,
          -0.00551589, 0.0140779, 0.00212647;

      // Check
      for (unsigned i = 0; i < dim; i++)
        REQUIRE(x_eigen(i) == Approx(check(i)).epsilon(Tolerance));

      SECTION("Eigen CG 50x50 sub solver type") {
        // Change CG to lscg
        REQUIRE_NOTHROW(eigen_matrix_solver->set_sub_solver_type("lscg"));

        // Solve
        x_eigen.setZero();
        x_eigen = eigen_matrix_solver->solve(A, b);

        // Check
        for (unsigned i = 0; i < dim; i++)
          REQUIRE(x_eigen(i) == Approx(check(i)).epsilon(Tolerance));

        // Change to bicgstab
        REQUIRE_NOTHROW(eigen_matrix_solver->set_sub_solver_type("bicgstab"));

        // Check for bi-cg-stab (expecting the results converge with more number
        // of iterations (max_iter is 100))
        REQUIRE_NOTHROW(eigen_matrix_solver->set_max_iteration(120));

        // Solve
        x_eigen = eigen_matrix_solver->solve(A, b);

        for (unsigned i = 0; i < dim; i++)
          REQUIRE(x_eigen(i) == Approx(check(i)).epsilon(Tolerance));

        // Change to error - unavailable subsolver type
        REQUIRE_NOTHROW(eigen_matrix_solver->set_sub_solver_type("error"));

        // Solve but x is zero size
        x_eigen = eigen_matrix_solver->solve(A, b);
        REQUIRE(x_eigen.size() == 0);
      }

      SECTION("Eigen CG convergence settings") {
        // Set tolerances and max iteration
        REQUIRE_NOTHROW(eigen_matrix_solver->set_tolerance(1.e-7));
        REQUIRE_NOTHROW(eigen_matrix_solver->set_abs_tolerance(1.e-7));
        REQUIRE_NOTHROW(eigen_matrix_solver->set_div_tolerance(1.e-7));

        // Reducing max iteration
        REQUIRE_NOTHROW(eigen_matrix_solver->set_max_iteration(10));

        // Solve
        x_eigen.setZero();
        x_eigen = eigen_matrix_solver->solve(A, b);

        // Check
        for (unsigned i = 0; i < dim; i++)
          REQUIRE(x_eigen(i) != Approx(check(i)).epsilon(Tolerance));

        // Increasing max iteration
        REQUIRE_NOTHROW(eigen_matrix_solver->set_max_iteration(100));

        // Solve
        x_eigen = eigen_matrix_solver->solve(A, b);

        // Check
        for (unsigned i = 0; i < dim; i++)
          REQUIRE(x_eigen(i) == Approx(check(i)).epsilon(Tolerance));

        // Set tolerances and max iteration to its limit
        REQUIRE_NOTHROW(eigen_matrix_solver->set_tolerance(
            std::numeric_limits<double>::epsilon()));
        REQUIRE_NOTHROW(eigen_matrix_solver->set_max_iteration(150));

        // Solve
        x_eigen = eigen_matrix_solver->solve(A, b);

        // Check
        for (unsigned i = 0; i < dim; i++)
          REQUIRE(x_eigen(i) == Approx(check(i)).epsilon(1.e-7));
      }

      SECTION("Eigen CG failure to converge") {
        // Check fail
        REQUIRE_NOTHROW(eigen_matrix_solver->set_tolerance(1.e-15));
        REQUIRE_NOTHROW(eigen_matrix_solver->set_max_iteration(10));

        // Solve cg
        REQUIRE_NOTHROW(eigen_matrix_solver->set_sub_solver_type("cg"));
        x_eigen.setZero();
        x_eigen = eigen_matrix_solver->solve(A, b);

        // Solve lscg
        REQUIRE_NOTHROW(eigen_matrix_solver->set_sub_solver_type("lscg"));
        x_eigen = eigen_matrix_solver->solve(A, b);

        // Solve bicgstab
        REQUIRE_NOTHROW(eigen_matrix_solver->set_sub_solver_type("bicgstab"));
        x_eigen = eigen_matrix_solver->solve(A, b);
      }
    }
  }

#if USE_PETSC
  SECTION("PETSC solver") {
    std::shared_ptr<mpm::SolverBase<Eigen::SparseMatrix<double>>>
        petsc_matrix_solver =
            Factory<mpm::SolverBase<Eigen::SparseMatrix<double>>, unsigned,
                    double>::instance()
                ->create("KrylovPETSC", std::move(max_iter),
                         std::move(rel_tolerance));

    // Initiate sub_solver_type
    std::string sub_solver_type = "cg";
    REQUIRE_NOTHROW(petsc_matrix_solver->set_sub_solver_type(sub_solver_type));

    SECTION("PETSC CG 3x3 solver") {
      std::vector<int> mapper{0, 1, 2};
      petsc_matrix_solver->assign_global_active_dof(3);
      petsc_matrix_solver->assign_rank_global_mapper(mapper);

      // Construct matrix and vector
      const auto& A = CreateSymmetricTestMatrix3x3();
      Eigen::VectorXd b(3);
      b << 6.0, 3.0, -9.0;

      // Solve
      auto x_eigen = petsc_matrix_solver->solve(A, b);

      // Check
      REQUIRE(x_eigen(0) == Approx(2.).epsilon(Tolerance));
      REQUIRE(x_eigen(1) == Approx(1.).epsilon(Tolerance));
      REQUIRE(x_eigen(2) == Approx(-3.).epsilon(Tolerance));

      SECTION("PETSC CG verbosity") {
        // Verbosity 3
        REQUIRE_NOTHROW(petsc_matrix_solver->set_verbosity(3));

        // Solve
        x_eigen.setZero();
        x_eigen = petsc_matrix_solver->solve(A, b);
      }
    }

    SECTION("PETSC CG 5x5 solver") {
      std::vector<int> mapper{0, 1, 2, 3, 4};
      petsc_matrix_solver->assign_global_active_dof(5);
      petsc_matrix_solver->assign_rank_global_mapper(mapper);

      // Construct matrix and vector
      const auto& A = CreateRandomSymmetricTestMatrix(5);
      const auto& b = CreateRandomRHSVector(5);

      // Solve
      const auto& x_eigen = petsc_matrix_solver->solve(A, b);

      // Check
      REQUIRE(x_eigen(0) == Approx(0.0120451).epsilon(Tolerance));
      REQUIRE(x_eigen(1) == Approx(0.0059998).epsilon(Tolerance));
      REQUIRE(x_eigen(2) == Approx(-0.0002562).epsilon(Tolerance));
      REQUIRE(x_eigen(3) == Approx(-0.00243321).epsilon(Tolerance));
      REQUIRE(x_eigen(4) == Approx(-0.0137986).epsilon(Tolerance));
    }

    SECTION("PETSC CG 50x50 solver") {
      unsigned dim = 50;

      std::vector<int> mapper;
      for (unsigned i = 0; i < dim; i++) mapper.push_back(i);
      petsc_matrix_solver->assign_global_active_dof(dim);
      petsc_matrix_solver->assign_rank_global_mapper(mapper);

      // Construct matrix and vector
      const auto& A = CreateRandomSymmetricTestMatrix(dim, 1);
      const auto& b = CreateRandomRHSVector(dim, 1);

      // Solve
      auto x_eigen = petsc_matrix_solver->solve(A, b);

      // Check vector
      Eigen::VectorXd check(dim);
      check << 0.00953087, 0.00267145, 0.00090515, -0.0051681, -0.0118201,
          -0.00524066, 0.00351178, -0.00226483, 0.00881954, 0.00402078,
          -0.00374931, 0.00389581, 0.0278946, 0.00086471, 0.0108373, 0.00994061,
          -0.0152397, -0.0131817, 0.00237895, 0.0107676, 0.0020594, 0.0167911,
          -0.00424347, 0.000822635, 0.00284699, -0.000252711, 0.004012,
          -0.0193392, -0.00497733, 0.00414727, -0.0191148, -0.00433217,
          -0.00845599, 0.00625809, 0.00608814, 0.00172608, -0.0319379,
          -0.0120693, -0.0199266, -0.0136787, 0.011165, 0.000247912,
          -0.00396258, -0.00169115, 0.00979188, 0.0168016, 0.0106938,
          -0.00551589, 0.0140779, 0.00212647;

      // Check
      for (unsigned i = 0; i < dim; i++)
        REQUIRE(x_eigen(i) == Approx(check(i)).epsilon(Tolerance));

      SECTION("PETSC CG 50x50 sub solver type") {
        // Change to gmres
        REQUIRE_NOTHROW(petsc_matrix_solver->set_sub_solver_type("gmres"));

        // Solve
        x_eigen.setZero();
        x_eigen = petsc_matrix_solver->solve(A, b);

        // Check
        for (unsigned i = 0; i < dim; i++)
          REQUIRE(x_eigen(i) == Approx(check(i)).epsilon(Tolerance));

        // Change to bicgstab
        REQUIRE_NOTHROW(petsc_matrix_solver->set_sub_solver_type("bicgstab"));

        // Solve
        x_eigen = petsc_matrix_solver->solve(A, b);

        for (unsigned i = 0; i < dim; i++)
          REQUIRE(x_eigen(i) == Approx(check(i)).epsilon(Tolerance));

        // Change to lsqr
        REQUIRE_NOTHROW(petsc_matrix_solver->set_sub_solver_type("lsqr"));

        // Solve - expected to diverge
        x_eigen = petsc_matrix_solver->solve(A, b);

        // Change to error - unavailable subsolver type (set to gmres)
        REQUIRE_NOTHROW(petsc_matrix_solver->set_sub_solver_type("error"));
        REQUIRE_NOTHROW(petsc_matrix_solver->set_verbosity(1));

        // Solve
        x_eigen = petsc_matrix_solver->solve(A, b);

        for (unsigned i = 0; i < dim; i++)
          REQUIRE(x_eigen(i) == Approx(check(i)).epsilon(Tolerance));
      }

      SECTION("PETSC CG 50x50 preconditioner type") {
        // Change to jacobi
        REQUIRE_NOTHROW(petsc_matrix_solver->set_preconditioner_type("jacobi"));

        // Solve - expected diverge
        x_eigen.setZero();
        x_eigen = petsc_matrix_solver->solve(A, b);

        // Change to block jacobi
        REQUIRE_NOTHROW(
            petsc_matrix_solver->set_preconditioner_type("bjacobi"));

        // Solve
        x_eigen = petsc_matrix_solver->solve(A, b);

        // Check
        for (unsigned i = 0; i < dim; i++)
          REQUIRE(x_eigen(i) == Approx(check(i)).epsilon(Tolerance));

        // Change to point block jacobi
        REQUIRE_NOTHROW(
            petsc_matrix_solver->set_preconditioner_type("pbjacobi"));

        // Solve - expected diverge
        x_eigen = petsc_matrix_solver->solve(A, b);

        // Change to additive Schwarz method
        REQUIRE_NOTHROW(petsc_matrix_solver->set_preconditioner_type("asm"));

        // Solve
        x_eigen = petsc_matrix_solver->solve(A, b);

        // Check
        for (unsigned i = 0; i < dim; i++)
          REQUIRE(x_eigen(i) == Approx(check(i)).epsilon(Tolerance));

        // Change to eisenstat
        REQUIRE_NOTHROW(
            petsc_matrix_solver->set_preconditioner_type("eisenstat"));

        // Solve - expected diverge
        x_eigen = petsc_matrix_solver->solve(A, b);

        // Change to icc
        REQUIRE_NOTHROW(petsc_matrix_solver->set_preconditioner_type("icc"));

        // Solve - expected diverge
        x_eigen = petsc_matrix_solver->solve(A, b);
      }

      SECTION("PETSC CG convergence settings") {
        // Set tolerance
        REQUIRE_NOTHROW(petsc_matrix_solver->set_tolerance(1.e-7));

        // Solve
        x_eigen.setZero();
        x_eigen = petsc_matrix_solver->solve(A, b);

        // Check
        for (unsigned i = 0; i < dim; i++)
          REQUIRE(x_eigen(i) == Approx(check(i)).epsilon(Tolerance));

        // Set tolerances and max iteration to its limit
        REQUIRE_NOTHROW(petsc_matrix_solver->set_tolerance(
            std::numeric_limits<double>::epsilon()));

        // Solve
        x_eigen = petsc_matrix_solver->solve(A, b);

        // Check
        for (unsigned i = 0; i < dim; i++)
          REQUIRE(x_eigen(i) == Approx(check(i)).epsilon(Tolerance));

        // Set abs_tolerance
        REQUIRE_NOTHROW(petsc_matrix_solver->set_abs_tolerance(
            std::numeric_limits<double>::epsilon()));

        // Solve
        x_eigen = petsc_matrix_solver->solve(A, b);

        // Check
        for (unsigned i = 0; i < dim; i++)
          REQUIRE(x_eigen(i) == Approx(check(i)).epsilon(Tolerance));

        // Set abs_tolerance - expecting faster convergence
        REQUIRE_NOTHROW(petsc_matrix_solver->set_abs_tolerance(1.e-7));

        // Solve
        x_eigen = petsc_matrix_solver->solve(A, b);

        // Check
        for (unsigned i = 0; i < dim; i++)
          REQUIRE(x_eigen(i) == Approx(check(i)).epsilon(Tolerance));

        // Set abs_tolerance - expecting to diverge
        REQUIRE_NOTHROW(petsc_matrix_solver->set_div_tolerance(
            std::numeric_limits<double>::epsilon()));

        // Solve
        x_eigen = petsc_matrix_solver->solve(A, b);

        // Set abs_tolerance - expecting faster convergence (default div_tol
        // is 1.e5)
        REQUIRE_NOTHROW(petsc_matrix_solver->set_div_tolerance(1.e1));

        // Solve
        x_eigen = petsc_matrix_solver->solve(A, b);

        // Check
        for (unsigned i = 0; i < dim; i++)
          REQUIRE(x_eigen(i) == Approx(check(i)).epsilon(Tolerance));
      }
    }
  }
#endif
}
#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"
#include <Eigen/Sparse>

#include "cg_eigen.h"
#include "solver_base.h"

// Generate 3x3 test matrix
Eigen::MatrixXd CreateSymmetricTestMatrix3x3() {
  Eigen::MatrixXd matrix(3, 3);
  matrix.setZero();
  matrix(0, 0) = 1.0;
  matrix(1, 1) = 0.0;
  matrix(2, 2) = 1.0;
  matrix(0, 1) = 1.0;
  matrix(0, 2) = 0.0;
  matrix(1, 2) = 1.0;
  matrix(1, 0) = matrix(0, 1);
  matrix(2, 0) = matrix(0, 2);
  matrix(2, 1) = matrix(1, 2);
  return matrix;
}

TEST_CASE("Linear solver conjugate gradient", "[linear_solver]") {

  // Maximum iteration possible
  unsigned max_iter = 100;
  // Allowed tolerance
  double tolerance = 10e-5;
  // Eigen::MatrixXd m = Eigen::MatrixXd::Random(testDim, testDim);
  // m = (m + Eigen::MatrixXd::Constant(testDim, testDim, 1.2)) * 50;
  // Eigen::SparseMatrix<double> A = m.sparseView();
  // Eigen::SparseMatrix<double> AT = A.transpose();
  // A = A + AT;
  // Eigen::VectorXd b = Eigen::VectorXd::Random(testDim);
  // double testTol = 0.000001;

  SECTION("Eigen solver") {

    std::shared_ptr<mpm::SolverBase<Eigen::SparseMatrix<double>>>
        eigen_matrix_solver =
            Factory<mpm::SolverBase<Eigen::SparseMatrix<double>>, unsigned,
                    double>::instance()
                ->create("EigenCG", std::move(max_iter), std::move(tolerance));
    // auto& x_eigen = eigen_matrix_solver->solve(A, b, "cg");
    // bool pass_compare = true;
    // pass_compare = (((A * x) - b).cwiseAbs().maxCoeff() < testTol);
    // REQUIRE(pass_compare);
  }

#if USE_PETSC
  // SECTION("PETSC solver") {
  //   auto petsc_matrix_solver =
  //       Factory<mpm::SolverBase<Eigen::SparseMatrix<double>>, unsigned,
  //               double>::instance()
  //           ->create("KrylovPETSC", std::move(&&max_iter),
  //                    std::move(&&tolerance));
  //   auto& x_petsc = petsc_matrix_solver->solve(A, b, "cg");
  //   bool pass_compare = true;
  //   pass_compare = (((A * x) - b).cwiseAbs().maxCoeff() < testTol);
  //   REQUIRE(pass_compare);
  // }
#endif
}
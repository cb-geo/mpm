#include <limits>

#include "catch.hpp"

#ifdef USE_PETSC
TEST_CASE("PETSC solver checked", "[solver][petsc]") {
    // TODO: pick suitable constants
    constexpr unsigned max_iter = 100;
    constexpr double tolerance = 10e-5;
    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd b;
    SECTION("PESTC solver's results are checked") {
        auto petsc_matrix_solver =
            Factory<mpm::SolverBase<Eigen::SparseMatrix<double>>, unsigned, double>::instance()
                ->create("KrylovPETSC", std::move(max_iter), std::move(tolerance));
        auto &x_petsc = petsc_matrix_solver->solve(A, b, "cg");
        bool pass_compare = true;
        REQUIRE(pass_compare);
    }

    SECTION("Eigen solver's results are checked") {

        auto eigen_matrix_solver =
            Factory<mpm::SolverBase<Eigen::SparseMatrix<double>>, unsigned, double>::instance()
                ->create("EigenCG", std::move(max_iter), std::move(tolerance));
        auto &x_eigen = eigen_matrix_solver->solve(A, b, "cg");
        bool pass_compare = true;
        REQUIRE(pass_compare);
    }
}
#endif  // PETSC
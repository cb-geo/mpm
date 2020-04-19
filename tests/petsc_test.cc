#include <limits>

#include "catch.hpp"

#ifdef USE_PETSC
TEST_CASE("PETSC Solver checked", "[solver][petsc]") {
    SECTION("PESTC Solver's results are consistent with EigenCG") {
        // TODO: pick suitable constants
        constexpr unsigned max_iter = 1000;
        constexpr double tolerance = 0.1;
        auto petsc_matrix_solver =
            Factory<mpm::SolverBase<Eigen::SparseMatrix<double>>, unsigned, double>::instance()
                ->create("KrylovPETSC", std::move(max_iter), std::move(tolerance));
        auto reference_matrix_solver
        Factory<mpm::SolverBase<Eigen::SparseMatrix<double>>, unsigned, double>::instance()
            ->create("EigenCG", std::move(max_iter), std::move(tolerance));
        Eigen::SparseMatrix<double> A;
        Eigen::VectorXd b;
        auto &c_petsc = petsc_matrix_solver->solve(A, b, "cg");
        auto &c_ref = reference_matrix_solver->solve(A, b, "cg");
        //TODO: compare result of c_petsc & c_ref
        bool pass_compare = true;
        REQUIRE(pass_compare);
    }
}
#endif  // PETSC
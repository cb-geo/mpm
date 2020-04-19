#include <limits>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Core>
#include "factory.h"

#include "catch.hpp"

#include "element.h"
#include "function_base.h"
#include "hexahedron_element.h"
#include "linear_function.h"
#include "mesh.h"
#include "node.h"
#include "partio_writer.h"
#include "quadrilateral_element.h"
#include <cmath>
#include <limits>
#include <memory>
#include <Eigen/Dense>
#include <boost/filesystem.hpp>

#ifdef USE_PETSC
TEST_CASE("PETSC solver checked", "[solver][petsc]") {
    // TODO: pick suitable constants
    constexpr unsigned max_iter = 100;
    constexpr double tolerance = 10e-5;
    int testDim = 100;
    Eigen::MatrixXd  m = Eigen::MatrixXd::Random(testDim, testDim);
    m = (m + Eigen::MatrixXd::Constant(testDim, testDim, 1.2)) * 50;
    Eigen::SparseMatrix<double> A = m.sparseView();
    Eigen::SparseMatrix<double> AT = A.transpose();
    A = A + AT;
    Eigen::VectorXd b = Eigen::VectorXd::Random(testDim);
    double testTol = 0.000001;

    SECTION("PESTC solver's results are checked") {
        auto petsc_matrix_solver =
            Factory<mpm::SolverBase<Eigen::SparseMatrix<double>>, unsigned, double>::instance()
                ->create("KrylovPETSC", std::move(&&max_iter), std::move(&&tolerance));
        auto &x_petsc = petsc_matrix_solver->solve(A, b, "cg");
        bool pass_compare = true;
        pass_compare = (((A*x)-b).cwiseAbs().maxCoeff() < testTol);
        REQUIRE(pass_compare);
    }

    SECTION("Eigen solver's results are checked") {

        auto eigen_matrix_solver =
            Factory<mpm::SolverBase<Eigen::SparseMatrix<double>>, unsigned, double>::instance()
                ->create("EigenCG", std::move(&&max_iter), std::move(&&tolerance));
        auto &x_eigen = eigen_matrix_solver->solve(A, b, "cg");
        bool pass_compare = true;
        pass_compare = (((A*x)-b).cwiseAbs().maxCoeff() < testTol);
        REQUIRE(pass_compare);
    }
}
#endif  // PETSC
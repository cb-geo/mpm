#include "assembler_base.h"
#include "assembler_eigen_semi_implicit_navierstokes.h"

#include "cg_eigen.h"
#include "solver_base.h"

// Assembler collections
// Asssembler 2D
static Register<mpm::AssemblerBase<2>,
                mpm::AssemblerEigenSemiImplicitNavierStokes<2>>
    assembler_eigen_semi_implicit_navierstokes_2d(
        "EigenSemiImplicitNavierStokes2D");
// Asssembler 3D
static Register<mpm::AssemblerBase<3>,
                mpm::AssemblerEigenSemiImplicitNavierStokes<3>>
    assembler_eigen_semi_implicit_navierstokes_3d(
        "EigenSemiImplicitNavierStokes3D");

// Linear Solver collections
// Eigen Conjugate Gradient
static Register<mpm::SolverBase<Eigen::SparseMatrix<double>>,
                mpm::CGEigen<Eigen::SparseMatrix<double>>, unsigned, double>
    solver_eigen_cg("EigenCG");
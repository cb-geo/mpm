#include "assembler_base.h"
#include "assembler_eigen_semi_implicit_navierstokes.h"
#include "assembler_parallel_semi_implicit_navierstokes.h"

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
// Asssembler parallel 2D
static Register<mpm::AssemblerBase<2>,
                mpm::AssemblerParallelSemiImplicitNavierStokes<2>>
    assembler_parallel_semi_implicit_navierstokes_2d(
        "ParallelSemiImplicitNavierStokes2D");
// Asssembler parallel 3D
static Register<mpm::AssemblerBase<3>,
                mpm::AssemblerParallelSemiImplicitNavierStokes<3>>
    assembler_parallel_semi_implicit_navierstokes_3d(
        "ParallelSemiImplicitNavierStokes3D");

// Linear Solver collections
// Solver 2D
static Register<mpm::SolverBase<2>, mpm::CGEigen<2>, unsigned, double>
    solver_eigen_cg_2d("EigenCG2D");
// Solver 3D
static Register<mpm::SolverBase<3>, mpm::CGEigen<3>, unsigned, double>
    solver_eigen_cg_3d("EigenCG3D");
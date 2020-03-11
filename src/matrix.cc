#include "matrix/assembler_base.h"
#include "matrix/assembler_eigen_semi_implicit_navierstokes.h"

#include "matrix/cg_eigen.h"
#include "matrix/solver_base.h"

// Assembler collections
// Asssembler 2D
static Register<mpm::AssemblerBase<2>,
                mpm::AssemblerEigenSemiImplicitNavierStokes<2>>
    AssemblerEigenSemiImplicitNavierStokes2d("EigenSemiImplicitNavierStokes2D");
// Asssembler 3D
static Register<mpm::AssemblerBase<3>,
                mpm::AssemblerEigenSemiImplicitNavierStokes<3>>
    AssemblerEigenSemiImplicitNavierStokes3d("EigenSemiImplicitNavierStokes3D");

// Linear Solver collections
// Solver 2D
static Register<mpm::SolverBase<2>, mpm::CGEigen<2>, unsigned, double>
    solvereigencg2d("EigenCG2D");
// Solver 3D
static Register<mpm::SolverBase<3>, mpm::CGEigen<3>, unsigned, double>
    solvereigencg3d("EigenCG3D");
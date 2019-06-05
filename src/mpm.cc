#include <memory>

#include "factory.h"
#include "io.h"
#include "mpm.h"
#include "mpm_explicit.h"
#include "mpm_explicit_twophase.h"

// 2D Explicit MPM USF
static Register<mpm::MPM, mpm::MPMExplicit<2>, std::unique_ptr<mpm::IO>&&>
    mpm_explicit_usf_2d("MPMExplicitUSF2D");

// 3D Explicit MPM USF
static Register<mpm::MPM, mpm::MPMExplicit<3>, std::unique_ptr<mpm::IO>&&>
    mpm_explicit_usf_3d("MPMExplicitUSF3D");

// 2D Explicit MPM USL
static Register<mpm::MPM, mpm::MPMExplicit<2>, std::unique_ptr<mpm::IO>&&>
    mpm_explicit_usl_2d("MPMExplicitUSL2D");

// 3D Explicit MPM USL
static Register<mpm::MPM, mpm::MPMExplicit<3>, std::unique_ptr<mpm::IO>&&>
    mpm_explicit_usl_3d("MPMExplicitUSL3D");

// 2D Explicit Two Phase MPM
static Register<mpm::MPM, mpm::MPMExplicitTwoPhase<2>,
                std::unique_ptr<mpm::IO>&&>
    mpm_explicit_twophase_2d("MPMExplicitTwoPhase2D");

// 3D Explicit Two Phase MPM
static Register<mpm::MPM, mpm::MPMExplicitTwoPhase<3>,
                std::unique_ptr<mpm::IO>&&>
    mpm_explicit_twophase_3d("MPMExplicitTwoPhase3D");

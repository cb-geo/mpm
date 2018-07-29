#include <memory>

#include "mpm.h"
#include "factory.h"
#include "io.h"
#include "mpm_explicit.h"

// 2D MPM
static Register<mpm::MPM<2>, mpm::MPMExplicit<2>, std::unique_ptr<mpm::IO>&&>
    mpm_explicit_2d("MPMExplicit2D");

// MPMExplicit3D (3 DoF, 1 Phase)
static Register<mpm::MPM<3>, mpm::MPMExplicit<3>, std::unique_ptr<mpm::IO>&&>
    mpm_explicit_3d("MPMExplicit3D");

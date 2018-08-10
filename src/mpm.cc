#include <memory>

#include "factory.h"
#include "io.h"
#include "mpm.h"
#include "mpm_explicit.h"
#include "mpm_explicit_usl.h"

// 2D Explicit MPM
static Register<mpm::MPM, mpm::MPMExplicit<2>, std::unique_ptr<mpm::IO>&&>
    mpm_explicit_2d("MPMExplicit2D");

// 3D Explicit MPM
static Register<mpm::MPM, mpm::MPMExplicit<3>, std::unique_ptr<mpm::IO>&&>
    mpm_explicit_3d("MPMExplicit3D");

// 2D Explicit MPM USL
static Register<mpm::MPM, mpm::MPMExplicitUSL<2>, std::unique_ptr<mpm::IO>&&>
    mpm_explicit_usl_2d("MPMExplicitUSL2D");

// 3D Explicit MPM USL
static Register<mpm::MPM, mpm::MPMExplicitUSL<3>, std::unique_ptr<mpm::IO>&&>
    mpm_explicit_usl_3d("MPMExplicitUSL3D");

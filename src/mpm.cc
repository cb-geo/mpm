#include <memory>

#include "factory.h"
#include "io.h"
#include "mpm.h"
#include "mpm_explicit.h"
#include "mpm_explicit_2p.h"

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

// 2D Explicit MPM Two phase USF
static Register<mpm::MPM, mpm::MPMExplicit2P<2>, std::unique_ptr<mpm::IO>&&>
    mpm_explicit_2p_usf_2d("MPMExplicit2PUSF2D");

// 3D Explicit MPM Two phase USF
static Register<mpm::MPM, mpm::MPMExplicit2P<3>, std::unique_ptr<mpm::IO>&&>
    mpm_explicit_2p_usf_3d("MPMExplicit2PUSF3D");

// 2D Explicit MPM Two phase USL
static Register<mpm::MPM, mpm::MPMExplicit2P<2>, std::unique_ptr<mpm::IO>&&>
    mpm_explicit_2p_usl_2d("MPMExplicit2PUSL2D");

// 3D Explicit MPM Two phase USL
static Register<mpm::MPM, mpm::MPMExplicit2P<3>, std::unique_ptr<mpm::IO>&&>
    mpm_explicit_2p_usl_3d("MPMExplicit2PUSL3D");

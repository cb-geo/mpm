#include <memory>

#include "factory.h"
#include "io.h"
#include "mpm.h"
#include "mpm_explicit.h"

namespace mpm {
// Stress update method
std::map<std::string, StressUpdate> stress_update = {
    {"usf", StressUpdate::USF},
    {"usl", StressUpdate::USL},
    {"musl", StressUpdate::MUSL}};
}  // namespace mpm

// 2D Explicit MPM USF
static Register<mpm::MPM, mpm::MPMExplicit<2>, std::shared_ptr<mpm::IO>&>
    mpm_explicit_usf_2d("MPMExplicitUSF2D");

// 3D Explicit MPM USF
static Register<mpm::MPM, mpm::MPMExplicit<3>, std::shared_ptr<mpm::IO>&>
    mpm_explicit_usf_3d("MPMExplicitUSF3D");

// 2D Explicit MPM USL
static Register<mpm::MPM, mpm::MPMExplicit<2>, std::shared_ptr<mpm::IO>&>
    mpm_explicit_usl_2d("MPMExplicitUSL2D");

// 3D Explicit MPM USL
static Register<mpm::MPM, mpm::MPMExplicit<3>, std::shared_ptr<mpm::IO>&>
    mpm_explicit_usl_3d("MPMExplicitUSL3D");

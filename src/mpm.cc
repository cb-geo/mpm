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

// 2D Explicit MPM
static Register<mpm::MPM, mpm::MPMExplicit<2>, std::unique_ptr<mpm::IO>&&>
    mpm_explicit_2d("MPMExplicit2D");

// 3D Explicit MPM
static Register<mpm::MPM, mpm::MPMExplicit<3>, std::unique_ptr<mpm::IO>&&>
    mpm_explicit_3d("MPMExplicit3D");

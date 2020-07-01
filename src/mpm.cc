#include <memory>

#include "factory.h"
#include "io.h"
#include "mpm.h"
#include "mpm_explicit.h"
#include "mpm_explicit_twophase.h"

namespace mpm {
// Stress update method
std::map<std::string, StressUpdate> stress_update = {
    {"usf", StressUpdate::USF},
    {"usl", StressUpdate::USL},
    {"musl", StressUpdate::MUSL}};
}  // namespace mpm

// 2D Explicit MPM
static Register<mpm::MPM, mpm::MPMExplicit<2>, const std::shared_ptr<mpm::IO>&>
    mpm_explicit_2d("MPMExplicit2D");

// 3D Explicit MPM
static Register<mpm::MPM, mpm::MPMExplicit<3>, const std::shared_ptr<mpm::IO>&>
    mpm_explicit_3d("MPMExplicit3D");

// 2D Explicit Two Phase MPM
static Register<mpm::MPM, mpm::MPMExplicitTwoPhase<2>,
                const std::shared_ptr<mpm::IO>&>
    mpm_explicit_twophase_2d("MPMExplicitTwoPhase2D");

// 3D Explicit Two Phase MPM
static Register<mpm::MPM, mpm::MPMExplicitTwoPhase<3>,
                const std::shared_ptr<mpm::IO>&>
    mpm_explicit_twophase_3d("MPMExplicitTwoPhase3D");
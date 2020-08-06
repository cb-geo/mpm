#include <memory>

#include "factory.h"
#include "io.h"
#include "mpm.h"
#include "mpm_explicit.h"

namespace mpm {
// Variable list
tsl::robin_map<std::string, VariableType> variables = {
    // Scalar variables
    {"mass", VariableType::Scalar},
    {"volume", VariableType::Scalar},
    // Vector variables
    {"displacements", VariableType::Vector},
    {"velocities", VariableType::Vector},
    // Tensor variables
    {"strains", VariableType::Tensor},
    {"stresses", VariableType::Tensor}};

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

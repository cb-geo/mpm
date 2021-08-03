#include <memory>

#include "factory.h"
#include "io.h"
#include "mpm.h"
#include "mpm_explicit.h"
#include "mpm_implicit_linear.h"

namespace mpm {
// 2D Explicit MPM
static Register<mpm::MPM, mpm::MPMExplicit<2>, const std::shared_ptr<mpm::IO>&>
    mpm_explicit_2d("MPMExplicit2D");

// 3D Explicit MPM
static Register<mpm::MPM, mpm::MPMExplicit<3>, const std::shared_ptr<mpm::IO>&>
    mpm_explicit_3d("MPMExplicit3D");

// 2D Implicit Linear MPM
static Register<mpm::MPM, mpm::MPMImplicitLinear<2>, const std::shared_ptr<mpm::IO>&>
    mpm_implicit_linear_2d("MPMImplicitLinear2D");

// 3D Implicit Linear MPM
static Register<mpm::MPM, mpm::MPMImplicitLinear<3>, const std::shared_ptr<mpm::IO>&>
    mpm_implicit_linear_3d("MPMImplicitLinear3D");

}  // namespace mpm

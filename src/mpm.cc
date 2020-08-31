#include <memory>

#include "factory.h"
#include "io.h"
#include "mpm.h"
#include "mpm_explicit.h"
#include "xmpm_explicit.h"

namespace mpm {
// 2D Explicit MPM
static Register<mpm::MPM, mpm::MPMExplicit<2>, const std::shared_ptr<mpm::IO>&>
    mpm_explicit_2d("MPMExplicit2D");

// 3D Explicit MPM
static Register<mpm::MPM, mpm::MPMExplicit<3>, const std::shared_ptr<mpm::IO>&>
    mpm_explicit_3d("MPMExplicit3D");

// 3D Explicit XMPM
static Register<mpm::MPM, mpm::XMPMExplicit<3>, const std::shared_ptr<mpm::IO>&>
    xmpm_explicit_3d("XMPMExplicit3D");
}  // namespace mpm

#include "particle.h"
#include "factory.h"
#include "fluid_particle.h"
#include "particle_base.h"

namespace mpm {
// ParticleType
std::map<std::string, int> ParticleType = {
    {"P2D", 0}, {"P3D", 1}, {"P2DFLUID", 2}, {"P3DFLUID", 3}};
std::map<int, std::string> ParticleTypeName = {
    {0, "P2D"}, {1, "P3D"}, {2, "P2DFLUID"}, {3, "P3DFLUID"}};
}  // namespace mpm

// Particle2D (2 Dim)
static Register<mpm::ParticleBase<2>, mpm::Particle<2>, mpm::Index,
                const Eigen::Matrix<double, 2, 1>&>
    particle2d("P2D");

// Particle3D (3 Dim)
static Register<mpm::ParticleBase<3>, mpm::Particle<3>, mpm::Index,
                const Eigen::Matrix<double, 3, 1>&>
    particle3d("P3D");

// Single phase (fluid) particle2D (2 Dim)
static Register<mpm::ParticleBase<2>, mpm::FluidParticle<2>, mpm::Index,
                const Eigen::Matrix<double, 2, 1>&>
    particle2dfluid("P2DFLUID");

// Single phase (fluid) particle3D (3 Dim)
static Register<mpm::ParticleBase<3>, mpm::FluidParticle<3>, mpm::Index,
                const Eigen::Matrix<double, 3, 1>&>
    particle3dfluid("P3DFLUID");
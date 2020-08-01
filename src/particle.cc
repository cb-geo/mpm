#include "particle.h"
#include "factory.h"
#include "particle_base.h"
#include "particle_xmpm.h"

// Particle2D (2 Dim)
static Register<mpm::ParticleBase<2>, mpm::Particle<2>, mpm::Index,
                const Eigen::Matrix<double, 2, 1>&>
    particle2d("P2D");

// Particle3D (3 Dim)
static Register<mpm::ParticleBase<3>, mpm::Particle<3>, mpm::Index,
                const Eigen::Matrix<double, 3, 1>&>
    particle3d("P3D");

// Particle3D_XMPM (3 Dim)
static Register<mpm::ParticleBase<3>, mpm::ParticleXMPM<3>, mpm::Index,
                const Eigen::Matrix<double, 3, 1>&>
    particle3dxmpm("P3DXMPM");

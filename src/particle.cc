#include "particle.h"
#include "factory.h"
#include "particle_base.h"

// Particle2D (2 Dim, 1 Phase)
static Register<mpm::ParticleBase<2>, mpm::Particle<2, 1>, mpm::Index,
                const Eigen::Matrix<double, 2, 1>&>
    particle2d("Particle2D");

// Particle3D (3 DoF, 1 Phase)
static Register<mpm::ParticleBase<3>, mpm::Particle<3, 1>, mpm::Index,
                const Eigen::Matrix<double, 3, 1>&>
    particle3d("Particle3D");

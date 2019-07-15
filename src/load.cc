#include "factory.h"
#include "loads/load_base.h"
#include "loads/step_load.h"
#include "loads/sine_load.h"

// 
static Register<mpm::LoadBase, mpm::StepLoad, unsigned,
                const std::map<double, double>&>
    stepload("STEP");

// Particle3D (3 DoF, 1 Phase)
static Register<mpm::LoadBase, mpm::SineLoad, unsigned> sineload("SINE");

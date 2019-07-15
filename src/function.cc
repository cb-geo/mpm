#include "factory.h"
#include "functions/function_base.h"
#include "functions/step_function.h"
#include "functions/sine_function.h"

// Step function
static Register<mpm::FunctionBase, mpm::StepFunction, unsigned,
                const std::vector<double>&, const std::vector<double>&>
    stepfunction("STEP");

// Sine function
static Register<mpm::FunctionBase, mpm::SineFunction, unsigned> sinefunction(
    "SINE");

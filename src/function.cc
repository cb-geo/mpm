#include "factory.h"
#include "functions/function_base.h"
#include "functions/step_function.h"

// Step function
static Register<mpm::FunctionBase, mpm::StepFunction, unsigned, const Json&>
    stepfunction("STEP");

#include "factory.h"
#include "function_base.h"
#include "linear_function.h"
#include "sin_function.h"

// Linear function
static Register<mpm::FunctionBase, mpm::LinearFunction, unsigned, const Json&>
    linearfn("Linear");

// Sin function
static Register<mpm::FunctionBase, mpm::SinFunction, unsigned, const Json&>
    sinfn("Sin");

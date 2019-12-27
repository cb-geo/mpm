#include "factory.h"
#include "function_base.h"
#include "linear_function.h"

// Linear function
static Register<mpm::FunctionBase, mpm::LinearFunction, unsigned, const Json&>
    linearfn("Linear");

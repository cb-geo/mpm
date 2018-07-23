#include "material/material.h"
#include "material/linear_elastic.h"

// LinearElastic
static Register<mpm::Material, mpm::LinearElastic, unsigned> linear_elastic(
    "LinearElastic");

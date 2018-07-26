#include "material/bingham.h"
#include "material/material.h"
#include "material/linear_elastic.h"

// LinearElastic 2D
static Register<mpm::Material<2>, mpm::LinearElastic<2>, unsigned>
    linear_elastic_2d("LinearElastic2D");

// LinearElastic 3D
static Register<mpm::Material<3>, mpm::LinearElastic<3>, unsigned>
    linear_elastic_3d("LinearElastic3D");

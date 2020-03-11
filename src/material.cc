#include "material.h"
#include "bingham.h"
#include "linear_elastic.h"
#include "mohr_coulomb.h"
#include "newtonian.h"
#include "norsand.h"

// Bingham 2D
static Register<mpm::Material<2>, mpm::Bingham<2>, unsigned, const Json&>
    bingham_2d("Bingham2D");

// Bingham 3D
static Register<mpm::Material<3>, mpm::Bingham<3>, unsigned, const Json&>
    bingham_3d("Bingham3D");

// LinearElastic 2D
static Register<mpm::Material<2>, mpm::LinearElastic<2>, unsigned, const Json&>
    linear_elastic_2d("LinearElastic2D");

// LinearElastic 3D
static Register<mpm::Material<3>, mpm::LinearElastic<3>, unsigned, const Json&>
    linear_elastic_3d("LinearElastic3D");

// MohrCoulomb 2D
static Register<mpm::Material<2>, mpm::MohrCoulomb<2>, unsigned, const Json&>
    mohr_coulomb_2d("MohrCoulomb2D");

// MohrCoulomb 3D
static Register<mpm::Material<3>, mpm::MohrCoulomb<3>, unsigned, const Json&>
    mohr_coulomb_3d("MohrCoulomb3D");

// Newtonian 2D
static Register<mpm::Material<2>, mpm::Newtonian<2>, unsigned, const Json&>
    newtonian_2d("Newtonian2D");

// Newtonian 3D
static Register<mpm::Material<3>, mpm::Newtonian<3>, unsigned, const Json&>
    newtonian_3d("Newtonian3D");

// Norsand 2D
static Register<mpm::Material<2>, mpm::NorSand<2>, unsigned, const Json&>
    nor_sand_2d("NorSand2D");

// Norsand 3D
static Register<mpm::Material<3>, mpm::NorSand<3>, unsigned, const Json&>
    nor_sand_3d("NorSand3D");

#include "shapefn.h"
#include "factory.h"
#include "hex_shapefn.h"
#include "quad_shapefn.h"

// Quadrilateral 4-noded shape function
static Register<mpm::ShapeFn<2>, mpm::QuadrilateralShapeFn<2, 4>> quad4("SFQ4");

// Quadrilateral 8-noded shape function
static Register<mpm::ShapeFn<2>, mpm::QuadrilateralShapeFn<2, 8>> quad8("SFQ8");

// Quadrilateral 9-noded shape function
static Register<mpm::ShapeFn<2>, mpm::QuadrilateralShapeFn<2, 9>> quad9("SFQ9");

// Hexahedron 8-noded shape function
static Register<mpm::ShapeFn<3>, mpm::HexahedronShapeFn<3, 8>> hex8("SFH8");

// Hexahedron 20-noded shape function
static Register<mpm::ShapeFn<3>, mpm::HexahedronShapeFn<3, 20>> hex20("SFH20");

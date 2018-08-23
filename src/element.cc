#include "element.h"
#include "factory.h"
#include "gimp_element.h"
#include "hexahedron_element.h"
#include "quadrilateral_element.h"

// Quadrilateral 4-noded element
static Register<mpm::Element<2>, mpm::QuadrilateralElement<2, 4>> quad4(
    "ED2Q4");

// Quadrilateral 8-noded element
static Register<mpm::Element<2>, mpm::QuadrilateralElement<2, 8>> quad8(
    "ED2Q8");

// Quadrilateral 9-noded element
static Register<mpm::Element<2>, mpm::QuadrilateralElement<2, 9>> quad9(
    "ED2Q9");

// Quadrilateral uGIMP element
static Register<mpm::Element<2>, mpm::GimpElement<2, 16>> GIMP("ED2GIMP");

// Hexahedron 8-noded element
static Register<mpm::Element<3>, mpm::HexahedronElement<3, 8>> hex8("ED3H8");

// Hexahedron 20-noded element
static Register<mpm::Element<3>, mpm::HexahedronElement<3, 20>> hex20("ED3H20");

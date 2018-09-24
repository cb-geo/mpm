#include "element.h"
#include "factory.h"
#include "hexahedron_element.h"
#include "quadrilateral_element.h"

// Quadrilateral 4-noded element
static Register<mpm::Element<2>, mpm::QuadrilateralElement<2, 4>> quad4(
    "E2DQ4");

// Quadrilateral 8-noded element
static Register<mpm::Element<2>, mpm::QuadrilateralElement<2, 8>> quad8(
    "E2DQ8");

// Quadrilateral 9-noded element
static Register<mpm::Element<2>, mpm::QuadrilateralElement<2, 9>> quad9(
    "E2DQ9");

// Hexahedron 8-noded element
static Register<mpm::Element<3>, mpm::HexahedronElement<3, 8>> hex8("E3DH8");

// Hexahedron 20-noded element
static Register<mpm::Element<3>, mpm::HexahedronElement<3, 20>> hex20("E3DH20");

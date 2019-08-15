#include "quadrature.h"
#include "factory.h"
#include "hexahedron_quadrature.h"
#include "quadrilateral_quadrature.h"
#include "triangle_quadrature.h"

// Triangle 1
static Register<mpm::Quadrature<2>, mpm::TriangleQuadrature<2, 1>>triangle1(
    "QT1");

// Triangle 3
static Register<mpm::Quadrature<2>, mpm::TriangleQuadrature<2, 3>>triangle3(
    "QT2");

// Quadrilateral 1
static Register<mpm::Quadrature<2>, mpm::QuadrilateralQuadrature<2, 1>>
    quadrilateral1("QQ1");

// Quadrilateral 4
static Register<mpm::Quadrature<2>, mpm::QuadrilateralQuadrature<2, 4>>
    quadrilateral4("QQ2");

// Quadrilateral 9
static Register<mpm::Quadrature<2>, mpm::QuadrilateralQuadrature<2, 9>>
    quadrilateral9("QQ3");

// Hexahedron 1
static Register<mpm::Quadrature<3>, mpm::HexahedronQuadrature<3, 1>>
    hexahedron1("QH1");

// Hexahedron 8
static Register<mpm::Quadrature<3>, mpm::HexahedronQuadrature<3, 8>>
    hexahedron8("QH2");

// Hexahedron 27
static Register<mpm::Quadrature<3>, mpm::HexahedronQuadrature<3, 27>>
    hexahedron("QH3");

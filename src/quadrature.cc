#include "quadrature.h"
#include "factory.h"
#include "hexahedron_quadrature.h"
#include "quadrilateral_quadrature.h"
#include "quadrilateral_random_quadrature.h"
#include "triangle_quadrature.h"

// Triangle 1
static Register<mpm::Quadrature<2>, mpm::TriangleQuadrature<2, 1>> triangle1(
    "QT1");

// Triangle 3
static Register<mpm::Quadrature<2>, mpm::TriangleQuadrature<2, 3>> triangle3(
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

// Quadrilateral 16
static Register<mpm::Quadrature<2>, mpm::QuadrilateralQuadrature<2, 16>>
    quadrilateral16("QQ4");

// Quadrilateral Random 1
static Register<mpm::Quadrature<2>, mpm::QuadrilateralRandomQuadrature<2, 1>>
    quadrilateralrandom1("QRQ1");

// Quadrilateral Random 4
static Register<mpm::Quadrature<2>, mpm::QuadrilateralRandomQuadrature<2, 4>>
    quadrilateralrandom4("QRQ4");

// Quadrilateral Random 9
static Register<mpm::Quadrature<2>, mpm::QuadrilateralRandomQuadrature<2, 9>>
    quadrilateralrandom9("QRQ9");

// Quadrilateral Random 16
static Register<mpm::Quadrature<2>, mpm::QuadrilateralRandomQuadrature<2, 16>>
    quadrilateralrandom16("QRQ16");

// Quadrilateral Random 25
static Register<mpm::Quadrature<2>, mpm::QuadrilateralRandomQuadrature<2, 25>>
    quadrilateralrandom25("QRQ25");

// Hexahedron 1
static Register<mpm::Quadrature<3>, mpm::HexahedronQuadrature<3, 1>>
    hexahedron1("QH1");

// Hexahedron 8
static Register<mpm::Quadrature<3>, mpm::HexahedronQuadrature<3, 8>>
    hexahedron8("QH2");

// Hexahedron 27
static Register<mpm::Quadrature<3>, mpm::HexahedronQuadrature<3, 27>>
    hexahedron27("QH3");

// Hexahedron 64
static Register<mpm::Quadrature<3>, mpm::HexahedronQuadrature<3, 64>>
    hexahedron64("QH4");

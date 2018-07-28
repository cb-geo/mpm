#include "cell.h"
#include "factory.h"
#include "shapefn.h"

// Cell2D
static Register<mpm::Cell<2>, mpm::Cell<2>, mpm::Index, unsigned,
                const std::shared_ptr<mpm::ShapeFn<2>>&>
    cell2d("C2D");

// Cell3D
static Register<mpm::Cell<3>, mpm::Cell<3>, mpm::Index, unsigned,
                const std::shared_ptr<mpm::ShapeFn<3>>&>
    cell3d("C3D");

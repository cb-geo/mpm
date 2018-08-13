#include "cell.h"
#include "element.h"
#include "factory.h"

// Cell1D
static Register<mpm::Cell<1>, mpm::Cell<1>, mpm::Index, unsigned,
                const std::shared_ptr<mpm::Element<1>>&>
    cell1d("C1D");

// Cell2D
static Register<mpm::Cell<2>, mpm::Cell<2>, mpm::Index, unsigned,
                const std::shared_ptr<mpm::Element<2>>&>
    cell2d("C2D");

// Cell3D
static Register<mpm::Cell<3>, mpm::Cell<3>, mpm::Index, unsigned,
                const std::shared_ptr<mpm::Element<3>>&>
    cell3d("C3D");

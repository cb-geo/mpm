#include "io_mesh.h"
#include "factory.h"
#include "io_mesh_ascii.h"

// IOMeshAscii
static Register<mpm::IOMesh<2>, mpm::IOMeshAscii<2>> iomesh_ascii_2d("Ascii2D");

// IOMeshAscii
static Register<mpm::IOMesh<3>, mpm::IOMeshAscii<3>> iomesh_ascii_3d("Ascii3D");

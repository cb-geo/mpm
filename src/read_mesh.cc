#include "read_mesh.h"
#include "factory.h"
#include "read_mesh_ascii.h"

// ReadMeshAscii
static Register<mpm::ReadMesh<2>, mpm::ReadMeshAscii<2>> readmesh_ascii_2d(
    "MeshAscii2D");

// ReadMeshAscii
static Register<mpm::ReadMesh<3>, mpm::ReadMeshAscii<3>> readmesh_ascii_3d(
    "MeshAscii3D");

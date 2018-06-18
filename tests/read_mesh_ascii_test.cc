#include "catch.hpp"

//! Alias for JSON
#include "json.hpp"
using Json = nlohmann::json;

#include "read_mesh_ascii.h"

// Check IO for input string
TEST_CASE("ReadMeshAscii is checked", "[ReadMesh][ReadMeshAscii]") {

  // Check input JSON
  SECTION("Check input JSON object") {
    const unsigned dim = 3;
    // Create a read_mesh object
    auto read_mesh = std::make_unique<mpm::ReadMeshAscii<dim>>();

    auto nodes = read_mesh->read_mesh_nodes("./mesh.dat");
  }
}

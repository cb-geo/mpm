#include "catch.hpp"

//! Alias for JSON
#include "json.hpp"
using Json = nlohmann::json;

#include "read_mesh_ascii.h"

// Check ReadMeshAscii
TEST_CASE("ReadMeshAscii is checked", "[ReadMesh][ReadMeshAscii]") {

  // Check read mesh nodes
  SECTION("Check read mesh nodes") {
    const unsigned dim = 3;

    // Vector of nodal coordinates
    std::vector<Eigen::Matrix<double, dim, 1>> coordinates;

    // Nodal coordinates
    Eigen::Matrix<double, dim, 1> node;

    // Cell 0
    // Node 0
    node << 0., 0., 0.;
    coordinates.emplace_back(node);
    // Node 1
    node << 0.5, 0., 0.;
    coordinates.emplace_back(node);
    // Node 2
    node << 0.5, 0.5, 0.;
    coordinates.emplace_back(node);
    // Node 3
    node << 0., 0.5, 0.;
    coordinates.emplace_back(node);
    // Node 4
    node << 0., 0., 0.5;
    coordinates.emplace_back(node);
    // Node 5
    node << 0.5, 0., 0.5;
    coordinates.emplace_back(node);
    // Node 6
    node << 0.5, 0.5, 0.5;
    coordinates.emplace_back(node);
    // Node 7
    node << 0., 0.5, 0.5;
    coordinates.emplace_back(node);

    // Cell 1
    // Node 8
    node << 0.5, 0., 0.;
    coordinates.emplace_back(node);
    // Node 9
    node << 1.0, 0., 0.;
    coordinates.emplace_back(node);
    // Node 10
    node << 1.0, 0.5, 0.;
    coordinates.emplace_back(node);
    // Node 11
    node << 0.5, 0.5, 0.;
    coordinates.emplace_back(node);
    // Node 12
    node << 0.5, 0., 0.5;
    coordinates.emplace_back(node);
    // Node 13
    node << 1.0, 0., 0.5;
    coordinates.emplace_back(node);
    // Node 14
    node << 1.0, 0.5, 0.5;
    coordinates.emplace_back(node);
    // Node 15
    node << 0.5, 0.5, 0.5;
    coordinates.emplace_back(node);

    // Cell with node ids
    std::vector<std::vector<unsigned>> cells{// cell #0
                                             {0, 1, 2, 3, 4, 5, 6, 7},
                                             // cell #1
                                             {8, 9, 10, 11, 12, 13, 14, 15}};

    // Dump JSON as an input file to be read
    std::ofstream file;
    file.open("mesh.txt");
    file.close();

    // Create a read_mesh object
    auto read_mesh = std::make_unique<mpm::ReadMeshAscii<dim>>();

    auto nodes = read_mesh->read_mesh_cells("./mesh.dat");
  }
}

#include "catch.hpp"

//! Alias for JSON
#include "json.hpp"
using Json = nlohmann::json;

#include "read_mesh_ascii.h"

// Check ReadMeshAscii
TEST_CASE("ReadMeshAscii is checked", "[ReadMesh][ReadMeshAscii]") {

  SECTION("ReadMeshAscii in 3D") {
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
    file << "! elementShape hexahedron\n";
    file << "! elementNumPoints 8\n";
    file << coordinates.size() << "\t" << cells.size() << "\n";

    // Write nodal coordinates
    for (const auto& coord : coordinates) {
      for (unsigned i = 0; i < coord.size(); ++i) {
        file << coord[i] << "\t";
      }
      file << "\n";
    }

    // Write cell node ids
    for (const auto& cell : cells) {
      for (const auto nid : cell) {
        file << nid << "\t";
      }
      file << "\n";
    }

    file.close();

    // Check read mesh nodes
    SECTION("Check read mesh nodes and cell ids") {
      const unsigned dim = 3;

      const double Tolerance = 1.E-7;
      // Create a read_mesh object
      auto read_mesh = std::make_unique<mpm::ReadMeshAscii<dim>>();

      // Check node ids in cell
      auto check_coords = read_mesh->read_mesh_nodes("./mesh.txt");
      // Check number of nodal coordinates
      REQUIRE(check_coords.size() == coordinates.size());

      // Check coordinates of nodes
      for (unsigned i = 0; i < coordinates.size(); ++i) {
        for (unsigned j = 0; j < dim; ++j) {
          REQUIRE(check_coords[i][j] ==
                  Approx(coordinates[i][j]).epsilon(Tolerance));
        }
      }
    }

    // Check read mesh cells
    SECTION("Check read mesh nodes and cell ids") {
      const unsigned dim = 3;

      const double Tolerance = 1.E-7;
      // Create a read_mesh object
      auto read_mesh = std::make_unique<mpm::ReadMeshAscii<dim>>();

      // Check node ids in cell
      auto check_node_ids = read_mesh->read_mesh_cells("./mesh.txt");
      // Check number of cells
      REQUIRE(check_node_ids.size() == cells.size());
      // Check node ids of cells
      for (unsigned i = 0; i < cells.size(); ++i) {
        for (unsigned j = 0; j < cells[i].size(); ++j) {
          REQUIRE(check_node_ids[i][j] == cells[i][j]);
        }
      }
    }
  }
}

#include "catch.hpp"

//! Alias for JSON
#include "json.hpp"
using Json = nlohmann::json;

#include "read_mesh_ascii.h"

// Check ReadMeshAscii
TEST_CASE("ReadMeshAscii is checked for 2D", "[ReadMesh][ReadMeshAscii][2D]") {

  // Dimension
  const unsigned dim = 2;
  // Tolerance
  const double Tolerance = 1.E-7;

  SECTION("Check mesh file ") {
    // Vector of nodal coordinates
    std::vector<Eigen::Matrix<double, dim, 1>> coordinates;

    // Nodal coordinates
    Eigen::Matrix<double, dim, 1> node;

    // Cell 0
    // Node 0
    node << 0., 0.;
    coordinates.emplace_back(node);
    // Node 1
    node << 0.5, 0.;
    coordinates.emplace_back(node);
    // Node 2
    node << 0.5, 0.5;
    coordinates.emplace_back(node);
    // Node 3
    node << 0., 0.5;
    coordinates.emplace_back(node);

    // Cell 1
    // Node 4
    node << 1.0, 0.;
    coordinates.emplace_back(node);
    // Node 5
    node << 1.0, 0.5;
    coordinates.emplace_back(node);

    // Cell with node ids
    std::vector<std::vector<unsigned>> cells{// cell #0
                                             {0, 1, 2, 3},
                                             // cell #1
                                             {1, 4, 5, 2}};

    // Dump mesh file as an input file to be read
    std::ofstream file;
    file.open("mesh-2d.txt");
    file << "! elementShape hexahedron\n";
    file << "! elementNumPoints 8\n";
    file << coordinates.size() << "\t" << cells.size() << "\n";

    // Write nodal coordinates
    for (const auto& coord : coordinates) {
      for (unsigned i = 0; i < coord.size(); ++i) file << coord[i] << "\t";
      file << "\n";
    }

    // Write cell node ids
    for (const auto& cell : cells) {
      for (auto nid : cell) file << nid << "\t";
      file << "\n";
    }

    file.close();

    // Check read mesh nodes
    SECTION("Check read mesh nodes") {
      // Create a read_mesh object
      auto read_mesh = std::make_unique<mpm::ReadMeshAscii<dim>>();

      // Try to read nodes from a non-existant file
      auto check_coords = read_mesh->read_mesh_nodes("mesh-missing.txt");
      // Check number of nodal coordinates
      REQUIRE(check_coords.size() == 0);

      check_coords = read_mesh->read_mesh_nodes("mesh-2d.txt");
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
    SECTION("Check read mesh cell ids") {
      // Create a read_mesh object
      auto read_mesh = std::make_unique<mpm::ReadMeshAscii<dim>>();

      // Try to read cells from a non-existant file
      auto check_node_ids = read_mesh->read_mesh_cells("mesh-missing.txt");
      // Check number of cells
      REQUIRE(check_node_ids.size() == 0);

      // Check node ids in cell
      check_node_ids = read_mesh->read_mesh_cells("mesh-2d.txt");
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

  SECTION("Check particles file") {
    // Vector of particle coordinates
    std::vector<Eigen::Matrix<double, dim, 1>> coordinates;
    coordinates.clear();

    // Particle coordinates
    Eigen::Matrix<double, dim, 1> particle;

    // Cell 0
    // Particle 0
    particle << 0.125, 0.125;
    coordinates.emplace_back(particle);
    // Particle 1
    particle << 0.25, 0.125;
    coordinates.emplace_back(particle);
    // Particle 2
    particle << 0.25, 0.25;
    coordinates.emplace_back(particle);
    // Particle 3
    particle << 0.125, 0.25;
    coordinates.emplace_back(particle);

    // Cell 1
    // Particle 4
    particle << 0.675, 0.125;
    coordinates.emplace_back(particle);
    // Particle 5
    particle << 0.85, 0.125;
    coordinates.emplace_back(particle);
    // Particle 6
    particle << 0.85, 0.25;
    coordinates.emplace_back(particle);
    // Particle 7
    particle << 0.675, 0.25;
    coordinates.emplace_back(particle);

    // Dump particles coordinates as an input file to be read
    std::ofstream file;
    file.open("particles-2d.txt");
    // Write particle coordinates
    for (const auto& coord : coordinates) {
      for (unsigned i = 0; i < coord.size(); ++i) {
        file << coord[i] << "\t";
      }
      file << "\n";
    }

    file.close();

    // Check read particle coordinates
    SECTION("Check particle coordinates") {
      // Create a read_mesh object
      auto read_mesh = std::make_unique<mpm::ReadMeshAscii<dim>>();

      // Try to read particle from a non-existant file
      auto particles = read_mesh->read_particles("particles-missing.txt");
      // Check number of particles
      REQUIRE(particles.size() == 0);

      // Check particle coordinates
      particles = read_mesh->read_particles("particles-2d.txt");
      // Check number of particles
      REQUIRE(particles.size() == coordinates.size());

      // Check coordinates of nodes
      for (unsigned i = 0; i < coordinates.size(); ++i) {
        for (unsigned j = 0; j < dim; ++j) {
          REQUIRE(particles[i][j] ==
                  Approx(coordinates[i][j]).epsilon(Tolerance));
        }
      }
    }
  }
  SECTION("Check velocity constraints file") {
    // Vector of particle coordinates
    std::vector<std::tuple<mpm::Index, unsigned, double>> velocity_constraints;

    // Constraint
    velocity_constraints.emplace_back(std::make_tuple(0, 0, 10.5));
    velocity_constraints.emplace_back(std::make_tuple(1, 1, -10.5));
    velocity_constraints.emplace_back(std::make_tuple(2, 0, -12.5));
    velocity_constraints.emplace_back(std::make_tuple(3, 1, 0.0));

    // Dump constratints as an input file to be read
    std::ofstream file;
    file.open("velocity-constraints-2d.txt");
    // Write particle coordinates
    for (const auto& velocity_constraint : velocity_constraints) {
      file << std::get<0>(velocity_constraint) << "\t";
      file << std::get<1>(velocity_constraint) << "\t";
      file << std::get<2>(velocity_constraint) << "\t";

      file << "\n";
    }

    file.close();

    // Check read velocity constraints
    SECTION("Check velocity constraints") {
      // Create a read_mesh object
      auto read_mesh = std::make_unique<mpm::ReadMeshAscii<dim>>();

      // Try to read constrtaints from a non-existant file
      auto constraints = read_mesh->read_velocity_constraints(
          "velocity-constraints-missing.txt");
      // Check number of constraints
      REQUIRE(constraints.size() == 0);

      // Check constraints
      constraints =
          read_mesh->read_velocity_constraints("velocity-constraints-2d.txt");
      // Check number of particles
      REQUIRE(constraints.size() == velocity_constraints.size());

      // Check coordinates of nodes
      for (unsigned i = 0; i < velocity_constraints.size(); ++i) {
        REQUIRE(
            std::get<0>(constraints.at(i)) ==
            Approx(std::get<0>(velocity_constraints.at(i))).epsilon(Tolerance));
        REQUIRE(
            std::get<1>(constraints.at(i)) ==
            Approx(std::get<1>(velocity_constraints.at(i))).epsilon(Tolerance));
        REQUIRE(
            std::get<2>(constraints.at(i)) ==
            Approx(std::get<2>(velocity_constraints.at(i))).epsilon(Tolerance));
      }
    }
  }
}

// Check ReadMeshAscii
TEST_CASE("ReadMeshAscii is checked for 3D", "[ReadMesh][ReadMeshAscii][3D]") {

  // Dimension
  const unsigned dim = 3;
  // Tolerance
  const double Tolerance = 1.E-7;

  SECTION("Check mesh file ") {
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
    node << 1.0, 0., 0.;
    coordinates.emplace_back(node);
    // Node 9
    node << 1.0, 0.5, 0.;
    coordinates.emplace_back(node);
    // Node 10
    node << 1.0, 0., 0.5;
    coordinates.emplace_back(node);
    // Node 11
    node << 1.0, 0.5, 0.5;
    coordinates.emplace_back(node);

    // Cell with node ids
    std::vector<std::vector<unsigned>> cells{// cell #0
                                             {0, 1, 2, 3, 4, 5, 6, 7},
                                             // cell #1
                                             {1, 8, 9, 2, 5, 10, 11, 6}};

    // Dump mesh file as an input file to be read
    std::ofstream file;
    file.open("mesh-3d.txt");
    file << "! elementShape hexahedron\n";
    file << "! elementNumPoints 8\n";
    file << coordinates.size() << "\t" << cells.size() << "\n";

    // Write nodal coordinates
    for (const auto& coord : coordinates) {
      for (unsigned i = 0; i < coord.size(); ++i) file << coord[i] << "\t";
      file << "\n";
    }

    // Write cell node ids
    for (const auto& cell : cells) {
      for (auto nid : cell) file << nid << "\t";
      file << "\n";
    }

    file.close();

    // Check read mesh nodes
    SECTION("Check read mesh nodes") {
      // Create a read_mesh object
      auto read_mesh = std::make_unique<mpm::ReadMeshAscii<dim>>();

      // Check nodal coordinates
      auto check_coords = read_mesh->read_mesh_nodes("mesh-missing.txt");
      // Check number of nodal coordinates
      REQUIRE(check_coords.size() == 0);

      // Check nodal coordinates
      check_coords = read_mesh->read_mesh_nodes("mesh-3d.txt");
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
    SECTION("Check read mesh cell ids") {
      // Create a read_mesh object
      auto read_mesh = std::make_unique<mpm::ReadMeshAscii<dim>>();

      // Try to read cells from a non-existant file
      auto check_node_ids = read_mesh->read_mesh_cells("mesh-missing.txt");
      // Check number of cells
      REQUIRE(check_node_ids.size() == 0);

      // Check node ids in cell
      check_node_ids = read_mesh->read_mesh_cells("mesh-3d.txt");
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

  SECTION("Check particles file") {
    // Vector of particle coordinates
    std::vector<Eigen::Matrix<double, dim, 1>> coordinates;

    // Particle coordinates
    Eigen::Matrix<double, dim, 1> particle;

    // Cell 0
    // Particle 0
    particle << 0.125, 0.125, 0.125;
    coordinates.emplace_back(particle);
    // Particle 1
    particle << 0.25, 0.125, 0.125;
    coordinates.emplace_back(particle);
    // Particle 2
    particle << 0.25, 0.25, 0.125;
    coordinates.emplace_back(particle);
    // Particle 3
    particle << 0.125, 0.25, 0.125;
    coordinates.emplace_back(particle);
    // Particle 4
    particle << 0.125, 0.125, 0.25;
    coordinates.emplace_back(particle);
    // Particle 5
    particle << 0.25, 0.125, 0.25;
    coordinates.emplace_back(particle);
    // Particle 6
    particle << 0.25, 0.25, 0.25;
    coordinates.emplace_back(particle);
    // Particle 7
    particle << 0.125, 0.25, 0.25;
    coordinates.emplace_back(particle);

    // Cell 1
    // Particle 8
    particle << 0.675, 0.125, 0.125;
    coordinates.emplace_back(particle);
    // Particle 9
    particle << 0.85, 0.125, 0.125;
    coordinates.emplace_back(particle);
    // Particle 10
    particle << 0.85, 0.25, 0.125;
    coordinates.emplace_back(particle);
    // Particle 11
    particle << 0.675, 0.25, 0.125;
    coordinates.emplace_back(particle);
    // Particle 12
    particle << 0.675, 0.125, 0.25;
    coordinates.emplace_back(particle);
    // Particle 13
    particle << 0.85, 0.125, 0.25;
    coordinates.emplace_back(particle);
    // Particle 14
    particle << 0.85, 0.25, 0.25;
    coordinates.emplace_back(particle);
    // Particle 15
    particle << 0.675, 0.25, 0.25;
    coordinates.emplace_back(particle);

    // Dump particles coordinates as an input file to be read
    std::ofstream file;
    file.open("particles-3d.txt");
    // Write particle coordinates
    for (const auto& coord : coordinates) {
      for (unsigned i = 0; i < coord.size(); ++i) {
        file << coord[i] << "\t";
      }
      file << "\n";
    }

    file.close();

    // Check read particles
    SECTION("Check particle coordinates") {
      // Create a read_mesh object
      auto read_mesh = std::make_unique<mpm::ReadMeshAscii<dim>>();

      // Try to read particle from a non-existant file
      auto particles = read_mesh->read_particles("particles-missing.txt");
      // Check number of particles
      REQUIRE(particles.size() == 0);

      // Check particle coordinates
      particles = read_mesh->read_particles("particles-3d.txt");
      // Check number of particles
      REQUIRE(particles.size() == coordinates.size());

      // Check coordinates of nodes
      for (unsigned i = 0; i < coordinates.size(); ++i) {
        for (unsigned j = 0; j < dim; ++j) {
          REQUIRE(particles[i][j] ==
                  Approx(coordinates[i][j]).epsilon(Tolerance));
        }
      }
    }
  }

  SECTION("Check velocity constraints file") {
    // Vector of particle coordinates
    std::vector<std::tuple<mpm::Index, unsigned, double>> velocity_constraints;

    // Constraint
    velocity_constraints.emplace_back(std::make_tuple(0, 0, 10.5));
    velocity_constraints.emplace_back(std::make_tuple(1, 1, -10.5));
    velocity_constraints.emplace_back(std::make_tuple(2, 2, -12.5));
    velocity_constraints.emplace_back(std::make_tuple(3, 0, 0.0));

    // Dump constratints as an input file to be read
    std::ofstream file;
    file.open("velocity-constraints-3d.txt");
    // Write particle coordinates
    for (const auto& velocity_constraint : velocity_constraints) {
      file << std::get<0>(velocity_constraint) << "\t";
      file << std::get<1>(velocity_constraint) << "\t";
      file << std::get<2>(velocity_constraint) << "\t";

      file << "\n";
    }

    file.close();

    // Check read velocity constraints
    SECTION("Check velocity constraints") {
      // Create a read_mesh object
      auto read_mesh = std::make_unique<mpm::ReadMeshAscii<dim>>();

      // Try to read constrtaints from a non-existant file
      auto constraints = read_mesh->read_velocity_constraints(
          "velocity-constraints-missing.txt");
      // Check number of constraints
      REQUIRE(constraints.size() == 0);

      // Check constraints
      constraints =
          read_mesh->read_velocity_constraints("velocity-constraints-3d.txt");
      // Check number of particles
      REQUIRE(constraints.size() == velocity_constraints.size());

      // Check coordinates of nodes
      for (unsigned i = 0; i < velocity_constraints.size(); ++i) {
        REQUIRE(
            std::get<0>(constraints.at(i)) ==
            Approx(std::get<0>(velocity_constraints.at(i))).epsilon(Tolerance));
        REQUIRE(
            std::get<1>(constraints.at(i)) ==
            Approx(std::get<1>(velocity_constraints.at(i))).epsilon(Tolerance));
        REQUIRE(
            std::get<2>(constraints.at(i)) ==
            Approx(std::get<2>(velocity_constraints.at(i))).epsilon(Tolerance));
      }
    }
  }
}

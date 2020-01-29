#include "catch.hpp"

//! Alias for JSON
#include "json.hpp"
using Json = nlohmann::json;

#include "io_mesh_ascii.h"

// Check IOMeshAscii
TEST_CASE("IOMeshAscii is checked for 2D", "[IOMesh][IOMeshAscii][2D]") {

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
      // Create a io_mesh object
      auto io_mesh = std::make_unique<mpm::IOMeshAscii<dim>>();

      // Try to read nodes from a non-existant file
      auto check_coords = io_mesh->read_mesh_nodes("mesh-missing.txt");
      // Check number of nodal coordinates
      REQUIRE(check_coords.size() == 0);

      check_coords = io_mesh->read_mesh_nodes("mesh-2d.txt");
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
      // Create a io_mesh object
      auto io_mesh = std::make_unique<mpm::IOMeshAscii<dim>>();

      // Try to read cells from a non-existant file
      auto check_node_ids = io_mesh->read_mesh_cells("mesh-missing.txt");
      // Check number of cells
      REQUIRE(check_node_ids.size() == 0);

      // Check node ids in cell
      check_node_ids = io_mesh->read_mesh_cells("mesh-2d.txt");
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
    file << coordinates.size() << "\n";
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
      // Create a io_mesh object
      auto io_mesh = std::make_unique<mpm::IOMeshAscii<dim>>();

      // Try to read particle from a non-existant file
      auto particles = io_mesh->read_particles("particles-missing.txt");
      // Check number of particles
      REQUIRE(particles.size() == 0);

      // Check particle coordinates
      particles = io_mesh->read_particles("particles-2d.txt");
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
    // Vector of velocity constraints
    std::vector<std::tuple<mpm::Index, unsigned, double>> velocity_constraints;

    // Constraint
    velocity_constraints.emplace_back(std::make_tuple(0, 0, 10.5));
    velocity_constraints.emplace_back(std::make_tuple(1, 1, -10.5));
    velocity_constraints.emplace_back(std::make_tuple(2, 0, -12.5));
    velocity_constraints.emplace_back(std::make_tuple(3, 1, 0.0));

    // Dump constraints as an input file to be read
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
      auto read_mesh = std::make_unique<mpm::IOMeshAscii<dim>>();

      // Try to read constraints from a non-existant file
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

  SECTION("Check friction constraints file") {
    // Vector of friction constraints
    std::vector<std::tuple<mpm::Index, unsigned, int, double>>
        friction_constraints;

    // Constraint
    friction_constraints.emplace_back(std::make_tuple(0, 0, 1, 10.5));
    friction_constraints.emplace_back(std::make_tuple(1, 1, -1, -10.5));
    friction_constraints.emplace_back(std::make_tuple(2, 0, 1, -12.5));
    friction_constraints.emplace_back(std::make_tuple(3, 1, -1, 0.0));

    // Dump constraints as an input file to be read
    std::ofstream file;
    file.open("friction-constraints-2d.txt");
    // Write particle coordinates
    for (const auto& friction_constraint : friction_constraints) {
      file << std::get<0>(friction_constraint) << "\t";
      file << std::get<1>(friction_constraint) << "\t";
      file << std::get<2>(friction_constraint) << "\t";
      file << std::get<3>(friction_constraint) << "\t";

      file << "\n";
    }

    file.close();

    // Check read friction constraints
    SECTION("Check friction constraints") {
      // Create a read_mesh object
      auto read_mesh = std::make_unique<mpm::IOMeshAscii<dim>>();

      // Try to read constraints from a non-existant file
      auto constraints = read_mesh->read_friction_constraints(
          "friction-constraints-missing.txt");
      // Check number of constraints
      REQUIRE(constraints.size() == 0);

      // Check constraints
      constraints =
          read_mesh->read_friction_constraints("friction-constraints-2d.txt");
      // Check number of particles
      REQUIRE(constraints.size() == friction_constraints.size());

      // Check coordinates of nodes
      for (unsigned i = 0; i < friction_constraints.size(); ++i) {
        REQUIRE(
            std::get<0>(constraints.at(i)) ==
            Approx(std::get<0>(friction_constraints.at(i))).epsilon(Tolerance));
        REQUIRE(
            std::get<1>(constraints.at(i)) ==
            Approx(std::get<1>(friction_constraints.at(i))).epsilon(Tolerance));
        REQUIRE(
            std::get<2>(constraints.at(i)) ==
            Approx(std::get<2>(friction_constraints.at(i))).epsilon(Tolerance));
        REQUIRE(
            std::get<3>(constraints.at(i)) ==
            Approx(std::get<3>(friction_constraints.at(i))).epsilon(Tolerance));
      }
    }
  }

  SECTION("Check forces file") {
    // Vector of particle forces
    std::vector<std::tuple<mpm::Index, unsigned, double>> nodal_forces;

    // Constraint
    nodal_forces.emplace_back(std::make_tuple(0, 0, 10.5));
    nodal_forces.emplace_back(std::make_tuple(1, 1, -10.5));
    nodal_forces.emplace_back(std::make_tuple(2, 0, -12.5));
    nodal_forces.emplace_back(std::make_tuple(3, 1, 0.0));

    // Dump constraints as an input file to be read
    std::ofstream file;
    file.open("forces-2d.txt");
    // Write particle coordinates
    for (const auto& force : nodal_forces) {
      file << std::get<0>(force) << "\t";
      file << std::get<1>(force) << "\t";
      file << std::get<2>(force) << "\t";

      file << "\n";
    }

    file.close();

    // Check read forces
    SECTION("Check forces") {
      // Create a read_mesh object
      auto read_mesh = std::make_unique<mpm::IOMeshAscii<dim>>();

      // Try to read constrtaints from a non-existant file
      auto forces = read_mesh->read_forces("forces-missing.txt");
      // Check number of forces
      REQUIRE(forces.size() == 0);

      // Check forces
      forces = read_mesh->read_forces("forces-2d.txt");
      // Check number of nodes
      REQUIRE(forces.size() == nodal_forces.size());

      // Check forces
      for (unsigned i = 0; i < nodal_forces.size(); ++i) {
        REQUIRE(std::get<0>(forces.at(i)) ==
                Approx(std::get<0>(nodal_forces.at(i))).epsilon(Tolerance));
        REQUIRE(std::get<1>(forces.at(i)) ==
                Approx(std::get<1>(nodal_forces.at(i))).epsilon(Tolerance));
        REQUIRE(std::get<2>(forces.at(i)) ==
                Approx(std::get<2>(nodal_forces.at(i))).epsilon(Tolerance));
      }
    }
  }

  SECTION("Check nodal euler angle file") {
    // Map of euler angles
    std::map<mpm::Index, Eigen::Matrix<double, dim, 1>> euler_angles;
    euler_angles.emplace(
        std::make_pair(0, (Eigen::Matrix<double, dim, 1>(10.5, 20.5))));
    euler_angles.emplace(
        std::make_pair(1, (Eigen::Matrix<double, dim, 1>(30.5, -40.5))));
    euler_angles.emplace(
        std::make_pair(2, (Eigen::Matrix<double, dim, 1>(-50.5, -60.5))));
    euler_angles.emplace(
        std::make_pair(3, (Eigen::Matrix<double, dim, 1>(-70.5, 80.5))));

    // Dump euler angles as an input file to be read
    std::ofstream file;
    file.open("nodal-euler-angles-2d.txt");
    // Write particle coordinates
    for (const auto& angles : euler_angles) {
      file << angles.first << "\t";
      for (unsigned i = 0; i < dim; ++i) file << (angles.second)(i) << "\t";
      file << "\n";
    }

    file.close();

    // Check read euler angles
    SECTION("Check euler angles") {
      // Create a io_mesh object
      auto io_mesh = std::make_unique<mpm::IOMeshAscii<dim>>();

      // Try to read euler angles from a non-existant file
      auto read_euler_angles =
          io_mesh->read_euler_angles("nodal-euler-angles-missing.txt");
      // Check number of euler angles
      REQUIRE(read_euler_angles.size() == 0);

      // Check euler angles
      read_euler_angles =
          io_mesh->read_euler_angles("nodal-euler-angles-2d.txt");

      // Check number of elements
      REQUIRE(read_euler_angles.size() == euler_angles.size());

      // Check euler angles
      for (unsigned i = 0; i < euler_angles.size(); ++i) {
        for (unsigned j = 0; j < dim; ++j) {
          REQUIRE(read_euler_angles.at(i)(j) ==
                  Approx(euler_angles.at(i)(j)).epsilon(Tolerance));
        }
      }
    }
  }

  SECTION("Check tractions file") {
    // Vector of particle tractions
    std::vector<std::tuple<mpm::Index, unsigned, double>> particles_tractions;

    // Constraint
    particles_tractions.emplace_back(std::make_tuple(0, 0, 10.5));
    particles_tractions.emplace_back(std::make_tuple(1, 1, -10.5));
    particles_tractions.emplace_back(std::make_tuple(2, 0, -12.5));
    particles_tractions.emplace_back(std::make_tuple(3, 1, 0.0));

    // Dump constraints as an input file to be read
    std::ofstream file;
    file.open("tractions-2d.txt");
    // Write particle coordinates
    for (const auto& traction : particles_tractions) {
      file << std::get<0>(traction) << "\t";
      file << std::get<1>(traction) << "\t";
      file << std::get<2>(traction) << "\t";

      file << "\n";
    }

    file.close();
  }

  SECTION("Check stresses file") {
    // Vector of particle stresses
    std::vector<Eigen::Matrix<double, 6, 1>> particles_stresses;

    // Stresses
    particles_stresses.emplace_back(
        Eigen::Matrix<double, 6, 1>::Constant(10.5));
    particles_stresses.emplace_back(
        Eigen::Matrix<double, 6, 1>::Constant(-12.5));
    particles_stresses.emplace_back(Eigen::Matrix<double, 6, 1>::Constant(0.4));

    // Dump initial stresses as an input file to be read
    std::ofstream file;
    file.open("particle-stresses-2d.txt");
    file << particles_stresses.size() << "\n";
    // Write particle coordinates
    for (const auto& stress : particles_stresses) {
      for (unsigned i = 0; i < stress.size(); ++i) file << stress[i] << "\t";
      file << "\n";
    }

    file.close();

    // Check read stresses
    SECTION("Check stresses") {
      // Create a io_mesh object
      auto io_mesh = std::make_unique<mpm::IOMeshAscii<dim>>();

      // Try to read stresses from a non-existant file
      auto stresses = io_mesh->read_particles_stresses("stresses-missing.txt");
      // Check number of stresses
      REQUIRE(stresses.size() == 0);

      // Check stresses
      stresses = io_mesh->read_particles_stresses("particle-stresses-2d.txt");
      // Check number of particles
      REQUIRE(stresses.size() == particles_stresses.size());

      // Check stresses
      for (unsigned i = 0; i < particles_stresses.size(); ++i) {
        for (unsigned j = 0; j < particles_stresses.at(0).size(); ++j)
          REQUIRE(stresses.at(i)[j] ==
                  Approx(particles_stresses.at(i)[j]).epsilon(Tolerance));
      }
    }
  }
}

// Check IOMeshAscii
TEST_CASE("IOMeshAscii is checked for 3D", "[IOMesh][IOMeshAscii][3D]") {

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
      // Create a io_mesh object
      auto io_mesh = std::make_unique<mpm::IOMeshAscii<dim>>();

      // Check nodal coordinates
      auto check_coords = io_mesh->read_mesh_nodes("mesh-missing.txt");
      // Check number of nodal coordinates
      REQUIRE(check_coords.size() == 0);

      // Check nodal coordinates
      check_coords = io_mesh->read_mesh_nodes("mesh-3d.txt");
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
      // Create a io_mesh object
      auto io_mesh = std::make_unique<mpm::IOMeshAscii<dim>>();

      // Try to read cells from a non-existant file
      auto check_node_ids = io_mesh->read_mesh_cells("mesh-missing.txt");
      // Check number of cells
      REQUIRE(check_node_ids.size() == 0);

      // Check node ids in cell
      check_node_ids = io_mesh->read_mesh_cells("mesh-3d.txt");
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
    file << coordinates.size() << "\n";
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
      // Create a io_mesh object
      auto io_mesh = std::make_unique<mpm::IOMeshAscii<dim>>();

      // Try to read particle from a non-existant file
      auto particles = io_mesh->read_particles("particles-missing.txt");
      // Check number of particles
      REQUIRE(particles.size() == 0);

      // Check particle coordinates
      particles = io_mesh->read_particles("particles-3d.txt");
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
      auto read_mesh = std::make_unique<mpm::IOMeshAscii<dim>>();

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

  SECTION("Check friction constraints file") {
    // Vector of friction constraints
    std::vector<std::tuple<mpm::Index, unsigned, int, double>>
        friction_constraints;

    // Constraint
    friction_constraints.emplace_back(std::make_tuple(0, 0, 1, 10.5));
    friction_constraints.emplace_back(std::make_tuple(1, 1, -1, -10.5));
    friction_constraints.emplace_back(std::make_tuple(2, 2, 1, -12.5));
    friction_constraints.emplace_back(std::make_tuple(3, 2, -1, 0.0));

    // Dump constraints as an input file to be read
    std::ofstream file;
    file.open("friction-constraints-3d.txt");
    // Write particle coordinates
    for (const auto& friction_constraint : friction_constraints) {
      file << std::get<0>(friction_constraint) << "\t";
      file << std::get<1>(friction_constraint) << "\t";
      file << std::get<2>(friction_constraint) << "\t";
      file << std::get<3>(friction_constraint) << "\t";

      file << "\n";
    }

    file.close();

    // Check read friction constraints
    SECTION("Check friction constraints") {
      // Create a read_mesh object
      auto read_mesh = std::make_unique<mpm::IOMeshAscii<dim>>();

      // Try to read constraints from a non-existant file
      auto constraints = read_mesh->read_friction_constraints(
          "friction-constraints-missing.txt");
      // Check number of constraints
      REQUIRE(constraints.size() == 0);

      // Check constraints
      constraints =
          read_mesh->read_friction_constraints("friction-constraints-3d.txt");
      // Check number of particles
      REQUIRE(constraints.size() == friction_constraints.size());

      // Check coordinates of nodes
      for (unsigned i = 0; i < friction_constraints.size(); ++i) {
        REQUIRE(
            std::get<0>(constraints.at(i)) ==
            Approx(std::get<0>(friction_constraints.at(i))).epsilon(Tolerance));
        REQUIRE(
            std::get<1>(constraints.at(i)) ==
            Approx(std::get<1>(friction_constraints.at(i))).epsilon(Tolerance));
        REQUIRE(
            std::get<2>(constraints.at(i)) ==
            Approx(std::get<2>(friction_constraints.at(i))).epsilon(Tolerance));
        REQUIRE(
            std::get<3>(constraints.at(i)) ==
            Approx(std::get<3>(friction_constraints.at(i))).epsilon(Tolerance));
      }
    }
  }

  SECTION("Check forces file") {
    // Vector of particle forces
    std::vector<std::tuple<mpm::Index, unsigned, double>> nodal_forces;

    // Constraint
    nodal_forces.emplace_back(std::make_tuple(0, 0, 10.5));
    nodal_forces.emplace_back(std::make_tuple(1, 1, -10.5));
    nodal_forces.emplace_back(std::make_tuple(2, 0, -12.5));
    nodal_forces.emplace_back(std::make_tuple(3, 2, 0.0));

    // Dump constraints as an input file to be read
    std::ofstream file;
    file.open("forces-3d.txt");
    // Write particle coordinates
    for (const auto& force : nodal_forces) {
      file << std::get<0>(force) << "\t";
      file << std::get<1>(force) << "\t";
      file << std::get<2>(force) << "\t";

      file << "\n";
    }

    file.close();

    // Check read forces
    SECTION("Check forces") {
      // Create a read_mesh object
      auto read_mesh = std::make_unique<mpm::IOMeshAscii<dim>>();

      // Try to read constrtaints from a non-existant file
      auto forces = read_mesh->read_forces("forces-missing.txt");
      // Check number of forces
      REQUIRE(forces.size() == 0);

      // Check forces
      forces = read_mesh->read_forces("forces-3d.txt");
      // Check number of nodes
      REQUIRE(forces.size() == nodal_forces.size());

      // Check forces
      for (unsigned i = 0; i < nodal_forces.size(); ++i) {
        REQUIRE(std::get<0>(forces.at(i)) ==
                Approx(std::get<0>(nodal_forces.at(i))).epsilon(Tolerance));
        REQUIRE(std::get<1>(forces.at(i)) ==
                Approx(std::get<1>(nodal_forces.at(i))).epsilon(Tolerance));
        REQUIRE(std::get<2>(forces.at(i)) ==
                Approx(std::get<2>(nodal_forces.at(i))).epsilon(Tolerance));
      }
    }
  }

  SECTION("Check nodal euler angle file") {
    // Map of euler angles
    std::map<mpm::Index, Eigen::Matrix<double, dim, 1>> euler_angles;
    euler_angles.emplace(
        std::make_pair(0, (Eigen::Matrix<double, dim, 1>(10.5, 20.5, 30.5))));
    euler_angles.emplace(
        std::make_pair(1, (Eigen::Matrix<double, dim, 1>(40.5, -50.5, -60.5))));
    euler_angles.emplace(std::make_pair(
        2, (Eigen::Matrix<double, dim, 1>(-70.5, -80.5, -90.5))));
    euler_angles.emplace(std::make_pair(
        3, (Eigen::Matrix<double, dim, 1>(-100.5, 110.5, 120.5))));

    // Dump euler angles as an input file to be read
    std::ofstream file;
    file.open("nodal-euler-angles-3d.txt");
    // Write particle coordinates
    for (const auto& angles : euler_angles) {
      file << angles.first << "\t";
      for (unsigned i = 0; i < dim; ++i) file << (angles.second)(i) << "\t";
      file << "\n";
    }

    file.close();

    // Check read euler angles
    SECTION("Check euler angles") {
      // Create a io_mesh object
      auto io_mesh = std::make_unique<mpm::IOMeshAscii<dim>>();

      // Try to read euler angles from a non-existant file
      auto read_euler_angles =
          io_mesh->read_euler_angles("nodal-euler-angles-missing.txt");
      // Check number of euler angles
      REQUIRE(read_euler_angles.size() == 0);

      // Check euler angles
      read_euler_angles =
          io_mesh->read_euler_angles("nodal-euler-angles-3d.txt");

      // Check number of elements
      REQUIRE(read_euler_angles.size() == euler_angles.size());

      // Check euler angles
      for (unsigned i = 0; i < euler_angles.size(); ++i) {
        for (unsigned j = 0; j < dim; ++j) {
          REQUIRE(read_euler_angles.at(i)(j) ==
                  Approx(euler_angles.at(i)(j)).epsilon(Tolerance));
        }
      }
    }
  }

  SECTION("Check tractions file") {
    // Vector of particle tractions
    std::vector<std::tuple<mpm::Index, unsigned, double>> particles_tractions;

    // Constraint
    particles_tractions.emplace_back(std::make_tuple(0, 0, 10.5));
    particles_tractions.emplace_back(std::make_tuple(1, 1, -10.5));
    particles_tractions.emplace_back(std::make_tuple(2, 0, -12.5));
    particles_tractions.emplace_back(std::make_tuple(3, 1, 0.0));

    // Dump constraints as an input file to be read
    std::ofstream file;
    file.open("tractions-3d.txt");
    // Write particle coordinates
    for (const auto& traction : particles_tractions) {
      file << std::get<0>(traction) << "\t";
      file << std::get<1>(traction) << "\t";
      file << std::get<2>(traction) << "\t";

      file << "\n";
    }

    file.close();
  }

  SECTION("Check stresses file") {
    // Vector of particle stresses
    std::vector<Eigen::Matrix<double, 6, 1>> particles_stresses;

    // Stresses
    particles_stresses.emplace_back(
        Eigen::Matrix<double, 6, 1>::Constant(100.5));
    particles_stresses.emplace_back(
        Eigen::Matrix<double, 6, 1>::Constant(-112.5));
    particles_stresses.emplace_back(
        Eigen::Matrix<double, 6, 1>::Constant(0.46));

    // Dump initial stresses as an input file to be read
    std::ofstream file;
    file.open("particle-stresses-3d.txt");
    file << particles_stresses.size() << "\n";
    // Write particle coordinates
    for (const auto& stress : particles_stresses) {
      for (unsigned i = 0; i < stress.size(); ++i) file << stress[i] << "\t";
      file << "\n";
    }

    file.close();

    // Check read stresses
    SECTION("Check stresses") {
      // Create a io_mesh object
      auto io_mesh = std::make_unique<mpm::IOMeshAscii<dim>>();

      // Try to read stresses from a non-existant file
      auto stresses = io_mesh->read_particles_stresses("stresses-missing.txt");
      // Check number of stresses
      REQUIRE(stresses.size() == 0);

      // Check stresses
      stresses = io_mesh->read_particles_stresses("particle-stresses-3d.txt");
      // Check number of particles
      REQUIRE(stresses.size() == particles_stresses.size());

      // Check stresses
      for (unsigned i = 0; i < particles_stresses.size(); ++i) {
        for (unsigned j = 0; j < particles_stresses.at(0).size(); ++j)
          REQUIRE(stresses.at(i)[j] ==
                  Approx(particles_stresses.at(i)[j]).epsilon(Tolerance));
      }
    }
  }
}

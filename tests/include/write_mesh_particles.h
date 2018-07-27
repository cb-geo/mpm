#include <fstream>

#include "mpm.h"

namespace mpm_test {

// Write JSON Configuration file
bool write_json(unsigned dim, const std::string& file_name) {
  // Make json object with input files
  // 2D
  std::string dimension = "2d";
  auto particle_type = "P2D";
  auto node_type = "N2D";
  auto cell_type = "SFQ4";
  auto mesh_reader = "Ascii2D";
  std::string material = "LinearElastic2D";
  std::vector<double> gravity{{0., -9.81}};

  // 3D
  if (dim == 3) {
    dimension = "3d";
    particle_type = "P3D";
    node_type = "N3D";
    cell_type = "SFH8";
    mesh_reader = "Ascii3D";
    material = "LinearElastic";
    gravity.clear();
    gravity = {0., 0., -9.81};
  }

  Json json_file = {
      {"title", "Example JSON Input for MPM"},
      {"input_files",
       {{"mesh", "mesh-" + dimension + ".txt"},
        {"constraints", "mesh_constraints.txt"},
        {"particles", "particles-" + dimension + ".txt"},
        {"initial_stresses", "initial_soil_stress.txt"},
        {"materials", "materials.txt"},
        {"traction", "traction.txt"}}},
      {"mesh",
       {{"mesh_reader", mesh_reader},
        {"node_type", node_type},
        {"material_id", 1},
        {"cell_type", cell_type},
        {"particle_type", particle_type}}},
      {"materials",
       {{{"id", 0},
         {"type", material},
         {"density", 1000.},
         {"youngs_modulus", 1.0E+8},
         {"poisson_ratio", 0.495}},
        {{"id", 1},
         {"type", material},
         {"density", 2300.},
         {"youngs_modulus", 1.5E+6},
         {"poisson_ratio", 0.25}}}},
      {"analysis",
       {{"dt", 0.001},
        {"nsteps", 10},
        {"gravity", gravity},
        {"boundary_friction", 0.5},
        {"damping", {{"damping", true}, {"damping_ratio", 0.02}}},
        {"newmark", {{"newmark", true}, {"gamma", 0.5}, {"beta", 0.25}}}}},
      {"post_processing", {{"path", "results/"}, {"output_steps", 10}}}};

  // Dump JSON as an input file to be read
  std::ofstream file;
  file.open((file_name + "-" + dimension + ".json").c_str());
  file << json_file.dump(2);
  file.close();

  return true;
}

// Write Mesh file in 2D
bool write_mesh_2d() {
  // Dimension
  const unsigned dim = 2;

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
  return true;
}

// Write particles file in 2D
bool write_particles_2d() {
  const unsigned dim = 2;
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
  return true;
}

// Write mesh file in 3D
bool write_mesh_3d() {

  // Dimension
  const unsigned dim = 3;
  // Tolerance
  const double Tolerance = 1.E-7;

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
  return true;
}

// Write particles file in 3D
bool write_particles_3d() {
  const unsigned dim = 3;
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
  return true;
}

}  // namespace mpm_test

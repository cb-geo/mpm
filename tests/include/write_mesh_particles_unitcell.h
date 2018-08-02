#include <fstream>

#include "mpm.h"

namespace mpm_test {

// Write JSON Configuration file
bool write_json_unitcell(unsigned dim, const std::string& file_name);

// Write Mesh file in 2D
bool write_mesh_2d_unitcell();
// Write particles file in 2D
bool write_particles_2d_unitcell();

// Write mesh file in 3D
bool write_mesh_3d_unitcell();
// Write particles file in 3D
<<<<<<< HEAD
bool write_particles_3d_unitcell();
=======
bool write_particles_3d_unitcell() {
  const unsigned dim = 3;
  // Vector of particle coordinates
  std::vector<Eigen::Matrix<double, dim, 1>> coordinates;

  // Particle coordinates
  Eigen::Matrix<double, dim, 1> particle;

  // Cell 0
  // Particle 0
  particle << 0.25, 0.25, 0.25;
  coordinates.emplace_back(particle);
  // Particle 1
  particle << 0.75, 0.25, 0.25;
  coordinates.emplace_back(particle);
  // Particle 2
  particle << 0.25, 0.75, 0.25;
  coordinates.emplace_back(particle);
  // Particle 3
  particle << 0.75, 0.75, 0.25;
  coordinates.emplace_back(particle);
  // Particle 4
  particle << 0.25, 0.25, 0.75;
  coordinates.emplace_back(particle);
  // Particle 5
  particle << 0.25, 0.25, 0.75;
  coordinates.emplace_back(particle);
  // Particle 6
  particle << 0.75, 0.75, 0.75;
  coordinates.emplace_back(particle);
  // Particle 7
  particle << 0.75, 0.75, 0.75;
  coordinates.emplace_back(particle);
  /*
    // Cell 1
    // Particle 8
    particle << 1.25, 0.25, 0.25;
    coordinates.emplace_back(particle);
    // Particle 9
    particle << 1.75, 0.25, 0.25;
    coordinates.emplace_back(particle);
    // Particle 10
    particle << 1.25, 0.75, 0.25;
    coordinates.emplace_back(particle);
    // Particle 11
    particle << 1.75, 0.75, 0.25;
    coordinates.emplace_back(particle);
    // Particle 12
    particle << 1.25, 0.25, 0.75;
    coordinates.emplace_back(particle);
    // Particle 13
    particle << 1.25, 0.25, 0.75;
    coordinates.emplace_back(particle);
    // Particle 14
    particle << 1.75, 0.75, 0.75;
    coordinates.emplace_back(particle);
    // Particle 15
    particle << 1.75, 0.75, 0.75;
    coordinates.emplace_back(particle);
    */
  // Dump particles coordinates as an input file to be read
  std::ofstream file;
  file.open("particles-3d-unitcell.txt");
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
>>>>>>> :construction: refactor bingham tcc

}  // namespace mpm_test

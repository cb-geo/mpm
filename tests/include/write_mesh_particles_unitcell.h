#include <fstream>

#include "mpm.h"

namespace mpm_test {

// Write JSON Configuration file
bool write_json_unitcell(unsigned dim, const std::string& analysis,
                         const std::string& stress_update,
                         const std::string& file_name);

// Write Mesh file in 2D
bool write_mesh_2d_unitcell();
// Write particles file in 2D
bool write_particles_2d_unitcell();

// Write mesh file in 3D
bool write_mesh_3d_unitcell();
// Write particles file in 3D
bool write_particles_3d_unitcell();

}  // namespace mpm_test

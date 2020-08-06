#include <fstream>

#include "mpm.h"

namespace mpm_test {

// Write JSON Configuration file
bool write_json(unsigned dim, bool resume, const std::string& analysis,
                const std::string& stress_update, const std::string& file_name);

// Write JSON Configuration file for twophase
bool write_json_twophase(unsigned dim, bool resume, const std::string& analysis,
                         const std::string& stress_update,
                         const std::string& file_name);

// Write JSON Entity Set
bool write_entity_set();

// Write Mesh file in 2D
bool write_mesh_2d();
// Write particles file in 2D
bool write_particles_2d();

// Write mesh file in 3D
bool write_mesh_3d();
// Write particles file in 3D
bool write_particles_3d();

}  // namespace mpm_test

#ifndef MPM_HDF5_H_
#define MPM_HDF5_H_

// HDF5
#include "hdf5.h"
#include "hdf5_hl.h"

#include "data_types.h"

namespace mpm {
// Define a struct of particle
typedef struct HDF5Particle {
  // Index
  mpm::Index id;
  // Mass
  double mass;
  // Volume
  double volume;
  // Pressure
  double pressure;
  // Coordinates
  double coord_x, coord_y, coord_z;
  // Displacement
  double displacement_x, displacement_y, displacement_z;
  // Natural particle size
  double nsize_x, nsize_y, nsize_z;
  // Velocity
  double velocity_x, velocity_y, velocity_z;
  // Stresses
  double stress_xx, stress_yy, stress_zz;
  double tau_xy, tau_yz, tau_xz;
  // Strains
  double strain_xx, strain_yy, strain_zz;
  double gamma_xy, gamma_yz, gamma_xz;
  // Volumetric strain centroid
  double epsilon_v;
  // Index
  mpm::Index cell_id;
  // Status
  bool status;
  // Material id
  unsigned material_id;
  // Number of state variables
  unsigned nstate_vars;
  // State variables (init to zero)
  double svars[20] = {0};
} HDF5Particle;

namespace hdf5 {
namespace particle {
const hsize_t NFIELDS = 53;

const size_t dst_size = sizeof(HDF5Particle);

// Destination offset
extern const size_t dst_offset[NFIELDS];

// Destination size
extern const size_t dst_sizes[NFIELDS];

// Define particle field information
extern const char* field_names[NFIELDS];

// Initialize field types
extern const hid_t field_type[NFIELDS];

}  // namespace particle
}  // namespace hdf5

}  // namespace mpm

#endif  // MPM_HDF5_H_

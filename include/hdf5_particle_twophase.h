#ifndef MPM_HDF5_TWOPHASE_H_
#define MPM_HDF5_TWOPHASE_H_

// HDF5 Particle
#include "hdf5_particle.h"

namespace mpm {
// Define a struct of particle
typedef struct HDF5ParticleTwoPhase : HDF5Particle {
  // Liquid Mass
  double liquid_mass;
  // Liquid Velocity
  double liquid_velocity_x, liquid_velocity_y, liquid_velocity_z;
  // Porosity
  double porosity;
  // Liquid Saturation
  double liquid_saturation;
  // Material id
  unsigned liquid_material_id;
  // Number of state variables
  unsigned nliquid_state_vars;
  // State variables (init to zero)
  double liquid_svars[5] = {0};
  // Destructor
  virtual ~HDF5ParticleTwoPhase() = default;
} HDF5ParticleTwoPhase;

namespace hdf5 {
namespace particletwophase {
const hsize_t NFIELDS = 66;

const size_t dst_size = sizeof(HDF5ParticleTwoPhase);

// Destination offset
extern const size_t dst_offset[NFIELDS];

// Destination size
extern const size_t dst_sizes[NFIELDS];

// Define particle field information
extern const char* field_names[NFIELDS];

// Initialize field types
extern const hid_t field_type[NFIELDS];

}  // namespace particletwophase
}  // namespace hdf5

}  // namespace mpm

#endif  // MPM_HDF5_TWOPHASE_H_

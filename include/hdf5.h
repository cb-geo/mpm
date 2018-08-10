#ifndef MPM_HDF5_H_
#define MPM_HDF5_H_

// HDF5
#include "hdf5.h"
#include "hdf5_hl.h"

namespace mpm {

// Define a struct of particle
typedef struct HDF5Particle {
  mpm::Index id;
  bool status;
  double coord_x, coord_y, coord_z;
  double velocity_x, velocity_y, velocity_z;
  double stress_xx, stress_yy, stress_zz;
  double tau_xy, tau_yz, tau_xz;
  double strain_xx, strain_yy, strain_zz;
  double gamma_xy, gamma_yz, gamma_xz;
} HDF5Particle;

}  // namespace mpm

#endif  // MPM_HDF5_H_

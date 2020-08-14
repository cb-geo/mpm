#ifndef MPM_MPI_HDF5_PARTICLE_TWOPHASE_H_
#define MPM_MPI_HDF5_PARTICLE_TWOPHASE_H_

// Offset macro
#include <cstddef>

// MPI
#ifdef USE_MPI
#include "mpi.h"

#include "hdf5_particle_twophase.h"

namespace mpm {
//! Initialize MPI particle data types
inline MPI_Datatype register_mpi_particle_type(
    const HDF5ParticleTwoPhase& particle) {
  // Number of blocks to create
  const unsigned nblocks = 43;
  // Array containing the length of each block
  int lengths[nblocks] = {1, 1, 1, 1,  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                          1, 1, 1, 1,  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                          1, 1, 1, 20, 1, 1, 1, 1, 1, 1, 1, 1, 5};
  // Array containing the displacement for each block, expressed in bytes. The
  // displacement is the distance between the start of the MPI datatype created
  // and the start of the block
  const MPI_Aint displacements[nblocks] = {
      // Solid phase
      offsetof(HDF5ParticleTwoPhase, id),              // id
      offsetof(HDF5ParticleTwoPhase, mass),            // mass
      offsetof(HDF5ParticleTwoPhase, volume),          // volume
      offsetof(HDF5ParticleTwoPhase, pressure),        // pressure
      offsetof(HDF5ParticleTwoPhase, coord_x),         // coord_x
      offsetof(HDF5ParticleTwoPhase, coord_y),         // coord_y
      offsetof(HDF5ParticleTwoPhase, coord_z),         // coord_z
      offsetof(HDF5ParticleTwoPhase, displacement_x),  // disp_x
      offsetof(HDF5ParticleTwoPhase, displacement_y),  // disp_y
      offsetof(HDF5ParticleTwoPhase, displacement_z),  // disp_z
      offsetof(HDF5ParticleTwoPhase, nsize_x),         // nsize_x
      offsetof(HDF5ParticleTwoPhase, nsize_y),         // nsize_y
      offsetof(HDF5ParticleTwoPhase, nsize_z),         // nsize_z
      offsetof(HDF5ParticleTwoPhase, velocity_x),      // vel_x
      offsetof(HDF5ParticleTwoPhase, velocity_y),      // vel_y
      offsetof(HDF5ParticleTwoPhase, velocity_z),      // vel_z
      offsetof(HDF5ParticleTwoPhase, stress_xx),       // stress_xx
      offsetof(HDF5ParticleTwoPhase, stress_yy),       // stress_yy
      offsetof(HDF5ParticleTwoPhase, stress_zz),       // stress_zz
      offsetof(HDF5ParticleTwoPhase, tau_xy),          // tau_xy
      offsetof(HDF5ParticleTwoPhase, tau_yz),          // tau_yz
      offsetof(HDF5ParticleTwoPhase, tau_xz),          // tau_xz
      offsetof(HDF5ParticleTwoPhase, strain_xx),       // strain_xx
      offsetof(HDF5ParticleTwoPhase, strain_yy),       // strain_yy
      offsetof(HDF5ParticleTwoPhase, strain_zz),       // strain_zz
      offsetof(HDF5ParticleTwoPhase, gamma_xy),        // gamma_xy
      offsetof(HDF5ParticleTwoPhase, gamma_yz),        // gamma_yz
      offsetof(HDF5ParticleTwoPhase, gamma_xz),        // gamma_xz
      offsetof(HDF5ParticleTwoPhase, epsilon_v),       // epsv
      offsetof(HDF5ParticleTwoPhase, cell_id),         // cellid
      offsetof(HDF5ParticleTwoPhase, status),          // status
      offsetof(HDF5ParticleTwoPhase, material_id),     // material_id
      offsetof(HDF5ParticleTwoPhase, nstate_vars),     // nstate_vars
      offsetof(HDF5ParticleTwoPhase, svars),           // state_vars
      // Fluid phase
      offsetof(HDF5ParticleTwoPhase, liquid_mass),         // liquid_mass
      offsetof(HDF5ParticleTwoPhase, liquid_velocity_x),   // liquid_vel_x
      offsetof(HDF5ParticleTwoPhase, liquid_velocity_y),   // liquid_vel_y
      offsetof(HDF5ParticleTwoPhase, liquid_velocity_z),   // liquid_vel_z
      offsetof(HDF5ParticleTwoPhase, porosity),            // porosity
      offsetof(HDF5ParticleTwoPhase, liquid_saturation),   // liquid_saturation
      offsetof(HDF5ParticleTwoPhase, liquid_material_id),  // liquid_material_id
      offsetof(HDF5ParticleTwoPhase, nliquid_state_vars),  // nliquid_state_vars
      offsetof(HDF5ParticleTwoPhase, liquid_svars)         // liquid_state_vars
  };

  // Array containing the MPI datatypes to replicate to make each block.
  MPI_Datatype types[nblocks] = {
      // Solid phase
      MPI_UNSIGNED_LONG_LONG,  // id
      MPI_DOUBLE,              // mass
      MPI_DOUBLE,              // volume
      MPI_DOUBLE,              // pressure
      MPI_DOUBLE,              // coord_x
      MPI_DOUBLE,              // coord_y
      MPI_DOUBLE,              // coord_z
      MPI_DOUBLE,              // disp_x
      MPI_DOUBLE,              // disp_y
      MPI_DOUBLE,              // disp_z
      MPI_DOUBLE,              // nsize_x
      MPI_DOUBLE,              // nsize_y
      MPI_DOUBLE,              // nsize_z
      MPI_DOUBLE,              // vel_x
      MPI_DOUBLE,              // vel_y
      MPI_DOUBLE,              // vel_z
      MPI_DOUBLE,              // stress_xx
      MPI_DOUBLE,              // stress_yy
      MPI_DOUBLE,              // stress_zz
      MPI_DOUBLE,              // tau_xx
      MPI_DOUBLE,              // tau_yy
      MPI_DOUBLE,              // tau_zz
      MPI_DOUBLE,              // strain_xx
      MPI_DOUBLE,              // strain_yy
      MPI_DOUBLE,              // strain_zz
      MPI_DOUBLE,              // gamma_xy
      MPI_DOUBLE,              // gamma_yz
      MPI_DOUBLE,              // gamma_zx
      MPI_DOUBLE,              // epsv
      MPI_UNSIGNED_LONG_LONG,  // cell_id
      MPI_C_BOOL,              // status
      MPI_UNSIGNED,            // material_id
      MPI_UNSIGNED,            // nstate_vars
      MPI_DOUBLE,              // state variables
      // Fluid phase
      MPI_DOUBLE,    // liquid_mass
      MPI_DOUBLE,    // liquid_velocity_x
      MPI_DOUBLE,    // liquid_velocity_y
      MPI_DOUBLE,    // liquid_velocity_z
      MPI_DOUBLE,    // porosity
      MPI_DOUBLE,    // liquid_saturation
      MPI_UNSIGNED,  // liquid_material_id
      MPI_UNSIGNED,  // nliquid_state_vars
      MPI_DOUBLE,    // liquid_svars
  };

  // Create the Particle datatype
  MPI_Datatype MPIParticle;

  // Create particle data types
  MPI_Type_create_struct(nblocks, lengths, displacements, types, &MPIParticle);
  MPI_Type_commit(&MPIParticle);
  return MPIParticle;
}
}  // namespace mpm

#endif  // MPI

#endif  // MPM_MPI_HDF5_PARTICLE_TWOPHASE_H_

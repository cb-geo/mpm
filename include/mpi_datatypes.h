#ifndef MPM_MPI_HDF5_PARTICLE_H_
#define MPM_MPI_HDF5_PARTICLE_H_

// Offset macro
#include <cstddef>

// MPI
#ifdef USE_MPI
#include "mpi.h"

#include "hdf5_particle.h"

namespace mpm {
//! Initialize MPI particle data types
inline MPI_Datatype register_mpi_particle_type(const HDF5Particle& particle) {
  // Number of blocks to create
  const unsigned nblocks = 34;
  // Array containing the length of each block
  int lengths[nblocks] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 20};
  // Array containing the displacement for each block, expressed in bytes. The
  // displacement is the distance between the start of the MPI datatype created
  // and the start of the block
  const std::size_t size_ull = sizeof(unsigned long long);

  const MPI_Aint displacements[nblocks] = {
      offsetof(HDF5Particle, id),              // id
      offsetof(HDF5Particle, mass),            // mass
      offsetof(HDF5Particle, volume),          // volume
      offsetof(HDF5Particle, pressure),        // pressure
      offsetof(HDF5Particle, coord_x),         // coord_x
      offsetof(HDF5Particle, coord_y),         // coord_y
      offsetof(HDF5Particle, coord_z),         // coord_z
      offsetof(HDF5Particle, displacement_x),  // disp_x
      offsetof(HDF5Particle, displacement_y),  // disp_y
      offsetof(HDF5Particle, displacement_z),  // disp_z
      offsetof(HDF5Particle, nsize_x),         // nsize_x
      offsetof(HDF5Particle, nsize_y),         // nsize_y
      offsetof(HDF5Particle, nsize_z),         // nsize_z
      offsetof(HDF5Particle, velocity_x),      // vel_x
      offsetof(HDF5Particle, velocity_y),      // vel_y
      offsetof(HDF5Particle, velocity_z),      // vel_z
      offsetof(HDF5Particle, stress_xx),       // stress_xx
      offsetof(HDF5Particle, stress_yy),       // stress_yy
      offsetof(HDF5Particle, stress_zz),       // stress_zz
      offsetof(HDF5Particle, tau_xy),          // tau_xy
      offsetof(HDF5Particle, tau_yz),          // tau_yz
      offsetof(HDF5Particle, tau_xz),          // tau_xz
      offsetof(HDF5Particle, strain_xx),       // strain_xx
      offsetof(HDF5Particle, strain_yy),       // strain_yy
      offsetof(HDF5Particle, strain_zz),       // strain_zz
      offsetof(HDF5Particle, gamma_xy),        // gamma_xy
      offsetof(HDF5Particle, gamma_yz),        // gamma_yz
      offsetof(HDF5Particle, gamma_xz),        // gamma_xz
      offsetof(HDF5Particle, epsilon_v),       // epsv
      offsetof(HDF5Particle, cell_id),         // cellid
      offsetof(HDF5Particle, status),          // status
      offsetof(HDF5Particle, material_id),     // material_id
      offsetof(HDF5Particle, nstate_vars),     // nstate_vars
      offsetof(HDF5Particle, svars)            // state_vars
  };

  // Array containing the MPI datatypes to replicate to make each block.
  MPI_Datatype types[nblocks] = {
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
  };

  // Create the Particle datatype
  MPI_Datatype MPIParticle;

  // Create particle data types
  MPI_Type_create_struct(nblocks, lengths, displacements, types, &MPIParticle);
  MPI_Type_commit(&MPIParticle);
  return MPIParticle;
}

//! Deregister MPI particle data type
inline void deregister_mpi_particle_type(MPI_Datatype& particle_type) {
  MPI_Type_free(&particle_type);
}
}  // namespace mpm

#endif  // MPI

#endif  // MPM_MPI_HDF5_PARTICLE_H_

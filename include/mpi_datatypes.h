#ifndef MPM_MPI_HDF5_PARTICLE_H_
#define MPM_MPI_HDF5_PARTICLE_H_

// MPI
#ifdef USE_MPI
#include "mpi.h"

#include "hdf5_particle.h"

namespace mpm {
// Create the Particle datatype
MPI_Datatype MPIParticle;

void init_mpi_particle_datatypes() {
  // Number of blocks to create
  const unsigned nblocks = 31;
  // Array containing the length of each block
  int lengths[nblocks] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  // Array containing the displacement for each block, expressed in bytes. The
  // displacement is the distance between the start of the MPI datatype created
  // and the start of the block
  // clang-format off
  const MPI_Aint displacements[nblocks] = {
      0,                                                 // id
      sizeof(unsigned long long),                        // mass
      sizeof(unsigned long long) +  1 * sizeof(double),  // volume
      sizeof(unsigned long long) +  2 * sizeof(double),  // pressure
      sizeof(unsigned long long) +  3 * sizeof(double),  // coord_x
      sizeof(unsigned long long) +  4 * sizeof(double),  // coord_y
      sizeof(unsigned long long) +  5 * sizeof(double),  // coord_z
      sizeof(unsigned long long) +  6 * sizeof(double),  // disp_x
      sizeof(unsigned long long) +  7 * sizeof(double),  // disp_y
      sizeof(unsigned long long) +  8 * sizeof(double),  // disp_z
      sizeof(unsigned long long) +  9 * sizeof(double),  // nsize_x
      sizeof(unsigned long long) + 10 * sizeof(double),  // nsize_y
      sizeof(unsigned long long) + 11 * sizeof(double),  // nsize_z
      sizeof(unsigned long long) + 12 * sizeof(double),  // vel_x
      sizeof(unsigned long long) + 13 * sizeof(double),  // vel_y
      sizeof(unsigned long long) + 14 * sizeof(double),  // vel_z
      sizeof(unsigned long long) + 15 * sizeof(double),  // stress_xx
      sizeof(unsigned long long) + 16 * sizeof(double),  // stress_yy
      sizeof(unsigned long long) + 17 * sizeof(double),  // stress_zz
      sizeof(unsigned long long) + 18 * sizeof(double),  // tau_xx
      sizeof(unsigned long long) + 19 * sizeof(double),  // tau_yy
      sizeof(unsigned long long) + 20 * sizeof(double),  // tau_zz
      sizeof(unsigned long long) + 21 * sizeof(double),  // strain_xx
      sizeof(unsigned long long) + 22 * sizeof(double),  // strain_yy
      sizeof(unsigned long long) + 23 * sizeof(double),  // strain_zz
      sizeof(unsigned long long) + 24 * sizeof(double),  // gamma_xy
      sizeof(unsigned long long) + 25 * sizeof(double),  // gamma_yz
      sizeof(unsigned long long) + 26 * sizeof(double),  // gamma_zx
      sizeof(unsigned long long) + 27 * sizeof(double),  // epsv
      sizeof(unsigned long long) + 28 * sizeof(double),  // status
      sizeof(unsigned long long) + 28 * sizeof(double) +
      sizeof(int)                                        // cellid
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
    MPI_INT,                 // status
    MPI_UNSIGNED_LONG_LONG   // cell_id
  };
  // clang-format on
  // Create particle data types
  MPI_Type_create_struct(nblocks, lengths, displacements, types, &MPIParticle);
  MPI_Type_commit(&MPIParticle);
}
}  // namespace mpm

#endif  // MPI

#endif  // MPM_MPI_HDF5_PARTICLE_H_

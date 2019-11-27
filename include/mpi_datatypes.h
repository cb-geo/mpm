#ifndef MPM_MPI_HDF5_PARTICLE_H_
#define MPM_MPI_HDF5_PARTICLE_H_

#include <cstddef>

// MPI
#ifdef USE_MPI
#include "mpi.h"

#include "hdf5_particle.h"

namespace mpm {
// Create the Particle datatype
MPI_Datatype MPIParticle;

//! Initialize MPI particle data types
void init_mpi_particle_datatypes() {
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
      0,                                                 // id
      size_ull,                                          // mass
      size_ull + 1 * sizeof(double),                     // volume
      size_ull + 2 * sizeof(double),                     // pressure
      size_ull + 3 * sizeof(double),                     // coord_x
      size_ull + 4 * sizeof(double),                     // coord_y
      size_ull + 5 * sizeof(double),                     // coord_z
      size_ull + 6 * sizeof(double),                     // disp_x
      size_ull + 7 * sizeof(double),                     // disp_y
      size_ull + 8 * sizeof(double),                     // disp_z
      size_ull + 9 * sizeof(double),                     // nsize_x
      size_ull + 10 * sizeof(double),                    // nsize_y
      size_ull + 11 * sizeof(double),                    // nsize_z
      size_ull + 12 * sizeof(double),                    // vel_x
      size_ull + 13 * sizeof(double),                    // vel_y
      size_ull + 14 * sizeof(double),                    // vel_z
      size_ull + 15 * sizeof(double),                    // stress_xx
      size_ull + 16 * sizeof(double),                    // stress_yy
      size_ull + 17 * sizeof(double),                    // stress_zz
      size_ull + 18 * sizeof(double),                    // tau_xx
      size_ull + 19 * sizeof(double),                    // tau_yy
      size_ull + 20 * sizeof(double),                    // tau_zz
      size_ull + 21 * sizeof(double),                    // strain_xx
      size_ull + 22 * sizeof(double),                    // strain_yy
      size_ull + 23 * sizeof(double),                    // strain_zz
      size_ull + 24 * sizeof(double),                    // gamma_xy
      size_ull + 25 * sizeof(double),                    // gamma_yz
      size_ull + 26 * sizeof(double),                    // gamma_zx
      size_ull + 27 * sizeof(double),                    // epsv
      size_ull + 28 * sizeof(double),                    // cellid
      2 * size_ull + 28 * sizeof(double),                // status
      2 * size_ull + 28 * sizeof(double) + sizeof(int),  // material_id
      2 * size_ull + 28 * sizeof(double) + sizeof(int) +
          sizeof(unsigned),  // nstate_vars
      2 * size_ull + 28 * sizeof(double) + sizeof(int) +
          2 * sizeof(unsigned)  // state_vars
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

  // Create particle data types
  MPI_Type_create_struct(nblocks, lengths, displacements, types, &MPIParticle);
  MPI_Type_commit(&MPIParticle);
}

//! Free MPI particle data types
void free_mpi_particle_datatypes() { MPI_Type_free(&MPIParticle); }
}  // namespace mpm

#endif  // MPI

#endif  // MPM_MPI_HDF5_PARTICLE_H_

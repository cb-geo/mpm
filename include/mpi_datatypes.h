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
  const unsigned nblocks = 3;
  // Array containing the length of each block
  int lengths[nblocks] = {1, 1, 1};
  // Array containing the displacement for each block, expressed in bytes. The
  // displacement is the distance between the start of the MPI datatype created
  // and the start of the block
  const MPI_Aint displacements[nblocks] = {
      0,                           // id
      sizeof(unsigned long long),  // mass
      sizeof(unsigned long long) + sizeof(double)};
  // Array containing the MPI datatypes to replicate to make each block.
  MPI_Datatype types[nblocks] = {MPI_UNSIGNED_LONG_LONG,  // id
                                 MPI_DOUBLE,              // mass
                                 MPI_DOUBLE};
  // Create particle data types
  MPI_Type_create_struct(nblocks, lengths, displacements, types, &MPIParticle);
  MPI_Type_commit(&MPIParticle);
}
}  // namespace mpm

#endif  // MPI

#endif  // MPM_MPI_HDF5_PARTICLE_H_

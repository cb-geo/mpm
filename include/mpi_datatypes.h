#ifndef MPM_MPI_HDF5_PARTICLE_H_
#define MPM_MPI_HDF5_PARTICLE_H_

// MPI
#ifdef USE_MPI
#include "mpi.h"

#include "hdf5_particle.h"

namespace mpm {
// Create the datatype
MPI_Datatype hdf5particle_type;

void init_mpi_datatypes() {
  int lengths[3] = {1, 1, 1};
  const MPI_Aint displacements[3] = {
      0, sizeof(unsigned long long),
      sizeof(unsigned long long) + sizeof(double)};
  MPI_Datatype types[3] = {MPI_UNSIGNED_LONG_LONG, MPI_DOUBLE, MPI_DOUBLE};
  MPI_Type_create_struct(3, lengths, displacements, types, &hdf5particle_type);
  MPI_Type_commit(&hdf5particle_type);
}
}  // namespace mpm

#endif  // MPI

#endif  // MPM_MPI_HDF5_PARTICLE_H_

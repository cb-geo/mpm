#ifndef MPM_MPI_WRAPPER_H_
#define MPM_MPI_WRAPPER_H_

#include <vector>

#include <Eigen/Dense>

#ifdef USE_MPI
#include "mpi.h"
#endif

namespace mpm {

//! Chunk quantities based on MPI ranks
//! \tparam Tnsize Size of quantity vector
//! \param[in] all_quantities Quantities to be chunked
//! \retval quantities Chunked quantities
template <int Tnsize>
void chunk_quantities(
    const std::vector<Eigen::Matrix<double, Tnsize, 1>>& all_quantities,
    std::vector<Eigen::Matrix<double, Tnsize, 1>>& quantities) {

#ifdef USE_MPI
  // Initialise MPI ranks and size
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  int mpi_size;
  // Get number of MPI ranks
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  // Create MPI array<double, Tnsize> type
  MPI_Datatype array_t;
  MPI_Type_vector(Tnsize, 1, 1, MPI_DOUBLE, &array_t);
  MPI_Type_commit(&array_t);

  // Calculate chunk size to split
  int chunk_size = all_quantities.size() / mpi_size;
  MPI_Bcast(&chunk_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  quantities.resize(chunk_size);

  // Send chunked quantities to different compute nodes
  MPI_Scatter(all_quantities.data(), chunk_size, array_t, quantities.data(),
              quantities.size(), array_t, 0, MPI_COMM_WORLD);

  // Calculate the remaining chunk of all_quantities to the last rank
  int chunk_remainder = all_quantities.size() % mpi_size;
  if (mpi_rank == (mpi_size - 1))
    quantities.insert(quantities.begin(),
                      all_quantities.end() - chunk_remainder,
                      all_quantities.end());

  MPI_Type_free(&array_t);
#else
  quantities = all_quantities;
#endif
}

}  // namespace mpm

#endif  // MPM_MPI_WRAPPER_H_

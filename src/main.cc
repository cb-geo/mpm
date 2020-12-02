#include <iostream>
#include <memory>

#ifdef USE_MPI
#include "mpi.h"
#endif
#ifdef USE_PETSC
#include <petscksp.h>
#endif
#include "spdlog/spdlog.h"

#include "io.h"
#include "mpm.h"

int main(int argc, char** argv) {

#ifdef USE_MPI
  // Initialise MPI
  MPI_Init(&argc, &argv);
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  // Get number of MPI ranks
  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  // Allocate enough space to issue the buffered send
  int mpi_buffer_size = 2000000000;
  void* mpi_buffer = malloc(mpi_buffer_size);
  // Pass the buffer allocated to MPI so it uses it when we issue MPI_Bsend
  MPI_Buffer_attach(mpi_buffer, mpi_buffer_size);

#endif

#ifdef USE_PETSC
  // Initialize PETSc
  PetscInitialize(&argc, &argv, 0, 0);
#endif

  try {
    // Logger level (trace, debug, info, warn, error, critical, off)
    spdlog::set_level(spdlog::level::trace);

    // Initialise logger
    auto console = spdlog::stdout_color_mt("main");

    // Create an IO object
    auto io = std::make_shared<mpm::IO>(argc, argv);

    // If number of threads are positive set to nthreads
    unsigned nthreads = io->nthreads();
#ifdef _OPENMP
    omp_set_num_threads(nthreads > 0 ? nthreads : omp_get_max_threads());
#endif

    // Get analysis type
    const std::string analysis = io->analysis_type();

    // Create an MPM analysis
    auto mpm =
        Factory<mpm::MPM, const std::shared_ptr<mpm::IO>&>::instance()->create(
            analysis, std::move(io));
    // Solve
    mpm->solve();

  } catch (std::exception& exception) {
    std::cerr << "MPM main: " << exception.what() << std::endl;

#ifdef USE_PETSC
    // Finalize PETSc
    PetscFinalize();
#endif

#ifdef USE_MPI
    free(mpi_buffer);
    MPI_Buffer_detach(&mpi_buffer, &mpi_buffer_size);
    MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    std::terminate();
  }

#ifdef USE_MPI
  free(mpi_buffer);
  MPI_Buffer_detach(&mpi_buffer, &mpi_buffer_size);
  MPI_Finalize();
#endif
}

#include <memory>

#ifdef USE_MPI
#include "mpi.h"
#endif
#include "spdlog/spdlog.h"
#include "tbb/task_scheduler_init.h"

#include "io.h"
#include "mpm.h"

int main(int argc, char** argv) {
  // Logger level (trace, debug, info, warn, error, critical, off)
  spdlog::set_level(spdlog::level::trace);

  // Initialise logger
  auto console = spdlog::stdout_color_mt("main");

#ifdef USE_MPI
  // Initialise MPI
  MPI_Init(&argc, &argv);
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  // Get number of MPI ranks
  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif

  try {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);

    // Set TBB threads
    unsigned nthreads = tbb::task_scheduler_init::default_num_threads();
    // If number of TBB threads are positive set to nthreads
    if (io->nthreads() > 0) nthreads = io->nthreads();
    tbb::task_scheduler_init init(nthreads);

    // Get analysis type
    const std::string analysis = io->analysis_type();

    // Create an MPM analysis
    auto mpm =
        Factory<mpm::MPM, std::unique_ptr<mpm::IO>&&>::instance()->create(
            analysis, std::move(io));
    // Solve
    mpm->solve();

  } catch (std::exception& exception) {
    console->error("MPM main: {}", exception.what());
  }

#ifdef USE_MPI
  MPI_Finalize();
#endif
}

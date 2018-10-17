#include <memory>

#include "spdlog/spdlog.h"
#include "mpi.h"

#include "io.h"
#include "mpm.h"

int main(int argc, char** argv) {
  // Logger level (trace, debug, info, warn, error, critical, off)
  spdlog::set_level(spdlog::level::trace);

  // Initialise logger
  auto console = spdlog::stdout_color_mt("main");

  // Initialise MPI
  MPI_Init(&argc, &argv);
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  // Get number of MPI ranks
  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  try {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);

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

  MPI_Finalize();
}

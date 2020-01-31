#define CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_RUNNER

#include <iostream>

#include "catch.hpp"
// MPI
#ifdef USE_MPI
#include "mpi.h"
#endif

int main(int argc, char** argv) {
  try {
    Catch::Session session;
    // Let Catch (using Clara) parse the command line
    int returnCode = session.applyCommandLine(argc, argv);
    if (returnCode != 0)  // Indicates a command line error
      return returnCode;

#ifdef USE_MPI
    // Initialise MPI
    MPI_Init(&argc, &argv);
#endif

    int result = session.run();

#ifdef USE_MPI
    MPI_Finalize();
#endif

    return result;
  } catch (std::exception& exception) {
    std::cerr << "Test: " << exception.what() << std::endl;
#ifdef USE_MPI
    MPI_Abort(MPI_COMM_WORLD, 1);
#endif
  }
}

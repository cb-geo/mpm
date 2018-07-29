#include <memory>

#include "spdlog/spdlog.h"

#include "io.h"
#include "mpm.h"
#include "vtk_writer.h"

int main(int argc, char** argv) {
  // Logger level (trace, debug, info, warn, error, critical, off)
  spdlog::set_level(spdlog::level::trace);

  // Initialise logger
  auto console = spdlog::stdout_color_mt("main");

  try {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);

    // Get analysis
    const std::string analysis = io->analysis_type();

    // Create an MPM analysis
    auto mpm =
        Factory<mpm::MPM, std::unique_ptr<mpm::IO>&&>::instance()->create(
            analysis, std::move(io));

    // Initialise mesh
    mpm->initialise_mesh_particles();

    // Initialise materials
    mpm->initialise_materials();

    // Solve
    mpm->solve();

  } catch (std::exception& exception) {
    console->error("MPM main: {}", exception.what());
  }
}

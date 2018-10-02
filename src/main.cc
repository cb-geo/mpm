#include <memory>

#include "spdlog/spdlog.h"

#include "io.h"
#include "mpm.h"

int main(int argc, char** argv) {
  // Logger level (trace, debug, info, warn, error, critical, off)
  spdlog::set_level(spdlog::level::trace);

  // Initialise logger
  auto console = spdlog::stdout_color_mt("main");

  try {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);

    // Get analysis type
    const std::string analysis = io->analysis_type();

    // Check if the specified analysis type is valid
    auto analyses =
        Factory<mpm::MPM, std::unique_ptr<mpm::IO>&&>::instance()->list();
    if (std::find(analyses.begin(), analyses.end(), analysis) == analyses.end())
      throw std::runtime_error("Specified MPM analysis type is invalid.");

    // Create an MPM analysis
    auto mpm =
        Factory<mpm::MPM, std::unique_ptr<mpm::IO>&&>::instance()->create(
            analysis, std::move(io));

    // Solve
    mpm->solve();

  } catch (std::exception& exception) {
    console->error("MPM main: {}", exception.what());
  }
}

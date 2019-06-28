#ifndef MPM_GENERATOR_FACTORY_H_
#define MPM_GENERATOR_FACTORY_H_

#include "point_generator.h"

#include "file_point_generator.h"
#include "gauss_point_generator.h"

// JSON
using Json = nlohmann::json;

namespace mpm {
//! Generator factory

template <unsigned Tdim>
std::vector<Eigen::Matrix<double, Tdim, 1>> generator_factory(
    const std::shared_ptr<mpm::Mesh<Tdim>>& mesh,
    const std::shared_ptr<mpm::IO>& io, const Json& generator) {
  // Logger
  auto console_ = spdlog::get("PointGenerator");
  // particles coordinates
  std::vector<Eigen::Matrix<double, Tdim, 1>> coordinates;
  try {
    // Particle generator
    const auto generator_type = generator["type"].template get<std::string>();

    // Generate particles from file
    if (generator_type == "file") {
      auto gen =
          std::make_shared<mpm::FilePointGenerator<Tdim>>(mesh, io, generator);
      coordinates = gen->generate_points();
    }
    // Generate material points at the Gauss location in all cells
    else if (generator_type == "gauss") {
      auto gen =
             std::make_shared<mpm::GaussPointGenerator<Tdim>>(mesh, io, generator);
      coordinates = gen->generate_points();
    } else
      throw std::runtime_error(
          "Particle generator type is not properly specified");

  } catch (std::exception& exception) {
    console_->error("{}: #{} Generating particle failed", __FILE__, __LINE__);
  }
  return coordinates;
}
}  // namespace mpm
#endif  // MPM_GENERATOR_FACTORY_H_

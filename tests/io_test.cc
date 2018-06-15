#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <utility>

#include "catch.hpp"

//! Alias for JSON
#include "json.hpp"
using Json = nlohmann::json;

#include "io.h"

//! \brief Check IO for input string
TEST_CASE("IO is checked for input parsing", "[IO][JSON]") {

  // Check input JSON
  SECTION("Check input JSON object") {
    //! Make json object with input files
    Json json_file = {
        {"title", "Example JSON Input for MPM"},
        {"input_files",
         {{"path", "../"},
          {"cmake", "CMakeLists.txt"},
          {"mesh", "mesh.dat"},
          {"submesh", "submesh.dat"},
          {"constraints", "mesh_constraints.dat"},
          {"soil_particles", "soil_particles.dat"},
          {"initial_stress", "initial_soil_stress.dat"},
          {"material", "material.dat"},
          {"traction", "traction.dat"}}},
        {"analysis",
         {{"dt", 0.001},
          {"number_steps", 1000},
          {"gravity", true},
          {"soil_particle_spacing", 0.01},
          {"boundary_friction", 0.5},
          {"damping", {{"damping", true}, {"damping_ratio", 0.02}}},
          {"newmark", {{"newmark", true}, {"gamma", 0.5}, {"beta", 0.25}}}}},
        {"post_processing", {{"path", "results/"}, {"output_steps", 10}}}};

    //! Create an IO object
    auto io = std::make_unique<IO>(json_file);

    //! Check if CMake and README files are present
    REQUIRE(io->check_file("../CMakeLists.txt") == true);
    REQUIRE(io->check_file("../README.md") == true);

    // Check if a non-existant file is present
    REQUIRE(io->check_file("../fail.txt") == false);
  }
}

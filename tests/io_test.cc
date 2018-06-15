#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <utility>

#include "catch.hpp"

//! Alias for JSON
#include "json.hpp"
using Json = nlohmann::json;

#include "io.h"

//! \brief Check IO for input string
TEST_CASE("IO is checked for input parsing", "[IO][JSON]"){

    // Check input JSON
    SECTION("Check input JSON object"){
        //! Make json object with input files
        Json json_file = {
          {"title", "Example JSON Input for MPM"},
          {"input_files",
           {{"path", "input/"},
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

          //! Make pointers to io
          // auto io = std::unique_ptr<IO>(new IO(json_file));

          //! Check for only string inputs
          // REQUIRE(io->inputMeshFileName() == "folder/inputFiles/mesh.smf");
        }
}

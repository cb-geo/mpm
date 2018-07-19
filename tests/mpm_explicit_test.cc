#include "catch.hpp"

//! Alias for JSON
#include "json.hpp"
using Json = nlohmann::json;

#include "mpm_explicit.h"

bool write_json(unsigned dim, const std::string& file_name) {
  // Make json object with input files
  std::string dimension = "2d";
  if (dim == 3) dimension = "3d";

  Json json_file = {
      {"title", "Example JSON Input for MPM"},
      {"input_files",
       {{"input", "mpm.json"},
        {"mesh", "mesh.dat"},
        {"constraints", "mesh_constraints.dat"},
        {"particles", "particles.dat"},
        {"initial_stresses", "initial_soil_stress.dat"},
        {"materials", "materials.dat"},
        {"traction", "traction.dat"}}},
      {"analysis",
       {{"dt", 0.001},
        {"nsteps", 1000},
        {"gravity", true},
        {"soil_particle_spacing", 0.01},
        {"boundary_friction", 0.5},
        {"damping", {{"damping", true}, {"damping_ratio", 0.02}}},
        {"newmark", {{"newmark", true}, {"gamma", 0.5}, {"beta", 0.25}}}}},
      {"post_processing", {{"path", "results/"}, {"output_steps", 10}}}};

  // Dump JSON as an input file to be read
  std::ofstream file;
  file.open((file_name + "-" + dimension + ".json").c_str());
  file << json_file.dump(2);
  file.close();

  return true;
}

// Check MPM Explicit
TEST_CASE("MPM 2D Explicit implementation is checked",
          "[MPM][2D][Explicit][1Phase]") {
  // Dimension
  const unsigned Dim = 2;

  // Write JSON file
  const std::string fname = "mpm-explicit";
  bool status = write_json(2, fname);
  REQUIRE(status == true);

  // Assign argc and argv to input arguments of MPM
  int argc = 5;
  char* argv[] = {(char*)"./mpm", (char*)"-f", (char*)"./", (char*)"-i",
                  (char*)"mpm-explicit-2d.json"};

  // Create an IO object
  auto io = std::make_unique<mpm::IO>(argc, argv);

  auto mpm = std::make_unique<mpm::MPMExplicit<Dim>>(std::move(io));
  REQUIRE(mpm->initialise() == true);
}

// Check MPM Explicit
TEST_CASE("MPM 3D Explicit implementation is checked",
          "[MPM][3D][Explicit][1Phase]") {
  // Dimension
  const unsigned Dim = 3;

  // Write JSON file
  const std::string fname = "mpm-explicit";
  bool status = write_json(3, fname);
  REQUIRE(status == true);

  // Assign argc and argv to input arguments of MPM
  int argc = 5;
  char* argv[] = {(char*)"./mpm", (char*)"-f", (char*)"./", (char*)"-i",
                  (char*)"mpm-explicit-3d.json"};

  // Create an IO object
  auto io = std::make_unique<mpm::IO>(argc, argv);

  auto mpm = std::make_unique<mpm::MPMExplicit<Dim>>(std::move(io));
  REQUIRE(mpm->initialise() == true);
}

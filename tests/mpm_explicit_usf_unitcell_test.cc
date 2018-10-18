#include "catch.hpp"

//! Alias for JSON
#include "json.hpp"
using Json = nlohmann::json;

#include "mpm_explicit_usf.h"
#include "write_mesh_particles_unitcell.h"

// Check MPM Explicit USF
TEST_CASE("MPM 2D Explicit USF implementation is checked in unitcells",
          "[MPM][2D][Explicit][1Phase][unitcell]") {
  // Dimension
  const unsigned Dim = 2;

  // Write JSON file
  const std::string fname = "mpm-explicit-usf";
  const std::string analysis = "MPMExplicitUSF2D";
  REQUIRE(mpm_test::write_json_unitcell(2, analysis, fname) == true);

  // Write Mesh
  REQUIRE(mpm_test::write_mesh_2d_unitcell() == true);

  // Write Particles
  REQUIRE(mpm_test::write_particles_2d_unitcell() == true);

  // Assign argc and argv to input arguments of MPM
  int argc = 5;
  // clang-format off
  char* argv[] = {(char*)"./mpm",
                  (char*)"-f",  (char*)"./",
                  (char*)"-i",  (char*)"mpm-explicit-usf-2d-unitcell.json"};
  // clang-format on

  SECTION("Check initialisation") {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run explicit MPM
    auto mpm = std::make_unique<mpm::MPMExplicitUSF<Dim>>(std::move(io));

    // Initialise mesh
    REQUIRE(mpm->initialise_mesh_particles() == true);

    // Initialise materials
    REQUIRE(mpm->initialise_materials() == true);

    // Reinitialise mesh
    REQUIRE(mpm->initialise_mesh_particles() == false);

    // Renitialise materials
    REQUIRE(mpm->initialise_materials() == false);
  }

  SECTION("Check solver") {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run explicit MPM
    auto mpm = std::make_unique<mpm::MPMExplicitUSF<Dim>>(std::move(io));
    // Solve
    REQUIRE(mpm->solve() == true);
  }
}

// Check MPM Explicit
TEST_CASE("MPM 3D Explicit USF implementation is checked in unitcells",
          "[MPM][3D][Explicit][USF][1Phase][unitcell]") {
  // Dimension
  const unsigned Dim = 3;

  // Write JSON file
  const std::string fname = "mpm-explicit-usf";
  const std::string analysis = "MPMExplicitUSF3D";
  REQUIRE(mpm_test::write_json_unitcell(3, analysis, fname) == true);

  // Write Mesh
  REQUIRE(mpm_test::write_mesh_3d_unitcell() == true);

  // Write Particles
  REQUIRE(mpm_test::write_particles_3d_unitcell() == true);

  // Assign argc and argv to input arguments of MPM
  int argc = 5;
  // clang-format off
  char* argv[] = {(char*)"./mpm",
                  (char*)"-f",  (char*)"./",
                  (char*)"-i",  (char*)"mpm-explicit-usf-3d-unitcell.json"};
  // clang-format on

  SECTION("Check initialisation") {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run explicit MPM
    auto mpm = std::make_unique<mpm::MPMExplicitUSF<Dim>>(std::move(io));

    // Initialise mesh
    REQUIRE(mpm->initialise_mesh_particles() == true);

    // Initialise materials
    REQUIRE(mpm->initialise_materials() == true);

    // Reinitialise mesh
    REQUIRE(mpm->initialise_mesh_particles() == false);

    // Renitialise materials
    REQUIRE(mpm->initialise_materials() == false);
  }

  SECTION("Check solver") {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run explicit MPM
    auto mpm = std::make_unique<mpm::MPMExplicitUSF<Dim>>(std::move(io));
    // Solve
    REQUIRE(mpm->solve() == true);
  }
}

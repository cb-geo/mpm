#include "catch.hpp"

//! Alias for JSON
#include "json.hpp"
using Json = nlohmann::json;

#include "mpm_explicit_usl.h"
#include "write_mesh_particles.h"

// Check MPM Explicit USL
TEST_CASE("MPM 2D Explicit USL implementation is checked",
          "[MPM][2D][Explicit][USL][1Phase]") {
  // Dimension
  const unsigned Dim = 2;

  // Write JSON file
  const std::string fname = "mpm-explicit-usl";
  const std::string analysis = "MPMExplicitUSL2D";
  const bool resume = false;
  REQUIRE(mpm_test::write_json(2, resume, analysis, fname) == true);

  // Write Mesh
  REQUIRE(mpm_test::write_mesh_2d() == true);

  // Write Particles
  REQUIRE(mpm_test::write_particles_2d() == true);

  // Assign argc and argv to input arguments of MPM
  int argc = 5;
  // clang-format off
  char* argv[] = {(char*)"./mpm",
                  (char*)"-f",  (char*)"./",
                  (char*)"-i",  (char*)"mpm-explicit-usl-2d.json"};
  // clang-format on

  SECTION("Check initialisation") {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run explicit MPM
    auto mpm = std::make_unique<mpm::MPMExplicitUSL<Dim>>(std::move(io));

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
    auto mpm = std::make_unique<mpm::MPMExplicitUSL<Dim>>(std::move(io));
    // Solve
    REQUIRE(mpm->solve() == true);
    // Test check point restart
    REQUIRE(mpm->checkpoint_resume() == false);
  }

  SECTION("Check resume") {
    // Write JSON file
    const std::string fname = "mpm-explicit-usl";
    const std::string analysis = "MPMExplicitUSL2D";
    bool resume = true;
    REQUIRE(mpm_test::write_json(2, resume, analysis, fname) == true);

    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run explicit MPM
    auto mpm = std::make_unique<mpm::MPMExplicitUSL<Dim>>(std::move(io));

    // Test check point restart
    REQUIRE(mpm->checkpoint_resume() == true);
    // Solve
    REQUIRE(mpm->solve() == true);
  }
}

// Check MPM Explicit USL
TEST_CASE("MPM 3D Explicit USL implementation is checked",
          "[MPM][3D][Explicit][USL][1Phase]") {
  // Dimension
  const unsigned Dim = 3;

  // Write JSON file
  const std::string fname = "mpm-explicit-usl";
  const std::string analysis = "MPMExplicitUSL3D";
  const bool resume = false;
  REQUIRE(mpm_test::write_json(3, resume, analysis, fname) == true);

  // Write Mesh
  REQUIRE(mpm_test::write_mesh_3d() == true);

  // Write Particles
  REQUIRE(mpm_test::write_particles_3d() == true);

  // Assign argc and argv to input arguments of MPM
  int argc = 5;
  // clang-format off
  char* argv[] = {(char*)"./mpm",
                  (char*)"-f",  (char*)"./",
                  (char*)"-i",  (char*)"mpm-explicit-usl-3d.json"};
  // clang-format on

  SECTION("Check initialisation") {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run explicit MPM
    auto mpm = std::make_unique<mpm::MPMExplicitUSL<Dim>>(std::move(io));

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
    auto mpm = std::make_unique<mpm::MPMExplicitUSL<Dim>>(std::move(io));
    // Solve
    REQUIRE(mpm->solve() == true);
    // Test check point restart
    REQUIRE(mpm->checkpoint_resume() == false);
  }

  SECTION("Check resume") {
    // Write JSON file
    const std::string fname = "mpm-explicit-usl";
    const std::string analysis = "MPMExplicitUSL3D";
    bool resume = true;
    REQUIRE(mpm_test::write_json(3, resume, analysis, fname) == true);

    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run explicit MPM
    auto mpm = std::make_unique<mpm::MPMExplicitUSL<Dim>>(std::move(io));

    // Test check point restart
    REQUIRE(mpm->checkpoint_resume() == true);
    // Solve
    REQUIRE(mpm->solve() == true);
  }
}

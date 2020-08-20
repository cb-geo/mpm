#include "catch.hpp"

//! Alias for JSON
#include "json.hpp"
using Json = nlohmann::json;

#include "mpm_explicit_twophase.h"
#include "write_mesh_particles.h"

// Check MPM Explicit
TEST_CASE("MPM 2D Explicit USL TwoPhase implementation is checked",
          "[MPM][2D][Explicit][USL][2Phase]") {
  // Dimension
  const unsigned Dim = 2;

  // Write JSON file
  const std::string fname = "mpm-explicit-twophase-usl";
  const std::string analysis = "MPMExplicitTwoPhase2D";
  const std::string mpm_scheme = "usl";
  bool resume = false;
  REQUIRE(mpm_test::write_json_twophase(2, resume, analysis, mpm_scheme,
                                        fname) == true);

  // Write JSON Entity Sets file
  REQUIRE(mpm_test::write_entity_set() == true);

  // Write Mesh
  REQUIRE(mpm_test::write_mesh_2d() == true);

  // Write Particles
  REQUIRE(mpm_test::write_particles_2d() == true);

  // Assign argc and argv to input arguments of MPM
  int argc = 5;
  // clang-format off
  char* argv[] = {(char*)"./mpm",
                  (char*)"-f",  (char*)"./",
                  (char*)"-i",  (char*)"mpm-explicit-twophase-usl-2d.json"};
  // clang-format on

  SECTION("Check initialisation") {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run explicit MPM
    auto mpm = std::make_unique<mpm::MPMExplicitTwoPhase<Dim>>(std::move(io));

    // Initialise materials
    REQUIRE_NOTHROW(mpm->initialise_materials());
    // Initialise mesh
    REQUIRE_NOTHROW(mpm->initialise_mesh());
    // Initialise particles
    REQUIRE_NOTHROW(mpm->initialise_particles());

    // Initialise external loading
    REQUIRE_NOTHROW(mpm->initialise_loads());

    // Renitialise materials
    REQUIRE_THROWS(mpm->initialise_materials());
  }

  SECTION("Check solver") {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run explicit MPM
    auto mpm = std::make_unique<mpm::MPMExplicitTwoPhase<Dim>>(std::move(io));
    // Solve
    REQUIRE(mpm->solve() == true);
    // Test check point restart
    REQUIRE(mpm->checkpoint_resume() == false);
  }

  // SECTION("Check resume") {
  //   // Write JSON file
  //   const std::string fname = "mpm-explicit-twophase-usl";
  //   const std::string analysis = "MPMExplicitTwoPhase2D";
  //   const std::string mpm_scheme = "usl";
  //   bool resume = true;
  //   REQUIRE(mpm_test::write_json_twophase(2, resume, analysis, mpm_scheme,
  //                                         fname) == true);

  //   // Create an IO object
  //   auto io = std::make_unique<mpm::IO>(argc, argv);
  //   // Run explicit MPM
  //   auto mpm =
  //   std::make_unique<mpm::MPMExplicitTwoPhase<Dim>>(std::move(io));

  //   // Test check point restart
  //   REQUIRE(mpm->checkpoint_resume() == true);
  //   // Solve
  //   REQUIRE(mpm->solve() == true);
  // }

  SECTION("Check pressure smoothing") {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run explicit MPM
    auto mpm = std::make_unique<mpm::MPMExplicitTwoPhase<Dim>>(std::move(io));
    // Pressure smoothing
    REQUIRE_NOTHROW(mpm->pressure_smoothing(mpm::ParticlePhase::Solid));
    REQUIRE_NOTHROW(mpm->pressure_smoothing(mpm::ParticlePhase::Liquid));
  }
}

// Check MPM Explicit
TEST_CASE("MPM 3D Explicit USL TwoPhase implementation is checked",
          "[MPM][3D][Explicit][USL][2Phase]") {
  // Dimension
  const unsigned Dim = 3;

  // Write JSON file
  const std::string fname = "mpm-explicit-twophase-usl";
  const std::string analysis = "MPMExplicitTwoPhase3D";
  const std::string mpm_scheme = "usl";
  const bool resume = false;
  REQUIRE(mpm_test::write_json_twophase(3, resume, analysis, mpm_scheme,
                                        fname) == true);

  // Write JSON Entity Sets file
  REQUIRE(mpm_test::write_entity_set() == true);

  // Write Mesh
  REQUIRE(mpm_test::write_mesh_3d() == true);

  // Write Particles
  REQUIRE(mpm_test::write_particles_3d() == true);

  // Assign argc and argv to input arguments of MPM
  int argc = 5;
  // clang-format off
  char* argv[] = {(char*)"./mpm",
                  (char*)"-f",  (char*)"./",
                  (char*)"-i",  (char*)"mpm-explicit-twophase-usl-3d.json"};
  // clang-format on

  SECTION("Check initialisation") {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run explicit MPM
    auto mpm = std::make_unique<mpm::MPMExplicitTwoPhase<Dim>>(std::move(io));

    // Initialise materials
    REQUIRE_NOTHROW(mpm->initialise_materials());
    // Initialise mesh
    REQUIRE_NOTHROW(mpm->initialise_mesh());
    // Initialise particles
    REQUIRE_NOTHROW(mpm->initialise_particles());

    // Renitialise materials
    REQUIRE_THROWS(mpm->initialise_materials());
  }

  SECTION("Check solver") {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run explicit MPM
    auto mpm = std::make_unique<mpm::MPMExplicitTwoPhase<Dim>>(std::move(io));
    // Solve
    REQUIRE(mpm->solve() == true);
    // Test check point restart
    REQUIRE(mpm->checkpoint_resume() == false);
  }

  // SECTION("Check resume") {
  //   // Write JSON file
  //   const std::string fname = "mpm-explicit-twophase-usl";
  //   const std::string analysis = "MPMExplicitTwoPhase3D";
  //   const std::string mpm_scheme = "usl";
  //   bool resume = true;
  //   REQUIRE(mpm_test::write_json_twophase(3, resume, analysis, mpm_scheme,
  //                                         fname) == true);

  //   // Create an IO object
  //   auto io = std::make_unique<mpm::IO>(argc, argv);
  //   // Run explicit MPM
  //   auto mpm =
  //   std::make_unique<mpm::MPMExplicitTwoPhase<Dim>>(std::move(io));

  //   // Test check point restart
  //   REQUIRE(mpm->checkpoint_resume() == true);
  //   // Solve
  //   REQUIRE(mpm->solve() == true);
  // }

  SECTION("Check pressure smoothing") {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run explicit MPM
    auto mpm = std::make_unique<mpm::MPMExplicitTwoPhase<Dim>>(std::move(io));
    // Pressure smoothing
    REQUIRE_NOTHROW(mpm->pressure_smoothing(mpm::ParticlePhase::Solid));
    REQUIRE_NOTHROW(mpm->pressure_smoothing(mpm::ParticlePhase::Liquid));
  }
}

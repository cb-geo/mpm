#include "catch.hpp"

//! Alias for JSON
#include "json.hpp"
using Json = nlohmann::json;

#include "mpm_explicit.h"
#include "write_mesh_particles.h"

// Check MPM Explicit
TEST_CASE("MPM 2D Explicit implementation is checked",
          "[MPM][2D][Explicit][1Phase]") {
  // Dimension
  const unsigned Dim = 2;

  // Write JSON file
  const std::string fname = "mpm-explicit";
  bool status = mpm_test::write_json(2, fname);
  REQUIRE(status == true);

  // Write Mesh
  bool mesh_status = mpm_test::write_mesh_2d();
  REQUIRE(mesh_status == true);

  // Write Particles
  bool particle_status = mpm_test::write_particles_2d();
  REQUIRE(particle_status == true);

  // Assign argc and argv to input arguments of MPM
  int argc = 7;
  // clang-format off
  char* argv[] = {(char*)"./mpm",
                  (char*)"-a",  (char*)"MPMExplicit2D",
                  (char*)"-f",  (char*)"./",
                  (char*)"-i",  (char*)"mpm-explicit-2d.json"};
  // clang-format on

  SECTION("Check initialisation") {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run explicit MPM
    auto mpm = std::make_unique<mpm::MPMExplicit<Dim>>(std::move(io));

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
    auto mpm = std::make_unique<mpm::MPMExplicit<Dim>>(std::move(io));
    // Solve
    REQUIRE(mpm->solve() == true);
    // Test check point restart
    REQUIRE(mpm->checkpoint_resume() == false);
  }
}

// Check MPM Explicit
TEST_CASE("MPM 3D Explicit implementation is checked",
          "[MPM][3D][Explicit][1Phase]") {
  // Dimension
  const unsigned Dim = 3;

  // Write JSON file
  const std::string fname = "mpm-explicit";
  bool status = mpm_test::write_json(3, fname);
  REQUIRE(status == true);

  // Write Mesh
  bool mesh_status = mpm_test::write_mesh_3d();
  REQUIRE(mesh_status == true);

  // Write Particles
  bool particle_status = mpm_test::write_particles_3d();
  REQUIRE(particle_status == true);

  // Assign argc and argv to input arguments of MPM
  int argc = 7;
  // clang-format off
  char* argv[] = {(char*)"./mpm",
                  (char*)"-a",  (char*)"MPMExplicit3D",
                  (char*)"-f",  (char*)"./",
                  (char*)"-i",  (char*)"mpm-explicit-3d.json"};
  // clang-format on

  SECTION("Check initialisation") {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run explicit MPM
    auto mpm = std::make_unique<mpm::MPMExplicit<Dim>>(std::move(io));

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
    auto mpm = std::make_unique<mpm::MPMExplicit<Dim>>(std::move(io));
    // Solve
    REQUIRE(mpm->solve() == true);
    // Test check point restart
    REQUIRE(mpm->checkpoint_resume() == false);
  }
}

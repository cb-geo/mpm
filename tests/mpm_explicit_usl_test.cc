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
  const bool resume = false;
  bool status = mpm_test::write_json(2, resume, fname);
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
                  (char*)"-a",  (char*)"MPMExplicitUSL2D",
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
    bool resume = true;
    bool status = mpm_test::write_json(2, resume, fname);

    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run explicit MPM
    auto mpm = std::make_unique<mpm::MPMExplicit<Dim>>(std::move(io));

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
  const bool resume = false;
  bool status = mpm_test::write_json(3, resume, fname);
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
                  (char*)"-a",  (char*)"MPMExplicitUSL3D",
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
    bool resume = true;
    bool status = mpm_test::write_json(3, resume, fname);

    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run explicit MPM
    auto mpm = std::make_unique<mpm::MPMExplicit<Dim>>(std::move(io));

    // Test check point restart
    REQUIRE(mpm->checkpoint_resume() == true);
    // Solve
    REQUIRE(mpm->solve() == true);
  }
}

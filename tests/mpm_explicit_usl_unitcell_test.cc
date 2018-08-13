#include "catch.hpp"

//! Alias for JSON
#include "json.hpp"
using Json = nlohmann::json;

#include "mpm_explicit_usl.h"
#include "write_mesh_particles_unitcell.h"

// Check MPM Explicit USL
TEST_CASE("MPM 2D Explicit USL implementation is checked in unitcells",
          "[MPM][2D][Explicit][1Phase][unitcell]") {
  // Dimension
  const unsigned Dim = 2;

  // Write JSON file
  const std::string fname = "mpm-explicit-usl";
  bool status = mpm_test::write_json_unitcell(2, fname);
  REQUIRE(status == true);

  // Write Mesh
  bool mesh_status = mpm_test::write_mesh_2d_unitcell();
  REQUIRE(mesh_status == true);

  // Write Particles
  bool particle_status = mpm_test::write_particles_2d_unitcell();
  REQUIRE(particle_status == true);

  // Assign argc and argv to input arguments of MPM
  int argc = 7;
  // clang-format off
  char* argv[] = {(char*)"./mpm",
                  (char*)"-a",  (char*)"MPMExplicitUSL2D",
                  (char*)"-f",  (char*)"./",
                  (char*)"-i",  (char*)"mpm-explicit-usl-2d-unitcell.json"};
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
  }
}

// Check MPM Explicit
TEST_CASE("MPM 3D Explicit USL implementation is checked in unitcells",
          "[MPM][3D][Explicit][USL][1Phase][unitcell]") {
  // Dimension
  const unsigned Dim = 3;

  // Write JSON file
  const std::string fname = "mpm-explicit-usl";
  bool status = mpm_test::write_json_unitcell(3, fname);
  REQUIRE(status == true);

  // Write Mesh
  bool mesh_status = mpm_test::write_mesh_3d_unitcell();
  REQUIRE(mesh_status == true);

  // Write Particles
  bool particle_status = mpm_test::write_particles_3d_unitcell();
  REQUIRE(particle_status == true);

  // Assign argc and argv to input arguments of MPM
  int argc = 7;
  // clang-format off
  char* argv[] = {(char*)"./mpm",
                  (char*)"-a",  (char*)"MPMExplicitUSL3D",
                  (char*)"-f",  (char*)"./",
                  (char*)"-i",  (char*)"mpm-explicit-usl-3d-unitcell.json"};
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
  }
}

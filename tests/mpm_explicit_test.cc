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
  int argc = 5;
  char* argv[] = {(char*)"./mpm", (char*)"-f", (char*)"./", (char*)"-i",
                  (char*)"mpm-explicit-2d.json"};

  // Create an IO object
  auto io = std::make_unique<mpm::IO>(argc, argv);

  // Run explicit MPM
  auto mpm = std::make_unique<mpm::MPMExplicit<Dim>>(std::move(io));
  // Initialise mesh
  REQUIRE(mpm->initialise_mesh_particles() == true);
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
  int argc = 5;
  char* argv[] = {(char*)"./mpm", (char*)"-f", (char*)"./", (char*)"-i",
                  (char*)"mpm-explicit-3d.json"};

  // Create an IO object
  auto io = std::make_unique<mpm::IO>(argc, argv);

  // Run explicit MPM
  auto mpm = std::make_unique<mpm::MPMExplicit<Dim>>(std::move(io));
  // Initialise mesh
  REQUIRE(mpm->initialise_mesh_particles() == true);
}

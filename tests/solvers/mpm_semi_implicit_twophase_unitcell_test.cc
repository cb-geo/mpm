#include "catch.hpp"

//! Alias for JSON
#include "json.hpp"
using Json = nlohmann::json;

#include "mpm_semi_implicit_twophase.h"
#include "write_mesh_particles_unitcell.h"

// Check MPM Semi-implicit TwoPhase
TEST_CASE(
    "MPM 2D Semi-implicit TwoPhase implementation is checked in unitcells",
    "[MPM][2D][Semi-implicit][2Phase][unitcell]") {
  // Dimension
  const unsigned Dim = 2;

  // Write JSON file
  const std::string fname = "mpm-explicit-twophase";
  const std::string analysis = "MPMSemiImplicitTwoPhase2D";
  const std::string mpm_scheme = "usf";
  const std::string fsd_type = "geometry";
  const std::string lin_solver_type = "IterativeEigen";
  REQUIRE(mpm_test::write_json_unitcell_twophase(2, analysis, mpm_scheme, fname,
                                                 fsd_type,
                                                 lin_solver_type) == true);

  // Write Mesh
  REQUIRE(mpm_test::write_mesh_2d_unitcell() == true);

  // Write Particles
  REQUIRE(mpm_test::write_particles_2d_unitcell() == true);

  // Assign argc and argv to input arguments of MPM
  int argc = 5;
  // clang-format off
  char* argv[] = {(char*)"./mpm",
                  (char*)"-f",  (char*)"./",
                  (char*)"-i",  (char*)"mpm-explicit-twophase-2d-unitcell.json"};
  // clang-format on

  SECTION("Check initialisation") {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run explicit MPM
    auto mpm =
        std::make_unique<mpm::MPMSemiImplicitTwoPhase<Dim>>(std::move(io));

    // Initialise materials
    REQUIRE_NOTHROW(mpm->initialise_materials());

    // Initialise mesh and particles
    REQUIRE_NOTHROW(mpm->initialise_mesh());
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
    auto mpm =
        std::make_unique<mpm::MPMSemiImplicitTwoPhase<Dim>>(std::move(io));
    // Solve
    REQUIRE(mpm->solve() == true);
  }
}

// Check MPM Semi Implicit
TEST_CASE(
    "MPM 3D Semi-implicit TwoPhase implementation is checked in unitcells",
    "[MPM][3D][Semi-implicit][2Phase][unitcell]") {
  // Dimension
  const unsigned Dim = 3;

  // Write JSON file
  const std::string fname = "mpm-explicit-twophase";
  const std::string analysis = "MPMSemiImplicitTwoPhase3D";
  const std::string mpm_scheme = "usf";
  const std::string fsd_type = "geometry";
  const std::string lin_solver_type = "IterativeEigen";
  REQUIRE(mpm_test::write_json_unitcell_twophase(3, analysis, mpm_scheme, fname,
                                                 fsd_type,
                                                 lin_solver_type) == true);

  // Write Mesh
  REQUIRE(mpm_test::write_mesh_3d_unitcell() == true);

  // Write Particles
  REQUIRE(mpm_test::write_particles_3d_unitcell() == true);

  // Assign argc and argv to input arguments of MPM
  int argc = 5;
  // clang-format off
  char* argv[] = {(char*)"./mpm",
                  (char*)"-f",  (char*)"./",
                  (char*)"-i",  (char*)"mpm-explicit-twophase-3d-unitcell.json"};
  // clang-format on

  SECTION("Check initialisation") {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run explicit MPM
    auto mpm =
        std::make_unique<mpm::MPMSemiImplicitTwoPhase<Dim>>(std::move(io));

    // Initialise materials
    REQUIRE_NOTHROW(mpm->initialise_materials());

    // Initialise mesh and particles
    REQUIRE_NOTHROW(mpm->initialise_mesh());
    REQUIRE_NOTHROW(mpm->initialise_particles());

    // Renitialise materials
    REQUIRE_THROWS(mpm->initialise_materials());
  }

  SECTION("Check solver") {
    // Create an IO object
    auto io = std::make_unique<mpm::IO>(argc, argv);
    // Run explicit MPM
    auto mpm =
        std::make_unique<mpm::MPMSemiImplicitTwoPhase<Dim>>(std::move(io));
    // Solve
    REQUIRE(mpm->solve() == true);
  }
}
